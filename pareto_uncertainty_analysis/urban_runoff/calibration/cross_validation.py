"""
Cross-validation calibration engine for the UrbanRunoffModeling pipeline.

Implements the full calibration loop from legacy notebook 01 (cell 12,
run_calibration) without global state:

  build_factor_combinations()
      Generates the 480-row factor grid using itertools.product on
      PARAMETER_GRID in the legacy iteration order.

  factors_row_to_dict(row, param_names)
      Converts one row of the factor array to a named dict.

  run_calibration(cv_dir, factor_combinations, param_names, base_inp_filename)
      Runs all combinations for all events.  Returns per-date DataFrame
      lists and the combined observed hydrograph, matching the legacy
      output structure consumed by compute_objective_functions().

  merge_calibration_results(df_list)
      Merges the list of per-date DataFrames on the shared 'Factors' columns,
      replicating the legacy reduce(merge_func, cali_results_df_l) step.

Global-state free: all accumulation is done in local variables.  No global
list anti-pattern.  The caller is responsible for merging per-date lists.

DataFrame column conventions (MultiIndex) match the legacy exactly:
  Factor labels : ('Factors', 'IMP'), ('Factors', 'Width'), ...
  Date results  : (date, 'Total Volume'), (date, 'Max Runoff'), ...
  Hydrograph    : (date, timestamp_string), ...
"""

import itertools
import logging
import os
from functools import reduce
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union  # Tuple used by filter_bounds

import numpy as np
import pandas as pd
from swmm_api import read_inp_file

from urban_runoff.config import PARAMETER_GRID, SWMM_BASE_INP_FILENAME, FILTER_BOUNDS
from urban_runoff.data.loaders import load_observed_runoff, load_swmm_output
from urban_runoff.data.preprocessors import align_obs_sim, resample_5min
from urban_runoff.swmm.inp_editor import apply_all_factors
from urban_runoff.swmm.output_parser import compute_storm_stats
from urban_runoff.swmm.runner import run_simulation

logger = logging.getLogger(__name__)

# Maps canonical config parameter names → legacy MultiIndex column labels
_PARAM_TO_LABEL: Dict[str, str] = {
    "imperviousness": "IMP",
    "storage":        "Storage",
    "width":          "Width",
    "n":              "N",
    "pct_zero":       "PCT_ZERO",
    "cn":             "CN",
    "pct_routed":     "PCT_ROUTED",
    "evaporation":    "EVAP",
}


def validate_calibration_cache(
    cached: dict,
    factor_combinations: np.ndarray,
) -> bool:
    """Return True if cached SWMM results match the current PARAMETER_GRID.

    Compares the factor_combinations array embedded in the cache pickle against
    the array produced by build_factor_combinations().  Returns False for old-
    format caches that pre-date this validation field.

    Parameters
    ----------
    cached : dict
        Pickle dict loaded from loocv_calibration_cache.pkl.
    factor_combinations : np.ndarray, shape (n_combos, n_params)
        Current grid from build_factor_combinations().

    Returns
    -------
    bool
        True  — cache is valid; safe to use as-is.
        False — cache is stale or old-format; must re-run SWMM.
    """
    cached_fc = cached.get("factor_combinations")
    if cached_fc is None:
        return False   # old-format cache without embedded grid
    return (
        cached_fc.shape == factor_combinations.shape
        and bool(np.allclose(cached_fc, factor_combinations, rtol=1e-5, atol=1e-9))
    )


def apply_filter_bounds(
    cali_merged: pd.DataFrame,
    filter_bounds: Optional[Dict[str, Optional[Tuple[float, float]]]] = None,
) -> pd.Index:
    """Return the row index of cali_merged that satisfies all filter bounds.

    Parameters
    ----------
    cali_merged : pd.DataFrame
        Merged calibration results with MultiIndex columns.
        Factor columns follow the pattern ('Factors', label).
    filter_bounds : dict, optional
        Maps canonical parameter names to (lo, hi) inclusive bounds, or None
        to skip that parameter.  Defaults to FILTER_BOUNDS from config.

    Returns
    -------
    pd.Index
        Subset of cali_merged.index whose factor values all fall within bounds.
    """
    if filter_bounds is None:
        filter_bounds = FILTER_BOUNDS

    mask = pd.Series(True, index=cali_merged.index)
    for param_name, bounds in filter_bounds.items():
        if bounds is None:
            continue
        lo, hi = bounds
        label = _PARAM_TO_LABEL.get(param_name, param_name)
        col   = ("Factors", label)
        if col not in cali_merged.columns:
            logger.warning("apply_filter_bounds: column %s not found, skipping", col)
            continue
        vals = pd.to_numeric(cali_merged[col], errors="coerce")
        mask &= (vals >= lo) & (vals <= hi)

    n_before = len(cali_merged)
    n_after  = int(mask.sum())
    if n_after < n_before:
        logger.info(
            "apply_filter_bounds: kept %d / %d rows after bounds filtering",
            n_after, n_before,
        )
    return cali_merged.index[mask]


def apply_filter_bounds_to_list(
    df_list: List[pd.DataFrame],
    filter_bounds: Optional[Dict[str, Optional[Tuple[float, float]]]] = None,
) -> List[pd.DataFrame]:
    """Filter every per-event DataFrame in df_list using the same factor bounds.

    All DataFrames in df_list share the same factor columns and values (only
    their date-result columns differ), so the mask derived from df_list[0]
    is applied uniformly to all elements.

    Parameters
    ----------
    df_list : list of pd.DataFrame
        Per-event calibration results from run_calibration() or hydrograph_df_l.
    filter_bounds : dict, optional
        Defaults to FILTER_BOUNDS from config.

    Returns
    -------
    list of pd.DataFrame
        Each element filtered to the same valid row set.
    """
    if filter_bounds is None:
        filter_bounds = FILTER_BOUNDS

    if not any(v is not None for v in filter_bounds.values()):
        return df_list

    valid_idx = apply_filter_bounds(df_list[0], filter_bounds)
    return [df.loc[valid_idx] for df in df_list]


def build_factor_combinations() -> np.ndarray:
    """Build the full calibration factor grid as a 2-D numpy array.

    Iterates ``PARAMETER_GRID`` in declaration order (which matches the legacy
    ``itertools.product`` argument order) and takes the Cartesian product of
    every parameter's values.

    Returns
    -------
    np.ndarray, shape (n_combinations, n_parameters)
        One row per unique factor combination.  Column order matches
        ``list(PARAMETER_GRID.keys())``.
        With the default config: 480 rows × 8 columns.
    """
    axes = [spec.values() for spec in PARAMETER_GRID.values()]
    combos = list(itertools.product(*axes))
    logger.debug("build_factor_combinations: %d combinations", len(combos))
    return np.array(combos)


def factors_row_to_dict(
    row: np.ndarray,
    param_names: Optional[List[str]] = None,
) -> Dict[str, float]:
    """Convert a single factor combination row to a named dict.

    Parameters
    ----------
    row : np.ndarray, shape (n_parameters,)
        One row from the array returned by ``build_factor_combinations()``.
    param_names : list of str, optional
        Parameter names in the same column order as ``row``.
        Defaults to ``list(PARAMETER_GRID.keys())``.

    Returns
    -------
    dict mapping parameter name → factor value.
    """
    if param_names is None:
        param_names = list(PARAMETER_GRID.keys())
    return dict(zip(param_names, row))


def run_calibration(
    cv_dir: Union[str, Path],
    factor_combinations: np.ndarray,
    param_names: Optional[List[str]] = None,
    base_inp_filename: str = SWMM_BASE_INP_FILENAME,
    work_dir: Optional[Union[str, Path]] = None,
) -> Tuple[List[pd.DataFrame], List[pd.DataFrame], pd.DataFrame]:
    """Run all factor combinations for all storm events in cv_dir.

    Replicates the legacy ``run_calibration(MAIN_PATH, FACTOR_COMBINATION)``
    from notebook 01 (cell 12) without global variables.

    Algorithm per date:
      1. Load observed discharge CSV and resample to 5-min grid.
      2. For each factor combination:
         a. Reload the original .inp from disk  (the reset mechanism).
         b. Apply all calibration factors via apply_all_factors().
         c. Write to calibration.inp and run SWMM.
         d. Load simulation output, build storm_df, apply float16 cast,
            clip negatives, fill NaN — exactly as in the legacy inner loop.
         e. Compute storm summary statistics via compute_storm_stats().
         f. Accumulate results in local lists.
      3. Build per-date MultiIndex DataFrames from the accumulated lists.
      4. Build observed hydrograph MultiIndex entry for this date.

    Parameters
    ----------
    cv_dir : str or Path
        Root directory containing one sub-folder per storm event.
        Each sub-folder is named ``YYYY_MM_DD`` and contains:
          - ``{date}.csv``                  observed discharge
          - ``{base_inp_filename}``         original SWMM .inp (never modified)
    factor_combinations : np.ndarray, shape (n_combinations, n_parameters)
        Factor grid from ``build_factor_combinations()``.
    param_names : list of str, optional
        Parameter names matching ``factor_combinations`` column order.
        Defaults to ``list(PARAMETER_GRID.keys())``.
    base_inp_filename : str
        Name of the original .inp file inside each event directory.
    work_dir : str or Path, optional
        Directory where SWMM working files (calibration.inp, calibration.out,
        calibration.rpt) are written.  A subdirectory per event date is created
        automatically.  If None, files are written into each event's cv_dir
        subdirectory (legacy behavior).
        MUST be set to a path outside the legacy read-only directories.

    Returns
    -------
    cali_results_df_l : list of pd.DataFrame
        One DataFrame per event date.  MultiIndex columns:
        ('Factors', label), (date, 'Total Volume'), (date, 'Max Runoff'),
        (date, 'Max Runoff Time').
    hydrograph_df_l : list of pd.DataFrame
        One DataFrame per event date.  MultiIndex columns:
        ('Factors', label), (date, timestamp_string), ...
    obs_hydrograph_df : pd.DataFrame
        Combined observed discharge for all events.
        MultiIndex (Date, Time); columns ['OBS runoff [CMS]', 'rainfall [mm/h]'].

    Notes
    -----
    The float64->float16 cast applied to storm_df inside the inner loop
    matches legacy notebook 01 behavior (memory optimization for the grid
    search).  The cast is intentional and preserved for methodological identity.
    """
    if param_names is None:
        param_names = list(PARAMETER_GRID.keys())

    cv_dir = Path(cv_dir)
    _work_dir = Path(work_dir) if work_dir is not None else None
    factor_labels = [_PARAM_TO_LABEL[p] for p in param_names]
    factors_mi = pd.MultiIndex.from_product([["Factors"], factor_labels])

    cali_results_df_l: List[pd.DataFrame] = []
    hydrograph_df_l:   List[pd.DataFrame] = []
    obs_hydrograph_df  = pd.DataFrame()

    for entry in sorted(os.listdir(str(cv_dir))):
        if not entry.split("_")[0].isdigit():
            continue

        date    = entry
        evt_dir = cv_dir / date
        logger.info("Processing event: %s", date)

        obs_csv  = evt_dir / f"{date}.csv"
        inp_path = evt_dir / base_inp_filename

        # run_dir: where calibration.inp and calibration.out are written.
        # Defaults to the event directory itself (legacy behavior).
        # When work_dir is provided, use an isolated subdirectory to preserve
        # the legacy source directory as strictly read-only.
        if _work_dir is not None:
            run_dir = _work_dir / date
            run_dir.mkdir(parents=True, exist_ok=True)
        else:
            run_dir = evt_dir

        obs_runoff_df = resample_5min(load_observed_runoff(obs_csv))

        # Per-date accumulators (local — no global state)
        factor_rows:            List[List[float]] = []
        swmm_total_runoff_l:    List[float]       = []
        swmm_max_runoff_l:      List[float]       = []
        swmm_max_runoff_time_l: List[object]      = []
        swmm_runoff_l:          List[List[float]] = []
        storm_df_last: Optional[pd.DataFrame] = None

        for factor_row in factor_combinations:
            factor_dict = factors_row_to_dict(factor_row, param_names)

            # Reset: reload original .inp for every combination
            inp = read_inp_file(str(inp_path))
            apply_all_factors(inp, factor_dict)

            cal_inp = run_dir / "calibration.inp"
            inp.write_file(str(cal_inp))
            run_simulation(cal_inp)

            sim_df   = load_swmm_output(run_dir / "calibration.out")
            storm_df = align_obs_sim(obs_runoff_df, sim_df)

            # Float16 cast — legacy notebook 01 memory optimization for grid search
            float64_cols = storm_df.select_dtypes(np.float64).columns
            storm_df[float64_cols] = storm_df[float64_cols].astype(np.float16)
            storm_df[storm_df < 0] = 0
            storm_df = storm_df.fillna(0)

            (swmm_total, _, swmm_max, _, swmm_max_time) = compute_storm_stats(storm_df)

            factor_rows.append(list(factor_row))
            swmm_total_runoff_l.append(swmm_total)
            swmm_max_runoff_l.append(swmm_max)
            swmm_max_runoff_time_l.append(swmm_max_time)
            swmm_runoff_l.append(list(storm_df["SWMM outflow [CMS]"].values))
            storm_df_last = storm_df

        # Build per-date results DataFrame (MultiIndex columns, one row per combo)
        results_cols = pd.MultiIndex.from_tuples(
            list(factors_mi)
            + [(date, "Total Volume"), (date, "Max Runoff"), (date, "Max Runoff Time")]
        )
        results_data = np.column_stack([
            np.array(factor_rows, dtype=object),
            swmm_total_runoff_l,
            swmm_max_runoff_l,
            swmm_max_runoff_time_l,
        ])
        cali_results_date_df = pd.DataFrame(results_data, columns=results_cols)
        for col in factors_mi:
            cali_results_date_df[col] = pd.to_numeric(
                cali_results_date_df[col], errors="coerce"
            )
        cali_results_df_l.append(cali_results_date_df)

        # Build per-date hydrograph DataFrame
        ts_labels = [str(t) for t in storm_df_last["SWMM outflow [CMS]"].index]
        hydro_cols = pd.MultiIndex.from_tuples(
            list(factors_mi) + [(date, ts) for ts in ts_labels]
        )
        hydro_data = np.column_stack([
            np.array(factor_rows, dtype=object),
            np.array(swmm_runoff_l),
        ])
        hydrograph_date_df = pd.DataFrame(hydro_data, columns=hydro_cols)
        hydrograph_df_l.append(hydrograph_date_df)

        # Accumulate observed hydrograph with MultiIndex (Date, Time)
        if storm_df_last is not None:
            mi = pd.MultiIndex.from_tuples(
                [(date, t) for t in storm_df_last.index], names=["Date", "Time"]
            )
            obs_entry = pd.DataFrame(
                storm_df_last[["OBS runoff [CMS]", "rainfall [mm/h]"]].values,
                index=mi,
                columns=["OBS runoff [CMS]", "rainfall [mm/h]"],
            )
            obs_hydrograph_df = pd.concat([obs_hydrograph_df, obs_entry], axis=0)

        logger.debug(
            "Event %s done: %d combos processed",
            date, len(factor_combinations)
        )

    return cali_results_df_l, hydrograph_df_l, obs_hydrograph_df


def merge_calibration_results(
    df_list: List[pd.DataFrame],
) -> pd.DataFrame:
    """Merge per-date calibration DataFrames into a single combined DataFrame.

    Replicates the legacy ``reduce(merge_func, cali_results_df_l)`` step that
    is called after ``run_calibration()`` in the legacy notebook.

    The merge is performed on all columns under the 'Factors' top-level label.
    Each date's result columns are appended as new columns, so the final
    DataFrame has one row per factor combination and one group of columns per
    event date.

    Parameters
    ----------
    df_list : list of pd.DataFrame
        Per-date DataFrames from ``run_calibration()``.

    Returns
    -------
    pd.DataFrame
        Combined DataFrame with all-events results merged on factor columns.
    """
    if not df_list:
        return pd.DataFrame()
    if len(df_list) == 1:
        return df_list[0]

    factor_cols = [col for col in df_list[0].columns if col[0] == "Factors"]

    def _merge(left: pd.DataFrame, right: pd.DataFrame) -> pd.DataFrame:
        return pd.merge(left, right, on=factor_cols)

    return reduce(_merge, df_list)
