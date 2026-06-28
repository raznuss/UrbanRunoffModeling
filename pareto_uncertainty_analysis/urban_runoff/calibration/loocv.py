"""
Leave-One-Out Cross-Validation (LOOCV) utilities for objective function selection.

Stage 1 of the two-stage calibration workflow:
  - Run calibration once on ALL N events.
  - For each fold k (k = 0..N-1):
      * Training set  : the N-1 events excluding event k.
      * Validation set: event k.
      * Pareto front is found on the training set using a CANDIDATE objective set.
      * The Euclidean optimum is then evaluated on the held-out event k.
  - Repeat for each CANDIDATE objective set.
  - The set that generalizes best across all held-out events is the "winning" set.

Public API:
  run_loocv_fold(...)       — single fold (one held-out event, one objective set)
  run_full_loocv(...)       — all folds x all candidate objective sets; returns summary
"""

import logging
import pickle
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd

from urban_runoff.calibration.cross_validation import merge_calibration_results
from urban_runoff.calibration.objectives import (
    compute_objective_functions,
    build_objectives_dataframe,
    compute_objectives,
)
from urban_runoff.optimization.pareto    import find_pareto_front
from urban_runoff.optimization.selection import closest_row_to_origin

logger = logging.getLogger(__name__)


def _get_event_date(df: pd.DataFrame) -> str:
    """Extract the event date string from a per-date calibration DataFrame.

    The date is the top-level column label for any column that is NOT 'Factors'.
    """
    for col in df.columns:
        if col[0] != "Factors":
            return str(col[0])
    raise ValueError("Could not determine event date from DataFrame columns.")


def run_loocv_fold(
    cali_results_df_l: List[pd.DataFrame],
    hydrograph_df_l: List[pd.DataFrame],
    obs_hydrograph_df: pd.DataFrame,
    held_out_idx: int,
    objective_names: List[str],
) -> Dict[str, Any]:
    """Run one LOOCV fold: train on N-1 events, evaluate the optimum on event k.

    Parameters
    ----------
    cali_results_df_l : list of pd.DataFrame
        Per-event calibration results from run_calibration().
        One element per storm event (N total).
    hydrograph_df_l : list of pd.DataFrame
        Per-event simulated hydrograph DataFrames from run_calibration().
    obs_hydrograph_df : pd.DataFrame
        Combined observed hydrograph for all N events.
        MultiIndex (Date, Time); column 'OBS runoff [CMS]'.
    held_out_idx : int
        Zero-based index of the event to hold out as validation.
        Must be in range [0, N-1].
    objective_names : list of str
        Flat metric names used for Pareto optimization on the training set.
        e.g. ['peak_1-kge', 'volume_1-kge', 'peak_rmsd', 'peak_bias']

    Returns
    -------
    dict with keys:
        held_out_date   : str            — event date of the held-out storm
        held_out_idx    : int            — 0-based fold index
        n_training      : int            — number of training events
        n_pareto        : int            — Pareto front size on training set
        optimum_row_idx : int            — row index of the selected optimum
        optimum_factors : dict           — factor values for the optimum
        validation      : dict[str,float]— objective values on the held-out storm
    """
    n = len(cali_results_df_l)
    if not (0 <= held_out_idx < n):
        raise IndexError(
            f"held_out_idx={held_out_idx} is out of range [0, {n - 1}]."
        )

    held_out_date = _get_event_date(cali_results_df_l[held_out_idx])

    # --- Build training set (exclude held-out event) ---
    train_results = [df for i, df in enumerate(cali_results_df_l) if i != held_out_idx]
    train_hydro   = [df for i, df in enumerate(hydrograph_df_l)   if i != held_out_idx]
    train_dates   = {_get_event_date(df) for df in train_results}

    train_obs = obs_hydrograph_df.loc[
        obs_hydrograph_df.index.get_level_values("Date").isin(train_dates)
    ]

    train_cali = merge_calibration_results(train_results)
    train_hydro_merged = merge_calibration_results(train_hydro)

    # --- Compute objectives on training set ---
    storm_obj, peak_obj, vol_obj = compute_objective_functions(
        train_obs, train_cali, train_hydro_merged
    )
    objectives_flat = build_objectives_dataframe(
        storm_obj, peak_obj, vol_obj, objective_names=objective_names
    )

    # --- Pareto front + Euclidean optimum ---
    pareto_df   = find_pareto_front(objectives_flat, normalize=True)
    optimum_row = closest_row_to_origin(pareto_df)

    # Extract factor values for the optimum (from training merged DataFrame)
    factor_cols = [c for c in train_cali.columns if c[0] == "Factors"]
    optimum_factors_series = train_cali.loc[optimum_row.name, factor_cols]
    optimum_factors = {
        col[1]: float(val) for col, val in optimum_factors_series.items()
    }

    # --- Evaluate optimum on the held-out event ---
    hold_out_results = cali_results_df_l[held_out_idx]
    hold_out_date_obs = obs_hydrograph_df.loc[
        obs_hydrograph_df.index.get_level_values("Date") == held_out_date
    ]

    # Extract this factor combination's simulated stats for the held-out event
    hold_out_row = hold_out_results.iloc[optimum_row.name]
    sim_peak  = float(hold_out_row[(held_out_date, "Max Runoff")])
    sim_vol   = float(hold_out_row[(held_out_date, "Total Volume")]) * 5 * 60
    obs_peak  = float(hold_out_date_obs["OBS runoff [CMS]"].max())
    obs_vol   = float(hold_out_date_obs["OBS runoff [CMS]"].sum()) * 5 * 60

    # Single-point comparison — compute scalar residuals as validation metrics
    validation = {
        "peak_error":  abs(sim_peak - obs_peak),
        "peak_rel_err": abs(sim_peak - obs_peak) / max(obs_peak, 1e-9),
        "vol_error":   abs(sim_vol  - obs_vol),
        "vol_rel_err": abs(sim_vol  - obs_vol)  / max(obs_vol, 1e-9),
        "obs_peak":    obs_peak,
        "sim_peak":    sim_peak,
        "obs_vol":     obs_vol,
        "sim_vol":     sim_vol,
    }

    logger.debug(
        "LOOCV fold %d (%s): n_pareto=%d, optimum=%d, "
        "peak_rel_err=%.3f, vol_rel_err=%.3f",
        held_out_idx, held_out_date, len(pareto_df),
        optimum_row.name,
        validation["peak_rel_err"], validation["vol_rel_err"],
    )

    return {
        "held_out_date":    held_out_date,
        "held_out_idx":     held_out_idx,
        "n_training":       len(train_results),
        "n_pareto":         len(pareto_df),
        "optimum_row_idx":  optimum_row.name,
        "optimum_factors":  optimum_factors,
        "validation":       validation,
    }


def run_full_loocv(
    cali_results_df_l: List[pd.DataFrame],
    hydrograph_df_l: List[pd.DataFrame],
    obs_hydrograph_df: pd.DataFrame,
    candidate_objective_sets: Dict[str, List[str]],
    results_dir: Optional[Union[str, Path]] = None,
) -> pd.DataFrame:
    """Run LOOCV across all folds and all candidate objective sets.

    For each (fold, objective_set) pair:
      1. Train on N-1 events with the candidate set.
      2. Find Pareto optimum.
      3. Evaluate on the held-out event.

    The result is a summary DataFrame showing average validation errors per
    objective set, allowing selection of the best-generalizing set.

    Parameters
    ----------
    cali_results_df_l : list of pd.DataFrame
        Per-event calibration results (N events).
    hydrograph_df_l : list of pd.DataFrame
        Per-event hydrograph DataFrames.
    obs_hydrograph_df : pd.DataFrame
        Combined observed hydrograph for all events.
    candidate_objective_sets : dict
        Keys = human-readable name; values = list of objective metric names.
        e.g. {'legacy_winner': ['peak_1-kge', 'volume_1-kge', 'peak_rmsd', 'peak_bias']}
    results_dir : str or Path, optional
        If given, per-fold results are saved as pickle files in this directory.
        Allows restarting an interrupted LOOCV without re-running completed folds.

    Returns
    -------
    pd.DataFrame
        Summary table, one row per objective set.  Columns:
        'obj_set_name', 'mean_peak_rel_err', 'mean_vol_rel_err',
        'std_peak_rel_err', 'std_vol_rel_err', 'n_folds'.
    """
    n_events   = len(cali_results_df_l)
    _res_dir   = Path(results_dir) if results_dir is not None else None
    if _res_dir is not None:
        _res_dir.mkdir(parents=True, exist_ok=True)

    summary_rows = []

    for set_name, obj_names in candidate_objective_sets.items():
        logger.info(
            "LOOCV: objective set '%s' (%d objectives, %d folds)",
            set_name, len(obj_names), n_events,
        )
        fold_results = []

        for fold_idx in range(n_events):
            # Checkpoint: skip if already computed
            cache_path = (
                _res_dir / f"fold_{fold_idx:02d}_{set_name}.pkl"
                if _res_dir is not None else None
            )
            if cache_path is not None and cache_path.exists():
                with open(cache_path, "rb") as fh:
                    fold_result = pickle.load(fh)
                logger.debug("Loaded cached fold %d for '%s'", fold_idx, set_name)
            else:
                fold_result = run_loocv_fold(
                    cali_results_df_l=cali_results_df_l,
                    hydrograph_df_l=hydrograph_df_l,
                    obs_hydrograph_df=obs_hydrograph_df,
                    held_out_idx=fold_idx,
                    objective_names=obj_names,
                )
                if cache_path is not None:
                    with open(cache_path, "wb") as fh:
                        pickle.dump(fold_result, fh)

            fold_results.append(fold_result)
            logger.info(
                "  Fold %2d/%d (%s): peak_rel_err=%.3f, vol_rel_err=%.3f",
                fold_idx + 1, n_events,
                fold_result["held_out_date"],
                fold_result["validation"]["peak_rel_err"],
                fold_result["validation"]["vol_rel_err"],
            )

        peak_errs = [r["validation"]["peak_rel_err"] for r in fold_results]
        vol_errs  = [r["validation"]["vol_rel_err"]  for r in fold_results]

        summary_rows.append({
            "obj_set_name":      set_name,
            "objectives":        ", ".join(obj_names),
            "mean_peak_rel_err": float(np.mean(peak_errs)),
            "std_peak_rel_err":  float(np.std(peak_errs)),
            "mean_vol_rel_err":  float(np.mean(vol_errs)),
            "std_vol_rel_err":   float(np.std(vol_errs)),
            "n_folds":           n_events,
        })

        if _res_dir is not None:
            all_folds_path = _res_dir / f"all_folds_{set_name}.pkl"
            with open(all_folds_path, "wb") as fh:
                pickle.dump(fold_results, fh)

    summary_df = pd.DataFrame(summary_rows).set_index("obj_set_name")
    summary_df = summary_df.sort_values("mean_peak_rel_err")
    return summary_df
