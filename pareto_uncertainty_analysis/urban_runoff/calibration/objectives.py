"""
Objective function computation for the UrbanRunoffModeling calibration pipeline.

  compute_objectives(observed, simulated)
      Exact replication of the legacy ``calculate_objectives(observed, simulated)``
      from notebook 01 (cell 10).  Returns a dict of all six metrics.

  compute_objective_functions(obs_hydrograph_df, cali_results_df, cali_hydrograph_df)
      Exact replication of the legacy ``compute_objective_functions()`` from
      notebook 01 (cell 11).  Computes storm-, peak-, and volume-domain
      objectives across all events and all calibration combinations.

Legacy behavior preserved exactly:
  - NSE uses the manual formula: 1 - ss_diff/ss_mean  (NOT hydroeval's nse)
  - KGE from hydroeval 0.1.0: he.evaluator(he.kge, sim, obs) returns shape (4,1);
    kge_result[0] is shape (1,); float(kge_result[0]) extracts the scalar.
  - bias = abs(mean(sim - obs))  — ABSOLUTE value, not signed
  - r_squared uses the same denominator as NSE (both equal 1 - ss_diff/ss_total
    with ss_total == ss_mean), so r_squared numerically equals NSE.
  - Pareto objectives are stored as minimize-form:
    '1-nse' = 1 - nse,  '1-kge' = 1 - kge,  rmsd/bias/mad kept as-is.
"""

import logging
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import hydroeval as he

logger = logging.getLogger(__name__)


def compute_objectives(
    observed: np.ndarray,
    simulated: np.ndarray,
) -> Dict[str, float]:
    """Compute all six calibration objectives between observed and simulated arrays.

    Exact replication of the legacy ``calculate_objectives(observed, simulated)``
    from notebook 01, cell 10.  Input arrays are converted to float64 numpy
    arrays if needed, matching the legacy isinstance / .to_numpy() guard.

    Parameters
    ----------
    observed : array-like
        Observed discharge or volume values.  Accepts numpy arrays, pandas Series,
        or any array-like.  Converted to float64 internally.
    simulated : array-like
        Simulated discharge or volume values.  Same shape and units as observed.

    Returns
    -------
    dict with keys:
        'nse'       : Nash-Sutcliffe Efficiency  (legacy manual formula)
        'rmsd'      : Root Mean Square Deviation
        'bias'      : abs(mean(sim - obs))  — absolute mean bias
        'mad'       : Mean Absolute Deviation
        'r_squared' : numerically equals NSE (same denominator as legacy)
        'kge'       : Kling-Gupta Efficiency  (hydroeval 0.1.0)
    """
    if not isinstance(observed, np.ndarray):
        observed = np.asarray(observed, dtype=np.float64)
    else:
        observed = observed.astype(np.float64)

    if not isinstance(simulated, np.ndarray):
        simulated = np.asarray(simulated, dtype=np.float64)
    else:
        simulated = simulated.astype(np.float64)

    observed_mean = np.mean(observed)
    ss_diff       = np.sum((observed - simulated) ** 2)
    ss_mean       = np.sum((observed - observed_mean) ** 2)

    nse       = 1.0 - (ss_diff / ss_mean)
    rmsd      = np.sqrt(np.mean((observed - simulated) ** 2))
    bias      = abs(np.mean(simulated - observed))
    mad       = np.mean(np.abs(simulated - observed))

    # hydroeval 0.1.0: he.evaluator(he.kge, sim, obs) → shape (4, 1)
    # Row 0 = KGE, shape (1,). Use .item() to extract scalar without deprecation.
    kge_result = he.evaluator(he.kge, simulated, observed)
    kge        = kge_result[0].item()

    # Legacy uses same ss_total == ss_mean, so r_squared == nse exactly
    ss_total  = np.sum((observed - observed_mean) ** 2)
    r_squared = 1.0 - (ss_diff / ss_total)

    return {
        "nse":       float(nse),
        "rmsd":      float(rmsd),
        "bias":      float(bias),
        "mad":       float(mad),
        "r_squared": float(r_squared),
        "kge":       float(kge),
    }


def compute_objective_functions(
    obs_hydrograph_df: pd.DataFrame,
    cali_results_df: pd.DataFrame,
    cali_swmm_hydrograph_df: pd.DataFrame,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Compute storm-, peak-, and volume-domain objectives for all calibration runs.

    Exact replication of the legacy ``compute_objective_functions()`` from
    notebook 01, cell 11.  Iterates over every calibration combination (row)
    and computes cross-event objectives.

    Parameters
    ----------
    obs_hydrograph_df : pd.DataFrame
        Observed discharge for all events.  MultiIndex (Date, Time).
        Must have column 'OBS runoff [CMS]'.
    cali_results_df : pd.DataFrame
        Combined calibration results DataFrame with MultiIndex columns.
        Must contain 'Max Runoff' and 'Total Volume' at level 1 for each date.
    cali_swmm_hydrograph_df : pd.DataFrame
        Combined simulated hydrograph DataFrame with MultiIndex columns.
        Top-level 'Factors' columns are dropped before objective computation.

    Returns
    -------
    storm_objective_function_df : pd.DataFrame
        One row per calibration combination; columns are the six storm metrics
        in minimize-form under MultiIndex ('storm objective functions', metric).
    peak_objective_function_df : pd.DataFrame
        Same structure for peak-flow domain.
    volume_objective_function_df : pd.DataFrame
        Same structure for total-volume domain.
    """
    storm_cols = [
        "storm_1-nse", "storm_rmsd", "storm_bias",
        "storm_mad",   "storm_r_squared", "storm_1-kge",
    ]
    peak_cols = [
        "peak_1-nse",  "peak_rmsd",  "peak_bias",
        "peak_mad",    "peak_r_squared",  "peak_1-kge",
    ]
    volume_cols = [
        "volume_1-nse", "volume_rmsd", "volume_bias",
        "volume_mad",   "volume_r_squared", "volume_1-kge",
    ]

    storm_obj_df  = pd.DataFrame(
        columns=pd.MultiIndex.from_tuples(
            [("storm objective functions",  c) for c in storm_cols]
        )
    )
    peak_obj_df   = pd.DataFrame(
        columns=pd.MultiIndex.from_tuples(
            [("peak objective functions",   c) for c in peak_cols]
        )
    )
    volume_obj_df = pd.DataFrame(
        columns=pd.MultiIndex.from_tuples(
            [("volume objective functions", c) for c in volume_cols]
        )
    )

    sim_hydrograph_only = cali_swmm_hydrograph_df.drop("Factors", level=0, axis=1)
    index_labels = list(sim_hydrograph_only.index)

    # --- Storm objectives (concatenated cross-event hydrograph) ---
    obs_arr = obs_hydrograph_df["OBS runoff [CMS]"].values.astype(np.float64)
    for i, idx in enumerate(index_labels):
        # Use enumerate+iloc so positional access is safe regardless of index label
        sim_arr = sim_hydrograph_only.iloc[i].astype(np.float64).values
        m = compute_objectives(obs_arr, sim_arr)
        storm_obj_df.loc[idx] = [
            (1 - m["nse"]), m["rmsd"], m["bias"],
            m["mad"],       m["r_squared"], (1 - m["kge"]),
        ]

    # --- Peak objectives (per-event max runoff vector) ---
    obs_max = obs_hydrograph_df.groupby(level=0)["OBS runoff [CMS]"].max().astype(np.float64)
    sim_peak_df = cali_results_df.xs("Max Runoff", level=1, axis=1, drop_level=False).astype(np.float64)
    for i, idx in enumerate(index_labels):
        sim_peak = sim_peak_df.iloc[i]
        m = compute_objectives(obs_max.values, sim_peak.values)
        peak_obj_df.loc[idx] = [
            (1 - m["nse"]), m["rmsd"], m["bias"],
            m["mad"],       m["r_squared"], (1 - m["kge"]),
        ]

    # --- Volume objectives (per-event total volume = sum x 5min x 60s) ---
    obs_vol = (
        obs_hydrograph_df.groupby(level=0)["OBS runoff [CMS]"].sum() * 5 * 60
    ).astype(np.float64)
    sim_vol_df = (
        cali_results_df.xs("Total Volume", level=1, axis=1, drop_level=False) * 5 * 60
    ).astype(np.float64)
    for i, idx in enumerate(index_labels):
        sim_vol = sim_vol_df.iloc[i]
        m = compute_objectives(obs_vol.values, sim_vol.values)
        volume_obj_df.loc[idx] = [
            (1 - m["nse"]), m["rmsd"], m["bias"],
            m["mad"],       m["r_squared"], (1 - m["kge"]),
        ]

    logger.debug(
        "compute_objective_functions: %d combinations, "
        "%d events in storm domain, %d events in peak/volume domain",
        len(sim_hydrograph_only),
        len(obs_arr) // max(1, len(obs_max)),
        len(obs_max),
    )

    return storm_obj_df, peak_obj_df, volume_obj_df


def build_objectives_dataframe(
    storm_obj_df: pd.DataFrame,
    peak_obj_df: pd.DataFrame,
    volume_obj_df: pd.DataFrame,
    objective_names: Optional[List[str]] = None,
) -> pd.DataFrame:
    """Flatten the three multi-domain objective DataFrames into one single-level DataFrame.

    Converts the MultiIndex column output of compute_objective_functions() into a
    flat DataFrame keyed by the short metric name (level 1 of the MultiIndex).
    This is what the legacy notebook does inline:
        objectives_df.columns = [col[1] for col in objectives_df.columns]

    Parameters
    ----------
    storm_obj_df, peak_obj_df, volume_obj_df : pd.DataFrame
        Output tuple from compute_objective_functions().
    objective_names : list of str, optional
        If given, return only these columns from the flattened DataFrame.
        Enables dynamic objective selection for Pareto and LOOCV.

    Returns
    -------
    pd.DataFrame
        Flat single-level DataFrame; one row per calibration combination,
        same integer index as the input DataFrames.
        Columns are the short metric names (e.g. 'peak_1-kge', 'peak_rmsd', ...).
    """
    flat = pd.concat(
        [
            storm_obj_df.droplevel(0, axis=1),
            peak_obj_df.droplevel(0, axis=1),
            volume_obj_df.droplevel(0, axis=1),
        ],
        axis=1,
    )
    if objective_names is not None:
        flat = flat[objective_names]
    return flat
