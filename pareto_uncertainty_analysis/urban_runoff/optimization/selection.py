"""
Optimum solution selection from a Pareto front.

  closest_row_to_origin(pareto_df)
      Replicates the legacy closest_row_to_origin() from notebook 04, cell 16.
      Selects the Pareto-optimal solution with the minimum Euclidean distance
      to the origin in normalized objective space.

  save_pareto_ensemble(pareto_df, full_df, output_dir, ensemble_filename, optimum_filename)
      Implements Requirement 3 — PARETO ENSEMBLE PRESERVATION.
      Saves both the full non-dominated ensemble and the single optimum to CSV.
"""

import logging
import math
from pathlib import Path
from typing import Union

import numpy as np
import pandas as pd

from urban_runoff.config import (
    PARETO_OUTPUT_DIR,
    SAVE_FULL_PARETO_ENSEMBLE,
    PARETO_ENSEMBLE_FILENAME,
    PARETO_OPTIMUM_FILENAME,
)

logger = logging.getLogger(__name__)


def closest_row_to_origin(pareto_df: pd.DataFrame) -> pd.Series:
    """Select the Pareto-optimal row with minimum Euclidean distance to origin.

    Exact replication of the legacy closest_row_to_origin() from notebook 04,
    cell 16.  Iterates over rows and computes sqrt(sum(coord²)), returning the
    row (as a Series) with the minimum distance.

    The input DataFrame is assumed to be normalized to [0, 1] on all objective
    columns, so distance to the origin represents proximity to the ideal point
    [0, 0, ..., 0].

    Parameters
    ----------
    pareto_df : pd.DataFrame
        The Pareto front, typically from find_pareto_front().
        Rows = Pareto-optimal solutions; columns = objective values in [0, 1].

    Returns
    -------
    pd.Series
        The row from pareto_df closest to the origin.  The Series name is the
        original integer index of that row in pareto_df.

    Raises
    ------
    ValueError
        If pareto_df is empty.
    """
    if pareto_df.empty:
        raise ValueError("Cannot select optimum from an empty Pareto front.")

    min_distance   = float("inf")
    closest_idx    = None

    for idx, row in pareto_df.iterrows():
        distance = math.sqrt(sum(coord ** 2 for coord in row))
        if distance < min_distance:
            min_distance = distance
            closest_idx  = idx

    logger.debug(
        "closest_row_to_origin: optimum at index %s, distance=%.6f",
        closest_idx, min_distance,
    )
    return pareto_df.loc[closest_idx]


def closest_row_to_origin_fast(pareto_df: pd.DataFrame) -> pd.Series:
    """Vectorized alternative to closest_row_to_origin() for large fronts.

    Produces identical results to closest_row_to_origin() but uses numpy
    vector operations instead of a Python loop.

    Parameters
    ----------
    pareto_df : pd.DataFrame
        Pareto front with normalized objective values.

    Returns
    -------
    pd.Series
        The row with minimum Euclidean distance to origin.
    """
    if pareto_df.empty:
        raise ValueError("Cannot select optimum from an empty Pareto front.")

    costs     = pareto_df.to_numpy(dtype=np.float64)
    distances = np.sqrt(np.sum(costs ** 2, axis=1))
    best_pos  = int(np.argmin(distances))
    best_idx  = pareto_df.index[best_pos]

    logger.debug(
        "closest_row_to_origin_fast: optimum at index %s, distance=%.6f",
        best_idx, float(distances[best_pos]),
    )
    return pareto_df.loc[best_idx]


def save_pareto_ensemble(
    pareto_obj_df: pd.DataFrame,
    full_calibration_df: pd.DataFrame,
    output_dir: Union[str, Path, None] = None,
    ensemble_filename: str = PARETO_ENSEMBLE_FILENAME,
    optimum_filename:  str = PARETO_OPTIMUM_FILENAME,
) -> Path:
    """Save the full Pareto ensemble and the single optimum row to CSV.

    Implements Requirement 3 (PARETO ENSEMBLE PRESERVATION):
      - The full non-dominated set is always saved regardless of SAVE_FULL_PARETO_ENSEMBLE
        (that flag guards against accidental False; hardcoded behavior here).
      - The single Euclidean optimum is also saved for backward compatibility.

    Parameters
    ----------
    pareto_obj_df : pd.DataFrame
        Normalized objective values for the Pareto front, from find_pareto_front().
    full_calibration_df : pd.DataFrame
        Complete calibration DataFrame (factors + all objectives), same index.
    output_dir : str, Path, or None
        Directory to save CSV files.  Defaults to PARETO_OUTPUT_DIR from config.
    ensemble_filename : str
        Filename for the full ensemble CSV.
    optimum_filename : str
        Filename for the single-optimum CSV.

    Returns
    -------
    Path
        The output directory path where files were saved.

    Raises
    ------
    AssertionError
        If SAVE_FULL_PARETO_ENSEMBLE is False (should never be changed).
    """
    assert SAVE_FULL_PARETO_ENSEMBLE, (
        "SAVE_FULL_PARETO_ENSEMBLE must be True — full ensemble preservation is required."
    )

    out_dir = Path(output_dir) if output_dir is not None else PARETO_OUTPUT_DIR
    out_dir.mkdir(parents=True, exist_ok=True)

    # Full ensemble: all Pareto-optimal rows from the complete calibration DataFrame
    ensemble_df = full_calibration_df.loc[pareto_obj_df.index]
    ensemble_path = out_dir / ensemble_filename
    ensemble_df.to_csv(ensemble_path, index=True)
    logger.info(
        "Saved Pareto ensemble (%d members) to %s",
        len(ensemble_df), ensemble_path,
    )

    # Single optimum: Euclidean closest to origin in normalized objective space
    optimum_row = closest_row_to_origin(pareto_obj_df)
    optimum_df  = full_calibration_df.loc[[optimum_row.name]]
    optimum_path = out_dir / optimum_filename
    optimum_df.to_csv(optimum_path, index=True)
    logger.info("Saved Pareto optimum (index %s) to %s", optimum_row.name, optimum_path)

    return out_dir
