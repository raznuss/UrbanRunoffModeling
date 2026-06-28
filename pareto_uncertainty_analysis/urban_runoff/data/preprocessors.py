"""
Data preprocessing functions for the UrbanRunoffModeling pipeline.

Each public function is a direct, cleaned-up replacement for the equivalent
inline code from the legacy calibration notebooks.  The mapping is:

  resample_5min(df)               <->  legacy data_5min_interpolate(df)
  filter_rainfall_columns(df)     <->  column-filtering logic in run_filter_38/87.py
  align_obs_sim(obs_df, sim_df)   <->  inline concat + cleanup inside run_calibration()

Design note on align_obs_sim dtype handling:
  The cross-validation notebook (01) applies a float64→float16 cast inside
  run_calibration() for memory efficiency during the large grid search.
  The final calibration notebook (04) does NOT apply this cast.
  align_obs_sim replicates notebook 04's behavior (no dtype cast) because
  it is used for the final validation pipeline.  The memory optimization
  cast is a concern for Phase 4 (calibration/grid_search.py) only.
"""

import logging
from typing import Tuple

import pandas as pd

logger = logging.getLogger(__name__)

# Keywords that identify rainfall-related columns (case-insensitive).
# Replicates the keyword list used in run_filter_38.py / run_filter_87.py.
_RAINFALL_KEYWORDS: Tuple[str, ...] = (
    "rainfall",
    "precip",
    "rain",
    "intensity",
    "outflow_time",
    "flow_time",
    "max_flow_time",
)


def resample_5min(df: pd.DataFrame) -> pd.DataFrame:
    """Resample a time-series DataFrame to a uniform 5-minute grid.

    Replicates legacy data_5min_interpolate() with identical output:
      1. resample('5min').mean()    — average any sub-5-minute samples into bins
      2. interpolate(method='linear') — fill NaN gaps created by the resampling

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame with a DatetimeIndex at any (potentially irregular) frequency.

    Returns
    -------
    pd.DataFrame
        DataFrame on a regular 5-minute DatetimeIndex, with NaN gaps filled by
        linear interpolation.

    Raises
    ------
    TypeError
        If df does not have a DatetimeIndex.
    """
    if not isinstance(df.index, pd.DatetimeIndex):
        raise TypeError(
            f"Input DataFrame must have a DatetimeIndex, got {type(df.index).__name__}."
        )

    resampled = df.resample("5min").mean()
    interpolated = resampled.interpolate(method="linear")

    logger.debug(
        "resample_5min: %d input rows → %d output rows", len(df), len(interpolated)
    )
    return interpolated


def _column_text(col) -> str:
    """Return a lowercase plain-string representation of a column label.

    Handles both flat string labels and MultiIndex tuple labels.
    """
    if isinstance(col, tuple):
        return " ".join(str(part) for part in col if part is not None).lower()
    return str(col).lower()


def _is_rainfall_column(col) -> bool:
    """Return True if the column label contains any rainfall-related keyword."""
    text = _column_text(col)
    return any(kw in text for kw in _RAINFALL_KEYWORDS)


def filter_rainfall_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Remove all rainfall-related columns from a DataFrame.

    Replicates the column-filtering logic from run_filter_38.py / run_filter_87.py.
    A column is dropped if its label (lowercased) contains any of the keywords:
    'rainfall', 'precip', 'rain', 'intensity', 'outflow_time', 'flow_time',
    'max_flow_time'.

    Works with both flat string column labels and MultiIndex column labels.

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame, possibly with MultiIndex columns.

    Returns
    -------
    pd.DataFrame
        Copy of df with all rainfall-related columns removed.
    """
    rainfall_cols = [col for col in df.columns if _is_rainfall_column(col)]
    result = df.drop(columns=rainfall_cols)

    logger.debug(
        "filter_rainfall_columns: removed %d column(s): %s",
        len(rainfall_cols),
        rainfall_cols,
    )
    return result


def align_obs_sim(
    obs_df: pd.DataFrame,
    sim_df: pd.DataFrame,
) -> pd.DataFrame:
    """Merge observed and simulated time-series into a single storm DataFrame.

    Replicates the inline concat + cleanup used in the final calibration
    notebook (04_swmm_calibration_after_cv_multobj.ipynb):
      1. pd.concat([obs_df, sim_df], axis=1)  — align on the shared datetime axis
      2. clip all values < 0 to 0             — physical constraint (no negative discharge)
      3. fillna(0)                            — fill any timesteps missing in either series

    Parameters
    ----------
    obs_df : pd.DataFrame
        Observed runoff, typically after resample_5min().
        Expected column: 'OBS runoff [CMS]'
    sim_df : pd.DataFrame
        SWMM simulation output from load_swmm_output().
        Expected columns: 'SWMM outflow [CMS]', 'rainfall [mm/h]'

    Returns
    -------
    pd.DataFrame
        Combined storm DataFrame with no negative values and no NaN values.
    """
    storm_df = pd.concat([obs_df, sim_df], axis=1)
    storm_df[storm_df < 0] = 0
    storm_df = storm_df.fillna(0)

    logger.debug(
        "align_obs_sim: %d timesteps, %d columns, "
        "negatives clipped, NaNs filled",
        len(storm_df),
        storm_df.shape[1],
    )
    return storm_df
