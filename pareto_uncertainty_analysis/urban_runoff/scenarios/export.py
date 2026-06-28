"""
Urbanization scenario filtering and export for the UrbanRunoffModeling pipeline.

Replaces the two nearly-identical hardcoded scripts:
  scripts/run_filter_38.py  (PREVIEW_SCENARIO = "38% urbanization")
  scripts/run_filter_87.py  (PREVIEW_SCENARIO = "87% urbanization")

The two scripts were identical except for PREVIEW_SCENARIO and the output
filename.  This module unifies them into parameterized functions.

Public API:
  PKL_MAP              : dict mapping scenario label → pickle filename
  export_scenario(...) : load, filter, save CSV + event stats for one scenario
  export_all_scenarios(...) : run export_scenario for every key in PKL_MAP
  compute_event_stats(df, event_col, exclude_cols) : group-by-event aggregation
"""

import logging
import os
from pathlib import Path
from typing import Dict, Optional, Tuple, Union

import pandas as pd

from urban_runoff.data.preprocessors import _column_text, _is_rainfall_column

logger = logging.getLogger(__name__)


# ──────────────────────────────────────────────────────────────────────────────
# Scenario registry
# ──────────────────────────────────────────────────────────────────────────────

PKL_MAP: Dict[str, str] = {
    "38% urbanization": "urbanization_1_0.pkl",
    "42% urbanization": "urbanization_1_1.pkl",
    "46% urbanization": "urbanization_1_2.pkl",
    "49% urbanization": "urbanization_1_3.pkl",
    "53% urbanization": "urbanization_1_4.pkl",
    "57% urbanization": "urbanization_1_5.pkl",
    "73% urbanization": "urbanization_2_0.pkl",
    "81% urbanization": "urbanization_3_0.pkl",
    "83% urbanization": "urbanization_4_0.pkl",
    "87% urbanization": "urbanization_5_0.pkl",
}

# Default source directory for urbanization scenario pickles (read-only legacy data)
_DEFAULT_SIM_DIR = (
    r"D:\Development\RESEARCH\Raanana\SWMM\from_radar"
    r"\Climate_Change\pickles\Urbanization_comparsion"
)

# Default output directory
_DEFAULT_OUTPUT_DIR = r"D:\Development\RESEARCH\Raanana\data"


def export_scenario(
    scenario_key: str,
    output_csv_filename: str,
    sim_dir: Union[str, Path, None] = None,
    output_dir: Union[str, Path, None] = None,
    pkl_map: Optional[Dict[str, str]] = None,
    save_event_stats: bool = True,
    event_stats_filename: Optional[str] = None,
) -> Tuple[pd.DataFrame, Optional[pd.DataFrame]]:
    """Load, filter, and export one urbanization scenario.

    Replaces the runtime behavior of run_filter_38.py / run_filter_87.py.
    The two legacy scripts were identical except for these parameters.

    Algorithm (verbatim from legacy):
      1. Load pickle for the given scenario key.
      2. Remove all rainfall-related columns using is_rainfall_column().
      3. Save the filtered DataFrame as a CSV.
      4. Optionally compute event-level stats (max/min/median grouped by
         the event-number column) and save as a second CSV.

    Parameters
    ----------
    scenario_key : str
        One of the keys in PKL_MAP, e.g. "38% urbanization".
    output_csv_filename : str
        Filename (not full path) for the filtered discharge CSV.
    sim_dir : str or Path, optional
        Directory containing the urbanization pickle files.
        Defaults to the legacy read-only path.
    output_dir : str or Path, optional
        Directory where output CSVs are written.
        Defaults to the legacy output path.
    pkl_map : dict, optional
        Mapping scenario_key → pickle filename.  Defaults to PKL_MAP.
    save_event_stats : bool
        If True, compute and save event-level stats CSV alongside the main CSV.
    event_stats_filename : str, optional
        Filename for the event stats CSV.  If None, defaults to
        ``output_csv_filename.replace('.csv', '_event_stats.csv')``.

    Returns
    -------
    clean_df : pd.DataFrame
        The filtered (rainfall columns removed) scenario DataFrame.
    event_stats : pd.DataFrame or None
        Event-level aggregation stats, or None if save_event_stats=False or
        no event column was found.

    Raises
    ------
    KeyError
        If scenario_key is not in pkl_map.
    FileNotFoundError
        If the pickle file does not exist at the expected path.
    """
    if pkl_map is None:
        pkl_map = PKL_MAP
    if scenario_key not in pkl_map:
        raise KeyError(
            f"Unknown scenario '{scenario_key}'. "
            f"Valid keys: {sorted(pkl_map.keys())}"
        )

    sim_dir    = Path(sim_dir)    if sim_dir    else Path(_DEFAULT_SIM_DIR)
    output_dir = Path(output_dir) if output_dir else Path(_DEFAULT_OUTPUT_DIR)

    pkl_path = sim_dir / pkl_map[scenario_key]
    if not pkl_path.exists():
        raise FileNotFoundError(f"Scenario pickle not found: {pkl_path}")

    logger.info("Loading scenario '%s' from %s", scenario_key, pkl_path.name)
    raw_df = pd.read_pickle(str(pkl_path))

    # Remove rainfall-related columns (reuses preprocessors._is_rainfall_column)
    rainfall_cols = [col for col in raw_df.columns if _is_rainfall_column(col)]
    clean_df = raw_df.drop(columns=rainfall_cols)

    output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = output_dir / output_csv_filename
    clean_df.to_csv(str(csv_path), index=False)

    logger.info(
        "Scenario '%s': raw shape %s → filtered shape %s, saved to %s",
        scenario_key, raw_df.shape, clean_df.shape, csv_path,
    )
    logger.debug("Removed %d rainfall column(s): %s", len(rainfall_cols), rainfall_cols)

    # Event-level statistics
    event_stats: Optional[pd.DataFrame] = None
    if save_event_stats:
        event_stats = compute_event_stats(clean_df)
        if event_stats is not None:
            if event_stats_filename is None:
                base = output_csv_filename.replace(".csv", "")
                event_stats_filename = f"{base}_event_stats.csv"
            stats_path = output_dir / event_stats_filename
            event_stats.to_csv(str(stats_path))
            logger.info("Saved event stats to %s, shape %s", stats_path, event_stats.shape)
        else:
            logger.warning("No event column found in '%s'; skipping event stats.", scenario_key)

    return clean_df, event_stats


def compute_event_stats(df: pd.DataFrame) -> Optional[pd.DataFrame]:
    """Compute per-event max/min/median statistics for all numeric columns.

    Replicates the event stats block from run_filter_38.py / run_filter_87.py.
    Detects the event number column by looking for a column whose label
    contains both 'event' and 'num' (case-insensitive).

    Parameters
    ----------
    df : pd.DataFrame
        Filtered scenario DataFrame with an event-number column.

    Returns
    -------
    pd.DataFrame or None
        MultiIndex-column DataFrame with (column, stat) for each numeric column,
        grouped by event number.  Returns None if no event column is found.
    """
    try:
        event_col = next(
            col for col in df.columns
            if "event" in _column_text(col) and "num" in _column_text(col)
        )
    except StopIteration:
        return None

    y_col = next(
        (col for col in df.columns if _column_text(col) == "y"), None
    )

    # Drop event, y, and single-letter spatial columns before aggregation
    exclude = {c for c in [event_col, y_col] if c is not None}
    spatial_labels = {"x", "y", "d"}

    summary_source = df.drop(columns=list(exclude))
    summary_source = summary_source.drop(
        columns=[
            col for col in summary_source.columns
            if (_column_text(col).split(" ")[0] if isinstance(col, str) else
                str(col[0] if isinstance(col, tuple) else col)).lower()
            in spatial_labels
        ],
        errors="ignore",
    )

    numeric_cols = summary_source.select_dtypes(include="number").columns
    event_stats  = summary_source.groupby(df[event_col])[numeric_cols].agg(
        ["max", "min", "median"]
    )
    return event_stats


def export_all_scenarios(
    output_dir: Union[str, Path, None] = None,
    sim_dir:    Union[str, Path, None] = None,
    pkl_map:    Optional[Dict[str, str]] = None,
) -> None:
    """Export every scenario in PKL_MAP to CSV.

    Generates output filenames automatically from the scenario key:
      "38% urbanization" → "scenario_38pct_urbanization_discharge.csv"

    Parameters
    ----------
    output_dir : str or Path, optional
        Output directory.  Defaults to the legacy output path.
    sim_dir : str or Path, optional
        Source directory for pickle files.  Defaults to the legacy path.
    pkl_map : dict, optional
        Scenario registry.  Defaults to PKL_MAP.
    """
    if pkl_map is None:
        pkl_map = PKL_MAP

    for key in pkl_map:
        # "38% urbanization" → "38pct_urbanization"
        safe_key = key.replace("%", "pct").replace(" ", "_")
        csv_name = f"scenario_{safe_key}_discharge.csv"
        try:
            export_scenario(
                scenario_key=key,
                output_csv_filename=csv_name,
                sim_dir=sim_dir,
                output_dir=output_dir,
                pkl_map=pkl_map,
            )
        except FileNotFoundError as exc:
            logger.warning("Skipping '%s': %s", key, exc)
