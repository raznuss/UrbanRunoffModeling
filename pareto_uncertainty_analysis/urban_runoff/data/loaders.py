"""
Data loading functions for the UrbanRunoffModeling pipeline.

Each public function is a direct, cleaned-up replacement for the equivalent
inline code from the legacy calibration notebooks.  The mapping is:

  load_observed_runoff(csv_path)  <->  legacy load_runoff_obs_to_df(path_no_ext)
  load_swmm_output(out_path)      <->  legacy load_sim_to_df(sim_output_path)
  load_pickle(path)               <->  inline pickle.load() calls throughout notebooks
  load_basins_rain(path)          <->  inline pickle.load() for Basin_radar_overlap_pkl

Interface change from legacy:
  load_observed_runoff now accepts the full path WITH the .csv extension.
  The legacy function appended ".csv" internally; callers no longer need to omit it.
"""

import logging
import pickle
from pathlib import Path
from typing import Any, Union

import pandas as pd
from swmm_api import read_out_file

logger = logging.getLogger(__name__)


def load_observed_runoff(csv_path: Union[str, Path]) -> pd.DataFrame:
    """Load an observed discharge CSV file and return a clean DataFrame.

    Replicates legacy load_runoff_obs_to_df() with identical output:
      - Drops rows with NaN values
      - Renames 'discharge_cms' to 'OBS runoff [CMS]'
      - Parses 'date_and_time' with format '%d/%m/%Y %H:%M:%S'
      - Sets 'date_and_time' as a DatetimeIndex
      - Drops the auto-generated integer column ('Unnamed: 0')

    Parameters
    ----------
    csv_path : str or Path
        Full path to the observed runoff CSV file, INCLUDING the .csv extension.

    Returns
    -------
    pd.DataFrame
        Single-column DataFrame with a DatetimeIndex and column 'OBS runoff [CMS]'.

    Raises
    ------
    FileNotFoundError
        If csv_path does not exist on disk.
    KeyError
        If expected columns ('discharge_cms', 'date_and_time', 'Unnamed: 0') are absent.
    """
    csv_path = Path(csv_path)
    if not csv_path.exists():
        raise FileNotFoundError(f"Observed runoff CSV not found: {csv_path}")

    logger.debug("Loading observed runoff from %s", csv_path)

    df = pd.read_csv(csv_path)
    df.dropna(inplace=True)
    df.rename(columns={"discharge_cms": "OBS runoff [CMS]"}, inplace=True, errors="raise")
    df["date_and_time"] = pd.to_datetime(
        df["date_and_time"], format="%d/%m/%Y %H:%M:%S"
    )
    df = df.set_index("date_and_time")
    df.drop(["Unnamed: 0"], axis=1, inplace=True)

    logger.debug("Loaded %d rows of observed runoff from %s", len(df), csv_path.name)
    return df


def load_swmm_output(out_path: Union[str, Path]) -> pd.DataFrame:
    """Load a SWMM binary output (.out) file and return a clean DataFrame.

    Replicates legacy load_sim_to_df() with identical output:
      - Reads system-level outflow and rainfall via swmm_api.read_out_file()
      - Renames columns to 'SWMM outflow [CMS]' and 'rainfall [mm/h]'
      - Ensures the index is a DatetimeIndex

    Parameters
    ----------
    out_path : str or Path
        Full path to the SWMM binary output file (.out).

    Returns
    -------
    pd.DataFrame
        Two-column DataFrame with a DatetimeIndex:
        ['SWMM outflow [CMS]', 'rainfall [mm/h]']

    Raises
    ------
    FileNotFoundError
        If out_path does not exist on disk.
    """
    out_path = Path(out_path)
    if not out_path.exists():
        raise FileNotFoundError(f"SWMM output file not found: {out_path}")

    logger.debug("Loading SWMM output from %s", out_path)

    df = read_out_file(str(out_path)).to_frame()["system"][""][["outflow", "rainfall"]]
    df.rename(
        columns={"outflow": "SWMM outflow [CMS]", "rainfall": "rainfall [mm/h]"},
        inplace=True,
        errors="raise",
    )
    df.index = pd.to_datetime(df.index)

    logger.debug("Loaded %d timesteps of SWMM output from %s", len(df), out_path.name)
    return df


def load_pickle(path: Union[str, Path]) -> Any:
    """Deserialize and return the contents of a pickle file.

    Parameters
    ----------
    path : str or Path
        Full path to the pickle file.

    Returns
    -------
    Any
        The deserialized Python object.

    Raises
    ------
    FileNotFoundError
        If path does not exist on disk.
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Pickle file not found: {path}")

    logger.debug("Loading pickle from %s", path)
    with open(path, "rb") as fh:
        obj = pickle.load(fh)
    return obj


def load_basins_rain(path: Union[str, Path]) -> dict:
    """Load the basin radar-overlap pickle and return the event-to-rain mapping.

    The pickle stores a dict where each key is an event date string
    (e.g. '2012_01_13') and each value is a DataFrame in which rows
    are basins and columns are timesteps with rainfall depth values.

    Parameters
    ----------
    path : str or Path
        Full path to the Basin_radar_overlap_pkl.pkl file.

    Returns
    -------
    dict
        {date_string: pd.DataFrame} for all available events.

    Raises
    ------
    FileNotFoundError
        If path does not exist on disk.
    TypeError
        If the deserialized object is not a dict.
    """
    obj = load_pickle(path)
    if not isinstance(obj, dict):
        raise TypeError(
            f"Expected a dict from {path}, got {type(obj).__name__}"
        )
    logger.debug("Loaded basin rain data: %d events available", len(obj))
    return obj
