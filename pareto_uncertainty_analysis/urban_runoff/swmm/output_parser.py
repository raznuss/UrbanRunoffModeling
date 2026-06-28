"""
SWMM output parsing and storm statistics for the UrbanRunoffModeling pipeline.

  parse_system_outflow(out_path) — thin wrapper around load_swmm_output()
                                   kept here so the swmm subpackage is self-
                                   contained for callers that only import swmm.
  compute_storm_stats(storm_df)  — replicates the legacy stat() function from
                                   notebook 01 (cell 8).
"""

import logging
from pathlib import Path
from typing import Tuple, Union

import pandas as pd

from urban_runoff.data.loaders import load_swmm_output

logger = logging.getLogger(__name__)


def parse_system_outflow(out_path: Union[str, Path]) -> pd.DataFrame:
    """Load a SWMM .out file and return a two-column system DataFrame.

    Delegates to ``data.loaders.load_swmm_output()``; provided here so callers
    that work exclusively with the ``swmm`` subpackage do not need to import
    from ``data``.

    Parameters
    ----------
    out_path : str or Path
        Path to the SWMM binary output (.out) file.

    Returns
    -------
    pd.DataFrame
        Two-column DataFrame with a DatetimeIndex:
        ['SWMM outflow [CMS]', 'rainfall [mm/h]']
    """
    return load_swmm_output(out_path)


def compute_storm_stats(
    storm_df: pd.DataFrame,
) -> Tuple[float, float, float, float, object]:
    """Compute summary statistics for a single storm event.

    Exact replication of the legacy ``stat(storm_df)`` function from notebook
    01, cell 8.  Returns the same five-tuple in the same order.

    Parameters
    ----------
    storm_df : pd.DataFrame
        Combined observed/simulated storm DataFrame produced by
        ``data.preprocessors.align_obs_sim()``.  Must contain columns:
        'SWMM outflow [CMS]' and 'OBS runoff [CMS]'.

    Returns
    -------
    swmm_total_runoff : float
        Sum of SWMM simulated outflow [CMS × 5-min timestep = CMS·step].
        Multiply by (5 × 60) to convert to cubic meters.
    obs_total_runoff : float
        Sum of observed outflow (same units).
    swmm_max_runoff : float
        Peak simulated outflow [CMS].
    obs_max_runoff : float
        Peak observed outflow [CMS].
    swmm_max_runoff_time : index label
        Timestamp of the simulated peak (first occurrence if tied).
    """
    sim_col = "SWMM outflow [CMS]"
    obs_col = "OBS runoff [CMS]"

    swmm_total_runoff = float(storm_df[sim_col].sum())
    obs_total_runoff  = float(storm_df[obs_col].sum())

    swmm_max_runoff      = float(storm_df[sim_col].max())
    swmm_max_runoff_time = storm_df[storm_df[sim_col] == swmm_max_runoff].index[0]

    obs_max_runoff = float(storm_df[obs_col].max())

    logger.debug(
        "compute_storm_stats: swmm_peak=%.4f CMS at %s, obs_peak=%.4f CMS",
        swmm_max_runoff,
        swmm_max_runoff_time,
        obs_max_runoff,
    )

    return (
        swmm_total_runoff,
        obs_total_runoff,
        swmm_max_runoff,
        obs_max_runoff,
        swmm_max_runoff_time,
    )
