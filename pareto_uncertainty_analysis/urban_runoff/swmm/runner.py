"""
SWMM simulation runner for the UrbanRunoffModeling calibration pipeline.

Wraps swmm_api.swmm5_run() with error detection and logging so that calling
code can treat any SWMM failure as a Python exception rather than a silent
bad output file.
"""

import logging
from pathlib import Path
from typing import Union

from swmm_api import swmm5_run

logger = logging.getLogger(__name__)


class SwmmRunError(RuntimeError):
    """Raised when a SWMM simulation fails (non-zero exit or .rpt ERROR lines)."""


def run_simulation(inp_path: Union[str, Path], progress_size: int = 1) -> None:
    """Run a SWMM simulation and raise SwmmRunError on any failure.

    Replicates the legacy bare ``swmm5_run(SIM_PATH + file_name + '.inp',
    progress_size=1)`` call but adds RPT-based error detection so that a
    successful Python return truly means a successful simulation.

    Parameters
    ----------
    inp_path : str or Path
        Path to the SWMM .inp file.  The companion .out and .rpt files are
        written to the same directory by SWMM.
    progress_size : int
        Passed directly to ``swmm5_run()``; controls internal progress
        reporting granularity.  Legacy value is 1.

    Raises
    ------
    FileNotFoundError
        If ``inp_path`` does not exist before the simulation starts.
    SwmmRunError
        If ``swmm5_run`` raises, or if the resulting .rpt file contains any
        lines that begin with ``ERROR``.
    """
    inp_path = Path(inp_path)
    if not inp_path.exists():
        raise FileNotFoundError(f"SWMM .inp file not found: {inp_path}")

    logger.debug("Running SWMM simulation: %s", inp_path.name)
    try:
        swmm5_run(str(inp_path), progress_size=progress_size)
    except Exception as exc:
        raise SwmmRunError(
            f"swmm5_run raised an exception for {inp_path}: {exc}"
        ) from exc

    rpt_path = inp_path.with_suffix(".rpt")
    if rpt_path.exists():
        _check_rpt_for_errors(rpt_path)

    logger.debug("SWMM simulation complete: %s", inp_path.name)


def _check_rpt_for_errors(rpt_path: Path) -> None:
    """Scan a SWMM .rpt report file and raise on any ERROR lines.

    Parameters
    ----------
    rpt_path : Path
        Path to the .rpt file written by SWMM after a simulation run.

    Raises
    ------
    SwmmRunError
        If any line in the report starts with the word ``ERROR``.
    """
    with open(rpt_path, "r", errors="replace") as fh:
        for line in fh:
            stripped = line.strip()
            if stripped.startswith("ERROR"):
                raise SwmmRunError(
                    f"SWMM reported an error in {rpt_path.name}: {stripped}"
                )
