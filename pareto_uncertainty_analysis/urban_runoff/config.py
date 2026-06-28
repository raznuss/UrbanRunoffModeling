"""
Central configuration for the UrbanRunoffModeling calibration pipeline.

This module implements three architectural requirements:

  REQUIREMENT 1 — EASY PARAMETER CONTROL
    Each SWMM calibration parameter is declared as either:
      Fixed(value)             : held at a constant; excluded from the grid search.
      Sweep(start, stop, step) : swept via np.arange; becomes one grid-search axis.
    To toggle a parameter, change Fixed(...) <-> Sweep(...) in PARAMETER_GRID.
    No other code needs to change.

  REQUIREMENT 2 — DYNAMIC OBJECTIVE SELECTION
    ALL_OBJECTIVES lists every computable metric (the full registry).
    PARETO_OBJECTIVES lists the currently active subset driving Pareto optimization.
    To swap objectives, edit PARETO_OBJECTIVES only. All objectives must use the
    "minimize = better" convention (1-NSE, 1-KGE, RMSD, BIAS, not raw NSE/KGE).

  REQUIREMENT 3 — PARETO ENSEMBLE PRESERVATION
    SAVE_FULL_PARETO_ENSEMBLE = True guarantees the complete non-dominated set is
    always written to disk. The single Euclidean optimum is saved separately but
    is NOT the primary output. The full ensemble is required for uncertainty analysis.
"""

import os
import numpy as np
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple, Union
from pathlib import Path


# ──────────────────────────────────────────────────────────────────────────────
# REQUIREMENT 1 — Parameter specification primitives
# ──────────────────────────────────────────────────────────────────────────────

@dataclass(frozen=True)
class Fixed:
    """Hold a SWMM parameter at a single constant value (not grid-searched).

    Example
    -------
    n = Fixed(1.0)   # Manning's n factor held at 1.0 for all calibration runs
    """
    value: float

    def values(self) -> np.ndarray:
        """Return a 1-element array so all grid-building code is uniform."""
        return np.array([self.value])

    def __repr__(self) -> str:
        return f"Fixed({self.value})"


@dataclass(frozen=True)
class Sweep:
    """Sweep a SWMM parameter over a range via np.arange(start, stop, step).

    Example
    -------
    width = Sweep(0.4, 1.11, 0.3)   # three candidate width factors: [0.4, 0.7, 1.0]
    """
    start: float
    stop:  float
    step:  float

    def values(self) -> np.ndarray:
        """Return the full array of candidate values for this parameter."""
        return np.arange(self.start, self.stop, self.step)

    def __repr__(self) -> str:
        return f"Sweep({self.start}, {self.stop}, {self.step})"


@dataclass(frozen=True)
class Values:
    """Specify exact grid points for a SWMM parameter (non-uniform spacing).

    Example
    -------
    width = Values((0.1, 0.5, 1.0, 3.0, 5.0, 8.0))   # six explicit width factors
    """
    points: tuple

    def values(self) -> np.ndarray:
        """Return the explicit array of candidate values for this parameter."""
        return np.array(self.points, dtype=float)

    def __repr__(self) -> str:
        return f"Values({list(self.points)})"


ParameterSpec = Union[Fixed, Sweep, Values]


# ──────────────────────────────────────────────────────────────────────────────
# PATH CONFIGURATION
# Override any path via the corresponding environment variable.
# ──────────────────────────────────────────────────────────────────────────────

# Root directory containing one sub-folder per storm event (folder name = YYYY_MM_DD)
CROSS_VALIDATION_DIR: Path = Path(
    os.environ.get(
        "URB_RUNOFF_CV_DIR",
        r"D:\Development\RESEARCH\Raanana\SWMM\from_radar\Cross_validation",
    )
)

# Name of the original SWMM input file present in every event directory.
# This file is NEVER modified — all parameter updates write to a separate copy.
SWMM_BASE_INP_FILENAME: str = "raanana_28subcatchments.inp"

# Root directory for calibration result pickles
RESULTS_DIR: Path = Path(
    os.environ.get(
        "URB_RUNOFF_RESULTS_DIR",
        r"D:\Development\RESEARCH\Raanana\SWMM",
    )
)

# Directory where Pareto ensemble outputs are written
PARETO_OUTPUT_DIR: Path = Path(
    os.environ.get(
        "URB_RUNOFF_PARETO_DIR",
        r"D:\MY_CODES\UrbanRunoffModeling_Refactored\outputs\pareto",
    )
)

# Directory for saved figures
FIGURES_DIR: Path = Path(
    os.environ.get(
        "URB_RUNOFF_FIGURES_DIR",
        r"D:\Development\RESEARCH\Raanana\figures\28_catchments_cross_validation\FINAL_MODEL",
    )
)

# Path to the basin radar-overlap rainfall pickle (read-only source data)
BASINS_RAIN_PKL: Path = Path(
    os.environ.get(
        "URB_RUNOFF_RAIN_PKL",
        r"D:\Development\RESEARCH\Raanana\data\rain_radar"
        r"\Basin_radar_overlap_pkl\Basin_radar_overlap_pkl.pkl",
    )
)


# ──────────────────────────────────────────────────────────────────────────────
# SWMM MODEL CONSTANTS
# ──────────────────────────────────────────────────────────────────────────────

# Total basin area in hectares — used for runoff coefficient calculations
BASIN_AREA_HA: float = 1325.6

# Number of subcatchments in the Raanana model
N_SUBCATCHMENTS: int = 28

# Routing destination for the pct_routed parameter (always 'PERVIOUS' in this model)
PCT_ROUTED_DESTINATION: str = "PERVIOUS"


# ──────────────────────────────────────────────────────────────────────────────
# INITIAL PARAMETER VALUE ARRAYS  (Type A parameters only)
#
# These arrays define the BASELINE values for width, CN, and pct_zero —
# the three parameters whose initial state is NOT read from the loaded .inp
# object at runtime, but instead supplied here as external constants.
#
# Source: W = 2*sqrt(Area) method for width; calibration baseline for CN and pct_zero.
# Order: must match the 28-subcatchment order in raanana_28subcatchments.inp.
# ──────────────────────────────────────────────────────────────────────────────

WIDTH_INITIAL_VALUES: np.ndarray = np.array([
    2470, 2800, 2800, 1890, 1160, 2840, 2660, 2000,
    1140, 1400,  540, 2200, 4160, 1070, 1310, 1310,
    1200, 1440, 1400, 1400, 1520, 1850, 1740, 2504,
     980, 1900, 2140,  950,
], dtype=float)

CN_INITIAL_VALUES: np.ndarray = np.full(N_SUBCATCHMENTS, 39, dtype=float)

PCT_ZERO_INITIAL_VALUES: np.ndarray = np.array([
    30, 50, 50, 50, 50, 50, 50, 50, 30, 50, 50, 30,
    30, 50, 30, 50, 50, 50, 50, 50, 50, 60, 50, 50,
    50, 60, 60, 30,
], dtype=float)


# ──────────────────────────────────────────────────────────────────────────────
# REQUIREMENT 1 — PARAMETER GRID
#
# Keys are the canonical parameter names used throughout the pipeline.
# Values are ParameterSpec instances (Fixed, Sweep, or Values).
#
# KEY ORDER defines the column order of the generated factor-combination matrix
# and must match the tuple-unpacking order in the calibration engine.
# The order below replicates the legacy itertools.product argument order exactly:
#   (imperviousness, storage, width, n, pct_zero, cn, pct_routed, evaporation)
#
# Parameter application modes (for inp_editor.py, Phase 3):
#   Type A — multiplicative against INITIAL VALUES arrays above:
#             width, pct_zero, cn
#   Type B — multiplicative against the value currently stored in the .inp object:
#             imperviousness, storage, n
#   Type C — set as an ABSOLUTE value (not a multiplier):
#             pct_routed, evaporation
#
# To switch a parameter from "fixed" to "grid-search":
#   Change Fixed(x) → Sweep(start, stop, step)
# To switch back:
#   Change Sweep(...) → Fixed(x)
# ──────────────────────────────────────────────────────────────────────────────

PARAMETER_GRID: Dict[str, ParameterSpec] = {
    # Type B — factor multiplied against .inp current imperviousness values
    # Range 0.5–1.5; step=0.25 hits calibrated value 1.0 exactly (0.5 + 2×0.25).
    "imperviousness": Values((0.5, 0.75, 0.9, 1.0, 1.1, 1.25, 1.5)),  #            7 values

    # Type B — factor multiplied against .inp current storage values
    # Non-uniform; covers low/dry, baseline, and high depression-storage regimes.
    "storage":        Values((0.1, 0.5, 1.0, 3.0, 5.0, 7.0, 8.0)),  # 7 values

    # Type A — factor multiplied against WIDTH_INITIAL_VALUES
    # Range 0.7–1.3; step=0.3, calibrated value 0.7 is the start.
    "width":          Sweep(0.7,  1.31, 0.3),   # [0.7, 1.0, 1.3]    3 values

    # Type B — factor multiplied against .inp current Manning's n values
    # Fixed at 1.0 per calibration table (not a calibrated parameter).
    "n":              Fixed(1.0),                # [1.0]                                     1 value

    # Type A — factor multiplied against PCT_ZERO_INITIAL_VALUES
    # Fixed at 1.0 per calibration table (not a calibrated parameter).
    "pct_zero":       Fixed(1.0),                # [1.0]                                     1 value

    # Type A — factor multiplied against CN_INITIAL_VALUES
    # Fixed at 1.0 per calibration table (not a calibrated parameter).
    "cn":             Fixed(1.0),                # [1.0]                                     1 value

    # Type C — ABSOLUTE percentage routed to pervious surface [%]
    # Non-uniform; finer resolution around calibrated region (30–45 %); includes 35 and 45.
    "pct_routed":     Values((0, 10, 15, 20, 25, 30, 35, 40, 45)),  # 9 values

    # Type C — ABSOLUTE daily evaporation constant [mm/day]
    # Fixed at 4.0 per calibration table (not a calibrated parameter).
    "evaporation":    Fixed(4.0),                # [4.0]                                     1 value
}

# Total grid combinations = 7 × 7 × 3 × 1 × 1 × 1 × 9 × 1 = 1,323


# ──────────────────────────────────────────────────────────────────────────────
# POST-RUN PARAMETER FILTER BOUNDS
#
# Applied AFTER loading cached SWMM results, BEFORE Pareto evaluation.
# Rows outside any bound are dropped from the search space without re-running
# any simulations.
#
# Usage:
#   Set a (lo, hi) tuple to keep only that inclusive factor range.
#   Set None to disable filtering for that parameter.
#
# Example — restrict imperviousness to physically plausible range:
#   "imperviousness": (1.0, 1.5)
#
# The filter is applied in the notebooks via apply_filter_bounds() from
# urban_runoff.calibration.cross_validation.
# ──────────────────────────────────────────────────────────────────────────────

FILTER_BOUNDS: Dict[str, Optional[Tuple[float, float]]] = {
    "imperviousness": (1.0, 1.5),   # keeps factors 1.0, 1.1, 1.25, 1.5
    "storage":        (0.1, 7.0),   # keeps factors 0.1–7.0; excludes 8.0
    "width":          None,
    "n":              None,
    "pct_zero":       None,
    "cn":             None,
    "pct_routed":     None,
    "evaporation":    None,
}

# Application mode metadata consumed by inp_editor.py (Phase 3)
PARAMETER_MODES: Dict[str, str] = {
    "imperviousness": "multiplicative_from_inp",
    "storage":        "multiplicative_from_inp",
    "width":          "multiplicative_from_array",
    "n":              "multiplicative_from_inp",
    "pct_zero":       "multiplicative_from_array",
    "cn":             "multiplicative_from_array",
    "pct_routed":     "absolute",
    "evaporation":    "absolute",
}

# Array lookup for Type A parameters (None = value read from .inp at runtime)
PARAMETER_INITIAL_VALUES: Dict[str, Optional[np.ndarray]] = {
    "imperviousness": None,
    "storage":        None,
    "width":          WIDTH_INITIAL_VALUES,
    "n":              None,
    "pct_zero":       PCT_ZERO_INITIAL_VALUES,
    "cn":             CN_INITIAL_VALUES,
    "pct_routed":     None,
    "evaporation":    None,
}


# ──────────────────────────────────────────────────────────────────────────────
# REQUIREMENT 2 — OBJECTIVE FUNCTIONS REGISTRY
#
# ALL_OBJECTIVES: complete list of every computable metric.
# PARETO_OBJECTIVES: the active subset currently driving Pareto optimization.
#
# Naming convention — all objectives use "minimize = better" form:
#   '1-nse'  : 1 − NSE  (perfect = 0)
#   '1-kge'  : 1 − KGE  (perfect = 0)
#   'rmsd'   : root mean squared deviation  (perfect = 0)
#   'bias'   : absolute mean bias           (perfect = 0)
#   'mad'    : mean absolute deviation      (perfect = 0)
#   'r2'     : R² stored raw (higher = better; handled specially if used in Pareto)
#
# To change the active set: edit PARETO_OBJECTIVES only.
# All entries in PARETO_OBJECTIVES must be present in ALL_OBJECTIVES.
# ──────────────────────────────────────────────────────────────────────────────

ALL_OBJECTIVES: List[str] = [
    # Storm hydrograph domain
    "storm_1-nse",   "storm_rmsd",   "storm_bias",   "storm_mad",   "storm_r2",  "storm_1-kge",
    # Peak flow domain
    "peak_1-nse",    "peak_rmsd",    "peak_bias",    "peak_mad",    "peak_r2",   "peak_1-kge",
    # Total volume domain
    "volume_1-nse",  "volume_rmsd",  "volume_bias",  "volume_mad",  "volume_r2", "volume_1-kge",
]

# Active objectives for the current Pareto run.
# These four replicate the legacy notebook 04 selection exactly.
PARETO_OBJECTIVES: List[str] = [
    "peak_1-kge",     # 1 − KGE on peak flows       (legacy: objective1)
    "volume_1-kge",   # 1 − KGE on total volumes    (legacy: objective2)
    "peak_rmsd",      # RMSD on peak flows           (legacy: objective3)
    "peak_bias",      # absolute BIAS on peak flows  (legacy: objective4)
]


# ──────────────────────────────────────────────────────────────────────────────
# REQUIREMENT 3 — PARETO ENSEMBLE PRESERVATION
#
# The full non-dominated set (all Pareto members) is ALWAYS saved to disk.
# The single Euclidean-optimum row is also saved for backward compatibility,
# but the ensemble file is the primary artifact for uncertainty analysis.
#
# SAVE_FULL_PARETO_ENSEMBLE must remain True.
# ──────────────────────────────────────────────────────────────────────────────

# Must remain True — setting False would discard the non-dominated ensemble
SAVE_FULL_PARETO_ENSEMBLE: bool = True

# CSV filename for the complete Pareto front (all non-dominated parameter sets + objectives)
PARETO_ENSEMBLE_FILENAME: str = "pareto_ensemble_full.csv"

# CSV filename for the single Euclidean-optimum solution row
PARETO_OPTIMUM_FILENAME: str = "pareto_optimum.csv"
