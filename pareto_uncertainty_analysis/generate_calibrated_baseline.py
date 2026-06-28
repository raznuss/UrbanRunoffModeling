"""
Generate the final calibrated SWMM baseline model.

Reads the Pareto-optimal parameter combination from the saved CSV
(outputs/pareto/final/pareto_optimum.csv), applies those factors to the
original read-only base .inp file, and writes the result to:

    outputs/models/calibrated_baseline.inp

This file is the definitive calibrated structural model for Part 2
(WRF climate change and urbanization scenario analysis).
The TIMESERIES section will be replaced with WRF forcing data for each scenario run;
all structural parameters (imperviousness, width, storage, N, CN, etc.) are
permanently baked in at their calibrated values.

Usage:
    python generate_calibrated_baseline.py

Legacy directories are never modified.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

import pandas as pd
from swmm_api import read_inp_file

from urban_runoff.config import SWMM_BASE_INP_FILENAME
from urban_runoff.swmm.inp_editor import apply_all_factors

# ── Paths ─────────────────────────────────────────────────────────────────────

PROJECT_ROOT = Path(__file__).parent

# Source: Pareto optimum CSV written by 02_Stage2_Final_Calibration.ipynb
OPTIMUM_CSV = PROJECT_ROOT / "outputs" / "pareto" / "final" / "pareto_optimum.csv"

# Source: original (read-only) base model from the first CV event
# All events share the same structural .inp — only timeseries differ.
CV_DIR = Path(
    r"D:\Development\RESEARCH\Raanana\SWMM\from_radar\Cross_validation"
)
FIRST_EVENT = "2012_01_13"
BASE_INP_PATH = CV_DIR / FIRST_EVENT / SWMM_BASE_INP_FILENAME

# Destination
OUTPUT_DIR  = PROJECT_ROOT / "outputs" / "models"
OUTPUT_PATH = OUTPUT_DIR / "calibrated_baseline.inp"

# ── Label → canonical name mapping (inverse of cross_validation._PARAM_TO_LABEL) ──

_LABEL_TO_PARAM = {
    "IMP":        "imperviousness",
    "Storage":    "storage",
    "Width":      "width",
    "N":          "n",
    "PCT_ZERO":   "pct_zero",
    "CN":         "cn",
    "PCT_ROUTED": "pct_routed",
    "EVAP":       "evaporation",
}


def load_optimum_factors(csv_path: Path) -> dict:
    """
    Read the Pareto optimum factors from the 2-row-header CSV.

    Returns a dict with canonical parameter names as keys
    (e.g. 'imperviousness', 'width', ...) and float values.
    """
    # The CSV was written with a MultiIndex column header (2 rows)
    df = pd.read_csv(csv_path, header=[0, 1], index_col=0)

    factor_dict = {}
    for label, param_name in _LABEL_TO_PARAM.items():
        col_val = df[("Factors", label)].iloc[0]
        factor_dict[param_name] = float(col_val)

    return factor_dict


def generate_baseline(
    base_inp_path: Path,
    factor_dict: dict,
    output_path: Path,
) -> None:
    """
    Apply calibrated factors to the base model and write the result.

    The input file at base_inp_path is loaded fresh (legacy read-only) and
    never modified on disk.  The calibrated model is written to output_path.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    inp = read_inp_file(str(base_inp_path))
    apply_all_factors(inp, factor_dict)
    inp.write_file(str(output_path))


def main():
    if not OPTIMUM_CSV.exists():
        print(f"ERROR: Pareto optimum CSV not found:\n  {OPTIMUM_CSV}")
        print("Run 02_Stage2_Final_Calibration.ipynb first.")
        sys.exit(1)

    if not BASE_INP_PATH.exists():
        print(f"ERROR: Base .inp not found:\n  {BASE_INP_PATH}")
        sys.exit(1)

    print("Reading Pareto optimum factors...")
    factor_dict = load_optimum_factors(OPTIMUM_CSV)

    print("Calibrated parameter factors:")
    for param, value in factor_dict.items():
        print(f"  {param:<16} = {value}")

    print(f"\nLoading base model from:\n  {BASE_INP_PATH}")
    print(f"Applying factors and writing to:\n  {OUTPUT_PATH}")
    generate_baseline(BASE_INP_PATH, factor_dict, OUTPUT_PATH)

    size_kb = OUTPUT_PATH.stat().st_size / 1024
    print(f"\nDone. File size: {size_kb:.1f} KB")
    print(f"\nCalibratedbaseline .inp path:\n  {OUTPUT_PATH.resolve()}")


if __name__ == "__main__":
    main()
