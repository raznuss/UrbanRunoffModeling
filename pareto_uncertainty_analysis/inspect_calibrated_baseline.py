"""
Inspect absolute physical parameter values from the calibrated SWMM baseline model.

Reads outputs/models/calibrated_baseline.inp and prints:
  1. Per-subcatchment table (area, IMP, width, Dstore imperv/perv, n, pct_zero, pct_routed).
  2. Factors vs basin-wide absolute values table (reads pareto_optimum.csv if available).
  3. Histogram grid of all parameter distributions.

Usage:
    python inspect_calibrated_baseline.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from swmm_api import read_inp_file

PROJECT_ROOT        = Path(__file__).parent
INP_PATH            = PROJECT_ROOT / "outputs" / "models" / "calibrated_baseline.inp"
PARETO_OPTIMUM_PATH = PROJECT_ROOT / "outputs" / "pareto" / "final" / "pareto_optimum.csv"

# Maps pareto_optimum.csv column label -> canonical config parameter name
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


def extract_subcatchment_table(inp_path: Path) -> pd.DataFrame:
    """Read all per-subcatchment physical parameters from a calibrated .inp."""
    inp           = read_inp_file(str(inp_path))
    subcatchments = inp["SUBCATCHMENTS"]
    subareas      = inp["SUBAREAS"]

    rows = []
    for name, sc in subcatchments.items():
        sa = subareas.get(name)
        rows.append({
            "subcatchment":   name,
            "area_ha":        float(sc.area),
            "imperv_pct":     float(sc.imperviousness),
            "width_m":        float(sc.width),
            "slope_pct":      float(sc.slope),
            "n_imperv":       float(sa.n_imperv)       if sa is not None else float("nan"),
            "storage_imperv": float(sa.storage_imperv) if sa is not None else float("nan"),
            "storage_perv":   float(sa.storage_perv)   if sa is not None else float("nan"),
            "pct_zero":       float(sa.pct_zero)        if sa is not None else float("nan"),
            "pct_routed":     float(sa.pct_routed)      if sa is not None else float("nan"),
        })

    return pd.DataFrame(rows).set_index("subcatchment")


def _read_evaporation(inp_path: Path) -> float:
    """Return the CONSTANT evaporation rate [mm/day] from the .inp."""
    inp = read_inp_file(str(inp_path))
    return float(inp["EVAPORATION"]["CONSTANT"])


def basin_avg_imperviousness(df: pd.DataFrame) -> float:
    """Area-weighted average imperviousness across all subcatchments [%]."""
    return float((df["area_ha"] * df["imperv_pct"]).sum() / df["area_ha"].sum())


def _weighted_avg(df: pd.DataFrame, col: str) -> float:
    """Area-weighted average of df[col] using df['area_ha'] as weights."""
    return float((df["area_ha"] * df[col]).sum() / df["area_ha"].sum())


def print_table(df: pd.DataFrame, avg_imp: float) -> None:
    """Print a formatted per-subcatchment table to terminal."""
    fmt = pd.DataFrame({
        "Area [ha]":      df["area_ha"].round(2),
        "IMP [%]":        df["imperv_pct"].round(1),
        "Width [m]":      df["width_m"].round(1),
        "Dstore-I [mm]":  df["storage_imperv"].round(3),
        "Dstore-P [mm]":  df["storage_perv"].round(3),
        "n_imperv":       df["n_imperv"].round(4),
        "PctZero [%]":    df["pct_zero"].round(1),
        "PctRouted [%]":  df["pct_routed"].round(1),
    })

    print("=" * 105)
    print("Calibrated baseline  -  per-subcatchment physical parameters")
    print("=" * 105)
    print(fmt.to_string())
    print("=" * 105)
    print(f"  N subcatchments           : {len(df)}")
    print(f"  Total basin area          : {df['area_ha'].sum():.1f} ha")
    print(f"  Area-weighted avg IMP     : {avg_imp:.2f} %")
    print(f"  IMP range                 : {df['imperv_pct'].min():.1f} - {df['imperv_pct'].max():.1f} %")
    print(f"  Width range               : {df['width_m'].min():.0f} - {df['width_m'].max():.0f} m")
    print(f"  Dstore imperv range       : {df['storage_imperv'].min():.3f} - {df['storage_imperv'].max():.3f} mm")
    print(f"  Dstore perv  range        : {df['storage_perv'].min():.3f} - {df['storage_perv'].max():.3f} mm")
    print("=" * 105)


def print_factors_table(df: pd.DataFrame, evap_mmd: float) -> None:
    """Print optimal factors vs basin-wide absolute values.

    Reads the 'Factors' columns from PARETO_OPTIMUM_PATH (pareto_optimum.csv)
    if the file exists.  Absolute values are computed from the calibrated .inp.

    The 'storage' factor applies to both Dstore imperv and Dstore perv.
    """
    # (config_param, csv_label, absolute_value_str, unit_note)
    param_info = [
        ("imperviousness", "IMP",
         f"{_weighted_avg(df, 'imperv_pct'):.2f}",
         "%   area-weighted avg IMP"),

        ("storage", "Storage",
         f"Dstore-I={_weighted_avg(df, 'storage_imperv'):.3f}  "
         f"Dstore-P={_weighted_avg(df, 'storage_perv'):.3f}",
         "mm  area-weighted avg (imperv / perv)"),

        ("width", "Width",
         f"{_weighted_avg(df, 'width_m'):.0f}",
         "m   area-weighted avg width"),

        ("n", "N",
         f"{df['n_imperv'].mean():.4f}",
         "-   mean Manning n imperv"),

        ("pct_zero", "PCT_ZERO",
         f"{_weighted_avg(df, 'pct_zero'):.1f}",
         "%   area-weighted avg zero-storage area"),

        ("cn", "CN",
         "39.00",
         "-   fixed, uniform (CN_INITIAL_VALUES = 39)"),

        ("pct_routed", "PCT_ROUTED",
         f"{float(df['pct_routed'].iloc[0]):.1f}",
         "%   absolute constant (same for all subcatchments)"),

        ("evaporation", "EVAP",
         f"{evap_mmd:.1f}",
         "mm/day  absolute daily evaporation constant"),
    ]

    # Read optimal factors from pareto_optimum.csv
    factors = {}
    if PARETO_OPTIMUM_PATH.exists():
        try:
            df_opt = pd.read_csv(PARETO_OPTIMUM_PATH, header=[0, 1], index_col=0)
            for _, csv_label, _, _ in param_info:
                col = ("Factors", csv_label)
                if col in df_opt.columns:
                    factors[csv_label] = float(df_opt[col].iloc[0])
        except Exception as exc:
            print(f"  WARNING: could not read {PARETO_OPTIMUM_PATH}: {exc}")
    else:
        print(f"  (pareto_optimum.csv not found at {PARETO_OPTIMUM_PATH})")
        print("  Run 03_Final_Calibration.ipynb first to produce the optimum file.")

    print()
    print("=" * 100)
    print("Pareto Optimum  -  Factors vs Basin-Wide Absolute Values")
    print("=" * 100)
    print(f"  {'Parameter':<16}  {'Factor':>7}  {'Basin-Wide Absolute Value':>38}  Note")
    print("-" * 100)
    for pname, csv_label, abs_val, note in param_info:
        factor_str = f"{factors[csv_label]:7.4f}" if csv_label in factors else "    N/A"
        print(f"  {pname:<16}  {factor_str}  {abs_val:>38}  {note}")
    print("=" * 100)


def plot_histograms(df: pd.DataFrame, avg_imp: float) -> None:
    """Histogram grid of key physical parameters (3 x 3, 7 panels)."""
    fig, axes = plt.subplots(3, 3, figsize=(15, 11))
    fig.suptitle(
        "Calibrated Baseline - Subcatchment Parameter Distributions",
        fontsize=14, fontweight="bold",
    )

    params = [
        ("imperv_pct",     "Imperviousness [%]",      "steelblue"),
        ("storage_imperv", "Dstore imperv [mm]",       "darkorange"),
        ("storage_perv",   "Dstore perv [mm]",         "coral"),
        ("width_m",        "Width [m]",                "seagreen"),
        ("n_imperv",       "Manning n (impervious)",   "mediumpurple"),
        ("pct_zero",       "Zero-storage area [%]",    "firebrick"),
        ("pct_routed",     "% routed to pervious",     "teal"),
    ]

    for i, ax in enumerate(axes.flat):
        if i >= len(params):
            ax.set_visible(False)
            continue
        col, label, color = params[i]
        ax.hist(df[col].dropna(), bins=10, color=color, edgecolor="white", alpha=0.85)
        if col == "imperv_pct":
            ax.axvline(avg_imp, color="red", linestyle="--", linewidth=1.5,
                       label=f"Avg = {avg_imp:.1f}%")
            ax.legend(fontsize=10)
        ax.set_xlabel(label, fontsize=11)
        ax.set_ylabel("Count", fontsize=11)
        ax.grid(True, linestyle=":", alpha=0.5)

    plt.tight_layout()
    out_path = PROJECT_ROOT / "outputs" / "metrics" / "calibrated_baseline_histograms.png"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"\nHistograms saved -> {out_path}")
    plt.show()


def main() -> None:
    if not INP_PATH.exists():
        print(f"ERROR: calibrated baseline not found:\n  {INP_PATH}")
        print("Run 03_Final_Calibration.ipynb, then: python generate_calibrated_baseline.py")
        sys.exit(1)

    print(f"Reading: {INP_PATH}")
    df      = extract_subcatchment_table(INP_PATH)
    evap    = _read_evaporation(INP_PATH)
    avg_imp = basin_avg_imperviousness(df)

    print_table(df, avg_imp)
    print_factors_table(df, evap)
    plot_histograms(df, avg_imp)


if __name__ == "__main__":
    main()
