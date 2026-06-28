"""
Performance and LOOCV validation scatter plots for the UrbanRunoffModeling pipeline.

Exact aesthetic match to legacy notebooks 04 (cell 73657570) and 02 (cell 43):
  - Font family: Arial
  - Scatter color (calibration): 'black'
  - Scatter color (LOOCV validation): 'navy'
  - 1:1 line: green dashed, alpha=0.8
  - Stats text box: NSE, KGE, R², RB (upper-left, fontsize=18)
  - Spine: edgecolor='black', linewidth=1.2
  - Grid: linestyle=':', alpha=0.6
  - Subplot labels (a)/(b): fontsize=20, bold
  - Volume axes: ScalarFormatter with useMathText, MultipleLocator(1e5)

Public API
----------
calculate_validation_objectives(observed, simulated) -> dict
    Compute NSE / KGE / R² / relative bias using the same hydroeval + scipy approach
    as the legacy notebook.  Used internally by the plot functions.

format_statistics(objectives) -> str
    Build the 4-line stats text box string — exact legacy format.

plot_performance_scatter(obs_peaks, sim_peaks, obs_volumes, sim_volumes, ...)
    Final-calibration 1x2 scatter: peak discharge (a) and total volume (b).

plot_loocv_summary_bar(loocv_summary_df, ...)
    Horizontal bar chart: mean LOOCV error per candidate objective set.

plot_loocv_validation_scatter(fold_results, ...)
    LOOCV held-out scatter: each of the N folds contributes one point.
"""

from typing import Dict, List, Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats
import hydroeval as he
from matplotlib.ticker import ScalarFormatter, MultipleLocator

# --- Legacy font constants from notebook 04 cell 73657570 ---
TITLE_FONT         = {"fontsize": 18, "fontweight": "bold"}
LABEL_FONT         = {"fontsize": 18, "fontweight": "normal"}
TICK_FONT          = {"labelsize": 18}
STATS_FONT         = {"fontsize": 18, "fontweight": "normal"}
SUBPLOT_LABEL_FONT = {"fontsize": 20, "fontweight": "bold"}


def calculate_validation_objectives(
    observed: np.ndarray,
    simulated: np.ndarray,
) -> Dict:
    """
    Compute display metrics for the stats text box on scatter plots.

    Matches the legacy calculate_objectives() from notebook 04 cell 73657570:
      - NSE / KGE via hydroeval
      - RMSD / bias / mad / rel_bias / rel_rmsd computed manually
      - R² from scipy.stats.pearsonr (Pearson r²)

    Parameters
    ----------
    observed, simulated : array-like
        Must have the same length.  Converted to float64 internally.

    Returns
    -------
    dict with keys: 'kge' (array shape (1,)), 'nse', 'rmsd', 'bias',
                    'rel_bias', 'rel_rmsd', 'mad', 'r_squared'
    """
    observed  = np.asarray(observed,  dtype=np.float64).flatten()
    simulated = np.asarray(simulated, dtype=np.float64).flatten()

    kge_result = he.evaluator(he.kge, simulated, observed)   # shape (4, 1)
    nse_result = he.evaluator(he.nse, simulated, observed)   # shape (1, 1)

    obs_mean = max(np.mean(observed), 1e-9)
    rmsd     = float(np.sqrt(np.mean((observed - simulated) ** 2)))
    bias     = float(np.mean(simulated - observed))
    mad      = float(np.mean(np.abs(simulated - observed)))
    rel_bias = float((bias / obs_mean) * 100)
    rel_rmsd = float((rmsd / obs_mean) * 100)

    pearson_r, _ = scipy.stats.pearsonr(observed, simulated)
    r_squared    = float(pearson_r ** 2)

    return {
        "kge":       kge_result[0],          # shape (1,) — matches legacy dict structure
        "nse":       float(nse_result[0]),
        "rmsd":      rmsd,
        "bias":      bias,
        "rel_bias":  rel_bias,
        "rel_rmsd":  rel_rmsd,
        "mad":       mad,
        "r_squared": r_squared,
    }


def format_statistics(objectives: Dict) -> str:
    """
    Build the 4-line stats text box — exact legacy format from notebook 04 cell 73657570.

    Renders as:
        NSE = X.XX
        KGE = X.XX
        R²  = X.XX
        RB  = XX%
    """
    kge_val = objectives["kge"]
    if hasattr(kge_val, "__len__"):
        kge_val = float(kge_val[0])
    return "\n".join([
        f"NSE = {objectives['nse']:.2f}",
        f"KGE = {kge_val:.2f}",
        f"R² = {objectives['r_squared']:.2f}",
        f"RB = {objectives['rel_bias']:.0f}%",
    ])


def _apply_frame_style(ax: plt.Axes) -> None:
    """Darken axis spines — exact match to legacy notebook 04."""
    for spine in ax.spines.values():
        spine.set_edgecolor("black")
        spine.set_linewidth(1.2)


def plot_performance_scatter(
    obs_peaks: np.ndarray,
    sim_peaks: np.ndarray,
    obs_volumes: np.ndarray,
    sim_volumes: np.ndarray,
    title_suffix: str = "",
    save_path: Optional[str] = None,
) -> plt.Figure:
    """
    1x2 subplot: (a) peak discharge scatter, (b) total volume scatter.

    Exact aesthetic match to notebook 04 cell 73657570 (calibration scatter).
    - scatter color: 'black'
    - 1:1 line: green dashed, alpha=0.8
    - stats text: upper-left corner, 4 lines (NSE / KGE / R² / RB)
    - volume axis: ScalarFormatter + MultipleLocator(1e5)
    - subplot labels (a) / (b): fontsize=20, bold, outside top-left of each panel

    Parameters
    ----------
    obs_peaks, sim_peaks : array-like, length = n_events
        Observed and simulated peak discharge per event [m³ s⁻¹].
    obs_volumes, sim_volumes : array-like, length = n_events
        Observed and simulated total storm volume per event [m³].
    title_suffix : str
        If non-empty, added as figure suptitle.
    save_path : str, optional
        Save figure to this path at 150 dpi.

    Returns
    -------
    plt.Figure
    """
    plt.rcParams["font.family"]    = "Arial"
    plt.rcParams["axes.linewidth"] = 1.0

    fig, axs = plt.subplots(1, 2, figsize=(12, 6))

    # --- (a) Peak discharge ---
    x_lim = (0, 50)
    axs[0].scatter(obs_peaks, sim_peaks, color="black", label="Data")
    axs[0].plot(x_lim, x_lim, color="green", linestyle="--", alpha=0.8, label="1:1 line")

    obj_peak = calculate_validation_objectives(
        np.asarray(obs_peaks,  dtype=np.float64),
        np.asarray(sim_peaks,  dtype=np.float64),
    )
    axs[0].text(
        x_lim[1] * 0.01, x_lim[1] * 0.98,
        format_statistics(obj_peak),
        **STATS_FONT, ha="left", va="top",
    )
    axs[0].set_xlabel("Observed peak discharge [m$^3$ s$^{-1}$]",  **LABEL_FONT)
    axs[0].set_ylabel("Simulated peak discharge [m$^3$ s$^{-1}$]", **LABEL_FONT)
    axs[0].set_xlim(*x_lim)
    axs[0].set_ylim(*x_lim)
    axs[0].grid(True, linestyle=":", alpha=0.6)
    axs[0].tick_params(axis="both", **TICK_FONT, width=1.0)
    axs[0].set_aspect("equal")
    axs[0].legend(frameon=True, fontsize=18, loc="lower right")
    axs[0].text(x_lim[1] * -0.06, x_lim[1] * 1.07, "(a)", **SUBPLOT_LABEL_FONT)
    _apply_frame_style(axs[0])

    # --- (b) Total volume ---
    v_max = 4e5
    v_lim = (0, v_max)
    axs[1].scatter(obs_volumes, sim_volumes, color="black", label="Data")
    axs[1].plot(v_lim, v_lim, color="green", linestyle="--", alpha=0.8, label="1:1 line")

    obj_vol = calculate_validation_objectives(
        np.asarray(obs_volumes, dtype=np.float64),
        np.asarray(sim_volumes, dtype=np.float64),
    )
    axs[1].text(
        v_max * 0.01, v_max * 0.98,
        format_statistics(obj_vol),
        **STATS_FONT, ha="left", va="top",
    )
    axs[1].set_xlabel("Observed total volume [m$^3$]",  **LABEL_FONT)
    axs[1].set_ylabel("Simulated total volume [m$^3$]", **LABEL_FONT)
    axs[1].set_xlim(*v_lim)
    axs[1].set_ylim(*v_lim)
    axs[1].grid(True, linestyle=":", alpha=0.6)
    axs[1].tick_params(axis="both", **TICK_FONT, width=1.0)
    axs[1].set_aspect("equal")
    axs[1].text(v_max * -0.05, v_max * 1.07, "(b)", **SUBPLOT_LABEL_FONT)
    _apply_frame_style(axs[1])

    # Scientific notation on volume axes (exact legacy)
    for formatter_attr in ("xaxis", "yaxis"):
        fmt = ScalarFormatter(useMathText=True)
        fmt.set_powerlimits((0, 0))
        getattr(axs[1], formatter_attr).set_major_formatter(fmt)
        getattr(axs[1], formatter_attr).set_major_locator(MultipleLocator(1e5))
        getattr(axs[1], formatter_attr).offsetText.set_fontsize(18)

    if title_suffix:
        fig.suptitle(title_suffix, fontsize=16, fontweight="bold")

    plt.tight_layout()

    if save_path is not None:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")

    return fig


def plot_hydrograph(
    event_date: str,
    obs_hydrograph_df: pd.DataFrame,
    hydro_merged: pd.DataFrame,
    row_idx: int,
    title_suffix: str = "",
    save_path: Optional[str] = None,
) -> plt.Figure:
    """
    Plot observed vs simulated discharge hydrograph for one storm event.

    Dual-axis layout: discharge lines on primary y-axis (bottom), rainfall
    bars on inverted secondary y-axis (top).  Includes a 4-line stats box
    (NSE / KGE / R² / RB) matching the legacy notebook aesthetic.

    Parameters
    ----------
    event_date : str
        Event date label, e.g. '2012_01_13'.  Must exist in obs_hydrograph_df
        and hydro_merged.
    obs_hydrograph_df : pd.DataFrame
        Combined observed hydrograph.  MultiIndex (Date, Time); columns
        'OBS runoff [CMS]' and 'rainfall [mm/h]'.
    hydro_merged : pd.DataFrame
        Merged simulated hydrograph DataFrame from merge_calibration_results().
        Columns: ('Factors', label) + (event_date, timestamp_string).
    row_idx : int
        Row index of the parameter combination to plot (e.g. optimum_row.name).
    title_suffix : str
        Appended to the title after the event date.
    save_path : str, optional
        Save figure to this path at 150 dpi.

    Returns
    -------
    plt.Figure
    """
    plt.rcParams["font.family"]    = "Arial"
    plt.rcParams["axes.linewidth"] = 1.0

    # --- Observed data ---
    event_obs = obs_hydrograph_df.xs(event_date, level="Date")
    obs_flow  = event_obs["OBS runoff [CMS]"].astype(np.float64).values
    rainfall  = event_obs["rainfall [mm/h]"].astype(np.float64).values

    # --- Simulated data for this parameter combination ---
    sim_cols = [c for c in hydro_merged.columns if c[0] == event_date]
    sim_flow = hydro_merged.loc[row_idx, sim_cols].astype(np.float64).values

    # Trim to common length (float16 cast in calibration can occasionally shift length by 1)
    n        = min(len(obs_flow), len(sim_flow), len(rainfall))
    obs_flow = obs_flow[:n]
    sim_flow = sim_flow[:n]
    rainfall = rainfall[:n]
    x        = np.arange(n)

    # --- Statistics ---
    obj        = calculate_validation_objectives(obs_flow, sim_flow)
    stats_text = format_statistics(obj)

    # --- Figure ---
    fig, ax1 = plt.subplots(figsize=(13, 5))
    ax2 = ax1.twinx()

    # Rainfall bars on inverted top axis (legacy style)
    ax2.bar(x, rainfall, color="steelblue", alpha=0.35, width=1.0, label="Rainfall [mm/h]")
    ax2.invert_yaxis()
    rain_max = max(rainfall.max() * 3.5, 1.0)   # push bars to top third
    ax2.set_ylim(rain_max, 0)
    ax2.set_ylabel("Rainfall [mm/h]", **LABEL_FONT)
    ax2.tick_params(axis="y", **TICK_FONT)
    ax2.spines["right"].set_edgecolor("steelblue")
    ax2.spines["right"].set_linewidth(1.2)
    ax2.yaxis.label.set_color("steelblue")
    ax2.tick_params(axis="y", colors="steelblue")

    # Discharge lines
    ax1.plot(x, obs_flow, color="black", linewidth=2.0, label="Observed")
    ax1.plot(x, sim_flow, color="red",   linewidth=2.0, linestyle="--", label="Simulated")

    # Stats text box (upper-left, inside axes)
    ax1.text(
        0.02, 0.97, stats_text,
        transform=ax1.transAxes,
        **STATS_FONT, ha="left", va="top",
        bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor="gray", alpha=0.9),
    )

    ax1.set_xlabel("Time step [5 min]", **LABEL_FONT)
    ax1.set_ylabel("Discharge [m³ s⁻¹]", **LABEL_FONT)
    ax1.tick_params(axis="both", **TICK_FONT, width=1.0)
    ax1.set_xlim(0, n - 1)
    ax1.set_ylim(bottom=0)
    ax1.grid(True, linestyle=":", alpha=0.5)
    ax1.legend(fontsize=14, loc="upper right")
    _apply_frame_style(ax1)

    title = f"Hydrograph  {event_date.replace('_', '-')}"
    if title_suffix:
        title += f"  |  {title_suffix}"
    ax1.set_title(title, **TITLE_FONT)

    plt.tight_layout()

    if save_path is not None:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")

    return fig


def plot_loocv_summary_bar(
    loocv_summary_df: pd.DataFrame,
    save_path: Optional[str] = None,
) -> plt.Figure:
    """
    Horizontal bar chart comparing LOOCV validation errors across candidate objective sets.

    Two paired bars per row: mean peak relative error (steelblue) and
    mean volume relative error (coral), with ±1 std error bars.
    Rows are ordered by ascending mean peak error (matching loocv_summary index order).

    Parameters
    ----------
    loocv_summary_df : pd.DataFrame
        Output of run_full_loocv().  Must have columns:
        'mean_peak_rel_err', 'std_peak_rel_err', 'mean_vol_rel_err', 'std_vol_rel_err'.
        Index = objective set names.
    save_path : str, optional
        Save figure to this path at 150 dpi.

    Returns
    -------
    plt.Figure
    """
    plt.rcParams["font.family"]    = "Arial"
    plt.rcParams["axes.linewidth"] = 1.0

    n      = len(loocv_summary_df)
    y_pos  = np.arange(n)
    labels = list(loocv_summary_df.index)

    peak_means = loocv_summary_df["mean_peak_rel_err"].values
    vol_means  = loocv_summary_df["mean_vol_rel_err"].values
    peak_stds  = loocv_summary_df["std_peak_rel_err"].values
    vol_stds   = loocv_summary_df["std_vol_rel_err"].values

    fig, ax = plt.subplots(figsize=(10, max(4, n * 1.4)))

    bar_h = 0.35
    ax.barh(y_pos + bar_h / 2, peak_means, bar_h,
            xerr=peak_stds, color="steelblue", alpha=0.85,
            error_kw={"linewidth": 1.5, "capsize": 4},
            label="Mean peak relative error")
    ax.barh(y_pos - bar_h / 2, vol_means, bar_h,
            xerr=vol_stds, color="coral", alpha=0.85,
            error_kw={"linewidth": 1.5, "capsize": 4},
            label="Mean volume relative error")

    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=13)
    ax.set_xlabel("Mean relative error (fraction)", **LABEL_FONT)
    ax.set_title("LOOCV Objective Set Comparison", **TITLE_FONT)
    ax.grid(True, axis="x", linestyle=":", alpha=0.6)
    ax.legend(fontsize=14, loc="lower right")
    ax.tick_params(axis="x", **TICK_FONT)
    _apply_frame_style(ax)

    plt.tight_layout()

    if save_path is not None:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")

    return fig


def plot_loocv_validation_scatter(
    fold_results: List[Dict],
    save_path: Optional[str] = None,
) -> plt.Figure:
    """
    LOOCV held-out validation scatter: each fold contributes one point.

    Mirrors the combined validation scatter from notebook 02 cell 43:
      - scatter color: 'navy' (validation held-out points)
      - 1:1 line: green dashed, alpha=0.8
      - Stats text: NSE / KGE / R² / RB (upper-left)
      - subplot labels (a) peak, (b) volume

    Parameters
    ----------
    fold_results : list of dict
        List of fold_result dicts from run_loocv_fold().  Each must have
        a 'validation' key with sub-keys:
        'obs_peak', 'sim_peak', 'obs_vol', 'sim_vol'.
    save_path : str, optional
        Save figure to this path at 150 dpi.

    Returns
    -------
    plt.Figure
    """
    plt.rcParams["font.family"]    = "Arial"
    plt.rcParams["axes.linewidth"] = 1.0

    obs_peaks = np.array([r["validation"]["obs_peak"] for r in fold_results], dtype=np.float64)
    sim_peaks = np.array([r["validation"]["sim_peak"] for r in fold_results], dtype=np.float64)
    obs_vols  = np.array([r["validation"]["obs_vol"]  for r in fold_results], dtype=np.float64)
    sim_vols  = np.array([r["validation"]["sim_vol"]  for r in fold_results], dtype=np.float64)

    fig, axs = plt.subplots(1, 2, figsize=(12, 6))

    # --- (a) Peak discharge ---
    x_lim = (0, 50)
    axs[0].scatter(obs_peaks, sim_peaks, color="navy", label="Held-out storms", zorder=3)
    axs[0].plot(x_lim, x_lim, color="green", linestyle="--", alpha=0.8, label="1:1 line")

    obj_peak = calculate_validation_objectives(obs_peaks, sim_peaks)
    axs[0].text(
        x_lim[1] * 0.01, x_lim[1] * 0.98,
        format_statistics(obj_peak),
        **STATS_FONT, ha="left", va="top",
    )
    axs[0].set_xlabel("Observed peak discharge [m$^3$ s$^{-1}$]",  **LABEL_FONT)
    axs[0].set_ylabel("Simulated peak discharge [m$^3$ s$^{-1}$]", **LABEL_FONT)
    axs[0].set_xlim(*x_lim)
    axs[0].set_ylim(*x_lim)
    axs[0].grid(True, linestyle=":", alpha=0.6)
    axs[0].tick_params(axis="both", **TICK_FONT, width=1.0)
    axs[0].set_aspect("equal")
    axs[0].legend(frameon=True, fontsize=15, loc="lower right")
    axs[0].text(x_lim[1] * -0.06, x_lim[1] * 1.07, "(a)", **SUBPLOT_LABEL_FONT)
    _apply_frame_style(axs[0])

    # --- (b) Total volume ---
    v_max = float(max(obs_vols.max(), sim_vols.max()) * 1.15)
    v_max = max(v_max, 4e5)
    v_lim = (0, v_max)
    axs[1].scatter(obs_vols, sim_vols, color="navy", label="Held-out storms", zorder=3)
    axs[1].plot(v_lim, v_lim, color="green", linestyle="--", alpha=0.8, label="1:1 line")

    obj_vol = calculate_validation_objectives(obs_vols, sim_vols)
    axs[1].text(
        v_max * 0.01, v_max * 0.98,
        format_statistics(obj_vol),
        **STATS_FONT, ha="left", va="top",
    )
    axs[1].set_xlabel("Observed total volume [m$^3$]",  **LABEL_FONT)
    axs[1].set_ylabel("Simulated total volume [m$^3$]", **LABEL_FONT)
    axs[1].set_xlim(*v_lim)
    axs[1].set_ylim(*v_lim)
    axs[1].grid(True, linestyle=":", alpha=0.6)
    axs[1].tick_params(axis="both", **TICK_FONT, width=1.0)
    axs[1].set_aspect("equal")
    axs[1].text(v_max * -0.05, v_max * 1.07, "(b)", **SUBPLOT_LABEL_FONT)
    _apply_frame_style(axs[1])

    for formatter_attr in ("xaxis", "yaxis"):
        fmt = ScalarFormatter(useMathText=True)
        fmt.set_powerlimits((0, 0))
        getattr(axs[1], formatter_attr).set_major_formatter(fmt)
        getattr(axs[1], formatter_attr).offsetText.set_fontsize(18)

    plt.tight_layout()

    if save_path is not None:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")

    return fig
