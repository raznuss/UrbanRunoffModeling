"""
Pareto front scatter plots for the UrbanRunoffModeling calibration pipeline.

Exact aesthetic match to legacy notebook 04 cells 'dc6dbf87' / '899e971e':
  - All calibration combinations: cornflowerblue
  - Non-dominated (Pareto front) points: lightsalmon
  - Euclidean optimum: red
  - Normalized axes [0, 1] with grid and origin lines
  - Legend: 'Objectives', 'Pareto Front', 'Optimum'

Public API
----------
plot_pareto_scatter(objectives_df, pareto_df, optimum_row, obj_x, obj_y, ax, title)
    Single 2-objective panel — matches the legacy single scatter exactly.
plot_pareto_grid(objectives_df, pareto_df, optimum_row, objective_names, save_path)
    Grid of all unique objective pairs.
"""

import math
from typing import List, Optional

import matplotlib.pyplot as plt
import pandas as pd

# --- Legacy color scheme (from notebook 04) ---
_COLOR_ALL    = "cornflowerblue"   # every calibration combination
_COLOR_PARETO = "lightsalmon"      # non-dominated subset
_COLOR_OPT    = "red"              # Euclidean closest-to-origin


def plot_pareto_scatter(
    objectives_df: pd.DataFrame,
    pareto_df: pd.DataFrame,
    optimum_row: pd.Series,
    obj_x: str,
    obj_y: str,
    ax: Optional[plt.Axes] = None,
    title: str = "Normalized Objective Functions",
) -> plt.Figure:
    """
    Exact legacy 2-objective Pareto scatter from notebook 04.

    Parameters
    ----------
    objectives_df : pd.DataFrame
        All calibration combinations — normalized, columns include obj_x and obj_y.
    pareto_df : pd.DataFrame
        Non-dominated subset of objectives_df (same index subset, same columns).
    optimum_row : pd.Series
        Single-row Series from pareto_df (Euclidean optimum).
        Must have obj_x and obj_y as index labels.
    obj_x, obj_y : str
        Column names to plot on x- and y-axis respectively.
    ax : plt.Axes, optional
        Axes to draw into.  If None, a new 6x6 figure is created.
    title : str
        Plot title (default matches legacy).

    Returns
    -------
    plt.Figure
    """
    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=(6, 6))
    else:
        fig = ax.get_figure()

    s1 = ax.scatter(objectives_df[obj_x], objectives_df[obj_y], color=_COLOR_ALL,    s=10)
    s2 = ax.scatter(pareto_df[obj_x],     pareto_df[obj_y],     color=_COLOR_PARETO, s=15)
    s3 = ax.scatter(optimum_row[obj_x],   optimum_row[obj_y],   color=_COLOR_OPT,    s=60, zorder=5)

    ax.set_xlabel(obj_x, size=15)
    ax.set_ylabel(obj_y, size=15)
    ax.set_title(title, size=18)
    ax.set_xlim(-0.01, 1)
    ax.set_ylim(-0.01, 1)
    ax.grid(True)
    ax.axhline(0, color="black", linewidth=1)
    ax.axvline(0, color="black", linewidth=1)
    ax.tick_params(axis="both", labelsize=11, width=1.5)
    ax.legend([s1, s2, s3], ["Objectives", "Pareto Front", "Optimum"],
              scatterpoints=1, loc="upper right")

    return fig


def plot_pareto_grid(
    objectives_df: pd.DataFrame,
    pareto_df: pd.DataFrame,
    optimum_row: pd.Series,
    objective_names: List[str],
    save_path: Optional[str] = None,
) -> plt.Figure:
    """
    Panel grid of all unique 2-objective pairs from objective_names.

    Lays out C(n, 2) sub-panels in at most 3 columns.  Each sub-panel uses
    the exact same style as plot_pareto_scatter.

    Parameters
    ----------
    objectives_df : pd.DataFrame
        Normalized objectives for all calibration combinations.
    pareto_df : pd.DataFrame
        Non-dominated subset (same normalization).
    optimum_row : pd.Series
        Euclidean optimum (same normalization).
    objective_names : list of str
        Which columns to include.  All pairwise combinations are shown.
    save_path : str, optional
        If given, saves the figure to this path at 150 dpi.

    Returns
    -------
    plt.Figure
    """
    pairs = [
        (objective_names[i], objective_names[j])
        for i in range(len(objective_names))
        for j in range(i + 1, len(objective_names))
    ]
    n_pairs = len(pairs)
    n_cols  = min(3, n_pairs)
    n_rows  = math.ceil(n_pairs / n_cols)

    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(5 * n_cols, 5 * n_rows),
        squeeze=False,
    )
    axes_flat = axes.flatten()

    for k, (obj_x, obj_y) in enumerate(pairs):
        plot_pareto_scatter(
            objectives_df, pareto_df, optimum_row,
            obj_x, obj_y,
            ax=axes_flat[k],
            title="",
        )
        axes_flat[k].set_title(f"{obj_x}\nvs {obj_y}", fontsize=10, pad=4)

    # Hide any empty slots in the grid
    for k in range(n_pairs, len(axes_flat)):
        axes_flat[k].set_visible(False)

    fig.suptitle(
        "Normalized Objective Functions — All Pairs",
        fontsize=14, fontweight="bold", y=1.01,
    )
    plt.tight_layout()

    if save_path is not None:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")

    return fig
