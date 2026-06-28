"""
Pareto front computation for the UrbanRunoffModeling calibration pipeline.

Provides two public functions that together replace the legacy Pareto block in
notebook 04 (cells 14–16, 25):

  normalize_columns(df)
      Exact replication of the legacy normalize_columns() — min-max scale to
      [0, 1] across each column independently.

  find_pareto_front(objectives_df, objective_names)
      Vectorized, mathematically sound Pareto front detection using STANDARD
      PARETO DOMINANCE:
          A dominates B  iff  all(A ≤ B)  and  any(A < B)
      This replaces the legacy O(n³) while-loop that contained an
      incorrectly-placed escape valve (count >= 99999) in the inner loop.

Legacy algorithm recap (notebook 04, cell 25):
  The legacy loop used ALL-STRICTLY-BETTER domination (all A < B), which is
  more conservative than standard Pareto:
    - Standard Pareto keeps FEWER points (drops more when ties exist).
    - Legacy kept MORE points (only drops when ALL objectives strictly better).
  In practice, with normalized floating-point objectives, exact ties are
  extremely rare so both algorithms return identical results.  When they
  differ, the standard definition is scientifically correct.

  The count >= 99999 escape valve was checked inside the inner loop but
  incremented in the outer loop, so it would fire only after ~100K outer
  iterations total (across all while passes).  For the 480-combo grid it
  never triggers; it was a latent correctness risk for larger grids.
"""

import logging
from typing import List, Optional, Union

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Min-max normalize each column of df to [0, 1].

    Exact replication of the legacy normalize_columns() from notebook 04,
    cell 14.  Each column is scaled independently:
        normalized = (x - min) / (max - min)

    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame.  All columns must be numeric.

    Returns
    -------
    pd.DataFrame
        New DataFrame with the same shape and column names, values in [0, 1].
        Columns where min == max are left as 0 (avoids division by zero).
    """
    out = df.copy()
    for col in df.columns:
        col_min = df[col].min()
        col_max = df[col].max()
        rng = col_max - col_min
        if rng == 0:
            out[col] = 0.0
        else:
            out[col] = (df[col] - col_min) / rng
    return out


def find_pareto_front(
    objectives_df: pd.DataFrame,
    objective_names: Optional[List[str]] = None,
    normalize: bool = True,
) -> pd.DataFrame:
    """Find the Pareto-optimal rows in an objectives DataFrame.

    Uses STANDARD PARETO DOMINANCE (all ≤ + any <) with a fully vectorized
    implementation: no while loop, no escape valve, O(n² m) time.

    Parameters
    ----------
    objectives_df : pd.DataFrame
        Rows = candidate solutions; columns = objective function values.
        All objectives must be in MINIMIZE form (lower = better).
        The DataFrame may contain all objective columns or a subset.
    objective_names : list of str, optional
        Column names to use as objectives.  If None, all columns are used.
        Providing a subset enables dynamic objective selection (Stage 1 LOOCV).
    normalize : bool
        If True (default), apply min-max normalization before computing the
        Pareto front.  Required for scale-independent dominance checking.
        Set to False only if objectives are already normalized.

    Returns
    -------
    pd.DataFrame
        Subset of objectives_df containing only Pareto-optimal rows,
        preserving the original integer index (for back-joining to the full
        calibration DataFrame via .iloc or .loc).

    Notes
    -----
    Domination definition used here:
        Solution A dominates solution B  iff:
            all(A[k] ≤ B[k])  for all objectives k
            AND  any(A[k] < B[k])  for some objective k

    Legacy definition (notebook 04, cell 15):
        A dominates B  iff  all(A[k] < B[k])  (ALL strictly better)

    The legacy definition is stricter: a tie in any single objective prevents
    domination.  For floating-point normalized values, ties are extremely rare
    and the two approaches yield identical results in practice.
    """
    work_df = objectives_df if objective_names is None else objectives_df[objective_names]

    if normalize:
        work_df = normalize_columns(work_df)

    costs = work_df.to_numpy(dtype=np.float64)
    n, m  = costs.shape
    logger.debug("find_pareto_front: %d candidates, %d objectives", n, m)

    # Vectorized pairwise dominance check.
    # costs_i[i, j, k] = costs[i, k]   (row i viewed across all j competitors)
    # costs_j[i, j, k] = costs[j, k]   (row j viewed as potential dominator of i)
    costs_i = costs[:, np.newaxis, :]   # (n, 1, m) broadcast → (n, n, m)
    costs_j = costs[np.newaxis, :, :]   # (1, n, m) broadcast → (n, n, m)

    weakly_better  = np.all(costs_j <= costs_i, axis=2)   # [i, j] = j ≤ i on all objectives
    strictly_somewhere = np.any(costs_j < costs_i, axis=2)  # [i, j] = j < i on some objective
    j_dominates_i = weakly_better & strictly_somewhere        # (n, n)

    # A point dominates itself trivially — exclude self-comparison
    np.fill_diagonal(j_dominates_i, False)

    # Point i is dominated if ANY j dominates it
    is_dominated = np.any(j_dominates_i, axis=1)
    pareto_mask  = ~is_dominated

    pareto_df = work_df.loc[work_df.index[pareto_mask]]
    n_pareto  = pareto_mask.sum()

    logger.info(
        "find_pareto_front: %d/%d solutions are Pareto-optimal (%.1f%%)",
        n_pareto, n, 100.0 * n_pareto / n,
    )
    return pareto_df


def pareto_front_with_full_row(
    full_df: pd.DataFrame,
    objectives_df: pd.DataFrame,
    objective_names: Optional[List[str]] = None,
) -> pd.DataFrame:
    """Return the Pareto-optimal rows from full_df, using objectives_df for ranking.

    Convenience wrapper for the two-stage workflow where objectives_df contains
    only the objective columns and full_df contains factor columns + objectives.
    The returned DataFrame is the full_df subset indexed by the Pareto front.

    Parameters
    ----------
    full_df : pd.DataFrame
        Complete calibration results DataFrame (factors + all objectives).
    objectives_df : pd.DataFrame
        Objective-only sub-DataFrame, same integer index as full_df.
    objective_names : list of str, optional
        Objective columns to use.  Defaults to all columns of objectives_df.

    Returns
    -------
    pd.DataFrame
        Pareto-optimal rows of full_df.
    """
    pareto_obj = find_pareto_front(objectives_df, objective_names=objective_names)
    return full_df.loc[pareto_obj.index]
