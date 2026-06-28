"""
Step 1: Extract and analyze the Pareto front from the calibration pickle file.
This script is READ-ONLY with respect to the original project files.
All outputs are saved exclusively in the Pareto_Uncertainty_Analysis folder.
"""

import os
import math
import pickle
import numpy as np
import pandas as pd

# ─────────────────────────────────────────────────────────────────────────────
# PATHS  (original data - READ ONLY)
# ─────────────────────────────────────────────────────────────────────────────
PICKLE_PATH = r"D:\Development\RESEARCH\Raanana\SWMM\cross_validation_cali_solution_dict.pickle"
OUTPUT_DIR  = r"D:\MY_CODES\UrbanRunoffModeling\Pareto_Uncertainty_Analysis\outputs"

# ─────────────────────────────────────────────────────────────────────────────
# Objective functions that were used to build the Pareto front
# (replicating exactly from notebook 04_swmm_calibration_after_cv_multobj)
# ─────────────────────────────────────────────────────────────────────────────
SELECTED_OBJECTIVES = ['peak_1-kge', 'volume_1-kge', 'peak_rmsd', 'peak_bias']

OBJECTIVES_COLS = [
    ('storm objective functions',  'storm_1-nse'),
    ('storm objective functions',  'storm_rmsd'),
    ('storm objective functions',  'storm_bias'),
    ('storm objective functions',  'storm_mad'),
    ('storm objective functions',  'storm_r_squared'),
    ('storm objective functions',  'storm_1-kge'),
    ('peak objective functions',   'peak_1-nse'),
    ('peak objective functions',   'peak_rmsd'),
    ('peak objective functions',   'peak_bias'),
    ('peak objective functions',   'peak_mad'),
    ('peak objective functions',   'peak_r_squared'),
    ('peak objective functions',   'peak_1-kge'),
    ('volume objective functions', 'volume_1-nse'),
    ('volume objective functions', 'volume_rmsd'),
    ('volume objective functions', 'volume_bias'),
    ('volume objective functions', 'volume_mad'),
    ('volume objective functions', 'volume_r_squared'),
    ('volume objective functions', 'volume_1-kge'),
]

FACTOR_COLS = [
    ('Factors', 'Width'),
    ('Factors', 'IMP'),
    ('Factors', 'Storage'),
    ('Factors', 'N'),
    ('Factors', 'PCT_ZERO'),
    ('Factors', 'CN'),
    ('Factors', 'PCT_ROUTED'),
    ('Factors', 'EVAP'),
]


# ─────────────────────────────────────────────────────────────────────────────
# Helper functions (copied from notebook 04)
# ─────────────────────────────────────────────────────────────────────────────

def normalize_columns(df):
    """Min-max normalize every column independently."""
    normalized_df = df.copy()
    for col in df.columns:
        col_data = df[col]
        normalized_df[col] = (col_data - col_data.min()) / (col_data.max() - col_data.min())
    return normalized_df


def objective_fun_idx_to_drop(objectives_df, row, row_comp):
    """Return row_comp if row dominates it (all objectives strictly smaller)."""
    row_comp_results = abs(objectives_df.iloc[row]) < abs(objectives_df.iloc[row_comp])
    return row_comp if row_comp_results.all() else None


def closest_row_to_origin(df):
    """Return the row with the smallest Euclidean distance to the origin."""
    min_distance = float('inf')
    closest_row_idx = None
    for idx, row in df.iterrows():
        dist = math.sqrt(sum(coord ** 2 for coord in row))
        if dist < min_distance:
            min_distance = dist
            closest_row_idx = idx
    return df.loc[closest_row_idx]


def compute_pareto_front(objectives_df):
    """Iteratively remove dominated solutions to find the Pareto front."""
    pareto_front_df = objectives_df.copy(deep=True)
    clean_flag = False
    count = 0

    while not clean_flag:
        idx_to_drop_l = []
        for row in range(len(pareto_front_df)):
            count += 1
            for row_comp in range(len(pareto_front_df)):
                idx = objective_fun_idx_to_drop(pareto_front_df, row, row_comp)
                if idx is not None:
                    idx_to_drop_l.append(idx)
            if count >= 99999:
                break

        if not idx_to_drop_l:
            clean_flag = True
        else:
            pareto_front_df = pareto_front_df.drop(pareto_front_df.index[idx_to_drop_l])

    return pareto_front_df


# ─────────────────────────────────────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────────────────────────────────────

print("Loading calibration solution pickle...")
with open(PICKLE_PATH, 'rb') as f:
    cross_validation_dict = pickle.load(f)

cali_total_df_l = cross_validation_dict['calibration_dfs']
cali_total_df = cali_total_df_l[0]   # Full 480-row results DataFrame

print(f"Full calibration DataFrame shape: {cali_total_df.shape}")
print(f"Columns (top level): {cali_total_df.columns.get_level_values(0).unique().tolist()}")

# ── Check whether the pickle already has the Pareto front appended ──────────
if len(cali_total_df_l) >= 5:
    print("\nPickle already contains pre-computed Pareto front (index 3) and optimum (index 4).")
    pareto_front_df_normalized = cali_total_df_l[3]   # Normalized objective space
    closest_point_series      = cali_total_df_l[4]    # Full row of the optimum solution
    print(f"Pareto front size (from pickle): {len(pareto_front_df_normalized)} members")
    pareto_indices = pareto_front_df_normalized.index
else:
    # Re-compute from scratch using the same logic as the notebook
    print("\nPre-computed Pareto front not found in pickle. Re-computing...")
    objectives_df = cali_total_df[OBJECTIVES_COLS].copy()
    objectives_df.columns = [col[1] for col in objectives_df.columns]
    objectives_df = abs(objectives_df[SELECTED_OBJECTIVES])
    objectives_df = normalize_columns(objectives_df)

    pareto_front_df_normalized = compute_pareto_front(objectives_df)
    closest_point_series = cali_total_df.iloc[closest_row_to_origin(pareto_front_df_normalized).name]
    pareto_indices = pareto_front_df_normalized.index
    print(f"Pareto front size (re-computed): {len(pareto_front_df_normalized)} members")

# ── Pull the full factor + objective rows for each Pareto member ─────────────
pareto_full_df = cali_total_df.loc[pareto_indices].copy()
pareto_full_df[FACTOR_COLS] = pareto_full_df[FACTOR_COLS].round(3)

print("\n" + "=" * 64)
print(f"PARETO FRONT: {len(pareto_full_df)} non-dominated solutions")
print("=" * 64)

# Display factors for each Pareto member
factor_summary = pareto_full_df[FACTOR_COLS].copy()
factor_summary.columns = [c[1] for c in factor_summary.columns]
print("\nParameter values for all Pareto members:")
print(factor_summary.to_string())

# Display the objective function values for each Pareto member
obj_cols_present = [c for c in OBJECTIVES_COLS if c in cali_total_df.columns]
if obj_cols_present:
    obj_summary = pareto_full_df[obj_cols_present].copy()
    obj_summary.columns = [c[1] for c in obj_summary.columns]
    peak_kge_col    = 'peak_1-kge'    if 'peak_1-kge'    in obj_summary.columns else None
    volume_kge_col  = 'volume_1-kge'  if 'volume_1-kge'  in obj_summary.columns else None
    peak_bias_col   = 'peak_bias'     if 'peak_bias'      in obj_summary.columns else None
    peak_rmsd_col   = 'peak_rmsd'     if 'peak_rmsd'      in obj_summary.columns else None

    display_obj_cols = [c for c in [peak_kge_col, volume_kge_col, peak_bias_col, peak_rmsd_col] if c]
    print("\nObjective function values for Pareto members:")
    print(obj_summary[display_obj_cols].round(4).to_string())

# ── Identify the OPTIMUM (closest to origin) ─────────────────────────────────
print("\n" + "=" * 64)
print("OPTIMUM solution (closest to Pareto-front origin):")
print("=" * 64)
opt_factors = {c[1]: closest_point_series[c] for c in FACTOR_COLS}
print(pd.Series(opt_factors).round(3).to_string())

# ── Normalize the objective space for diversity scoring ──────────────────────
obj_subset = pareto_front_df_normalized.copy() if len(cali_total_df_l) >= 5 else (
    pareto_full_df[obj_cols_present].copy().rename(columns={c: c[1] for c in obj_cols_present})
)

# ── Select 4 diverse members ─────────────────────────────────────────────────
# Strategy:
#   P1 = closest to origin  (balanced – the "optimum")
#   P2 = farthest on peak_1-kge axis (optimized for peak flow KGE)
#   P3 = farthest on volume_1-kge axis (optimized for total volume KGE)
#   P4 = farthest on peak_bias axis (optimized for peak bias – conservative)

# Use the normalized Pareto-front DataFrame to rank
norm_df = pareto_front_df_normalized if len(cali_total_df_l) >= 5 else normalize_columns(
    abs(pareto_full_df[[c for c in OBJECTIVES_COLS if c in pareto_full_df.columns]])
)
norm_df.columns = [c[1] if isinstance(c, tuple) else c for c in norm_df.columns]

# Euclidean distance to origin
norm_df['dist_to_origin'] = norm_df[SELECTED_OBJECTIVES].apply(
    lambda row: math.sqrt(sum(v**2 for v in row)), axis=1
)

p1_idx = norm_df['dist_to_origin'].idxmin()   # balanced / optimum

# For the remaining 3, pick the ones that are most extreme along each objective axis
def farthest_from(df, col):
    """Return index of the row with the maximum value along 'col'."""
    return df[col].idxmax()

p2_idx = farthest_from(norm_df, 'peak_1-kge')    # peak-KGE-optimized
p3_idx = farthest_from(norm_df, 'volume_1-kge')  # volume-KGE-optimized
p4_idx = farthest_from(norm_df, 'peak_bias')      # bias-optimized (conservative)

# Guard against duplicates
selected_indices = list(dict.fromkeys([p1_idx, p2_idx, p3_idx, p4_idx]))
print(f"\nSelected {len(selected_indices)} diverse Pareto members: {selected_indices}")

selected_pareto_df = cali_total_df.loc[selected_indices, FACTOR_COLS].copy()
selected_pareto_df.columns = [c[1] for c in selected_pareto_df.columns]
selected_pareto_df.index = [
    'P1 – Balanced (closest to origin)',
    'P2 – Peak-KGE optimized',
    'P3 – Volume-KGE optimized',
    'P4 – Peak-bias optimized',
][:len(selected_indices)]

# Also pull the objective values
selected_obj_df = cali_total_df.loc[selected_indices, obj_cols_present].copy()
selected_obj_df.columns = [c[1] for c in selected_obj_df.columns]
selected_obj_df.index = selected_pareto_df.index
selected_obj_df = selected_obj_df[[c for c in display_obj_cols if c in selected_obj_df.columns]].round(4)

print("\n" + "#" * 64)
print("SELECTED PARETO PARAMETER SETS")
print("#" * 64)
print(selected_pareto_df.round(3).to_string())
print("\nCorresponding objective function values (lower = better):")
print(selected_obj_df.to_string())

# ── Save selected parameter sets to CSV for Step 2 ───────────────────────────
os.makedirs(OUTPUT_DIR, exist_ok=True)

selected_pareto_df['original_df_index'] = selected_indices
selected_pareto_df.to_csv(os.path.join(OUTPUT_DIR, 'selected_pareto_sets.csv'))

# Also save the original DataFrame row indices so Step 2 can look them up
pd.Series(selected_indices, name='pareto_original_index').to_csv(
    os.path.join(OUTPUT_DIR, 'pareto_selected_indices.csv'), index=False
)

print(f"\nSaved selected parameter sets to: {OUTPUT_DIR}")
print("Done.")
