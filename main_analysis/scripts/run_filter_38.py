import os
import pandas as pd

SIM_DIR = r"D:\Development\RESEARCH\Raanana\SWMM\from_radar\Climate_Change\pickles\Urbanization_comparsion"
PREVIEW_SCENARIO = "38% urbanization"
PKL_MAP = {
    '38% urbanization': 'urbanization_1_0.pkl',
    '42% urbanization': 'urbanization_1_1.pkl',
    '46% urbanization': 'urbanization_1_2.pkl',
    '49% urbanization': 'urbanization_1_3.pkl',
    '53% urbanization': 'urbanization_1_4.pkl',
    '57% urbanization': 'urbanization_1_5.pkl',
    '73% urbanization': 'urbanization_2_0.pkl',
    '81% urbanization': 'urbanization_3_0.pkl',
    '83% urbanization': 'urbanization_4_0.pkl',
    '87% urbanization': 'urbanization_5_0.pkl',
}


def column_text(col):
    if isinstance(col, tuple):
        return " ".join(str(part) for part in col if part is not None).lower()
    return str(col).lower()


def is_rainfall_column(col):
    text = column_text(col)
    return any(keyword in text for keyword in ("rainfall", "precip", "rain", "intensity", "outflow_time", "flow_time", "max_flow_time"))


pkl_path = os.path.join(SIM_DIR, PKL_MAP[PREVIEW_SCENARIO])
print(f"Reading pickle: {pkl_path}")
raw_df = pd.read_pickle(pkl_path)
rainfall_columns = [col for col in raw_df.columns if is_rainfall_column(col)]
clean_df = raw_df.drop(columns=rainfall_columns)

output_dir = r"D:\Development\RESEARCH\Raanana\data"
os.makedirs(output_dir, exist_ok=True)
output_path = os.path.join(output_dir, "efrat_row_discharge_results_lowurbanization_level.csv")
clean_df.to_csv(output_path, index=False)

print(f"Loaded scenario: {PREVIEW_SCENARIO}")
print(f"Raw shape: {raw_df.shape}")
print(f"Filtered shape: {clean_df.shape}")
print(f"Saved filtered dataframe to: {output_path}")
print("Removed rainfall-related columns:")
for col in rainfall_columns:
    print(col)

# --- compute and save event-level stats (max/min/median for numeric cols) ---
try:
    event_col = next(col for col in clean_df.columns if "event" in column_text(col) and "num" in column_text(col))
except StopIteration:
    event_col = None

y_col = next((col for col in clean_df.columns if column_text(col) == "y"), None)

if event_col is not None:
    summary_source = clean_df.drop(columns=[col for col in [event_col, y_col] if col is not None])
    summary_source = summary_source.drop(columns=[col for col in summary_source.columns if (col[0] if isinstance(col, tuple) else str(col)).lower() in {"x", "y", "d"}], errors="ignore")
    numeric_cols = summary_source.select_dtypes(include="number").columns
    event_stats = summary_source.groupby(clean_df[event_col])[numeric_cols].agg(["max", "min", "median"])
    summary_output_path = os.path.join(output_dir, "efrat_event_stats_lowurbanization.csv")
    event_stats.to_csv(summary_output_path)
    print(f"Saved event stats to: {summary_output_path}")
    print(f"Summary shape: {event_stats.shape}")
else:
    print("No event column found; skipping event stats generation.")
