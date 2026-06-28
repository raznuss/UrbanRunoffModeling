"""
Correct per-event peak-flow uncertainty visualization.

For each storm event the script plots:
  - A vertical line spanning the MIN to MAX simulated peak across the 4 Pareto members
  - 4 individual dots (one per Pareto set) along that vertical line
  - The observed peak as a prominent black star

This is methodologically correct: uncertainty is evaluated within each event,
not by mixing flows from different storms into a single distribution.

Output: plots/peak_uncertainty_per_event.png
        plots/peak_uncertainty_per_event_sorted.png  (sorted by observed magnitude)

READ-ONLY: this script only reads ensemble_results.pkl and writes to plots/.
No files outside Pareto_Uncertainty_Analysis/ are touched.
"""

import os
import pickle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.patheffects as pe

# ─────────────────────────────────────────────────────────────────────────────
# PATHS
# ─────────────────────────────────────────────────────────────────────────────
BASE_DIR    = os.path.dirname(os.path.abspath(__file__))
PICKLE_PATH = os.path.join(BASE_DIR, 'outputs', 'ensemble_results.pkl')
PLOTS_DIR   = os.path.join(BASE_DIR, 'plots')
os.makedirs(PLOTS_DIR, exist_ok=True)

# ─────────────────────────────────────────────────────────────────────────────
# LOAD PRE-COMPUTED ENSEMBLE RESULTS
# ─────────────────────────────────────────────────────────────────────────────
with open(PICKLE_PATH, 'rb') as f:
    data = pickle.load(f)

ensemble_results = data['ensemble_results']   # {label: {date: storm_df}}
obs_data         = data['obs_data']           # {date: obs_df}
PARETO_SETS      = data['pareto_sets']        # list of dicts with 'label', 'color'

# ─────────────────────────────────────────────────────────────────────────────
# EXTRACT PER-EVENT PEAK FLOWS
# ─────────────────────────────────────────────────────────────────────────────
# Sort events chronologically (dates are already in 'YYYY_MM_DD' format)
all_dates = sorted(obs_data.keys())

# Build arrays:
#   obs_peaks[i]        → observed peak for event i
#   sim_peaks[i][j]     → simulated peak for event i, Pareto member j
obs_peaks = []
sim_peaks = []   # shape: (n_events, n_pareto_sets)

for date in all_dates:
    # Observed peak for this event
    obs_peak = obs_data[date]['OBS runoff [CMS]'].max()
    obs_peaks.append(obs_peak)

    # Simulated peaks from each Pareto member for this event
    event_sim_peaks = []
    for ps in PARETO_SETS:
        storm_df = ensemble_results[ps['label']][date]
        sim_peak = storm_df['SWMM outflow [CMS]'].max()
        event_sim_peaks.append(sim_peak)
    sim_peaks.append(event_sim_peaks)

obs_peaks = np.array(obs_peaks)                   # shape: (n_events,)
sim_peaks = np.array(sim_peaks)                   # shape: (n_events, n_pareto_sets)

# Min and Max simulated peak across Pareto members for each event
sim_min = sim_peaks.min(axis=1)
sim_max = sim_peaks.max(axis=1)

n_events    = len(all_dates)
n_pareto    = len(PARETO_SETS)
x_positions = np.arange(n_events)

# Format date labels: '2016_01_08' -> '08/01/2016'
date_labels = [f"{d[8:10]}/{d[5:7]}/{d[:4]}" for d in all_dates]


# ─────────────────────────────────────────────────────────────────────────────
# HELPER: draw the per-event uncertainty figure on given axes
# ─────────────────────────────────────────────────────────────────────────────
def draw_peak_uncertainty(ax, x_pos, date_labels, obs_peaks, sim_peaks,
                          sim_min, sim_max, pareto_sets, title_suffix=''):
    """
    Draw the per-event peak-flow uncertainty chart onto ax.

    Parameters
    ----------
    ax           : matplotlib Axes
    x_pos        : array of x-axis positions (integers)
    date_labels  : list of formatted date strings for x-tick labels
    obs_peaks    : array, observed peak per event
    sim_peaks    : 2-D array (n_events x n_pareto_sets), simulated peaks
    sim_min/max  : arrays, min/max across Pareto members per event
    pareto_sets  : list of dicts with 'label' and 'color'
    title_suffix : appended to the subplot title
    """
    # ── Uncertainty band: vertical line from min to max per event ────────────
    for i, x in enumerate(x_pos):
        ax.vlines(x=x,
                  ymin=sim_min[i],
                  ymax=sim_max[i],
                  colors='#aaaaaa',
                  linewidths=3.5,
                  zorder=2,
                  label='_nolegend_')

    # ── Individual Pareto member peaks: colored dots ─────────────────────────
    # Use horizontal jitter within each event so overlapping dots are visible
    n_p = len(pareto_sets)
    jitter = np.linspace(-0.18, 0.18, n_p)   # small offset per member

    for j, ps in enumerate(pareto_sets):
        x_jittered = x_pos + jitter[j]
        ax.scatter(x_jittered, sim_peaks[:, j],
                   color=ps['color'],
                   s=70,
                   zorder=4,
                   edgecolors='white',
                   linewidths=0.6,
                   label=ps['label'])

    # ── Observed peak: large black star, highly visible ───────────────────────
    ax.scatter(x_pos, obs_peaks,
               marker='*',
               s=220,
               color='#000000',
               zorder=6,
               label='Observed peak',
               edgecolors='white',
               linewidths=0.4)

    # ── Axes formatting ───────────────────────────────────────────────────────
    ax.set_xticks(x_pos)
    ax.set_xticklabels(date_labels, rotation=45, ha='right', fontsize=11)
    ax.set_xlim(-0.6, len(x_pos) - 0.4)
    ax.set_ylim(bottom=0)
    ax.set_ylabel('Peak discharge [m$^3$ s$^{-1}$]', fontsize=14)
    ax.set_xlabel('Storm event', fontsize=14)
    ax.tick_params(axis='y', labelsize=13)
    ax.grid(True, axis='y', linestyle=':', alpha=0.55, color='#888888')
    ax.set_axisbelow(True)

    # Darken frame
    for spine in ax.spines.values():
        spine.set_edgecolor('#333333')
        spine.set_linewidth(1.0)

    # ── Subtitle showing which ordering is used ───────────────────────────────
    ax.set_title(title_suffix, fontsize=12, pad=6, style='italic', color='#444444')


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 1: Chronological order
# ─────────────────────────────────────────────────────────────────────────────
fig1, ax1 = plt.subplots(figsize=(16, 6))
plt.rcParams['font.family'] = 'Arial'

draw_peak_uncertainty(
    ax=ax1,
    x_pos=x_positions,
    date_labels=date_labels,
    obs_peaks=obs_peaks,
    sim_peaks=sim_peaks,
    sim_min=sim_min,
    sim_max=sim_max,
    pareto_sets=PARETO_SETS,
    title_suffix='Events ordered chronologically'
)

# Shared legend above the axes
legend_handles = []
for ps in PARETO_SETS:
    legend_handles.append(
        mlines.Line2D([], [], color=ps['color'],
                      marker='o', markersize=8, linestyle='None',
                      label=ps['label'])
    )
legend_handles.append(
    mlines.Line2D([], [], color='black',
                  marker='*', markersize=12, linestyle='None',
                  label='Observed peak')
)
# Uncertainty band element
legend_handles.append(
    mlines.Line2D([], [], color='#aaaaaa',
                  linewidth=5, label='Uncertainty band (min–max)')
)

ax1.legend(handles=legend_handles,
           loc='upper left',
           fontsize=10,
           framealpha=0.9,
           ncol=2)

# Main title
fig1.suptitle(
    'Per-Event Peak-Flow Uncertainty Across 4 Pareto-Optimal Parameter Sets\n'
    '(vertical bar = min/max simulated range; dots = individual members; '
    'star = observed)',
    fontsize=13, y=1.01, fontweight='normal'
)

plt.tight_layout()
chron_path = os.path.join(PLOTS_DIR, 'peak_uncertainty_per_event.png')
fig1.savefig(chron_path, dpi=300, bbox_inches='tight')
print(f'Saved: {chron_path}')
plt.close(fig1)


# ─────────────────────────────────────────────────────────────────────────────
# FIGURE 2: Same plot but events SORTED by observed peak magnitude (ascending)
#           This ordering reveals whether uncertainty grows with storm size.
# ─────────────────────────────────────────────────────────────────────────────
sort_idx     = np.argsort(obs_peaks)
obs_sorted   = obs_peaks[sort_idx]
sim_sorted   = sim_peaks[sort_idx, :]
sim_min_sort = sim_sorted.min(axis=1)
sim_max_sort = sim_sorted.max(axis=1)
labels_sorted = [date_labels[i] for i in sort_idx]

fig2, ax2 = plt.subplots(figsize=(16, 6))

draw_peak_uncertainty(
    ax=ax2,
    x_pos=x_positions,
    date_labels=labels_sorted,
    obs_peaks=obs_sorted,
    sim_peaks=sim_sorted,
    sim_min=sim_min_sort,
    sim_max=sim_max_sort,
    pareto_sets=PARETO_SETS,
    title_suffix='Events ordered by observed peak magnitude (ascending)'
)

ax2.legend(handles=legend_handles,
           loc='upper left',
           fontsize=10,
           framealpha=0.9,
           ncol=2)

fig2.suptitle(
    'Per-Event Peak-Flow Uncertainty Across 4 Pareto-Optimal Parameter Sets\n'
    '(vertical bar = min/max simulated range; dots = individual members; '
    'star = observed)',
    fontsize=13, y=1.01, fontweight='normal'
)

plt.tight_layout()
sorted_path = os.path.join(PLOTS_DIR, 'peak_uncertainty_per_event_sorted.png')
fig2.savefig(sorted_path, dpi=300, bbox_inches='tight')
print(f'Saved: {sorted_path}')
plt.close(fig2)


# ─────────────────────────────────────────────────────────────────────────────
# CONSOLE SUMMARY: per-event peak values for quick verification
# ─────────────────────────────────────────────────────────────────────────────
import pandas as pd

rows = []
for i, date in enumerate(all_dates):
    row = {'Date': date_labels[i], 'Obs peak': round(obs_peaks[i], 2)}
    for j, ps in enumerate(PARETO_SETS):
        short = ps['label'].split('-')[0].strip()   # 'P1', 'P2', etc.
        row[short] = round(sim_peaks[i, j], 2)
    row['Sim min'] = round(sim_min[i], 2)
    row['Sim max'] = round(sim_max[i], 2)
    row['Band / obs [%]'] = round((sim_max[i] - sim_min[i]) / max(obs_peaks[i], 0.1) * 100, 1)
    rows.append(row)

summary_df = pd.DataFrame(rows)
print('\nPer-event peak discharge summary (all values in m^3/s):')
print(summary_df.to_string(index=False))

# Save the numerical summary alongside the plots
summary_path = os.path.join(
    BASE_DIR, 'outputs', 'peak_uncertainty_per_event_summary.csv'
)
summary_df.to_csv(summary_path, index=False)
print(f'\nNumerical summary saved: {summary_path}')
print('\nDone.')
