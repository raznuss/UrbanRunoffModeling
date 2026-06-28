# Main Analysis

This folder contains the primary analysis pipeline for the paper, organized in sequential
processing stages from raw data ingestion through climate change impact visualization.

It corresponds to the `UrbanRunoffModeling` repository as developed during the original study.

---

## Folder Structure and Processing Order

```
main_analysis/
├── 01_rainfall_runoff_data/           Stage 1: data extraction & preprocessing
├── 02_swmm_calibration/               Stage 2: model calibration
├── 03_sensitivity_analysis/           Stage 3: parameter & rainfall sensitivity
├── 04_impact_of_climate_change_and_urbanization/   Stage 4: scenario analysis
├── cmip_projection/                   Stage 5: large-scale precipitation context
├── Pareto_Uncertainty_Analysis/       Stage 6: Pareto ensemble uncertainty
└── scripts/                           Utility: data export helpers
```

Stages 1–4 must be executed in order when starting from raw data. Stages 5 and 6 are
independent of each other and can be run once Stage 2 is complete.

---

## Stage 1 — Rainfall and Runoff Data (`01_rainfall_runoff_data/`)

Prepares the observed rainfall (radar) and discharge (stream gauge) time series used as
forcing and calibration targets for SWMM.

### `Radar_Files_Extraction/` (MATLAB + Python)

Extracts storm events from IMS weather-radar raw data, applies a gauge-based bias
correction, and writes gridded rainfall time series as `.mat` files and SWMM-compatible
ASCII text files.

**Run order:**
1. `RainGaugesCorrelationRaanana.m` — determines per-event radar–gauge bias correction factors
2. `BoundStormsInRadarData.m` — identifies event start/end times in the radar archive
3. `CreateRadar_MATandASCII.m` — extracts and saves per-event `.mat` and ASCII text files
4. `MissingRadarTimeCorrection.m` — fills missing timesteps by linear interpolation
5. `Basin_Radar_Overlap.ipynb` — computes the area-weighted basin-average rainfall per
   subcatchment and saves `Basin_radar_overlap_pkl.pkl` (main input to calibration)

Bias correction factors (per-event, from gauge–radar regression) are hardcoded in
`Basin_Radar_Overlap.ipynb` and documented in `README.txt`.

> **Restricted data:** The `.mat` files in this folder contain IMS radar data and SERS
> gauge records and are not redistributable:
> `gaugeData_updated_042023.mat`, `Gauges_Radar_raanana_byevent*.mat`,
> `event_structure*.mat`.

> **Hardcoded paths:** `MissingRadarTimeCorrection.m` and `RainGaugesCorrelationRaanana.m`
> contain absolute paths to the original radar archive. Update them before running.

### `Discharge_Files/` (Python notebooks)

Pre-processes stream discharge records from the SERS gauging station.

- `filter_observation_event .ipynb` — filters raw discharge time series to storm events
  above the 6.5 m³/s threshold used throughout the study
- `raanana_hydrodaya_encoding.ipynb` — fixes character-encoding issues in the raw CSV export
- `split_storms_by_cumulative_rain.ipynb` — separates overlapping discharge records into
  individual storm events using cumulative rainfall as the segmentation criterion

> **Restricted data:** Raw discharge CSV files (SERS data) are not included in this
> repository. They are loaded by these notebooks from the local data directory.

---

## Stage 2 — SWMM Calibration (`02_swmm_calibration/`)

Calibrates the 28-subcatchment SWMM model of Nahal Ra'anana using 23 storm events
(2012–2020) and a 480-combination grid search over 8 parameters (imperviousness,
depression storage, flow width, Manning's n, CN, percent zero-storage, percent routed,
and daily evaporation).

**Notebooks (run in order):**

| Notebook | Purpose |
|----------|---------|
| `01_swmm_cross_validation.ipynb` | Leave-one-out cross-validation (LOOCV): runs all 480 parameter combinations with each event withheld once |
| `01a_cross_validation_concat_runs.ipynb` | Concatenates per-chunk LOOCV pickle files into a single results DataFrame |
| `02_swmm_cross_validation_multobj.ipynb` | Evaluates multi-objective performance (NSE, KGE, RMSD, bias); selects the objective set used for Pareto calibration |
| `03_swmm_calibration_after_cv.ipynb` | Runs all 480 combinations on the full 23-event dataset |
| `04_swmm_calibration_after_cv_multobj.ipynb` | Builds the Pareto front from full-dataset results; identifies the calibrated parameter set |

**Other files:**
- `Ensemble_Experiments/Ensemble_Optimizaztion.ipynb` — exploratory ensemble approach, not used in the final paper
- `Radar_Hydrograph_gif.m` — generates an animated GIF of radar rainfall over the basin
- `storm_data.mat` — intermediate MATLAB data file derived from radar processing; input to calibration notebooks

> **Hardcoded paths:** Calibration notebooks load pickle files from the original data
> directory. Update path variables in the first cell of each notebook before running.

---

## Stage 3 — Sensitivity Analysis (`03_sensitivity_analysis/`)

Three complementary analyses assess how model outputs respond to parameter and rainfall uncertainty.

### `global_sa/` — Global response surface

A full-factorial grid search over SWMM parameters to visualize the response surface
of peak discharge and runoff volume.

- `01_global_sa_runs.ipynb` — runs SWMM for a wide parameter grid; saves results as pickles
- `02_global_sa_visualization.ipynb` — loads results and plots response surfaces

### `variance_based_sa/` — Variance-based sensitivity

Quantifies the contribution of individual parameters and rain-scaling factors to output
variance using a variance-decomposition approach.

- `sa_time_series_txtfile_creator.ipynb` — applies Box-Cox transformation and uniform
  rain-factor scaling to radar time series; writes SWMM input text files
- `vb_rainfactor.ipynb` / `vb_coxbox.ipynb` — runs SWMM with the scaled inputs; saves pickles
- `vb_rainfactor_pickle_to_df.ipynb` / `vb_coxbox_pickle_to_df.ipynb` — aggregates chunks,
  computes variance decomposition, displays final results
- `vb_norain.ipynb` / `vb_norain_pickle_to_df.ipynb` — same workflow for the no-rainfall
  (static parameter only) sensitivity baseline

### `rain_shifts_spatial_sa/` — Spatial rainfall sensitivity

Evaluates how spatial displacement of the storm cell (x/y shifts of the radar field)
changes simulated peak discharge across the basin.

- `01_rain_shift_txtfile_creator.ipynb` — shifts the gridded radar field in x and y;
  writes SWMM input text files for each shift combination
- `02_rain_shift_run_to_plot.ipynb` — runs SWMM for all shifted fields and generates the
  spatial sensitivity map of peak discharge

### `Archive/` — Draft notebooks

Early exploratory sensitivity scripts superseded by the above. Retained for reference only.

---

## Stage 4 — Climate Change and Urbanization Impact (`04_impact_of_climate_change_and_urbanization/`)

The core scenario analysis of the paper. Compares hydrological response under historical
vs. future (WRF-PGW) rainfall, crossed with two urbanization levels (current ~55% and
projected ~81% impervious).

### `01_data_processing/`

Converts WRF model output (MAT-files produced by MATLAB) into SWMM-compatible
time series. Run these before the scenario notebooks.

- `wrf_matfile_adjustment.m` — converts WRF MAT-file structure for Python compatibility
- `matfile_to_pickle.ipynb` — reads adjusted MAT-files and saves Python pickles
  (lat/lon grid, rainfall rate, timestamps)
- `wrf__rainshift_Timeseries_txtfile_creator__historical.ipynb` — applies area-weighting
  and bias correction; writes SWMM input text files for historical WRF events
- `wrf__rainshift_Timeseries_txtfile_creator__future.ipynb` — same for future PGW events

> `drafts/` contains exploratory pre-processing scripts not used in the final workflow.

### `02_scanrios_comparsion/` — Scenario runs and visualization

- `01_scanrios_rain_shift_run_to_pickle.ipynb` — runs SWMM for all historical and future
  WRF events across the 81-shift ensemble and 2 urbanization multipliers; saves pickles
- `01_scanrios_rain_shift_run_to_pickle_adding_infiltration.ipynb` — variant that also
  modifies the Green-Ampt infiltration parameters
- `02a_paired_shifted_events_load_and_visu.ipynb` — loads paired historical/future
  results and produces the main 3-panel CDF comparison figure (primary paper figure)
- `02b_visualization_discharge_param.ipynb` — discharge parameter statistics and CDFs
- `02c_wrf_visualization_rain_param.ipynb` — spatial and temporal WRF rainfall characteristics
- `02d_swmm_visualization_rain_param.ipynb` — precipitation statistics at the SWMM input level
- `filter_data_for_efrat.ipynb` — exports selected scenario results to CSV
- `areal_coverage_by_subcatchment/` — per-subcatchment areal rainfall coverage analysis

Published figure PDFs (`Fig_3_Split_Final_Robust_v2.pdf`, `Nature_Style_Figure_180mm.pdf`, etc.)
are included as reference outputs from the paper.

### `wrf_rainfall_analysis/` — Spatial rainfall structure

Analyses the areal coverage and wet-fraction of WRF rainfall within the basin and
the full WRF domain.

- `wrf_basin_analysis.m`, `wrf_full_domain_analysis.m` — compute per-event statistics over the basin/domain
- `compute_basin_wet_fraction.m`, `compute_domain_wet_fraction.m` — wet-area fraction per threshold
- `grid_overlap_raanana.m` — maps WRF grid cells to SWMM subcatchment boundaries
- `rainfall_coverage_per_pixle.ipynb` — Python visualization of areal coverage results

### `other_approaches/` — Exploratory analyses

Alternative approaches considered during the study (ensemble sizing, Wilcoxon tests,
visualization animations). Not part of the published workflow; retained for transparency.

---

## Stage 5 — CMIP Precipitation Context (`cmip_projection/`)

Provides large-scale precipitation trend context for the eastern Mediterranean region,
supporting the interpretation of the WRF-PGW forcing used in Stage 4.

- `CMIP5_CMIP6_combined_analysis.ipynb` — main analysis: loads CMIP5 and CMIP6
  multi-model ensembles, computes precipitation trend statistics for the eastern
  Mediterranean, and produces CDF comparison figures. Uses **publicly available**
  data from ESGF/KNMI Climate Explorer.
- `CMIP5_prec_analysis.ipynb` — CMIP5-only analysis (precursor to the combined notebook)
- `prec_analysis.ipynb` — early exploratory precipitation analysis

---

## Stage 6 — Pareto Ensemble Uncertainty (`Pareto_Uncertainty_Analysis/`)

An initial uncertainty analysis using 4 diverse Pareto-optimal parameter sets selected
from the Stage 2 calibration. This analysis was subsequently extended to 12 Pareto
members in the peer-review response; see `pareto_uncertainty_analysis/` for the full version.

- `step1_extract_pareto.py` — reads the calibration pickle; extracts the 4 most diverse
  Pareto members (balanced, peak-KGE-optimal, volume-KGE-optimal, bias-optimal);
  saves parameter sets and indices to `outputs/`
- `pareto_uncertainty_analysis.ipynb` — runs SWMM with each of the 4 parameter sets on
  all 23 storm events; saves the ensemble time series to `outputs/ensemble_results.pkl`
- `plot_peak_uncertainty_per_event.py` — visualizes per-event peak discharge uncertainty
  as a min–max band across the 4 Pareto members alongside observed peaks

> **Hardcoded paths:** `step1_extract_pareto.py` contains an absolute path to the
> calibration pickle (`PICKLE_PATH`). `plot_peak_uncertainty_per_event.py` contains an
> absolute path to the repo root (`BASE_DIR`). Update both before running.

Pre-computed outputs are available in `outputs/`:
- `ensemble_results.pkl` — full ensemble time series (4 members × 23 events)
- `selected_pareto_sets.csv` — parameter values for the 4 selected Pareto members
- `peak_uncertainty_per_event_summary.csv` — per-event min/max/observed peak summary

---

## Utility Scripts (`scripts/`)

One-time data export helpers used to share simulation results with collaborators.

- `run_filter_38.py` — loads the 38% urbanization scenario pickle, removes rainfall
  columns, and exports discharge and event-statistics CSVs
- `run_filter_87.py` — same for the 87% urbanization scenario

> **Hardcoded paths:** Both scripts contain absolute paths (`SIM_DIR`, `output_dir`)
> pointing to `D:\Development\RESEARCH\Raanana\...`. Update before running.

---

## Auxiliary Files

- `schematic_hydrograph.ipynb` — generates a schematic hydrograph figure for the paper methods section
- `ipynb_remove_output.ipynb` — strips embedded cell outputs from notebooks before committing
