# Pareto Uncertainty Analysis — Refactored Pipeline

This folder contains the refactored calibration pipeline and the 12-member Pareto ensemble
analysis added in response to peer review. It corresponds to the `UrbanRunoffModeling_Refactored`
repository.

The core scientific contribution here is the **propagation of parameter uncertainty** from the
full Pareto-optimal ensemble (12 non-dominated solutions) through the complete climate change
scenario analysis. This extends the main analysis (`main_analysis/`) which used a single
calibrated parameter set.

---

## What Is Here and Why

During peer review, reviewers requested that the climate change results be conditioned on
parameter uncertainty rather than a single optimal parameter set. In response, this refactored
pipeline was built to:

1. **Re-run calibration** in a clean, testable Python package (`urban_runoff/`)
2. **Preserve all 12 Pareto-front members** rather than selecting a single optimum
3. **Run the 41-event WRF scenario ensemble** through all 12 Pareto members (~79,700 SWMM runs)
4. **Quantify and visualize parametric uncertainty** in the CDF comparisons

Results produced here are the primary uncertainty figures in the peer-reviewed paper.

---

## Folder Structure

```
pareto_uncertainty_analysis/
├── urban_runoff/               Core Python package
│   ├── config.py               Parameter grid, paths, objective definitions
│   ├── calibration/            Cross-validation engine, objective functions
│   ├── data/                   Data loaders and preprocessors
│   ├── optimization/           Pareto front detection, ensemble selection
│   ├── scenarios/              Urbanization scenario export
│   ├── swmm/                   SWMM runner, inp editor, output parser
│   └── visualization/          Plot functions for Pareto and calibration figures
│
├── 01_Generate_LOOCV_Cache.ipynb        Stage 1: run LOOCV calibration
├── 02_Analyze_CV_and_Select_Objectives.ipynb  Stage 2: select objective set
├── 03_Final_Calibration.ipynb           Stage 3: full-dataset Pareto calibration
├── 04_Pareto_Uncertainty_Analysis.ipynb Stage 4: parameter uncertainty on observed events
├── 04_Prepare_WRF_Scenarios.ipynb       Stage 4b: urbanization scenario CSV export
│
├── 05_Climate_Change_Impact/            Stage 5: climate scenario ensemble
│   ├── 05_01_SWMM_Ensemble_Runner.ipynb   Run all 12 members × 41 events × 81 shifts
│   ├── 05_02_CDF_Uncertainty_Analysis.ipynb   CDF figures with Pareto uncertainty bands
│   ├── 05_03_Sensitivity_Analysis.ipynb   Rain threshold and shift threshold sensitivity
│   └── 06_Spatial_Rainfall_Analysis.ipynb   Spatial rainfall structure analysis
│
├── outputs/                    All generated outputs (pre-computed versions provided)
│   ├── pareto/final/           pareto_ensemble_full.csv, pareto_optimum.csv
│   ├── models/calibrated_baseline.inp   Calibrated SWMM model (baked-in parameters)
│   ├── models/scenarios/       Urbanization scenario discharge CSVs
│   ├── results/loocv/          LOOCV fold results (pickle files)
│   ├── results/final_calibration/   Full-dataset calibration cache
│   ├── climate/figures/        Paper figures (PDF)
│   ├── climate/raw/            Per-member simulation pickles (one per Pareto index)
│   ├── metrics/                Calibration diagnostic plots (PNG)
│   └── figures/spatial_rainfall/   Spatial rainfall plots per WRF event
│
├── generate_calibrated_baseline.py   One-time script: writes calibrated_baseline.inp
├── inspect_calibrated_baseline.py    Diagnostic: prints subcatchment parameter table
├── validate_phase1_2.py        Test suite: config, data loading (15 tests)
├── validate_phase3_4.py        Test suite: SWMM execution, objectives (287 tests)
└── validate_phase5_6.py        Test suite: Pareto, selection, scenarios (5 tests)
```

The `_deprecated/` folder contains early notebook drafts superseded by the above.

---

## The `urban_runoff` Python Package

A clean, modular Python package that re-implements the legacy calibration logic from
`main_analysis/02_swmm_calibration/` with full test coverage and no global state.

All paths are configurable via environment variables (see `config.py`). The package
never touches the legacy data directories.

### Key modules

**`config.py`** — central configuration. Defines:
- `PARAMETER_GRID`: 8 calibration parameters with their sweep ranges (`Fixed`, `Sweep`,
  or `Values` spec); currently 1,323 combinations (7×7×3×1×1×1×9×1)
- `PARETO_OBJECTIVES`: the 4 active objectives (`peak_1-kge`, `volume_1-kge`,
  `peak_rmsd`, `peak_bias`) — change here only to alter the active set
- `FILTER_BOUNDS`: post-run parameter filter applied before Pareto evaluation
- Data path constants (overridable via environment variables)

**`calibration/cross_validation.py`** — builds the factor-combination matrix and
runs all SWMM simulations with checkpoint-aware caching.

**`calibration/objectives.py`** — computes NSE, KGE (via `hydroeval 0.1.0`), RMSD,
bias, and MAD in minimize-form (1-NSE, 1-KGE, etc.).

**`calibration/loocv.py`** — leave-one-out cross-validation: holds out one event,
trains on the remaining 22, evaluates on the held-out event.

**`optimization/pareto.py`** — vectorized Pareto dominance (standard: all ≤ + any <)
with normalization before comparison.

**`swmm/inp_editor.py`** — applies calibrated factors to a base SWMM `.inp` file.
Three parameter types: multiplicative from array (Type A), multiplicative from `.inp`
current value (Type B), or absolute set (Type C).

**`swmm/runner.py`** — wraps `swmm5_run` to execute a single SWMM simulation and
return the binary output path.

---

## Calibration Workflow

### Stages 1–3 (require restricted data)

These notebooks require access to the SWMM event directories (restricted).
Pre-computed outputs are already in `outputs/` so re-running is optional.

**Stage 1 — LOOCV (`01_Generate_LOOCV_Cache.ipynb`)**
- Runs 1,323 parameter combinations × 23 events × 23 LOOCV folds
- Each fold withholds one event, calibrates on the remaining 22
- Saves per-fold results as `outputs/results/loocv/fold_NN_*.pkl`
- Saves combined results as `outputs/results/loocv/all_folds_*.pkl`

**Stage 2 — Objective Selection (`02_Analyze_CV_and_Select_Objectives.ipynb`)**
- Evaluates 4 candidate objective sets on the LOOCV validation events
- Selects the winning set: `['peak_1-kge', 'volume_1-kge', 'peak_rmsd', 'peak_bias']`
- Produces `outputs/metrics/loocv_*.png` diagnostic plots

**Stage 3 — Final Calibration (`03_Final_Calibration.ipynb`)**
- Runs all 1,323 combinations on the full 23-event dataset
- Builds the Pareto front (12 non-dominated solutions at the time of the peer-review run)
- Saves: `outputs/pareto/final/pareto_ensemble_full.csv` (all Pareto members)
  and `outputs/pareto/final/pareto_optimum.csv` (single Euclidean optimum)
- Run `python generate_calibrated_baseline.py` afterward to create `calibrated_baseline.inp`

### Stage 4 — Pareto Uncertainty on Observed Events (`04_Pareto_Uncertainty_Analysis.ipynb`)

Runs all 12 Pareto members on the 23 calibration events and visualizes the uncertainty
envelope (min–max band across members) against observed hydrographs.

Key statistics reported in the paper:
- P-factor (fraction of observed peaks within the Pareto band): 17.4% for peaks, 26.1% for volume
- Mean band width: ~0.50 m³/s

### Stage 5 — Climate Change Impact (requires restricted WRF data)

**`05_01_SWMM_Ensemble_Runner.ipynb`**
Runs 12 Pareto members × 41 WRF events × 2 conditions (historical/future) × 81 spatial
shifts × 2 urbanization multipliers. One pickle per Pareto member is written to
`outputs/climate/raw/pareto_{idx}_results.pkl`. Pre-computed pickles are provided
for all 12 members.

**`05_02_CDF_Uncertainty_Analysis.ipynb`**
Loads the 12-member ensemble pickles and produces:
- CDF of peak discharge and runoff volume: historical vs. future, with Pareto uncertainty band
- Percentile decomposition bar charts (climate change vs. urbanization contribution)

Active filters matching the main analysis:
- `RAIN_FILTER_MM = 10.0 mm` — removes events where all shifts yield < 10 mm
- `SHIFT_THRESHOLD = 20,000 m` — restricts to ±20 km y-shifts

**`05_03_Sensitivity_Analysis.ipynb`**
Sweeps the rain-threshold and shift-threshold filters to test robustness of the
historical-vs-future difference signal.

**`06_Spatial_Rainfall_Analysis.ipynb`**
Analyses spatial structure of WRF rainfall across the 41 events; produces per-event
spatial maps saved to `outputs/figures/spatial_rainfall/`.

---

## Pre-Computed Outputs (No Restricted Data Needed)

These files are included and can be used directly:

| File | Description |
|------|-------------|
| `outputs/pareto/final/pareto_ensemble_full.csv` | All 12 Pareto members: factors + objectives |
| `outputs/pareto/final/pareto_optimum.csv` | Single Euclidean-optimal solution |
| `outputs/models/calibrated_baseline.inp` | SWMM model with calibrated parameters baked in |
| `outputs/models/scenarios/*.csv` | Urbanization scenario discharge results |
| `outputs/climate/figures/*.pdf` | Paper figures (CDF, percentile bars, sensitivity) |
| `outputs/metrics/*.png` | Calibration diagnostic plots |

---

## Running the Test Suite

All 307 tests verify methodological identity with the legacy calibration notebooks:

```bash
cd pareto_uncertainty_analysis/
python validate_phase1_2.py   # 15 tests — config, data loading, preprocessing
python validate_phase3_4.py   # 287 tests — SWMM execution, objective functions
python validate_phase5_6.py   # 5 tests  — Pareto front, selection, scenario export
```

A test failure indicates a methodological regression and should be treated as a blocker.

---

## Hardcoded Paths

`generate_calibrated_baseline.py` contains a hardcoded path to the original SWMM
event directory that is also used as a fallback in `urban_runoff/config.py`. This path
(`D:\Development\RESEARCH\Raanana\SWMM\from_radar\Cross_validation`) must be overridden
via the `URB_RUNOFF_CV_DIR` environment variable or updated in `config.py` before running
any calibration notebook from scratch. The pre-computed outputs do not require this path.

---

## Relationship to Main Analysis

| Aspect | `main_analysis/` | `pareto_uncertainty_analysis/` |
|--------|-----------------|-------------------------------|
| Calibration code | Legacy notebooks (one large notebook per stage) | Modular `urban_runoff` package with test suite |
| Parameter grid | 480 combinations | 1,323 combinations |
| Pareto members used | 1 (optimum) or 4 (initial uncertainty) | 12 (full ensemble) |
| Climate scenarios | Single parameter set | 12-member Pareto band |
| Primary purpose | Original study analysis | Peer-review uncertainty quantification |
