# Urban Runoff Response to Climate Change and Urbanization

Code repository accompanying the peer-reviewed paper:

> **"Urban Runoff Response to Climate Change and Urbanization in the Eastern Mediterranean: A SWMM Ensemble Approach"**  
> Nahal Ra'anana catchment, Israel

## Overview

This repository provides the complete analysis pipeline for a study that quantifies how concurrent climate change and urban expansion alter flood peak discharge and runoff volume in a Mediterranean urban watershed. The study uses:

- **IMS weather-radar rainfall** (bias-corrected, gridded at 1-km resolution) as observed forcing
- **WRF-PGW downscaled rainfall** (pseudo-global-warming, ~4 km resolution) as future climate forcing
- **CMIP5 and CMIP6 projections** for large-scale precipitation trend context
- **SWMM** (EPA Storm Water Management Model 5.1) as the hydrological engine
- **Multi-objective Pareto calibration** (leave-one-out cross-validation over 23 storm events)
- **Ensemble uncertainty propagation** through 12 Pareto-optimal parameter sets

The analysis covers 23 historical storm events (2012–2020), 41 WRF-simulated events, 2 urbanization scenarios, and 81 spatial rainfall shifts per event — producing ~79,700 SWMM runs in total.

---

## Repository Structure

```
UrbanRunoffModeling_PUBLIC/
├── main_analysis/               ← Primary analysis (calibration, sensitivity, climate impact)
│   ├── 01_rainfall_runoff_data/ ← Data extraction and preprocessing
│   ├── 02_swmm_calibration/     ← Model calibration and cross-validation
│   ├── 03_sensitivity_analysis/ ← Variance-based and spatial sensitivity analyses
│   ├── 04_impact_of_climate_change_and_urbanization/  ← Main scenario analysis
│   ├── cmip_projection/         ← Large-scale CMIP5/CMIP6 precipitation context
│   ├── Pareto_Uncertainty_Analysis/   ← Initial Pareto ensemble uncertainty analysis
│   └── scripts/                 ← Utility scripts for data export
│
└── pareto_uncertainty_analysis/ ← Peer-review addition: refactored Pareto pipeline
    ├── urban_runoff/            ← Clean Python package (calibration, optimization, scenarios)
    ├── 05_Climate_Change_Impact/← Ensemble climate change analysis notebooks
    ├── outputs/                 ← Pre-computed results (Pareto CSV, calibrated model, figures)
    └── 01–04_*.ipynb            ← Refactored orchestration notebooks
```

See the `README.md` inside each subfolder for detailed descriptions.

---

## Data Availability

| Data type | Source | Availability |
|-----------|--------|--------------|
| IMS weather-radar rainfall | Israel Meteorological Service | **Restricted** — requires IMS data agreement |
| SERS stream discharge | Surface and Estuarial Research Station | **Restricted** — requires data agreement |
| Municipal GIS / catchment delineation | Ra'anana Municipality | **Restricted** — not redistributable |
| WRF-PGW rainfall (historical + future) | Regional climate modeling (in-house) | **Restricted** — contact authors |
| CMIP5 / CMIP6 precipitation | [ESGF / KNMI Climate Explorer](https://climexp.knmi.nl) | **Public** |
| GHSL urban extent data | [JRC Global Human Settlement Layer](https://ghsl.jrc.ec.europa.eu) | **Public** |
| Calibrated SWMM baseline model | This repository | `pareto_uncertainty_analysis/outputs/models/calibrated_baseline.inp` |
| Pareto ensemble parameter sets | This repository | `pareto_uncertainty_analysis/outputs/pareto/final/` |
| Processed simulation results | This repository | CSV files in `outputs/` subdirectories |

Scripts that require restricted data will fail at the data-loading step and print an informative error. All analysis code is fully reproducible given the inputs described above.

---

## Software Dependencies

### Python

Tested with **Python 3.9.13** on Windows 10/11.

```
numpy
pandas
scipy
matplotlib
swmm-api==0.3.2        # SWMM input/output file handling
pyswmm                 # SWMM simulation runner
hydroeval==0.1.0       # KGE / NSE evaluation
SALib                  # Variance-based sensitivity analysis
netCDF4                # CMIP NetCDF file reading
xarray
```

Install all dependencies:
```bash
pip install numpy pandas scipy matplotlib "swmm-api==0.3.2" pyswmm "hydroeval==0.1.0" SALib netCDF4 xarray
```

### MATLAB

MATLAB R2019b or later is required for the radar data extraction and WRF pre-processing
scripts in:
- `main_analysis/01_rainfall_runoff_data/Radar_Files_Extraction/`
- `main_analysis/04_impact_of_climate_change_and_urbanization/wrf_rainfall_analysis/`

MATLAB scripts are upstream of the Python pipeline and only need to be re-run if you have
access to the raw restricted data files.

### SWMM

EPA SWMM 5.1 must be installed. Download from:
https://www.epa.gov/water-research/storm-water-management-model-swmm

---

## How to Run

### Using pre-computed results (no restricted data needed)

The key calibration and simulation outputs are already provided:

- **Calibrated SWMM model:** `pareto_uncertainty_analysis/outputs/models/calibrated_baseline.inp`
- **Pareto ensemble parameters:** `pareto_uncertainty_analysis/outputs/pareto/final/pareto_ensemble_full.csv`
- **Climate scenario figures:** `pareto_uncertainty_analysis/outputs/climate/figures/`
- **CMIP precipitation analysis:** run `main_analysis/cmip_projection/CMIP5_CMIP6_combined_analysis.ipynb`
  (uses publicly available ESGF/KNMI data)

### Full reproduction from raw data (requires restricted data)

```
Stage 1 — Data preparation:
  Run MATLAB scripts in main_analysis/01_rainfall_runoff_data/Radar_Files_Extraction/
  Then run Basin_Radar_Overlap.ipynb

Stage 2 — Calibration:
  Run main_analysis/02_swmm_calibration/ notebooks 01 → 04 in order

Stage 3 — Sensitivity:
  Run main_analysis/03_sensitivity_analysis/ notebooks (independent of each other)

Stage 4 — Climate scenarios:
  Run main_analysis/04_impact_of_climate_change_and_urbanization/ notebooks 01 → 02

Stage 5 (peer-review ensemble):
  Run pareto_uncertainty_analysis/ notebooks 01 → 04
  Then run pareto_uncertainty_analysis/05_Climate_Change_Impact/ notebooks 01 → 03
```

---

## Notes on Hardcoded Paths

Several scripts in `main_analysis/` contain absolute paths pointing to the original
research data directory (`D:\Development\RESEARCH\Raanana\`). These must be updated
to your local data directory before running. Affected files are documented in each
subfolder's README.

The `pareto_uncertainty_analysis/urban_runoff/config.py` module uses environment variables
for all data paths, with the hardcoded paths as fallbacks. Override them with:

```bash
export URB_RUNOFF_CV_DIR="path/to/Cross_validation"
export URB_RUNOFF_RESULTS_DIR="path/to/SWMM"
export URB_RUNOFF_RAIN_PKL="path/to/Basin_radar_overlap_pkl.pkl"
```

---

## License

Code is released under the MIT License. See `LICENSE`.
Data files are subject to the data agreements described above and are not covered by this license.

## Contact

For questions about data access or methods, contact the corresponding author.
