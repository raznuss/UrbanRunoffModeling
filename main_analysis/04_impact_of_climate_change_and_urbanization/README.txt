# WRF to SWMM Data Processing and Rain Shift Analysis

This repository contains scripts for processing WRF (Weather Research and Forecasting) model outputs, preparing them for hydrological simulations in SWMM (Storm Water Management Model), and analyzing the impact of urbanization and climate change on hydrological responses.

## 1. WRF Output Data Preparation
This section includes scripts for converting WRF model outputs stored in MAT-files into a format suitable for SWMM simulations.

### Scripts:
1. **wrf_matfile_adjustment.m**
   - Converts MAT-file structures into a readable format for further processing in Python.

2. **matfile_to_pickle.ipynb**
   - Converts the MAT-files into pickle files for easier handling in Python.
   - The pickle files store the following data:
     - `lat`, `lon`: Latitude and longitude values.
     - `rainRatePGW`: Rainfall rate under the Pseudo Global Warming (PGW) scenario.
     - `timesList`: List of timestamps for rainfall data.
     - `totalRainPGW`: Total accumulated rainfall in the PGW scenario.

3. **wrf_rainshift_Timeseries_txtfile_creator__historical.ipynb**
   - Processes the extracted data to match SWMM input file requirements for historical rainfall scenarios.
   - Converts rainfall time series into text files compatible with SWMM.

4. **wrf_rainshift_Timeseries_txtfile_creator__future.ipynb**
   - Similar to the historical script, but processes future rainfall scenarios under climate change conditions.
   - Formats rainfall time series for SWMM input.

---

## 2. Rain Shift Analysis for SWMM Simulations
This section includes scripts for running multiple SWMM model simulations with an ensemble of spatially shifted rainfall scenarios. These simulations analyze hydrological responses under different levels of urbanization and climate change scenarios.

### Scripts:
1. **01_scenarios_rain_shift_run_to_pickle**
   - Runs multiple SWMM simulations, considering different urbanization levels and climate scenarios (historical vs. future), and stores results in pickle files.

2. **02a_paired_shifted_events_load_and_visu**
   - Loads and visualizes paired rainfall events with identical urbanization levels but different climate conditions to compare hydrological responses.

3. **02b_visualization_discharge_param**
   - Analyzes and visualizes discharge-related parameters to assess hydrological differences between scenarios.

4. **02c_visualization_rain_param**
   - Visualizes spatial and temporal variations in rainfall characteristics across different scenarios.

### Additional:
- **more_visu/**: Contains supplementary visualization scripts for further analysis.

---

## Summary
This repository provides a complete workflow from WRF model outputs to SWMM simulations, allowing for the assessment of hydrological impacts under urbanization and climate change. The ensemble of rain shift scenarios enables a comprehensive evaluation of spatial variability in storm events and their influence on watershed responses.

