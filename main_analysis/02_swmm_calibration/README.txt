# Overview

This repository contains Python scripts for SWMM model calibration and validation using cross-validation, multi-objective functions, and parameter optimization.

## Workflow

### 1. Cross-Validation (`01_swmm_cross_validation.ipynb`)
- Performs leave-one-out cross-validation on storm events.
- Explores parameter combinations scaled by predefined factors.
- Outputs a DataFrame with performance metrics for all iterations.

### 2. Cross-Validation Concatenation (`01a_cross_validation_concat_runs.ipynb`)
- Aggregates cross-validation results for further analysis.

### 3. Objective Function Selection (`02_swmm_cross_validation_multobj.ipynb`)
- Loads the DataFrame from Step 1.
- Evaluates parameter performance using multi-objective functions.
- Determines the best set of objective functions for calibration.

### 4. Full Dataset Calibration (`03_swmm_calibration_after_cv.ipynb`)
- Runs parameter combinations using all storm events.
- Outputs a DataFrame with results for all parameter combinations.

### 5. Final Parameter Optimization (`04_swmm_calibration_after_cv_multobj.ipynb`)
- Identifies calibrated parameters using the DataFrame from Step 3 and objective functions from Step 2.
- Outputs the final calibrated parameters.

## Notes
- Scripts must be executed in order.
- Parameters and objective functions can be customized as needed.
- Additional experimental scripts can be found in the `Ensemble_Experiments/` directory.

