# ================================================
# Rainfall Areal Coverage Analysis – MATLAB + Python
# ================================================

## MATLAB

### compute_basin_wet_fraction.m
- Computes wet-area fraction **inside the Raanana basin** for each event and rain-intensity threshold.
- Saves:
    • basin_wet_fraction_results.mat  
    • basin_wet_fraction_events.csv

### compute_domain_wet_fraction.m
- Same calculation but for the **entire WRF domain**.
- Saves:
    • domain_wet_fraction_results.mat  
    • domain_wet_fraction_events.csv

### run_basin_domain_fraction.m
- Main driver script.
- Loads shapefile and event directories, defines thresholds,
  runs the two functions above.
- Produces two tables with **identical column structure** for easy comparison.


## Python

- Loads either the basin or domain CSV.
- Builds paired rows per **(event × threshold)**.
- Normalizes areal-coverage values.
- Computes:
    • Mean Δ normalized coverage (future − historic)  
    • Fraction of events where historical > future
- Generates a two-panel figure:
    (a) Δ normalized coverage with bootstrap CIs  
    (b) Fraction of events where historical coverage is larger

- Used for **visualization and summary statistics only**.
