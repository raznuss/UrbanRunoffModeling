"""
Validation script for Phase 1 (config.py) and Phase 2 (data/loaders.py, data/preprocessors.py).

Strategy:
  1. Define each legacy function verbatim (copy-pasted from the legacy notebooks).
  2. Call both the legacy and new implementations on identical real-data inputs.
  3. Assert strict equality of all outputs using pandas/numpy testing utilities.

Real data path used: Cross_validation/2012_01_13/
Exit code: 0 if all assertions pass, 1 if any assertion fails.
"""

import sys
import itertools
import logging

import numpy as np
import pandas as pd
from swmm_api import read_out_file

# Add the refactored package root to sys.path
sys.path.insert(0, r"D:\MY_CODES\UrbanRunoffModeling_Refactored")

from urban_runoff import config
from urban_runoff.data.loaders import (
    load_observed_runoff,
    load_swmm_output,
    load_pickle,
    load_basins_rain,
)
from urban_runoff.data.preprocessors import (
    resample_5min,
    filter_rainfall_columns,
    align_obs_sim,
)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)-5s] %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger("validate")

# ──────────────────────────────────────────────────────────────────────────────
# Real-data paths (legacy layout, read-only)
# ──────────────────────────────────────────────────────────────────────────────

EVENT_DIR      = r"D:\Development\RESEARCH\Raanana\SWMM\from_radar\Cross_validation\2012_01_13"
LEGACY_OBS     = EVENT_DIR + r"\2012_01_13"        # legacy style: no extension
NEW_OBS_CSV    = EVENT_DIR + r"\2012_01_13.csv"    # new style: full path with extension
SWMM_OUT       = EVENT_DIR + r"\Final.out"
RAIN_PKL       = (
    r"D:\Development\RESEARCH\Raanana\data\rain_radar"
    r"\Basin_radar_overlap_pkl\Basin_radar_overlap_pkl.pkl"
)


# ──────────────────────────────────────────────────────────────────────────────
# LEGACY FUNCTIONS — verbatim copy from legacy notebooks (ground-truth reference)
# ──────────────────────────────────────────────────────────────────────────────

def _legacy_load_runoff_obs_to_df(OBS_RUNOFF_DATA_PATH):
    df_obs_runoff = pd.read_csv(OBS_RUNOFF_DATA_PATH + ".csv")
    df_obs_runoff.dropna(inplace=True)
    df_obs_runoff.rename(
        columns={"discharge_cms": "OBS runoff [CMS]"}, inplace=True, errors="raise"
    )
    df_obs_runoff["date_and_time"] = pd.to_datetime(
        df_obs_runoff["date_and_time"], format="%d/%m/%Y %H:%M:%S"
    )
    df_obs_runoff = df_obs_runoff.set_index("date_and_time")
    df_obs_runoff.drop(["Unnamed: 0"], axis=1, inplace=True)
    return df_obs_runoff


def _legacy_data_5min_interpolate(df_to_interpolate):
    df_interpol = df_to_interpolate.resample("5min").mean()
    df_interpol = df_interpol.interpolate(method="linear")
    return df_interpol


def _legacy_load_sim_to_df(sim_output_path):
    sim_df = read_out_file(sim_output_path).to_frame()["system"][""][
        ["outflow", "rainfall"]
    ]
    sim_df.rename(
        columns={"outflow": "SWMM outflow [CMS]", "rainfall": "rainfall [mm/h]"},
        inplace=True,
        errors="raise",
    )
    sim_df.index = pd.to_datetime(sim_df.index)
    return sim_df


def _legacy_is_rainfall_column(col):
    """Verbatim from run_filter_38.py / run_filter_87.py."""
    def column_text(c):
        if isinstance(c, tuple):
            return " ".join(str(part) for part in c if part is not None).lower()
        return str(c).lower()
    text = column_text(col)
    return any(
        keyword in text
        for keyword in (
            "rainfall", "precip", "rain", "intensity",
            "outflow_time", "flow_time", "max_flow_time",
        )
    )


# ──────────────────────────────────────────────────────────────────────────────
# Test runner
# ──────────────────────────────────────────────────────────────────────────────

_PASS = []
_FAIL = []


def _test(name: str):
    """Decorator that wraps a test function, catches failures, and logs results."""
    def decorator(fn):
        def wrapper():
            try:
                fn()
                log.info("  PASS  %s", name)
                _PASS.append(name)
            except Exception as exc:
                log.error("  FAIL  %s\n         %s: %s", name, type(exc).__name__, exc)
                _FAIL.append(name)
        return wrapper
    return decorator


# ──────────────────────────────────────────────────────────────────────────────
# CONFIG TESTS
# ──────────────────────────────────────────────────────────────────────────────

@_test("config — PARAMETER_GRID Cartesian product matches legacy itertools.product")
def test_config_grid_cartesian_product():
    # Legacy factor arrays — verbatim from the calibration setup cell
    imp_factors        = np.arange(1,    1.2,  0.05)
    storage_factors    = np.arange(3,    8,    1)
    width_factors      = np.arange(0.4,  1.11, 0.3)
    n_factors          = np.arange(1,    1.1,  100)   # produces [1.0]
    pct_zero_factors   = np.arange(1,    1.41, 0.7)   # produces [1.0]
    cn_factors         = np.arange(1,    1.1,  100)   # produces [1.0]
    pct_routed_factors = np.arange(30,   46,   5)
    evap_factors       = np.arange(3,    4.1,  1)

    legacy_combo = np.array(list(itertools.product(
        imp_factors, storage_factors, width_factors, n_factors,
        pct_zero_factors, cn_factors, pct_routed_factors, evap_factors,
    )))

    # Build the same combination from the new config
    new_values = [config.PARAMETER_GRID[k].values() for k in config.PARAMETER_GRID]
    new_combo  = np.array(list(itertools.product(*new_values)))

    assert legacy_combo.shape == new_combo.shape, (
        f"Shape mismatch: legacy {legacy_combo.shape} vs new {new_combo.shape}"
    )
    np.testing.assert_allclose(legacy_combo, new_combo, rtol=1e-6)
    log.info("         Grid shape: %s  (%d total combinations)",
             new_combo.shape, len(new_combo))


@_test("config — WIDTH_INITIAL_VALUES matches legacy hardcoded array")
def test_config_width_initial():
    legacy = np.array([
        2470, 2800, 2800, 1890, 1160, 2840, 2660, 2000,
        1140, 1400,  540, 2200, 4160, 1070, 1310, 1310,
        1200, 1440, 1400, 1400, 1520, 1850, 1740, 2504,
         980, 1900, 2140,  950,
    ], dtype=float)
    np.testing.assert_array_equal(config.WIDTH_INITIAL_VALUES, legacy)
    assert len(config.WIDTH_INITIAL_VALUES) == 28


@_test("config — CN_INITIAL_VALUES matches legacy hardcoded array")
def test_config_cn_initial():
    legacy = np.array([
        39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39,
        39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39,
    ], dtype=float)
    np.testing.assert_array_equal(config.CN_INITIAL_VALUES, legacy)
    assert len(config.CN_INITIAL_VALUES) == 28


@_test("config — PCT_ZERO_INITIAL_VALUES matches legacy hardcoded array")
def test_config_pctzero_initial():
    legacy = np.array([
        30, 50, 50, 50, 50, 50, 50, 50, 30, 50, 50, 30,
        30, 50, 30, 50, 50, 50, 50, 50, 50, 60, 50, 50,
        50, 60, 60, 30,
    ], dtype=float)
    np.testing.assert_array_equal(config.PCT_ZERO_INITIAL_VALUES, legacy)
    assert len(config.PCT_ZERO_INITIAL_VALUES) == 28


@_test("config — PARETO_OBJECTIVES matches legacy notebook-04 selection and is subset of ALL_OBJECTIVES")
def test_config_pareto_objectives():
    legacy_selection = ["peak_1-kge", "volume_1-kge", "peak_rmsd", "peak_bias"]
    assert config.PARETO_OBJECTIVES == legacy_selection, (
        f"Got: {config.PARETO_OBJECTIVES}"
    )
    for obj in config.PARETO_OBJECTIVES:
        assert obj in config.ALL_OBJECTIVES, f"'{obj}' missing from ALL_OBJECTIVES"
    assert config.SAVE_FULL_PARETO_ENSEMBLE is True
    log.info("         Active objectives: %s", config.PARETO_OBJECTIVES)


@_test("config — parameter key order matches legacy itertools.product argument order")
def test_config_parameter_key_order():
    # The legacy product call order was: imp, storage, width, n, pct_zero, cn, pct_routed, evap
    expected_order = [
        "imperviousness", "storage", "width", "n",
        "pct_zero", "cn", "pct_routed", "evaporation",
    ]
    actual_order = list(config.PARAMETER_GRID.keys())
    assert actual_order == expected_order, (
        f"Expected: {expected_order}\nGot:      {actual_order}"
    )


# ──────────────────────────────────────────────────────────────────────────────
# LOADER TESTS
# ──────────────────────────────────────────────────────────────────────────────

@_test("loaders — load_observed_runoff produces frame identical to legacy function")
def test_load_observed_runoff():
    legacy_df = _legacy_load_runoff_obs_to_df(LEGACY_OBS)
    new_df    = load_observed_runoff(NEW_OBS_CSV)

    pd.testing.assert_frame_equal(new_df, legacy_df, check_names=True)
    assert isinstance(new_df.index, pd.DatetimeIndex)
    assert list(new_df.columns) == ["OBS runoff [CMS]"]
    log.info(
        "         Rows: %d  |  index: %s → %s",
        len(new_df), new_df.index[0], new_df.index[-1],
    )


@_test("loaders — load_swmm_output produces frame identical to legacy function")
def test_load_swmm_output():
    legacy_df = _legacy_load_sim_to_df(SWMM_OUT)
    new_df    = load_swmm_output(SWMM_OUT)

    pd.testing.assert_frame_equal(new_df, legacy_df, check_names=True)
    assert isinstance(new_df.index, pd.DatetimeIndex)
    assert list(new_df.columns) == ["SWMM outflow [CMS]", "rainfall [mm/h]"]
    log.info(
        "         Rows: %d  |  columns: %s", len(new_df), new_df.columns.tolist()
    )


@_test("loaders — load_basins_rain returns a dict with expected event keys")
def test_load_basins_rain():
    data = load_basins_rain(RAIN_PKL)
    assert isinstance(data, dict)
    assert len(data) > 0
    # Verify that the test event is present
    assert "2012_01_13" in data, "Expected event '2012_01_13' not found in rain dict"
    log.info("         Events available: %d  (includes '2012_01_13')", len(data))


@_test("loaders — load_observed_runoff raises FileNotFoundError for missing file")
def test_load_observed_runoff_missing():
    try:
        load_observed_runoff(r"D:\nonexistent\path\file.csv")
        assert False, "Expected FileNotFoundError was not raised"
    except FileNotFoundError:
        pass  # expected


@_test("loaders — load_swmm_output raises FileNotFoundError for missing file")
def test_load_swmm_output_missing():
    try:
        load_swmm_output(r"D:\nonexistent\path\file.out")
        assert False, "Expected FileNotFoundError was not raised"
    except FileNotFoundError:
        pass  # expected


# ──────────────────────────────────────────────────────────────────────────────
# PREPROCESSOR TESTS
# ──────────────────────────────────────────────────────────────────────────────

@_test("preprocessors — resample_5min produces frame identical to legacy function")
def test_resample_5min():
    raw_legacy = _legacy_load_runoff_obs_to_df(LEGACY_OBS)
    legacy_df  = _legacy_data_5min_interpolate(raw_legacy)

    raw_new = load_observed_runoff(NEW_OBS_CSV)
    new_df  = resample_5min(raw_new)

    pd.testing.assert_frame_equal(new_df, legacy_df)
    log.info(
        "         Rows before resampling: %d  →  after: %d",
        len(raw_new), len(new_df),
    )


@_test("preprocessors — filter_rainfall_columns removes 'rainfall [mm/h]' and nothing else")
def test_filter_rainfall_columns():
    sim_df   = load_swmm_output(SWMM_OUT)
    filtered = filter_rainfall_columns(sim_df)

    assert "rainfall [mm/h]"    not in filtered.columns, "Rainfall column was not removed"
    assert "SWMM outflow [CMS]"     in filtered.columns, "Outflow column was wrongly removed"
    assert filtered.shape[1] == sim_df.shape[1] - 1

    # Verify behaviour matches the legacy keyword check
    for col in filtered.columns:
        assert not _legacy_is_rainfall_column(col), (
            f"Column '{col}' should have been removed by the legacy check too"
        )
    log.info(
        "         Columns: %s  →  %s",
        sim_df.columns.tolist(), filtered.columns.tolist(),
    )


@_test("preprocessors — align_obs_sim output is identical to legacy inline logic")
def test_align_obs_sim():
    obs_new = resample_5min(load_observed_runoff(NEW_OBS_CSV))
    sim_new = load_swmm_output(SWMM_OUT)

    # Legacy inline logic (verbatim from notebook 04 run_calibration loop)
    obs_leg  = _legacy_data_5min_interpolate(_legacy_load_runoff_obs_to_df(LEGACY_OBS))
    sim_leg  = _legacy_load_sim_to_df(SWMM_OUT)
    legacy_storm = pd.concat([obs_leg, sim_leg], axis=1)
    legacy_storm[legacy_storm < 0] = 0
    legacy_storm = legacy_storm.fillna(0)

    new_storm = align_obs_sim(obs_new, sim_new)

    pd.testing.assert_frame_equal(new_storm, legacy_storm)
    assert (new_storm < 0).sum().sum() == 0, "Negative values found after alignment"
    assert new_storm.isna().sum().sum() == 0, "NaN values found after alignment"
    log.info(
        "         storm_df shape: %s  |  no negatives, no NaNs", new_storm.shape
    )


@_test("preprocessors — resample_5min raises TypeError for non-DatetimeIndex input")
def test_resample_5min_type_error():
    bad_df = pd.DataFrame({"x": [1, 2, 3]})  # integer index
    try:
        resample_5min(bad_df)
        assert False, "Expected TypeError was not raised"
    except TypeError:
        pass  # expected


# ──────────────────────────────────────────────────────────────────────────────
# MAIN
# ──────────────────────────────────────────────────────────────────────────────

def _run_all():
    all_tests = [
        test_config_grid_cartesian_product,
        test_config_width_initial,
        test_config_cn_initial,
        test_config_pctzero_initial,
        test_config_pareto_objectives,
        test_config_parameter_key_order,
        test_load_observed_runoff,
        test_load_swmm_output,
        test_load_basins_rain,
        test_load_observed_runoff_missing,
        test_load_swmm_output_missing,
        test_resample_5min,
        test_filter_rainfall_columns,
        test_align_obs_sim,
        test_resample_5min_type_error,
    ]

    log.info("=" * 65)
    log.info("Phase 1 + 2 Validation — UrbanRunoffModeling_Refactored")
    log.info("=" * 65)

    log.info("--- CONFIG TESTS ---")
    for t in all_tests[:6]:
        t()

    log.info("--- LOADER TESTS ---")
    for t in all_tests[6:11]:
        t()

    log.info("--- PREPROCESSOR TESTS ---")
    for t in all_tests[11:]:
        t()

    log.info("=" * 65)
    total  = len(_PASS) + len(_FAIL)
    passed = len(_PASS)
    log.info("Results: %d / %d passed", passed, total)
    if _FAIL:
        log.error("Failed tests:")
        for name in _FAIL:
            log.error("  - %s", name)
        log.info("=" * 65)
        sys.exit(1)
    else:
        log.info("ALL %d TESTS PASSED", total)
        log.info("=" * 65)
        sys.exit(0)


if __name__ == "__main__":
    _run_all()
