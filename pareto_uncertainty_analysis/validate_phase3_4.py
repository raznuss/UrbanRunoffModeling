"""
Phase 3 & 4 validation script — SWMM execution layer and calibration objectives.

Test 1 — INP Mutation Identity
    Loads the same original .inp twice (read-only legacy directory).
    Applies calibration factors via the legacy 8-function approach (verbatim
    copies of notebook 01 cell 6) and via the refactored apply_all_factors().
    Asserts that every subcatchment/subarea/infiltration attribute matches
    exactly between the two approaches.

Test 2 — Objective Math Identity
    Generates synthetic observed/simulated arrays with a fixed NumPy seed.
    Computes all six objectives via the legacy calculate_objectives() (verbatim
    copy of notebook 01 cell 10) and via the refactored compute_objectives().
    Asserts float equality to 12 decimal places.

Both tests are read-only with respect to the legacy directories.
No files are written to D:\\Development\\RESEARCH\\Raanana or
D:\\MY_CODES\\UrbanRunoffModeling.
"""

import sys
import os
import traceback

import numpy as np
import hydroeval as he
from swmm_api import read_inp_file
from swmm_api.input_file import section_labels as sections

# ──────────────────────────────────────────────────────────────────────────────
# Paths — read-only legacy data
# ──────────────────────────────────────────────────────────────────────────────

LEGACY_INP = (
    r"D:\Development\RESEARCH\Raanana\SWMM\from_radar"
    r"\Cross_validation\2012_01_13\raanana_28subcatchments.inp"
)

# ──────────────────────────────────────────────────────────────────────────────
# Refactored module under test
# ──────────────────────────────────────────────────────────────────────────────

sys.path.insert(0, r"D:\MY_CODES\UrbanRunoffModeling_Refactored")
from urban_runoff.swmm.inp_editor import apply_all_factors
from urban_runoff.calibration.objectives import compute_objectives


# ──────────────────────────────────────────────────────────────────────────────
# Verbatim copies of legacy functions (notebook 01, cell 6)
# ──────────────────────────────────────────────────────────────────────────────

def _legacy_update_imperviousness(subcatchment_dict, factor):
    for _, subcatchment in subcatchment_dict.items():
        subcatchment.imperviousness = subcatchment.imperviousness * factor
    return subcatchment_dict


def _legacy_update_storage(subareas_dict, factor):
    for _, subarea in subareas_dict.items():
        subarea.storage_imperv = subarea.storage_imperv * factor
        subarea.storage_perv   = subarea.storage_perv   * factor
    return subareas_dict


def _legacy_update_n(subareas_dict, factor):
    for _, subarea in subareas_dict.items():
        subarea.n_imperv = subarea.n_imperv * factor
        subarea.n_perv   = subarea.n_perv   * factor
    return subareas_dict


def _legacy_update_width(subcatchment_dict, factor, width_initial_values):
    for index, (_, subcatchment) in enumerate(subcatchment_dict.items()):
        subcatchment.width = width_initial_values[index] * factor
    return subcatchment_dict


def _legacy_update_curve_num(infiltration_dict, factor, curve_no_initial_values):
    for index, (_, inf) in enumerate(infiltration_dict.items()):
        inf.curve_no = curve_no_initial_values[index] * factor
    return infiltration_dict


def _legacy_update_pct_zero(subareas_dict, pct_zero_factor, pct_zero_initial_values):
    for index, (_, subarea) in enumerate(subareas_dict.items()):
        subarea.pct_zero = pct_zero_initial_values[index] * pct_zero_factor
    return subareas_dict


def _legacy_update_pct_routed(subareas_dict, routed_to, percent):
    for _, subarea in subareas_dict.items():
        subarea.route_to   = routed_to
        subarea.pct_routed = percent
    return subareas_dict


# Verbatim copy of notebook 01, cell 10
def _legacy_calculate_objectives(observed, simulated):
    if not isinstance(observed, np.ndarray):
        observed = observed.to_numpy()
    if not isinstance(simulated, np.ndarray):
        simulated = simulated.to_numpy()
    observed_mean = np.mean(observed)
    ss_diff   = np.sum((observed - simulated) ** 2)
    ss_mean   = np.sum((observed - observed_mean) ** 2)
    nse       = 1 - (ss_diff / ss_mean)
    rmsd      = np.sqrt(np.mean((observed - simulated) ** 2))
    bias      = abs(np.mean(simulated - observed))
    mad       = np.mean(np.abs(simulated - observed))
    kge_tuple = he.evaluator(he.kge, simulated, observed)
    kge       = float(kge_tuple[0])
    ss_total  = np.sum((observed - observed_mean) ** 2)
    r_squared = 1 - (ss_diff / ss_total)
    return nse, rmsd, bias, mad, r_squared, kge


# ──────────────────────────────────────────────────────────────────────────────
# Initial value arrays (from config — kept inline for test self-containment)
# ──────────────────────────────────────────────────────────────────────────────

_WIDTH_INIT = np.array([
    2470, 2800, 2800, 1890, 1160, 2840, 2660, 2000,
    1140, 1400,  540, 2200, 4160, 1070, 1310, 1310,
    1200, 1440, 1400, 1400, 1520, 1850, 1740, 2504,
     980, 1900, 2140,  950,
], dtype=float)

_CN_INIT = np.full(28, 39, dtype=float)

_PCT_ZERO_INIT = np.array([
    30, 50, 50, 50, 50, 50, 50, 50, 30, 50, 50, 30,
    30, 50, 30, 50, 50, 50, 50, 50, 50, 60, 50, 50,
    50, 60, 60, 30,
], dtype=float)

# Factor combination used for both tests
_TEST_FACTORS = {
    "imperviousness": 1.05,
    "storage":        5.0,
    "width":          0.7,
    "n":              1.0,
    "pct_zero":       1.0,
    "cn":             1.0,
    "pct_routed":     35.0,
    "evaporation":    3.0,
}


# ──────────────────────────────────────────────────────────────────────────────
# Test 1 — INP Mutation Identity
# ──────────────────────────────────────────────────────────────────────────────

def test_inp_mutation():
    """Assert that apply_all_factors() produces identical attributes to the
    legacy 8-function sequence for all 28 subcatchments."""

    print("\n" + "=" * 60)
    print("TEST 1 — INP Mutation Identity")
    print("=" * 60)

    passed = 0
    failed = 0

    # --- Legacy approach: apply via the 8 separate functions ---
    inp_legacy = read_inp_file(LEGACY_INP)
    sc_d   = dict(inp_legacy[sections.SUBCATCHMENTS])
    sa_d   = dict(inp_legacy[sections.SUBAREAS])
    inf_d  = dict(inp_legacy[sections.INFILTRATION])

    _legacy_update_imperviousness(sc_d,  _TEST_FACTORS["imperviousness"])
    _legacy_update_storage(sa_d,         _TEST_FACTORS["storage"])
    _legacy_update_width(sc_d,           _TEST_FACTORS["width"], _WIDTH_INIT)
    _legacy_update_n(sa_d,               _TEST_FACTORS["n"])
    _legacy_update_pct_zero(sa_d,        _TEST_FACTORS["pct_zero"], _PCT_ZERO_INIT)
    _legacy_update_curve_num(inf_d,      _TEST_FACTORS["cn"], _CN_INIT)
    _legacy_update_pct_routed(sa_d,      "PERVIOUS", _TEST_FACTORS["pct_routed"])
    inp_legacy["EVAPORATION"]["CONSTANT"] = _TEST_FACTORS["evaporation"]

    # --- Refactored approach: apply via apply_all_factors() ---
    inp_new = read_inp_file(LEGACY_INP)
    apply_all_factors(inp_new, _TEST_FACTORS)

    # --- Compare subcatchment attributes ---
    sc_d_new  = dict(inp_new[sections.SUBCATCHMENTS])
    sa_d_new  = dict(inp_new[sections.SUBAREAS])
    inf_d_new = dict(inp_new[sections.INFILTRATION])

    sc_keys = list(sc_d.keys())
    assert len(sc_keys) == 28, f"Expected 28 subcatchments, got {len(sc_keys)}"

    for key in sc_keys:
        sc_leg = sc_d[key];       sc_ref = sc_d_new[key]
        sa_leg = sa_d[key];       sa_ref = sa_d_new[key]
        inf_leg = inf_d[key];     inf_ref = inf_d_new[key]

        checks = [
            ("imperviousness", sc_leg.imperviousness,  sc_ref.imperviousness),
            ("width",          sc_leg.width,            sc_ref.width),
            ("storage_imperv", sa_leg.storage_imperv,  sa_ref.storage_imperv),
            ("storage_perv",   sa_leg.storage_perv,    sa_ref.storage_perv),
            ("n_imperv",       sa_leg.n_imperv,        sa_ref.n_imperv),
            ("n_perv",         sa_leg.n_perv,          sa_ref.n_perv),
            ("pct_zero",       sa_leg.pct_zero,        sa_ref.pct_zero),
            ("pct_routed",     sa_leg.pct_routed,      sa_ref.pct_routed),
            ("route_to",       sa_leg.route_to,        sa_ref.route_to),
            ("curve_no",       inf_leg.curve_no,       inf_ref.curve_no),
        ]

        for attr_name, leg_val, ref_val in checks:
            if isinstance(leg_val, str):
                ok = (leg_val == ref_val)
            else:
                ok = np.isclose(float(leg_val), float(ref_val), rtol=1e-12)
            if ok:
                passed += 1
            else:
                failed += 1
                print(f"  FAIL [{key}] {attr_name}: legacy={leg_val!r}  refactored={ref_val!r}")

    # --- Compare evaporation constant ---
    leg_evap = inp_legacy["EVAPORATION"]["CONSTANT"]
    ref_evap = inp_new["EVAPORATION"]["CONSTANT"]
    if np.isclose(float(leg_evap), float(ref_evap), rtol=1e-12):
        passed += 1
        print(f"  PASS  evaporation: legacy={leg_evap}  refactored={ref_evap}")
    else:
        failed += 1
        print(f"  FAIL  evaporation: legacy={leg_evap}  refactored={ref_evap}")

    total = passed + failed
    print(f"\n  {passed}/{total} attribute checks passed")
    if failed == 0:
        print("  RESULT: PASS — all 28 × 10 attributes + evaporation match exactly")
    else:
        print(f"  RESULT: FAIL — {failed} attribute(s) did not match")
    return failed == 0


# ──────────────────────────────────────────────────────────────────────────────
# Test 2 — Objective Math Identity
# ──────────────────────────────────────────────────────────────────────────────

def test_objective_math():
    """Assert that compute_objectives() produces the same six values as the
    legacy calculate_objectives() on identical synthetic input arrays."""

    print("\n" + "=" * 60)
    print("TEST 2 — Objective Math Identity")
    print("=" * 60)

    rng = np.random.default_rng(42)
    obs = rng.uniform(0.01, 5.0, size=100)
    sim = obs + rng.normal(0, 0.3, size=100)
    sim = np.clip(sim, 0.0, None)

    # Legacy outputs (verbatim)
    leg_nse, leg_rmsd, leg_bias, leg_mad, leg_r2, leg_kge = \
        _legacy_calculate_objectives(obs, sim)

    # Refactored outputs
    ref = compute_objectives(obs, sim)

    pairs = [
        ("nse",       leg_nse,  ref["nse"]),
        ("rmsd",      leg_rmsd, ref["rmsd"]),
        ("bias",      leg_bias, ref["bias"]),
        ("mad",       leg_mad,  ref["mad"]),
        ("r_squared", leg_r2,   ref["r_squared"]),
        ("kge",       leg_kge,  ref["kge"]),
    ]

    passed = 0
    failed = 0
    for name, leg_val, ref_val in pairs:
        ok = np.isclose(leg_val, ref_val, rtol=1e-12, atol=1e-12)
        status = "PASS" if ok else "FAIL"
        if ok:
            passed += 1
        else:
            failed += 1
        print(f"  {status}  {name:12s}: legacy={leg_val:.12f}  refactored={ref_val:.12f}")

    # Additional check: minimize-form values
    print()
    print(f"  1-nse = {1 - ref['nse']:.6f},  1-kge = {1 - ref['kge']:.6f}  (Pareto minimize-form)")

    print(f"\n  {passed}/6 objective checks passed")
    if failed == 0:
        print("  RESULT: PASS — all six objectives match to 12 decimal places")
    else:
        print(f"  RESULT: FAIL — {failed} objective(s) did not match")
    return failed == 0


# ──────────────────────────────────────────────────────────────────────────────
# Runner
# ──────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    all_passed = True
    try:
        all_passed &= test_inp_mutation()
    except Exception:
        all_passed = False
        print("  ERROR in test_inp_mutation:")
        traceback.print_exc()

    try:
        all_passed &= test_objective_math()
    except Exception:
        all_passed = False
        print("  ERROR in test_objective_math:")
        traceback.print_exc()

    print()
    print("=" * 60)
    if all_passed:
        print("OVERALL: ALL TESTS PASSED (Phase 3 & 4 validation complete)")
    else:
        print("OVERALL: SOME TESTS FAILED")
    print("=" * 60)
    sys.exit(0 if all_passed else 1)
