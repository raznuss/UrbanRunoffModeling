"""
SWMM .inp file parameter editor for the UrbanRunoffModeling calibration pipeline.

Implements two public functions:

  apply_all_factors(inp, factor_dict)
      Replaces the eight legacy update_X_with_factor() calls in run_calibration().
      Applies all calibration factors to a loaded inp object in a single call.

  update_timeseries(inp, event_date, basins_rain_df)
      Replaces the legacy update_TimeSeriesData() called in the preprocessing step
      (cell 17 of notebook 01).  event_date is passed explicitly to fix the legacy
      closure bug where the outer-loop variable 'date' was used inside the function.

Parameter application modes (from config.PARAMETER_MODES):
  Type A — multiplicative_from_array : new_val = INITIAL_ARRAY[i] * factor
  Type B — multiplicative_from_inp   : new_val = current_inp_value * factor
  Type C — absolute                  : new_val = factor  (direct assignment)

swmm_api 0.3.2 field names confirmed against legacy notebook cell 6:
  SubCatchment : .imperviousness, .width
  SubArea      : .storage_imperv, .storage_perv, .n_imperv, .n_perv,
                 .pct_zero, .route_to, .pct_routed
  Infiltration : .curve_no
  Evaporation  : inp['EVAPORATION']['CONSTANT']  (direct string-key access)
"""

import logging
from typing import Any, Dict

import pandas as pd
from swmm_api.input_file import section_labels as sections

from urban_runoff.config import (
    WIDTH_INITIAL_VALUES,
    CN_INITIAL_VALUES,
    PCT_ZERO_INITIAL_VALUES,
)

logger = logging.getLogger(__name__)


def apply_all_factors(inp: Any, factor_dict: Dict[str, float]) -> None:
    """Apply a complete set of calibration factors to a loaded SWMM inp object.

    Replicates the combined effect of the eight legacy update_X_with_factor()
    functions called inside run_calibration():

      update_imperviousness_with_factor  → Type B on subcatchment.imperviousness
      update_width_with_factor           → Type A on subcatchment.width
      update_storage_with_factor         → Type B on subarea.storage_imperv/_perv
      update_n_with_factor               → Type B on subarea.n_imperv/n_perv
      update_pct_zero_with_factor        → Type A on subarea.pct_zero
      update_pct_routed_with_factor      → Type C on subarea.route_to + pct_routed
      update_curve_num_with_factor       → Type A on infiltration.curve_no
      inp['EVAPORATION']['CONSTANT'] =   → Type C on evaporation constant

    The caller is responsible for loading a fresh inp from disk before each call
    (the reset mechanism that prevents factor accumulation across iterations).

    Parameters
    ----------
    inp : SwmmInput
        A loaded SWMM input object returned by ``read_inp_file()``.
        The object is mutated in-place; no return value.
    factor_dict : dict
        Mapping from canonical parameter name to factor value.  All eight keys
        should be present:
        'imperviousness', 'storage', 'width', 'n', 'pct_zero', 'cn',
        'pct_routed', 'evaporation'.
        Missing keys fall back to neutral values (1.0 for factors, defaults for
        absolute parameters) so that partial dicts can be used in tests.
    """
    imp_f  = factor_dict.get("imperviousness", 1.0)
    sto_f  = factor_dict.get("storage",        1.0)
    wid_f  = factor_dict.get("width",          1.0)
    n_f    = factor_dict.get("n",              1.0)
    pz_f   = factor_dict.get("pct_zero",       1.0)
    cn_f   = factor_dict.get("cn",             1.0)
    pr_val = factor_dict.get("pct_routed",    30.0)
    ev_val = factor_dict.get("evaporation",    3.0)

    # Shallow dicts with live references — mutations on the value objects are
    # immediately visible inside inp without any explicit write-back.
    subcatchment_d = dict(inp[sections.SUBCATCHMENTS])
    subareas_d     = dict(inp[sections.SUBAREAS])
    infiltration_d = dict(inp[sections.INFILTRATION])

    for i, sc_key in enumerate(subcatchment_d):
        sc  = subcatchment_d[sc_key]
        sa  = subareas_d[sc_key]
        inf = infiltration_d[sc_key]

        # Type B: imperviousness — multiply against current .inp value
        sc.imperviousness = sc.imperviousness * imp_f

        # Type A: width — multiply against initial array (W=2L baseline)
        sc.width = float(WIDTH_INITIAL_VALUES[i]) * wid_f

        # Type B: depression storage — multiply both surfaces
        sa.storage_imperv = sa.storage_imperv * sto_f
        sa.storage_perv   = sa.storage_perv * sto_f

        # Type B: Manning's n — multiply both surfaces
        sa.n_imperv = sa.n_imperv * n_f
        sa.n_perv   = sa.n_perv * n_f

        # Type A: pct_zero — multiply against initial array
        sa.pct_zero = float(PCT_ZERO_INITIAL_VALUES[i]) * pz_f

        # Type C: pct_routed — set routing destination and absolute percentage
        sa.route_to   = "PERVIOUS"
        sa.pct_routed = pr_val

        # Type A: curve number — multiply against initial array
        inf.curve_no = float(CN_INITIAL_VALUES[i]) * cn_f

    # Type C: evaporation constant [mm/day] — direct string-key access
    inp["EVAPORATION"]["CONSTANT"] = ev_val

    logger.debug(
        "apply_all_factors: imp=%.3f sto=%.3f wid=%.3f n=%.3f "
        "pz=%.3f cn=%.3f pr=%.1f ev=%.1f",
        imp_f, sto_f, wid_f, n_f, pz_f, cn_f, pr_val, ev_val,
    )


def update_timeseries(
    inp: Any,
    event_date: str,
    basins_rain_df: pd.DataFrame,
) -> None:
    """Update the TIMESERIES section of a SWMM inp object with radar rainfall.

    Replaces the legacy ``update_TimeSeriesData(timeseries_dict, basins_rain_df)``
    used in notebook cell 17.  The key fix over the legacy is that ``event_date``
    is passed explicitly rather than captured from an outer loop variable.

    The generated timeseries names follow the legacy convention:
        ``{event_date}_S{basin_name}``  e.g. ``2012_01_13_S3``

    Parameters
    ----------
    inp : SwmmInput
        Loaded SWMM input object.  The TIMESERIES section is replaced in-place
        with the raw-text representation generated from basins_rain_df.
    event_date : str
        Event date string in ``YYYY_MM_DD`` format (e.g. ``'2012_01_13'``).
        Used as the prefix in timeseries entry names.
    basins_rain_df : pd.DataFrame
        Radar rainfall DataFrame for this event.  Rows are basins; columns are
        datetime timestamps except the first column ('Basin_name').  Values are
        rainfall depths in mm.

    Notes
    -----
    swmm_api 0.3.2 accepts a raw text string assignment to a section:
        ``inp[sections.TIMESERIES] = text``
    This replicates the exact pattern used in the legacy preprocessing cell.
    """
    header = (
        ";;Name                 Date          Time         Value     \n"
        ";;------------ ----------------- ------------- -------------\n"
    )
    lines = [header]
    separator = ";;------------ ----------------- ------------- -------------\n"

    timestamps = basins_rain_df.columns[1:]  # all columns except 'Basin_name'

    for _, row in basins_rain_df.iterrows():
        basin_name = int(row["Basin_name"])
        values     = row.values[1:]
        dates_str  = [ts.strftime("%m/%d/%Y") for ts in timestamps]
        times_str  = [ts.strftime("%H:%M:%S") for ts in timestamps]

        for i in range(len(timestamps)):
            lines.append(
                f"  {event_date}_S{basin_name}      "
                f"{dates_str[i]}    {times_str[i]}   {values[i]:.4f}\n"
            )
        lines.append(separator)

    inp[sections.TIMESERIES] = "".join(lines)

    logger.debug(
        "update_timeseries: event=%s, %d basins, %d timesteps each",
        event_date,
        len(basins_rain_df),
        len(timestamps),
    )
