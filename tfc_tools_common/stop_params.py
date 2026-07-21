"""tfc_tools_common.stop_params

Single source of truth for the spatial parameters controlling stop extraction
and trip-stop sequencing in the TfC SDI pipeline.

These parameters are exposed in two tools:
  - rl2sdi (RouteLab → SDI PostGIS, RouteLab → GeoPackage)
  - sdi_tools.refresh_sdi_derived (rebuild derived layers after GIS edits)

Keeping defaults and parameter keys identical across both tools ensures that
the "edit in GIS, refresh" workflow reproduces the same result as a fresh
export from RouteLab.

All distance values are in meters (EPSG:3857 / geography in PostGIS).
"""

from __future__ import annotations

import contextlib
from dataclasses import dataclass, asdict, field


# ---------------------------------------------------------------------------
# Defaults  —  historical hard-coded values preserved exactly.
# ---------------------------------------------------------------------------

# Stop clustering (DBSCAN on raw operator-tagged stops)
# v2 (direction-aware): clustering now runs PER direction, so eps/minpoints are
# tuned for a single carriageway rather than for both merged. Validated on the
# Sousse export at eps=40 m, minpoints=2.
DEFAULT_DBSCAN_EPS_M = 40.0
DEFAULT_DBSCAN_MINPOINTS = 2
DEFAULT_CELL_M = 120.0              # DEPRECATED (v1 grid dedup); unused by the
                                   # direction-aware SQL. Retained so the
                                   # geopandas mirror still imports until ported.

# Stop snapping / sequencing
DEFAULT_SNAP_MAX_M = 40.0           # max cluster-to-trip snap distance
DEFAULT_TERMINAL_M = 75.0           # distance from trip start/end for Terminal
DEFAULT_PAIR_MAX_M = 40.0           # max dir-0<->dir-1 separation to link as a pair
DEFAULT_RESCUE_GAP_M = 300.0        # coverage rescue: keep an isolated on-route raw
                                    # stop (DBSCAN noise) as a stop when it is farther
                                    # than this from any confident cluster stop, so
                                    # sparsely-surveyed corridors are not left without
                                    # any stop. 0 disables rescue.
DEFAULT_STOP_TRIP_BUFFER_M = 1.0    # stop-to-trip join radius for sequence
DEFAULT_MIN_STOP_SPACING_M = 100.0  # drop consecutive stops closer than this

# OD stats (travel-time computation)
DEFAULT_TRACKPOINT_BUFFER_M = 30.0  # trackpoint → stop snap radius
# When True, travel times are computed per (o_id, d_id, interval, vehicle_type).
# When False, durations are pooled across vehicle types for the same OD+interval,
# which gives a larger sample when trackpoints are sparse per vehicle type.
# (Output rows are still emitted per vehicle type so downstream GTFS joins work.)
DEFAULT_DISTINGUISH_SPEEDS_BY_VEHICLE = True


# ---------------------------------------------------------------------------
# Parameter keys (stable strings, used as QGIS parameter names)
# ---------------------------------------------------------------------------

class StopParamKeys:
    """QGIS parameter names. Kept identical across algorithms."""
    DBSCAN_EPS_M        = "dbscan_eps_m"
    DBSCAN_MINPOINTS    = "dbscan_minpoints"
    CELL_M              = "cell_m"
    SNAP_MAX_M          = "snap_max_m"
    TERMINAL_M          = "terminal_m"
    PAIR_MAX_M          = "pair_max_m"
    RESCUE_GAP_M        = "rescue_gap_m"
    STOP_TRIP_BUFFER_M  = "stop_trip_buffer_m"
    MIN_STOP_SPACING_M  = "min_stop_spacing_m"
    TRACKPOINT_BUFFER_M = "trackpoint_buffer_m"
    DISTINGUISH_SPEEDS_BY_VEHICLE = "distinguish_speeds_by_vehicle"


# ---------------------------------------------------------------------------
# Dataclass carrying the resolved values
# ---------------------------------------------------------------------------

@dataclass
class StopParams:
    """Resolved stop-extraction parameters for a single pipeline run."""

    dbscan_eps_m:        float = DEFAULT_DBSCAN_EPS_M
    dbscan_minpoints:    int   = DEFAULT_DBSCAN_MINPOINTS
    cell_m:              float = DEFAULT_CELL_M
    snap_max_m:          float = DEFAULT_SNAP_MAX_M
    terminal_m:          float = DEFAULT_TERMINAL_M
    pair_max_m:          float = DEFAULT_PAIR_MAX_M
    rescue_gap_m:        float = DEFAULT_RESCUE_GAP_M
    stop_trip_buffer_m:  float = DEFAULT_STOP_TRIP_BUFFER_M
    min_stop_spacing_m:  float = DEFAULT_MIN_STOP_SPACING_M
    trackpoint_buffer_m: float = DEFAULT_TRACKPOINT_BUFFER_M
    distinguish_speeds_by_vehicle: bool = DEFAULT_DISTINGUISH_SPEEDS_BY_VEHICLE

    # --------------------------------------------------------------
    # Convenience
    # --------------------------------------------------------------
    def as_dict(self) -> dict:
        return asdict(self)

    def is_default(self) -> bool:
        """True if every parameter matches the historical default.

        Used by refresh_sdi_derived's Postgres path to decide between a cheap
        `REFRESH MATERIALIZED VIEW` (if all defaults) and a full DROP+CREATE
        (if any parameter has been customized).
        """
        return (
            self.dbscan_eps_m        == DEFAULT_DBSCAN_EPS_M
            and self.dbscan_minpoints    == DEFAULT_DBSCAN_MINPOINTS
            and self.cell_m              == DEFAULT_CELL_M
            and self.snap_max_m          == DEFAULT_SNAP_MAX_M
            and self.terminal_m          == DEFAULT_TERMINAL_M
            and self.pair_max_m          == DEFAULT_PAIR_MAX_M
            and self.rescue_gap_m        == DEFAULT_RESCUE_GAP_M
            and self.stop_trip_buffer_m  == DEFAULT_STOP_TRIP_BUFFER_M
            and self.min_stop_spacing_m  == DEFAULT_MIN_STOP_SPACING_M
            and self.trackpoint_buffer_m == DEFAULT_TRACKPOINT_BUFFER_M
            and self.distinguish_speeds_by_vehicle == DEFAULT_DISTINGUISH_SPEEDS_BY_VEHICLE
        )

    def describe(self) -> str:
        """One-line summary for feedback.pushInfo()."""
        vehicle_mode = "per vehicle type" if self.distinguish_speeds_by_vehicle else "pooled across vehicle types"
        return (
            f"DBSCAN(eps={self.dbscan_eps_m}m, minpoints={self.dbscan_minpoints}), "
            f"snap<={self.snap_max_m}m, pair<={self.pair_max_m}m, "
            f"rescue-gap>{self.rescue_gap_m}m, "
            f"terminal<={self.terminal_m}m, stop-trip<={self.stop_trip_buffer_m}m, "
            f"min-spacing>={self.min_stop_spacing_m}m, "
            f"trackpoint<={self.trackpoint_buffer_m}m, "
            f"OD speeds {vehicle_mode}"
        )


# ---------------------------------------------------------------------------
# QGIS parameter-registration helper (imported only when QGIS is available)
# ---------------------------------------------------------------------------

def add_stop_params_to_algorithm(algorithm) -> None:
    """Register all stop-extraction parameters on a QgsProcessingAlgorithm.

    All parameters are flagged `Advanced` so they're hidden by default
    and won't clutter the main dialog.
    """
    from qgis.core import (
        QgsProcessingParameterNumber,
        QgsProcessingParameterBoolean,
        QgsProcessingParameterDefinition,
    )

    def _mark_advanced(p):
        with contextlib.suppress(Exception):
            p.setFlags(p.flags() | QgsProcessingParameterDefinition.Flag.FlagAdvanced)
        return p

    def _add_number(key, label, default, is_int=False, min_val=0.0):
        ptype = (QgsProcessingParameterNumber.Type.Integer if is_int
                 else QgsProcessingParameterNumber.Type.Double)
        p = QgsProcessingParameterNumber(
            key, label, type=ptype, defaultValue=default, minValue=min_val,
        )
        algorithm.addParameter(_mark_advanced(p))

    def _add_bool(key, label, default):
        p = QgsProcessingParameterBoolean(key, label, defaultValue=default)
        algorithm.addParameter(_mark_advanced(p))

    # --- Stop clustering (DBSCAN) ---
    _add_number(
        StopParamKeys.DBSCAN_EPS_M,
        "DBSCAN eps — stop clustering neighborhood radius (m)",
        DEFAULT_DBSCAN_EPS_M,
    )
    _add_number(
        StopParamKeys.DBSCAN_MINPOINTS,
        "DBSCAN minpoints — min neighbors to form a cluster core",
        DEFAULT_DBSCAN_MINPOINTS,
        is_int=True, min_val=1,
    )
    _add_number(
        StopParamKeys.CELL_M,
        "Grid cell size for cluster deduplication (m)",
        DEFAULT_CELL_M,
    )

    # --- Stop snapping / sequencing ---
    _add_number(
        StopParamKeys.SNAP_MAX_M,
        "Max cluster-to-trip snap distance (m)",
        DEFAULT_SNAP_MAX_M,
    )
    _add_number(
        StopParamKeys.TERMINAL_M,
        "Terminal classification radius from trip endpoints (m)",
        DEFAULT_TERMINAL_M,
    )
    _add_number(
        StopParamKeys.PAIR_MAX_M,
        "Max separation to pair opposite-direction stops (m)",
        DEFAULT_PAIR_MAX_M,
    )
    _add_number(
        StopParamKeys.RESCUE_GAP_M,
        "Coverage rescue gap — keep an isolated on-route stop when no other stop "
        "is within this distance (m); 0 disables",
        DEFAULT_RESCUE_GAP_M,
    )
    _add_number(
        StopParamKeys.STOP_TRIP_BUFFER_M,
        "Stop-to-trip buffer for sequencing (m)",
        DEFAULT_STOP_TRIP_BUFFER_M,
    )
    _add_number(
        StopParamKeys.MIN_STOP_SPACING_M,
        "Min spacing between consecutive stops on a trip (m)",
        DEFAULT_MIN_STOP_SPACING_M,
    )

    # --- OD stats (travel time) ---
    _add_number(
        StopParamKeys.TRACKPOINT_BUFFER_M,
        "Trackpoint-to-stop buffer for OD travel-time (m)",
        DEFAULT_TRACKPOINT_BUFFER_M,
    )
    _add_bool(
        StopParamKeys.DISTINGUISH_SPEEDS_BY_VEHICLE,
        "Distinguish speeds by vehicle type "
        "(uncheck to pool durations across vehicle types when per-vehicle "
        "trackpoint samples are too small)",
        DEFAULT_DISTINGUISH_SPEEDS_BY_VEHICLE,
    )


def read_stop_params_from_algorithm(algorithm, parameters, context) -> StopParams:
    """Pull the stop-extraction parameter values out of a QGIS parameter dict."""
    def _num(key):
        return algorithm.parameterAsDouble(parameters, key, context)

    def _int(key):
        return algorithm.parameterAsInt(parameters, key, context)

    def _bool(key):
        return algorithm.parameterAsBool(parameters, key, context)

    return StopParams(
        dbscan_eps_m        = _num(StopParamKeys.DBSCAN_EPS_M)        or DEFAULT_DBSCAN_EPS_M,
        dbscan_minpoints    = _int(StopParamKeys.DBSCAN_MINPOINTS)    or DEFAULT_DBSCAN_MINPOINTS,
        cell_m              = _num(StopParamKeys.CELL_M)              or DEFAULT_CELL_M,
        snap_max_m          = _num(StopParamKeys.SNAP_MAX_M)          or DEFAULT_SNAP_MAX_M,
        terminal_m          = _num(StopParamKeys.TERMINAL_M)          or DEFAULT_TERMINAL_M,
        pair_max_m          = _num(StopParamKeys.PAIR_MAX_M)          or DEFAULT_PAIR_MAX_M,
        rescue_gap_m        = _num(StopParamKeys.RESCUE_GAP_M)        or DEFAULT_RESCUE_GAP_M,
        stop_trip_buffer_m  = _num(StopParamKeys.STOP_TRIP_BUFFER_M)  or DEFAULT_STOP_TRIP_BUFFER_M,
        min_stop_spacing_m  = _num(StopParamKeys.MIN_STOP_SPACING_M)  or DEFAULT_MIN_STOP_SPACING_M,
        trackpoint_buffer_m = _num(StopParamKeys.TRACKPOINT_BUFFER_M) or DEFAULT_TRACKPOINT_BUFFER_M,
        # Booleans: QGIS always returns True/False, so no default-fallback needed.
        distinguish_speeds_by_vehicle = _bool(StopParamKeys.DISTINGUISH_SPEEDS_BY_VEHICLE),
    )
