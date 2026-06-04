import geopandas as gpd
from shapely.ops import substring
from shapely.geometry import LineString, Point
from shapely import STRtree
import numpy as np
import pandas as pd

# import similaritymeasures

import functools
# from loguru import logger 

# loguru is an external logging library that QGIS’s Python environment does not include it by default.
# so we did the following:
import logging
logger = logging.getLogger("vehicle_passenger_flow") 
logging.basicConfig(level=logging.INFO) # This ensures that .info(), .debug(), .error() will actually output something.

# THIS IS A NEW HELPER FUNCTION TO TURN MULTILINESTRING TO LINESTRING IN FRECHET DISTANCE CALCULATION
from shapely.geometry import MultiLineString, GeometryCollection
from shapely.ops import linemerge
def to_single_linestring(geom):
    """
    Return a LineString for frechet distance:
    - LineString => itself
    - MultiLineString => longest component (or linemerge if it helps)
    - GeometryCollection => try merging/extracting the longest LineString
    - Anything else / empty => None
    """
    if geom is None or geom.is_empty:
        return None
    if isinstance(geom, LineString):
        return geom
    if isinstance(geom, MultiLineString):
        # pick the longest linestring component if any
        if len(geom.geoms) == 0:
            return None
        longest = max(geom.geoms, key=lambda g: g.length)
        # try to merge first if it helps, else use longest
        merged = linemerge(geom)
        if isinstance(merged, LineString):
            return merged
        return longest
    if isinstance(geom, GeometryCollection):
        # extract all lines, pick the longest
        lines = [g for g in geom.geoms if isinstance(g, (LineString, MultiLineString))]
        if not lines:
            return None
        merged = linemerge(lines)
        if isinstance(merged, LineString):
            return merged
        if isinstance(merged, MultiLineString) and len(merged.geoms) > 0:
            return max(merged.geoms, key=lambda g: g.length)
        # fallback to longest single linestring among originals
        flat = []
        for g in lines:
            if isinstance(g, LineString):
                flat.append(g)
            elif isinstance(g, MultiLineString):
                flat.extend(list(g.geoms))
        return max(flat, key=lambda g: g.length) if flat else None
    return None


import math
def discrete_frechet(P, Q):
    if len(P) == 0 or len(Q) == 0:  # WE ADDED THE IF CONDITION AS A DEFENSIVE CHECK TO AVOID AN ERROR
        return float("inf")  # or return None
    """
    Discrete Frechet distance between curves P and Q (each a list of [x, y] coordinates).
    """
    ca = [[-1 for _ in range(len(Q))] for _ in range(len(P))]

    def c(i, j):
        if ca[i][j] > -1:
            return ca[i][j]
        elif i == 0 and j == 0:
            ca[i][j] = math.dist(P[0], Q[0])
        elif i > 0 and j == 0:
            ca[i][j] = max(c(i - 1, 0), math.dist(P[i], Q[0]))
        elif i == 0 and j > 0:
            ca[i][j] = max(c(0, j - 1), math.dist(P[0], Q[j]))
        elif i > 0 and j > 0:
            ca[i][j] = max(
                min(c(i - 1, j), c(i - 1, j - 1), c(i, j - 1)),
                math.dist(P[i], Q[j])
            )
        else:
            ca[i][j] = float("inf")
        return ca[i][j]

    return c(len(P) - 1, len(Q) - 1)


def _agency_vehicle_name_map(agencies_df: pd.DataFrame, vehicles_df: pd.DataFrame) -> pd.DataFrame:
    """
    Returns a df: agency_id, vehicle_name.
    Assumes agencies has vehicle_id referencing vehicles.gid (common in your SDI).
    """
    a = agencies_df.copy()
    v = vehicles_df.copy()

    if "vehicle_id" in a.columns and "gid" in v.columns:
        out = a.merge(v[["gid", "name"]], left_on="vehicle_id", right_on="gid", how="left")
        out = out.rename(columns={"name": "vehicle_name"})
        return out[["agency_id", "vehicle_name"]].drop_duplicates()

    if "vehicle_name" in a.columns:
        return a[["agency_id", "vehicle_name"]].drop_duplicates()

    raise ValueError("Could not build agency->vehicle_name mapping (need agencies.vehicle_id + vehicles.gid or agencies.vehicle_name).")


def build_trip_stop_pair_segments_from_sdi(
    trips_view_gdf: gpd.GeoDataFrame,
    trip_stops_sequence_df: pd.DataFrame,
    agencies_df: pd.DataFrame,
    vehicles_df: pd.DataFrame,
) -> gpd.GeoDataFrame:
    """
    Builds trip stop-to-stop segments using:
      - transit.trips_view geometry (LineString)
      - transit.trip_stops_sequence distance_frac (0..1) and stop_sequence

    IMPORTANT: In RL2SDI GeoPackage exports, trip_stops_sequence often contains:
      - t_id (numeric internal trip gid)
      - observer_trip_id (string, the RouteLab/Observer trip id)
    While trips_view uses observer_id (string). So we prefer observer_trip_id when present.

    Output columns:
      trip_id, segment_order, from_id, to_id, vehicle_name, geometry, length_m
    """
    tv = trips_view_gdf.copy()
    tss = trip_stops_sequence_df.copy()

    # Attach vehicle_name to trips_view via agency->vehicle mapping (if needed)
    if "vehicle_name" not in tv.columns:
        vn = _agency_vehicle_name_map(agencies_df, vehicles_df)
        if "agency_id" in tv.columns:
            tv = tv.merge(vn, on="agency_id", how="left")

    # Expect canonical trip_id already normalized upstream
    if "trip_id" not in tss.columns:
        raise ValueError("trip_stops_sequence_df must include canonical trip_id (observer trip id).")

    tss["trip_id"] = tss["trip_id"].astype(str).str.strip()
    tv["trip_id"] = tv["trip_id"].astype(str).str.strip()

    if "stop_sequence" not in tss.columns:
        raise ValueError("trip_stops_sequence must contain stop_sequence.")
    if "stop_id" not in tss.columns:
        raise ValueError("trip_stops_sequence must contain stop_id.")
    if "distance_frac" not in tss.columns:
        raise ValueError("trip_stops_sequence must contain distance_frac (0..1).")

    # Ensure distance_frac is numeric
    tss["distance_frac"] = pd.to_numeric(tss["distance_frac"], errors="coerce")

    # --- Prep stop pairs ---
    tss = tss.sort_values(["trip_id", "stop_sequence"]).copy()
    tss["from_id"] = tss["stop_id"]
    tss["to_id"] = tss.groupby("trip_id")["stop_id"].shift(-1)
    tss["from_frac"] = tss["distance_frac"]
    tss["to_frac"] = tss.groupby("trip_id")["distance_frac"].shift(-1)
    tss["segment_order"] = tss.groupby("trip_id").cumcount()

    seg = tss.dropna(subset=["to_id", "to_frac", "from_frac"]).copy()

    # --- Map trip geometry & vehicle_name ---
    tv["trip_id"] = tv["trip_id"].astype(str).str.strip()
    geom_map = tv.set_index("trip_id")["geometry"].to_dict()
    veh_map = tv.set_index("trip_id")["vehicle_name"].to_dict() if "vehicle_name" in tv.columns else {}

    seg["trip_id"] = seg["trip_id"]
    seg["trip_geometry"] = seg["trip_id"].map(geom_map)
    seg["vehicle_name"] = seg["trip_id"].map(veh_map)

    # Defensive diagnostics: if nothing maps, fail loudly with hints
    if seg["trip_geometry"].isna().all():
        raise ValueError(
            "No trip geometries could be mapped from trips_view to trip_stops_sequence. "
            "This usually means trip IDs don't match. "
            "Check trips_view.trip_id values vs trip_stops_sequence.observer_trip_id (preferred) / trip_id / t_id."
        )

    # --- Build segment geometry using shapely substring with normalized fractions ---
    seg["geometry"] = seg.apply(
        lambda r: substring(
            r["trip_geometry"],
            float(min(r["from_frac"], r["to_frac"])),
            float(max(r["from_frac"], r["to_frac"])),
            normalized=True,
        )
        if isinstance(r["trip_geometry"], LineString)
        else None,
        axis=1,
    )

    out = gpd.GeoDataFrame(
        seg[["trip_id", "segment_order", "from_id", "to_id", "vehicle_name", "geometry"]],
        geometry="geometry",
        crs=tv.crs,
    )
    out = out[~out.geometry.isna()].copy()
    out["length_m"] = out.geometry.length
    out = out[out["length_m"] > 0].copy()

    return out


def build_flow_stop_pair_segments(trip_segments_gdf: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
    """
    Deduplicate stop-pair segments across trips (like the spirit of od_stats),
    keeping vehicle_name in the key.
    """
    cols = ["from_id", "to_id", "vehicle_name", "geometry"]
    base = trip_segments_gdf[cols].drop_duplicates().copy()
    base["flow_seg_id"] = (
        base["from_id"].astype(str) + "__" +
        base["to_id"].astype(str) + "__" +
        base["vehicle_name"].astype(str)
    )
    return gpd.GeoDataFrame(base, geometry="geometry", crs=trip_segments_gdf.crs)





def get_segments(curve):
    return list(map(LineString, zip(curve.coords[:-1], curve.coords[1:])))


def get_avg_occupancy_per_segment_v3_sdi_timeproxy(
    trip_segments_gdf,
    intervals_df,
    raw_stops_gdf,
    onboard_instances_gdf,
    feedback=None
):
    """
    Per-(trip, segment, interval) median occupancy from SDI onboard surveys.

    Each onboard survey carries its RouteLab-assigned interval as a TEXT id
    (raw.onboard_instances.interval_id). We join that directly to
    transit_intervals.observer_id to obtain the canonical interval — every stop
    event from a given survey is attributed to exactly one canonical interval,
    regardless of where its UTC timestamp lands relative to the interval window.
    This replaces the earlier approach of cumulating od_stats.duration along
    each trip and filtering rows whose estimated wall-clock arrival fell inside
    the window, which produced mid-trip-segment gaps for boundary-straddling
    rides and was vulnerable to UTC-vs-local-time bookkeeping errors. ID-based
    joining is also robust to weekday/weekend windows that share clock times,
    which a (start, end) crosswalk would silently mis-attribute.

    Inputs:
      - trip_segments_gdf : per-(trip, segment_order) line geometries
      - intervals_df      : canonical intervals; must include observer_id and
                            interval_name (sourced from transit.intervals)
      - raw_stops_gdf     : per-stop board/alight events with onboard_instance_observer_id
      - onboard_instances_gdf : per-survey rows; must include id, trip_id, status,
                                valid, departed_at, interval_id

    Returns:
      avg_occupancy_per_trip_segment_per_interval,
      onboard_segments_with_occupancy,
      matched_stops,
      filtered_stops
    """
    
    # ------------------------------------------------------------
    # Option A: snap raw stop events directly to stop-pair segments
    # (no micro-segmentation)
    # ------------------------------------------------------------
    
    segments_targets = trip_segments_gdf[[
        "trip_id", "segment_order", "from_id", "to_id", "vehicle_name", "geometry"
    ]].copy()
    segments_targets = gpd.GeoDataFrame(segments_targets, geometry="geometry", crs=trip_segments_gdf.crs)
    segments_targets = segments_targets[~segments_targets.geometry.isna()].copy()
    
    if segments_targets.empty:
        raise ValueError("trip_segments_gdf produced no valid geometries (segments_targets is empty).")
    
    # Keep the RouteLab-assigned interval id on each survey so we can attribute the
    # entire ride to one canonical interval without re-bucketing by timestamps.
    _required_oi_cols = ["id", "trip_id", "status", "valid", "departed_at", "interval_id"]
    _missing = [c for c in _required_oi_cols if c not in onboard_instances_gdf.columns]
    if _missing:
        raise ValueError(
            "onboard_instances_gdf is missing required column(s): "
            f"{_missing}. raw.onboard_instances must include interval_id (TEXT) — the "
            "RouteLab interval id stamped on the survey at assignment time."
        )
    valid_onboard_instances = onboard_instances_gdf[_required_oi_cols].copy()
    valid_onboard_instances = valid_onboard_instances.query("status == 'finished'")
    valid_onboard_instances["valid"] = (
    valid_onboard_instances["valid"]
    .astype(str)
    .str.strip()
    .str.lower()
    .map({"true": True, "false": False}))

    valid_onboard_instances = valid_onboard_instances[
        valid_onboard_instances["valid"] == True
    ]
    
    valid_raw_stops = raw_stops_gdf.drop(columns="gid", errors="ignore").merge(
        valid_onboard_instances,
        left_on="onboard_instance_observer_id",
        right_on="id",
        how="inner",
    )
    
    filtered_stops = (
        valid_raw_stops[valid_raw_stops["trip_id"].isin(segments_targets["trip_id"].unique())]
        .reset_index(drop=True)
        .copy()
    )
    filtered_stops["orig_index"] = filtered_stops.index
    
    def _snap_points_to_stop_pair_segments(group: gpd.GeoDataFrame) -> pd.DataFrame:
        """Return a DataFrame with segment_order/from_id/to_id/vehicle_name + orig_index for this trip."""
        segs = segments_targets[segments_targets["trip_id"] == group.name]
        if segs.empty:
            out = group[["orig_index"]].copy()
            out["segment_order"] = np.nan
            out["from_id"] = np.nan
            out["to_id"] = np.nan
            out["vehicle_name"] = np.nan
            return out[["segment_order", "from_id", "to_id", "vehicle_name", "orig_index"]]

        segs = segs.reset_index(drop=True)
        geoms = list(segs.geometry.values)
        tree = STRtree(geoms)
        geom_pos = {id(g): i for i, g in enumerate(geoms)}

        pos = []
        for pt in group.geometry.values:
            if pt is None or getattr(pt, "is_empty", False):
                pos.append(None)
                continue

            nearest_result = tree.nearest(pt)

            if nearest_result is None:
                pos.append(None)
            elif isinstance(nearest_result, (int, np.integer)):
                pos.append(int(nearest_result))
            else:
                pos.append(geom_pos.get(id(nearest_result)))

        # Build result
        rows = []
        for p in pos:
            if p is None:
                rows.append((np.nan, np.nan, np.nan, np.nan))
            else:
                r = segs.loc[p, ["segment_order", "from_id", "to_id", "vehicle_name"]]
                rows.append(tuple(r.values.tolist()))

        out = pd.DataFrame(rows, columns=["segment_order", "from_id", "to_id", "vehicle_name"])
        out["orig_index"] = group["orig_index"].values
        return out[["segment_order", "from_id", "to_id", "vehicle_name", "orig_index"]]

    matched_pairs = (
        filtered_stops.groupby("trip_id", group_keys=False)
        .apply(_snap_points_to_stop_pair_segments)
        .reset_index(drop=True)
    )

    # if feedback:
    #     feedback.pushInfo(f"raw stops rows: {len(valid_raw_stops)}")
    #     feedback.pushInfo(f"raw stops trip id: {valid_raw_stops["trip_id"]}")
    #     feedback.pushInfo(f"segments targets id: {segments_targets["trip_id"]}")
    #     feedback.pushInfo(f"segments columns: {matched_pairs.columns}")
    
    matched_stops = (
        matched_pairs.merge(
            filtered_stops,
            on="orig_index",
            how="left",
        )
        .groupby(
            ["onboard_instance_observer_id", "departed_at", "trip_id", "segment_order"],
            as_index=False,
        )
        .agg({"board": "sum", "alight": "sum", "created_at": "first"})
    )

    
    onboard_segments_with_occupancy = (
        segments_targets
        .merge(
            matched_stops[["trip_id", "onboard_instance_observer_id"]].drop_duplicates(),
            on="trip_id",
            how="inner",
        )
        .merge(
            matched_stops,
            on=["trip_id", "segment_order", "onboard_instance_observer_id"],
            how="left",
        )
        .sort_values(["trip_id", "onboard_instance_observer_id", "segment_order"])
        .fillna({"board": 0, "alight": 0})
        .assign(
            survey_departed_at=lambda df: df.groupby(
                ["trip_id", "onboard_instance_observer_id"]
            )["departed_at"].ffill().bfill()
        )
        .assign(
            survey_departed_at=lambda df: pd.to_datetime(df.survey_departed_at, utc=True)
            .dt.tz_localize(None)
        )
        .assign(
            survey_departed_at=lambda df: (
                df.survey_departed_at - df.survey_departed_at.dt.floor("d")
            ).dt.total_seconds()
        )
        .assign(
            vehicle_occupancy=lambda df: (
                df.groupby(["trip_id", "onboard_instance_observer_id"])
                .apply(lambda g: (g["board"] - g["alight"]).cumsum())
                .reset_index(level=[0, 1], drop=True)
            )
        )
    )
    
# ---- Attribute every onboard ride to its RouteLab-assigned canonical interval ----
    # Each onboard survey row carries an interval_id (TEXT) referencing the RouteLab
    # interval that the surveyor selected at assignment time. We join it directly to
    # transit_intervals.observer_id to get the canonical (gid, name). This is robust to
    # weekday/weekend windows that share clock times — a (start_time, end_time) crosswalk
    # would silently mis-attribute those.
    #
    # If your RouteLab export pipeline currently writes a different agency-specific id to
    # raw.onboard_instances.interval_id (so the join below fails), fix the exporter to
    # write the canonical RouteLab interval id (= transit_intervals.observer_id) instead.

    # ---- 1) Normalize intervals_df: must carry observer_id and interval_name ----
    iv = intervals_df.copy()
    if "interval_name" not in iv.columns:
        raise ValueError("intervals_df must include interval_name.")
    if "observer_id" not in iv.columns:
        raise ValueError(
            "intervals_df must include observer_id (the RouteLab interval id) so it can be "
            "joined to raw.onboard_instances.interval_id."
        )
    iv = iv[["observer_id", "interval_name"]].drop_duplicates("observer_id").copy()
    iv["observer_id"] = iv["observer_id"].astype(str).str.strip()

    # ---- 2) Build per-survey interval lookup keyed on the survey's instance id ----
    if "interval_id" not in valid_onboard_instances.columns:
        raise ValueError(
            "valid_onboard_instances must include interval_id (TEXT) — the RouteLab interval "
            "id stamped on each onboard survey at assignment time."
        )
    survey_intervals = (
        valid_onboard_instances[["id", "interval_id"]]
        .rename(columns={"id": "onboard_instance_observer_id"})
        .copy()
    )
    survey_intervals["interval_id"] = survey_intervals["interval_id"].astype(str).str.strip()

    # ---- 3) Direct ID join: onboard.interval_id ↔ transit_intervals.observer_id ----
    survey_intervals = survey_intervals.merge(
        iv.rename(columns={"observer_id": "interval_id"}),
        on="interval_id",
        how="left",
    )
    _unmatched = survey_intervals["interval_name"].isna().sum()
    if _unmatched > 0:
        bad = (
            survey_intervals.loc[survey_intervals["interval_name"].isna(), ["interval_id"]]
            .drop_duplicates()
            .head(5)["interval_id"]
            .tolist()
        )
        raise ValueError(
            f"{_unmatched} valid onboard surveys reference interval_id values not present in "
            f"transit_intervals.observer_id. Examples: {bad}. "
            "This typically means the RouteLab→SDI exporter is writing an agency-specific "
            "interval id to raw.onboard_instances.interval_id instead of the canonical "
            "RouteLab interval id. Fix the exporter, or align the ids in the gpkg."
        )

    if feedback:
        feedback.pushInfo(
            f"Joined {len(survey_intervals)} valid surveys to "
            f"{survey_intervals['interval_name'].nunique()} canonical interval(s) via interval_id."
        )

    # ---- 4) Attach canonical interval_name to every (instance × segment) occupancy row ----
    occ = onboard_segments_with_occupancy.merge(
        survey_intervals[["onboard_instance_observer_id", "interval_name"]],
        on="onboard_instance_observer_id",
        how="left",
    )

    if occ["interval_name"].isna().any():
        n_missing = occ["interval_name"].isna().sum()
        raise ValueError(
            f"{n_missing} occupancy rows have no interval after crosswalk. "
            "Check that every onboard_instance_observer_id is present in valid_onboard_instances."
        )

    # ---- 5) Median occupancy per (trip, segment, interval) — no time-window filter ----
    avg_occupancy_per_trip_segment_per_interval = (
        occ.groupby(["trip_id", "segment_order", "interval_name"], as_index=False)
           .agg({"vehicle_occupancy": "median"})
    )

    return (
        avg_occupancy_per_trip_segment_per_interval,
        onboard_segments_with_occupancy,
        matched_stops,
        filtered_stops,
    )
