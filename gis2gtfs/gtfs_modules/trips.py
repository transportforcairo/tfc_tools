'''
Build trips.txt with one row per (base_trip, active_interval) pair.

Changes vs the original version
-------------------------------
Originally trips.txt had one row per base trip with service_id copied
from transit_trips (typically the hardcoded 'Ground_Daily'). Two things
broke:
  - frequencies.txt rows per interval collided on (trip_id, start_time)
    whenever two intervals shared a start time.
  - calendar.txt claimed the service ran every day even when intervals
    carried mon..sun flags saying otherwise.

This version:
  - Cross-joins base trips x active intervals so every frequencies row
    attaches to a distinct GTFS trip. trip_id = '<base>_<suffix>' where
    suffix comes from the interval name.
  - Takes each expanded trip's service_id from its interval's day-pattern
    (via the shared helper, which maps mon..sun to a canonical service
    name like 'svc_weekday' / 'svc_weekend' / 'svc_daily'). A base trip
    whose intervals span multiple day-patterns therefore produces GTFS
    trips on multiple services, which is the canonical GTFS way to
    express that.
  - Falls back to the legacy service_id from transit_trips only when the
    intervals.csv lacks day columns (older SDI exports).

shape_id continues to point at the BASE trip's shape so shapes.txt does
NOT need to be duplicated -- every interval-copy of a trip shares the
same shape, which is correct and keeps the feed small.
'''
import os
import pandas as pd
import geopandas as gpd

from ._interval_expansion import load_intervals_and_suffixes, expanded_trip_id


def generate(data_dir, data_raw_dir):
    '''Generate GTFS trips.txt by cross-joining base trips with active intervals.'''
    print("Generating trips.txt...")

    trips_path = os.path.join(data_raw_dir, "trips.geojson")

    # Load trips and remove duplicates on the BASE trip_id.
    trips_gdf = gpd.read_file(trips_path)
    trips_gdf["base_trip_id"] = trips_gdf["observer_id"].astype(str)
    trips_gdf = trips_gdf.drop_duplicates(subset="base_trip_id").copy()

    # trip_headsign from destination, shape_id keyed to the BASE trip so
    # we can reuse one shape across the interval-copies.
    trips_gdf["trip_headsign"] = trips_gdf["destination"]
    trips_gdf["shape_id"] = trips_gdf["base_trip_id"].astype(str) + "_Shape"

    # Load active intervals + per-interval suffix map.
    iv, _suffix_by_gid = load_intervals_and_suffixes(data_raw_dir)
    if iv.empty:
        raise ValueError(
            "intervals.csv has no active intervals; cannot expand trips."
        )

    # Cross-join base trips x active intervals.
    #
    # Each interval carries its own service_id (derived from mon..sun in
    # the raw intervals.csv by the shared helper). Each expanded trip
    # adopts the service_id of ITS interval -- so a base trip that has
    # frequency rows on both weekday and weekend intervals ends up as
    # multiple GTFS trips, one per service. That's choice (a) from the
    # design discussion: split by service_id rather than fold everything
    # onto a single 7-day calendar.
    #
    # If the intervals.csv didn't carry day columns (legacy export), the
    # 'service_id' column from the helper is None; we then fall back to
    # the service_id column carried by transit_trips (historical
    # 'Ground_Daily'), matching the legacy behaviour.
    base_cols = [
        "base_trip_id", "route_id", "service_id", "trip_headsign",
        "direction_id", "shape_id",
    ]
    base_df = trips_gdf[base_cols].copy().rename(
        columns={"service_id": "_service_id_legacy"}
    )
    intervals_keep = iv[["gid", "suffix", "service_id"]].rename(
        columns={"gid": "interval_id", "service_id": "_service_id_interval"}
    )
    intervals_keep["interval_id"] = intervals_keep["interval_id"].astype(str)
    expanded = base_df.merge(intervals_keep, how="cross")

    # Prefer the interval-derived service_id. Only fall back to the legacy
    # value when the intervals didn't carry day columns.
    expanded["service_id"] = expanded["_service_id_interval"].where(
        expanded["_service_id_interval"].notna(),
        expanded["_service_id_legacy"],
    )

    expanded["trip_id"] = [
        expanded_trip_id(bt, suf)
        for bt, suf in zip(expanded["base_trip_id"], expanded["suffix"])
    ]

    # Guard: every expanded trip_id must be unique.
    dup = expanded["trip_id"].duplicated()
    if dup.any():
        raise ValueError(
            f"Expanded trips.txt would contain duplicate trip_id values "
            f"(sample: {expanded.loc[dup, 'trip_id'].head().tolist()})"
        )

    trips_out = expanded[
        ["route_id", "service_id", "trip_headsign",
         "direction_id", "shape_id", "trip_id"]
    ].copy()

    trips_out.to_csv(os.path.join(data_dir, "trips.txt"),
                     index=False, encoding="utf-8")

    # Companion file: stop_times.py does NOT need this (it cross-joins
    # intervals itself), but frequencies.py and any downstream tooling can
    # use it to know the full (base_trip_id, interval_id, trip_id) mapping
    # that trips.txt represents.
    companion_cols = [
        "base_trip_id", "interval_id", "trip_id",
        "route_id", "service_id", "trip_headsign",
        "direction_id", "shape_id",
    ]
    expanded[companion_cols].to_csv(
        os.path.join(data_raw_dir, "trips_with_intervals.csv"),
        index=False, encoding="utf-8",
    )

    print(
        f"trips.txt written ({len(trips_out)} rows = "
        f"{len(trips_gdf)} base trips x {len(iv)} intervals)."
    )
