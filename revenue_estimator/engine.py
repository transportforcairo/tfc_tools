from __future__ import annotations

import math
import os
import sqlite3
from dataclasses import dataclass
from typing import Dict, List, Tuple

import geopandas as gpd
import pandas as pd
from shapely.ops import nearest_points

from ..tfc_tools_common.sdi_io import SDISource, read_df, read_gdf


@dataclass
class RevenueConfig:
    behavior_model: str = "proportional"  # fifo | proportional | distance_weighted
    qa_mode: str = "warn"                 # strict | warn | repair_minor
    drop_zero_stops: bool = True
    band1_max_fraction: float = 0.33
    band1_share: float = 0.50
    band2_max_fraction: float = 0.66
    band2_share: float = 0.75
    min_headway_secs: float = 60.0
    final_load_tolerance: float = 0.0
    max_snap_distance_m: float = 100.0
    write_od_matrix: bool = True
    write_stop_profile: bool = True
    write_qa_tables: bool = True


class RevenueEstimator:
    def __init__(self, source: SDISource, output_folder: str, config: RevenueConfig):
        self.source = source
        self.output_folder = output_folder
        self.config = config

    def run(self, feedback=None) -> Dict[str, str]:
        os.makedirs(self.output_folder, exist_ok=True)
        gpkg_out = os.path.join(self.output_folder, "revenue_estimation.gpkg")
        csv_out = os.path.join(self.output_folder, "revenue_route_direction_summary.csv")
        if os.path.exists(gpkg_out):
            os.remove(gpkg_out)

        if feedback:
            feedback.pushInfo("Reading SDI layers …")
        onboard = read_gdf(self.source, "raw.onboard_instances", None, "geometry")
        stops = read_gdf(self.source, "raw.stops")
        trips = read_gdf(self.source, "transit.trips_view")
        trips_intervals = read_df(self.source, "transit.trips_intervals")
        intervals = read_df(self.source, "transit.intervals")

        onboard, stops, trips, trips_intervals, intervals = self._standardize(onboard, stops, trips, trips_intervals, intervals)
        onboard = self._filter_onboard(onboard)
        stops = stops[stops["onboard_instance_observer_id"].isin(onboard["onboard_instance_id"])].copy()
        if self.config.drop_zero_stops:
            stops = stops[(stops["board"] != 0) | (stops["alight"] != 0)].copy()

        if feedback:
            feedback.pushInfo("Building stop profiles and inferred OD matrices …")
        stop_profiles: List[gpd.GeoDataFrame] = []
        od_rows: List[Dict] = []
        trip_rows: List[Dict] = []
        qa_rows: List[Dict] = []

        trips_idx = trips.set_index("trip_ref", drop=False)
        onboard = onboard.sort_values(["trip_ref", "onboard_instance_id"])

        for _, oi in onboard.iterrows():
            inst_id = oi["onboard_instance_id"]
            inst_stops = stops[stops["onboard_instance_observer_id"] == inst_id].copy().sort_values("stop_time")
            if inst_stops.empty:
                qa_rows.append(self._qa_row(oi, "NO_STOPS", "No raw stops found for onboard instance", "error"))
                continue
            if oi["trip_ref"] not in trips_idx.index:
                qa_rows.append(self._qa_row(oi, "MISSING_TRIP", "Trip not found in transit_trips_view", "error"))
                continue
            trip = trips_idx.loc[oi["trip_ref"]]
            try:
                prof, trip_qa = self._build_stop_profile(inst_stops, oi, trip)
                stop_profiles.append(prof)
                qa_rows.extend(trip_qa)
                od, od_qa = self._infer_od(prof, oi, trip)
                od_rows.extend(od)
                qa_rows.extend(od_qa)
                trip_rows.append(self._build_trip_trace_row(prof, od, oi, trip))
            except Exception as e:
                if self.config.qa_mode == "strict":
                    raise
                qa_rows.append(self._qa_row(oi, "PROCESSING_FAILED", str(e), "error"))

        stop_profile_gdf = gpd.GeoDataFrame(pd.concat(stop_profiles, ignore_index=True), geometry="geometry", crs=stops.crs) if stop_profiles else gpd.GeoDataFrame(columns=["geometry"], geometry="geometry", crs=stops.crs)
        od_df = pd.DataFrame(od_rows)
        trip_trace_df = pd.DataFrame(trip_rows)
        qa_df = pd.DataFrame(qa_rows)

        if feedback:
            feedback.pushInfo("Applying interval scaling and creating summaries …")
        trip_trace_df = self._apply_interval_scaling(trip_trace_df, trips_intervals, intervals)
        route_summary_df = self._build_route_summary(trip_trace_df, trips, trips_intervals, intervals)
        qa_summary_df = self._build_qa_summary(qa_df)

        self._write_outputs(gpkg_out, csv_out, stop_profile_gdf, od_df, trip_trace_df, route_summary_df, qa_df, qa_summary_df)
        if feedback:
            feedback.pushInfo(f"Wrote outputs to {gpkg_out}")
            feedback.pushInfo(f"Wrote route-direction summary to {csv_out}")
        return {
            "OUTPUT_GPKG": gpkg_out,
            "OUTPUT_SUMMARY_CSV": csv_out,
        }

    def _standardize(self, onboard, stops, trips, trips_intervals, intervals):
        onboard = onboard.rename(columns={
            "id": "onboard_instance_id",
            "trip_id": "trip_ref",
            "interval_id": "interval_observer_id",  # Postgres path: interval FK id
            "interval": "interval_name_survey",      # GeoPackage path: interval name text
            "created_at": "onboard_created_at",
        }).copy()
        if "geometry" not in onboard.columns:
            onboard = onboard.rename(columns={onboard.geometry.name: "geometry"}).set_geometry("geometry")
        onboard["valid_norm"] = onboard.get("valid", pd.Series(index=onboard.index)).astype(str).str.lower()
        onboard["status_norm"] = onboard.get("status", pd.Series(index=onboard.index)).astype(str).str.lower()

        stops = stops.rename(columns={
            "created_at": "stop_time",
        }).copy()
        if "geometry" not in stops.columns:
            stops = stops.rename(columns={stops.geometry.name: "geometry"}).set_geometry("geometry")
        for c in ["board", "alight"]:
            stops[c] = pd.to_numeric(stops[c], errors="coerce").fillna(0.0)
        stops["board"] = stops["board"].clip(lower=0)
        stops["alight"] = stops["alight"].clip(lower=0)

        trips = trips.rename(columns={
            "observer_id": "trip_ref",
            "fare": "full_fare",
        }).copy()
        if "geometry" not in trips.columns:
            trips = trips.rename(columns={trips.geometry.name: "geometry"}).set_geometry("geometry")
        trips["full_fare"] = pd.to_numeric(trips["full_fare"], errors="coerce")
        trips["trip_gid"] = pd.to_numeric(trips["gid"], errors="coerce")
        trips["trip_len_m"] = pd.to_numeric(trips.get("len_km"), errors="coerce") * 1000.0

        trips_intervals = trips_intervals.rename(columns={"trip_id": "trip_gid"}).copy()
        trips_intervals["trip_gid"] = pd.to_numeric(trips_intervals["trip_gid"], errors="coerce")
        trips_intervals["headway_secs"] = pd.to_numeric(trips_intervals["headway_secs"], errors="coerce")

        intervals = intervals.rename(columns={
            "gid": "interval_gid",
            "observer_id": "interval_observer_id",
            "name": "interval_name",
        }).copy()
        intervals["interval_gid"] = pd.to_numeric(intervals["interval_gid"], errors="coerce")
        intervals["active"] = pd.to_numeric(intervals["active"], errors="coerce").fillna(0).astype(int)
        return onboard, stops, trips, trips_intervals, intervals

    def _filter_onboard(self, onboard: gpd.GeoDataFrame) -> gpd.GeoDataFrame:
        mask = onboard["status_norm"].eq("finished")
        if "valid_norm" in onboard.columns:
            mask &= onboard["valid_norm"].isin(["true", "1", "yes"])
        return onboard[mask].copy()

    def _build_stop_profile(self, inst_stops: gpd.GeoDataFrame, oi: pd.Series, trip: pd.Series) -> Tuple[gpd.GeoDataFrame, List[Dict]]:
        qa = []
        inst_stops = inst_stops.copy().sort_values(["stop_time", "gid" if "gid" in inst_stops.columns else "stop_time"]).reset_index(drop=True)
        inst_stops["stop_sequence"] = range(1, len(inst_stops) + 1)
        line = trip.geometry
        if line is None or line.is_empty:
            raise ValueError("Missing trip geometry")
        line_len_m = trip.get("trip_len_m")
        if not pd.notna(line_len_m) or line_len_m <= 0:
            line_len_m = 1.0
        # normalized fraction along line; robust enough for route fractions
        fractions = []
        snap_dists = []
        for geom in inst_stops.geometry:
            frac = float(line.project(geom, normalized=True)) if geom is not None and not geom.is_empty else float("nan")
            fractions.append(frac)
            try:
                nearest = nearest_points(geom, line)[1]
                snap_d = float(geom.distance(nearest))
            except Exception:
                snap_d = float("nan")
            snap_dists.append(snap_d)
        inst_stops["distance_fraction"] = fractions
        inst_stops["distance_along_m"] = inst_stops["distance_fraction"].fillna(0) * line_len_m
        inst_stops["snap_distance_native"] = snap_dists

        prev = 0.0
        loads_before = []
        loads_after = []
        min_load = 0.0
        for _, r in inst_stops.iterrows():
            loads_before.append(prev)
            after = prev + float(r["board"]) - float(r["alight"])
            loads_after.append(after)
            prev = after
            min_load = min(min_load, after)
        inst_stops["load_before"] = loads_before
        inst_stops["load_after"] = loads_after
        inst_stops["trip_ref"] = oi["trip_ref"]
        inst_stops["route_id"] = trip.get("route_id")
        inst_stops["direction_id"] = trip.get("direction_id")

        if min_load < 0:
            qa.append(self._qa_row(oi, "NEGATIVE_INTERMEDIATE_LOAD", f"Minimum load {min_load}", "error"))
        final_load = float(inst_stops["load_after"].iloc[-1])
        if abs(final_load) > self.config.final_load_tolerance:
            qa.append(self._qa_row(oi, "NONZERO_FINAL_LOAD", f"Final load {final_load}", "warning"))
        if inst_stops["stop_time"].duplicated().any():
            qa.append(self._qa_row(oi, "DUPLICATE_STOP_TIMESTAMPS", "Duplicate stop timestamps detected", "warning"))
        return inst_stops, qa

    def _infer_od(self, prof: gpd.GeoDataFrame, oi: pd.Series, trip: pd.Series) -> Tuple[List[Dict], List[Dict]]:
        behavior = self.config.behavior_model
        qa = []
        cohorts: List[Dict] = []
        out: List[Dict] = []
        for _, row in prof.iterrows():
            board = float(row["board"])
            alight = float(row["alight"])
            if board > 0:
                cohorts.append({
                    "remaining": board,
                    "origin_stop_seq": int(row["stop_sequence"]),
                    "origin_distance_m": float(row["distance_along_m"]),
                    "origin_fraction": float(row["distance_fraction"]),
                })
            if alight > 0:
                allocations, shortfall = self._allocate_alightings(cohorts, alight, float(row["distance_along_m"]), int(row["stop_sequence"]), behavior)
                for a in allocations:
                    travelled_m = max(0.0, float(row["distance_along_m"]) - a["origin_distance_m"])
                    frac = max(0.0, min(1.0, float(row["distance_fraction"]) - a["origin_fraction"]))
                    fare_share, fare_band = self._fare_share(frac)
                    est_paid = float(trip["full_fare"]) * fare_share if pd.notna(trip.get("full_fare")) else math.nan
                    out.append({
                        "onboard_instance_id": oi["onboard_instance_id"],
                        "trip_ref": oi["trip_ref"],
                        "trip_gid": trip.get("trip_gid"),
                        "route_id": trip.get("route_id"),
                        "direction_id": trip.get("direction_id"),
                        "origin_stop_seq": a["origin_stop_seq"],
                        "destination_stop_seq": int(row["stop_sequence"]),
                        "origin_distance_m": a["origin_distance_m"],
                        "destination_distance_m": float(row["distance_along_m"]),
                        "traveled_m": travelled_m,
                        "traveled_fraction": frac,
                        "passenger_count": a["count"],
                        "fare_band": fare_band,
                        "fare_share": fare_share,
                        "full_fare": trip.get("full_fare"),
                        "estimated_paid_fare": est_paid,
                        "flow_revenue": est_paid * a["count"] if pd.notna(est_paid) else math.nan,
                        "behavior_model": behavior,
                        "inferred_terminal_alight": False,
                    })
                if shortfall > 1e-9:
                    qa.append(self._qa_row(oi, "ALIGHTING_SHORTFALL", f"Could not allocate {shortfall} alightings at stop {int(row['stop_sequence'])}", "warning"))
        # Alight any passengers still onboard at the terminal (last observed stop).
        # These are real passengers who rode to the end of the route but whose alighting
        # was not recorded (enumerator disembarked early, or survey ended at the terminal).
        # We assign them destination = terminal (fraction = 1.0, distance = trip_len_m).
        terminal_seq = int(prof["stop_sequence"].iloc[-1]) + 1
        terminal_dist_m = float(trip.get("trip_len_m") or 0.0)
        terminal_frac = 1.0
        remaining_cohorts = [c for c in cohorts if c["remaining"] > 1e-12]
        if remaining_cohorts:
            for c in remaining_cohorts:
                travelled_m = max(0.0, terminal_dist_m - c["origin_distance_m"])
                frac = max(0.0, min(1.0, terminal_frac - c["origin_fraction"]))
                fare_share, fare_band = self._fare_share(frac)
                est_paid = float(trip["full_fare"]) * fare_share if pd.notna(trip.get("full_fare")) else math.nan
                out.append({
                    "onboard_instance_id": oi["onboard_instance_id"],
                    "trip_ref": oi["trip_ref"],
                    "trip_gid": trip.get("trip_gid"),
                    "route_id": trip.get("route_id"),
                    "direction_id": trip.get("direction_id"),
                    "origin_stop_seq": c["origin_stop_seq"],
                    "destination_stop_seq": terminal_seq,
                    "origin_distance_m": c["origin_distance_m"],
                    "destination_distance_m": terminal_dist_m,
                    "traveled_m": travelled_m,
                    "traveled_fraction": frac,
                    "passenger_count": c["remaining"],
                    "fare_band": fare_band,
                    "fare_share": fare_share,
                    "full_fare": trip.get("full_fare"),
                    "estimated_paid_fare": est_paid,
                    "flow_revenue": est_paid * c["remaining"] if pd.notna(est_paid) else math.nan,
                    "behavior_model": behavior,
                    "inferred_terminal_alight": True,
                })
                c["remaining"] = 0.0

        # unmatched riders still onboard at trip end (should now always be zero)
        leftover = sum(max(0.0, c["remaining"]) for c in cohorts)
        if leftover > 1e-9:
            qa.append(self._qa_row(oi, "UNMATCHED_ONBOARD", f"{leftover} passengers remained unmatched at trip end", "warning"))
        return out, qa

    def _allocate_alightings(self, cohorts: List[Dict], alight: float, dest_distance_m: float, dest_seq: int, behavior: str):
        allocations = []
        remaining_to_allocate = alight
        active = [c for c in cohorts if c["remaining"] > 1e-12]
        if behavior == "fifo":
            for c in active:
                if remaining_to_allocate <= 1e-12:
                    break
                take = min(c["remaining"], remaining_to_allocate)
                if take > 0:
                    c["remaining"] -= take
                    remaining_to_allocate -= take
                    allocations.append({"origin_stop_seq": c["origin_stop_seq"], "origin_distance_m": c["origin_distance_m"], "origin_fraction": c["origin_fraction"], "count": take})
            return allocations, remaining_to_allocate

        while remaining_to_allocate > 1e-9:
            active = [c for c in cohorts if c["remaining"] > 1e-12]
            if not active:
                break
            if behavior == "proportional":
                weights = [c["remaining"] for c in active]
            else:  # distance_weighted
                weights = [c["remaining"] * max(dest_distance_m - c["origin_distance_m"], 1.0) for c in active]
            total_w = sum(weights)
            if total_w <= 0:
                break
            changed = False
            for c, w in zip(active, weights):
                share = remaining_to_allocate * (w / total_w)
                take = min(c["remaining"], share)
                if take > 1e-12:
                    c["remaining"] -= take
                    allocations.append({"origin_stop_seq": c["origin_stop_seq"], "origin_distance_m": c["origin_distance_m"], "origin_fraction": c["origin_fraction"], "count": take})
                    changed = True
            new_remaining = alight - sum(a["count"] for a in allocations)
            # cleanup tiny negatives
            for c in cohorts:
                if abs(c["remaining"]) < 1e-12:
                    c["remaining"] = 0.0
            if (not changed) or abs(new_remaining - remaining_to_allocate) < 1e-9:
                break
            remaining_to_allocate = new_remaining
        return allocations, max(0.0, alight - sum(a["count"] for a in allocations))

    def _fare_share(self, frac: float):
        frac = max(0.0, min(1.0, frac))
        if frac <= self.config.band1_max_fraction:
            return self.config.band1_share, f"<= {self.config.band1_max_fraction:.2f}"
        if frac <= self.config.band2_max_fraction:
            return self.config.band2_share, f"<= {self.config.band2_max_fraction:.2f}"
        return 1.0, "> upper_band"

    def _build_trip_trace_row(self, prof: gpd.GeoDataFrame, od_rows: List[Dict], oi: pd.Series, trip: pd.Series) -> Dict:
        od_df = pd.DataFrame([r for r in od_rows if r["onboard_instance_id"] == oi["onboard_instance_id"]])
        observed_trip_revenue = float(od_df["flow_revenue"].sum()) if not od_df.empty else 0.0
        total_pax = float(od_df["passenger_count"].sum()) if not od_df.empty else 0.0
        avg_frac = float((od_df["traveled_fraction"] * od_df["passenger_count"]).sum() / total_pax) if total_pax > 0 else math.nan
        avg_fare = float((od_df["estimated_paid_fare"] * od_df["passenger_count"]).sum() / total_pax) if total_pax > 0 else math.nan
        terminal_alight_pax = float(od_df.loc[od_df["inferred_terminal_alight"] == True, "passenger_count"].sum()) if not od_df.empty else 0.0
        return {
            "onboard_instance_id": oi["onboard_instance_id"],
            "trip_ref": oi["trip_ref"],
            "trip_gid": trip.get("trip_gid"),
            "route_id": trip.get("route_id"),
            "direction_id": trip.get("direction_id"),
            "vehicle_name": trip.get("vehicle_name"),
            "origin": trip.get("origin"),
            "destination": trip.get("destination"),
            "interval_observer_id": oi.get("interval_observer_id"),
            "interval_name": oi.get("interval_name"),
            "full_fare": trip.get("full_fare"),
            "total_board": float(prof["board"].sum()),
            "total_alight": float(prof["alight"].sum()),
            "terminal_alight": terminal_alight_pax,
            "final_load": float(prof["load_after"].iloc[-1]),
            "observed_trip_revenue": observed_trip_revenue,
            "avg_inferred_trip_fraction": avg_frac,
            "avg_inferred_fare": avg_fare,
            "behavior_model": self.config.behavior_model,
            "fare_model_type": "banded_fraction_of_full_fare",
        }

    def _apply_interval_scaling(self, trip_trace_df: pd.DataFrame, trips_intervals: pd.DataFrame, intervals: pd.DataFrame) -> pd.DataFrame:
        if trip_trace_df.empty:
            return trip_trace_df
        iv = intervals[["interval_gid", "interval_observer_id", "interval_name", "start_time", "end_time", "active"]].copy()
        ti = trips_intervals.copy()

        # Primary join: via interval_observer_id (Postgres path — interval FK on onboard instance)
        merged = trip_trace_df.merge(iv, on="interval_observer_id", how="left", suffixes=("", "_iv"))

        # Fallback join: via interval_name_survey (GeoPackage path — interval name text on onboard instance)
        # Used for any rows where the primary join produced no interval_gid (i.e. interval_id was absent).
        if "interval_name_survey" in merged.columns:
            missing = merged["interval_gid"].isna()
            if missing.any():
                iv_by_name = iv[["interval_gid", "interval_name", "start_time", "end_time", "active"]].rename(
                    columns={
                        "interval_gid": "interval_gid_fb",
                        "start_time": "start_time_fb",
                        "end_time": "end_time_fb",
                        "active": "active_fb",
                    }
                )
                fallback = merged.loc[missing, ["interval_name_survey"]].merge(
                    iv_by_name, left_on="interval_name_survey", right_on="interval_name", how="left"
                )
                merged.loc[missing, "interval_gid"] = fallback["interval_gid_fb"].values
                merged.loc[missing, "start_time"] = fallback["start_time_fb"].values
                merged.loc[missing, "end_time"] = fallback["end_time_fb"].values
                merged.loc[missing, "active"] = fallback["active_fb"].values
                if "interval_name_iv" not in merged.columns:
                    merged.loc[missing, "interval_name"] = fallback["interval_name"].values
                else:
                    merged.loc[missing, "interval_name_iv"] = fallback["interval_name"].values

        merged = merged.merge(ti[["trip_gid", "interval_id", "headway_secs", "headway_estimation_method"]], left_on=["trip_gid", "interval_gid"], right_on=["trip_gid", "interval_id"], how="left")
        merged["interval_duration_secs"] = merged.apply(lambda r: _interval_duration_secs(r.get("start_time"), r.get("end_time")), axis=1)
        merged["vehicle_trips_in_interval"] = merged.apply(lambda r: (r["interval_duration_secs"] / r["headway_secs"]) if pd.notna(r["headway_secs"]) and r["headway_secs"] >= self.config.min_headway_secs else math.nan, axis=1)
        merged["estimated_interval_revenue"] = merged["observed_trip_revenue"] * merged["vehicle_trips_in_interval"]
        merged["interval_name"] = merged["interval_name_iv"].combine_first(merged["interval_name"]) if "interval_name_iv" in merged.columns else merged["interval_name"]
        return merged.drop(columns=[c for c in ["interval_name_iv", "interval_id", "interval_name_survey"] if c in merged.columns])

    def _build_route_summary(self, trip_trace_df: pd.DataFrame, trips: gpd.GeoDataFrame, trips_intervals: pd.DataFrame, intervals: pd.DataFrame) -> pd.DataFrame:
        base_cols = ["route_id", "direction_id", "origin", "destination", "vehicle_name", "full_fare"]
        if trip_trace_df.empty:
            return pd.DataFrame(columns=base_cols + ["sample_trip_count", "sample_total_board", "sample_avg_trip_revenue", "daily_estimated_revenue", "active_interval_count"])
        base = trip_trace_df.groupby(base_cols, dropna=False).agg(
            sample_trip_count=("onboard_instance_id", "nunique"),
            sample_total_board=("total_board", "sum"),
            sample_avg_trip_revenue=("observed_trip_revenue", "mean"),
        ).reset_index()

        # avg observed trip revenue by route/direction/interval from sample
        sample_iv = trip_trace_df.groupby(base_cols + ["interval_name"], dropna=False).agg(
            sample_avg_trip_revenue_interval=("observed_trip_revenue", "mean")
        ).reset_index()

        # supply intervals from trip definitions
        supply = trips[["trip_gid", "route_id", "direction_id", "origin", "destination", "vehicle_name", "full_fare"]].merge(
            trips_intervals[["trip_gid", "interval_id", "headway_secs"]], on="trip_gid", how="left"
        ).merge(
            intervals[["interval_gid", "interval_name", "start_time", "end_time", "active"]], left_on="interval_id", right_on="interval_gid", how="left"
        )
        supply["interval_duration_secs"] = supply.apply(lambda r: _interval_duration_secs(r.get("start_time"), r.get("end_time")), axis=1)
        supply["vehicle_trips_in_interval"] = supply.apply(lambda r: (r["interval_duration_secs"] / r["headway_secs"]) if pd.notna(r["headway_secs"]) and r["headway_secs"] >= self.config.min_headway_secs else math.nan, axis=1)
        supply = supply.groupby(base_cols + ["interval_name"], dropna=False).agg(
            headway_secs=("headway_secs", "mean"),
            vehicle_trips_in_interval=("vehicle_trips_in_interval", "mean"),
            active=("active", "max"),
        ).reset_index()

        summary = supply.merge(sample_iv, on=base_cols + ["interval_name"], how="left")
        summary["interval_revenue"] = summary["sample_avg_trip_revenue_interval"] * summary["vehicle_trips_in_interval"]

        # pivot wide
        rows = []
        for keys, grp in summary.groupby(base_cols, dropna=False):
            row = {col: val for col, val in zip(base_cols, keys if isinstance(keys, tuple) else (keys,))}
            for _, g in grp.iterrows():
                iv = g.get("interval_name")
                if pd.isna(iv):
                    continue
                safe = str(iv)
                row[f"headway_secs__{safe}"] = g.get("headway_secs")
                row[f"vehicle_trips__{safe}"] = g.get("vehicle_trips_in_interval")
                row[f"sample_avg_trip_revenue__{safe}"] = g.get("sample_avg_trip_revenue_interval")
                row[f"interval_revenue__{safe}"] = g.get("interval_revenue")
                row[f"active__{safe}"] = g.get("active")
            rows.append(row)
        wide = pd.DataFrame(rows)
        out = base.merge(wide, on=base_cols, how="left")
        interval_rev_cols = [c for c in out.columns if c.startswith("interval_revenue__")]
        active_cols = [c for c in out.columns if c.startswith("active__")]
        out["daily_estimated_revenue"] = out[interval_rev_cols].sum(axis=1, skipna=True) if interval_rev_cols else math.nan
        out["active_interval_count"] = out[active_cols].fillna(0).astype(float).sum(axis=1) if active_cols else 0
        return out

    def _build_qa_summary(self, qa_df: pd.DataFrame) -> pd.DataFrame:
        if qa_df.empty:
            return pd.DataFrame(columns=["qa_issue_code", "trip_count", "record_count"])
        return qa_df.groupby("qa_issue_code", dropna=False).agg(
            trip_count=("onboard_instance_id", "nunique"),
            record_count=("onboard_instance_id", "size"),
        ).reset_index()

    def _write_outputs(self, gpkg_out, csv_out, stop_profile_gdf, od_df, trip_trace_df, route_summary_df, qa_df, qa_summary_df):
        if self.config.write_stop_profile and not stop_profile_gdf.empty:
            stop_profile_gdf.to_file(gpkg_out, layer="revenue_stop_profile", driver="GPKG", index=False)
        with sqlite3.connect(gpkg_out) as conn:
            if self.config.write_od_matrix:
                od_df.to_sql("revenue_od_matrix", conn, if_exists="replace", index=False)
            trip_trace_df.to_sql("revenue_trip_trace", conn, if_exists="replace", index=False)
            route_summary_df.to_sql("revenue_route_direction_summary", conn, if_exists="replace", index=False)
            if self.config.write_qa_tables:
                qa_df.to_sql("revenue_qa_trips", conn, if_exists="replace", index=False)
                qa_summary_df.to_sql("revenue_qa_summary", conn, if_exists="replace", index=False)
        route_summary_df.to_csv(csv_out, index=False)

    def _qa_row(self, oi: pd.Series, code: str, detail: str, severity: str) -> Dict:
        return {
            "onboard_instance_id": oi.get("onboard_instance_id"),
            "trip_ref": oi.get("trip_ref"),
            "route_id": oi.get("route_id"),
            "direction_id": oi.get("direction_id"),
            "qa_issue_code": code,
            "qa_issue_detail": detail,
            "severity": severity,
        }


def _interval_duration_secs(start, end) -> float:
    if pd.isna(start) or pd.isna(end):
        return math.nan
    start = pd.to_timedelta(str(start))
    end = pd.to_timedelta(str(end))
    secs = (end - start).total_seconds()
    if secs < 0:
        secs += 24 * 3600
    return secs
