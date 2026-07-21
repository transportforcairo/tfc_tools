#!/usr/bin/env python3
"""
key_repair.py — safe, unambiguous cross-source key repairs (spec section 9).

The Refresh tool relies on full rebuild + validation for correctness, and blocks
/ warns on broken foreign keys. This module handles the *narrow* subset of FK
breakages that can be repaired mechanically with zero guessing, because they are
representation mismatches rather than genuine semantic renames:

    A. agency_id gid -> code
       transit_trips_view.agency_id (and transit_trips, if present) sometimes
       stores the agency SURROGATE gid (e.g. "1") instead of the text code
       (e.g. "AMUGA"). GTFS requires the code, and the value "1" is an orphan
       against agencies.agency_id. Because "1" maps unambiguously to exactly one
       agencies.gid, we can rewrite it to that agency's code. (This is the same
       remap routes.py already does at export; doing it here makes the gpkg
       internally consistent — see commit 4c2408e, "foreign_key_violation".)

    B. numeric-key canonicalization
       Integer FK columns (trips_intervals.trip_id / interval_id,
       agencies.vehicle_id) can end up with REAL/float affinity — typically when
       a null forced the column to float64 — so a gid "1" is stored as "1.0" and
       every downstream string join misses. We coerce those columns back to
       INTEGER affinity, turning "1.0" into 1.

Genuine renames ("AMU" -> "AMUGA") and deletes have no history to reconstruct,
so they are left to the block-and-warn validator. This module never guesses.

Behaviour
---------
``plan_key_repairs`` is pure (DataFrames in, plan out) and unit-testable.
``apply_key_repairs`` writes the plan back to the GeoPackage. The Refresh tool
runs plan+report on every refresh and only calls apply when the user opts in.

pandas + sqlite3 only. No QGIS / shapely.
"""

from __future__ import annotations

import re
import sqlite3
from dataclasses import dataclass, field

import pandas as pd

# Layers whose agency_id may hold a gid instead of the text code.
_AGENCY_ID_LAYERS = ("transit_trips_view", "transit_trips")
# Non-spatial layers + their integer-key columns to coerce to INTEGER affinity.
_NUMERIC_KEY_COLS = {
    "transit_trips_intervals": ("trip_id", "interval_id"),
    "transit_agencies": ("vehicle_id",),
}

_FLOATINT_RE = re.compile(r"^\d+\.0+$")


@dataclass
class RepairPlan:
    # (layer, gid, old_value, new_value)
    agency_remap: list = field(default_factory=list)
    # (layer, column, n_cells)
    numeric_coerce: list = field(default_factory=list)

    @property
    def is_empty(self):
        return not self.agency_remap and not self.numeric_coerce

    def to_lines(self):
        lines = []
        by_layer = {}
        for layer, gid, old, new in self.agency_remap:
            by_layer.setdefault(layer, []).append((gid, old, new))
        for layer, items in by_layer.items():
            sample = [f"{o}->{n}" for _, o, n in items[:5]]
            lines.append(f"agency_id gid->code: {layer} — {len(items)} row(s), e.g. {sample}")
        for layer, col, n in self.numeric_coerce:
            lines.append(f"numeric key canonicalization: {layer}.{col} — {n} cell(s) '1.0'->1")
        return lines


# --------------------------------------------------------------------------- #
# Helpers
# --------------------------------------------------------------------------- #
def _canon(v):
    """String form, stripped, with a float-int like '1.0' collapsed to '1'."""
    if v is None or (isinstance(v, float) and pd.isna(v)):
        return ""
    s = str(v).strip()
    if _FLOATINT_RE.match(s):
        s = s.split(".", 1)[0]
    return s


def _to_int(v):
    try:
        return int(float(v))
    except (TypeError, ValueError):
        return None


def _agency_maps(agencies):
    """Return (gid_to_code, valid_codes) from the agencies table."""
    if agencies is None or "agency_id" not in agencies.columns or "gid" not in agencies.columns:
        return {}, set()
    gid_to_code, codes = {}, set()
    for gid, code in zip(agencies["gid"], agencies["agency_id"]):
        c = _canon(code)
        g = _canon(gid)
        if c:
            codes.add(c)
        if g and c:
            gid_to_code[g] = c
    return gid_to_code, codes


def _count_floatint(series):
    return int(series.map(lambda v: bool(_FLOATINT_RE.match(str(v).strip()))).sum())


# --------------------------------------------------------------------------- #
# Planning (pure)
# --------------------------------------------------------------------------- #
def plan_key_repairs(tables: dict) -> RepairPlan:
    """Detect the safe repairs. ``tables`` maps physical layer name -> DataFrame."""
    plan = RepairPlan()

    gid_to_code, valid_codes = _agency_maps(tables.get("transit_agencies"))
    if gid_to_code:
        for layer in _AGENCY_ID_LAYERS:
            df = tables.get(layer)
            if df is None or "agency_id" not in df.columns or "gid" not in df.columns:
                continue
            for gid, val in zip(df["gid"], df["agency_id"]):
                s = _canon(val)
                # Only remap a value that is NOT already a valid code but IS a
                # known gid. If codes are themselves numeric and collide with
                # gids, the "already a valid code" guard wins and we skip.
                if s and s not in valid_codes and s in gid_to_code:
                    plan.agency_remap.append((layer, _to_int(gid), s, gid_to_code[s]))

    for layer, cols in _NUMERIC_KEY_COLS.items():
        df = tables.get(layer)
        if df is None:
            continue
        for c in cols:
            if c in df.columns:
                n = _count_floatint(df[c])
                if n:
                    plan.numeric_coerce.append((layer, c, n))

    return plan


# --------------------------------------------------------------------------- #
# Apply (mutates the GeoPackage)
# --------------------------------------------------------------------------- #
def _sqlite_type(series):
    if pd.api.types.is_integer_dtype(series):
        return "INTEGER"
    if pd.api.types.is_float_dtype(series):
        return "REAL"
    if pd.api.types.is_bool_dtype(series):
        return "INTEGER"
    return "TEXT"


def _table_exists(conn, name) -> bool:
    row = conn.execute(
        "SELECT 1 FROM sqlite_master WHERE type IN ('table','view') AND name=?", (name,)
    ).fetchone()
    return row is not None


def _rewrite_nonspatial(conn, table, coerce_cols):
    """Recreate a NON-SPATIAL table, coercing coerce_cols to INTEGER affinity."""
    df = pd.read_sql_query(f'SELECT * FROM "{table}"', conn)  # nosec B608 - fixed name
    if "geom" in df.columns or "geometry" in df.columns:
        raise ValueError(f"refusing to rewrite spatial table {table} in place")
    for c in coerce_cols:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce").astype("Int64")
    if "gid" in df.columns:
        df["gid"] = pd.to_numeric(df["gid"], errors="coerce").astype("Int64")

    cols_sql = []
    for col in df.columns:
        if col == "gid":
            cols_sql.append('"gid" INTEGER PRIMARY KEY')
        elif col in coerce_cols:
            cols_sql.append(f'"{col}" INTEGER')
        else:
            cols_sql.append(f'"{col}" {_sqlite_type(df[col])}')
    conn.execute(f'DROP TABLE IF EXISTS "{table}"')
    conn.execute(f'CREATE TABLE "{table}" ({", ".join(cols_sql)})')
    df.to_sql(name=table, con=conn, if_exists="append", index=False)


def apply_key_repairs(gpkg_path: str, plan: RepairPlan) -> dict:
    """Write the plan back to the GeoPackage. Returns a summary dict of counts."""
    summary = {"agency_remapped": 0, "columns_coerced": 0}
    if plan.is_empty:
        return summary

    with sqlite3.connect(gpkg_path) as conn:
        # A) agency_id gid -> code, per-row UPDATE keyed by gid (preserves
        #    geometry / gpkg registration; agency_id is a TEXT column so the
        #    written code sticks).
        by_layer = {}
        for layer, gid, _old, new in plan.agency_remap:
            if gid is not None:
                by_layer.setdefault(layer, []).append((new, gid))
        for layer, params in by_layer.items():
            if not _table_exists(conn, layer):
                continue
            conn.executemany(
                f'UPDATE "{layer}" SET "agency_id"=? WHERE "gid"=?', params  # nosec B608
            )
            summary["agency_remapped"] += len(params)

        # B) numeric coercion — recreate the non-spatial table with INTEGER
        #    affinity so '1.0' can never come back.
        coerce_by_layer = {}
        for layer, col, _n in plan.numeric_coerce:
            coerce_by_layer.setdefault(layer, set()).add(col)
        for layer, cols in coerce_by_layer.items():
            if not _table_exists(conn, layer):
                continue
            _rewrite_nonspatial(conn, layer, cols)
            summary["columns_coerced"] += len(cols)

        conn.commit()
    return summary
