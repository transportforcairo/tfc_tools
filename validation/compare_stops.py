#!/usr/bin/env python3
"""
compare_stops.py — validation harness for TfC automatic stop creation.

Scores a *candidate* stops layer (produced by the automatic pipeline) against a
*reference* stops layer (a human-curated ground truth), on the dimensions that
matter for the GIS2GTFS redesign:

    * cardinality & physical-place count
    * directional pairing (are the two directions of a stop created as a pair?)
    * spacing (nearest-neighbour distribution)
    * coverage (does every reference stop have a candidate nearby, and vice-versa?)
    * naming quality (blank/"Unnamed" rate, stray whitespace, pair name-sharing)

It can also run in --parity mode to compare two *candidate* layers against each
other (e.g. the PostGIS SQL output vs the offline-GeoPackage geopandas mirror),
to prove the two implementations agree.

Design notes
------------
Deliberately depends on the Python standard library ONLY (sqlite3, json,
struct, math, argparse). No geopandas / shapely / fiona. This keeps the harness
runnable in any environment — including a bare RL-suite CI runner — and means it
never fails to install. A GeoPackage is just a SQLite DB; its geometry blobs are
parsed inline (points only, which is all a stops layer contains).

Inputs may be:
    * a GeoPackage file  ->  path.gpkg[:layer_name]   (default layer transit_stops)
    * a GeoJSON file     ->  path.geojson

Usage
-----
    # score auto vs human
    python compare_stops.py --candidate export.gpkg:transit_stops \
                            --reference human_stops.geojson

    # SQL-vs-geopandas parity check (no ground truth needed)
    python compare_stops.py --candidate sql_out.gpkg:transit_stops \
                            --reference gpkg_out.gpkg:transit_stops --parity

    # emit machine-readable JSON as well as the text scorecard
    python compare_stops.py -c a.gpkg -r human.geojson --json report.json
"""

from __future__ import annotations

import argparse
import json
import math
import re
import struct
import sqlite3
import sys
import unicodedata
from collections import Counter, defaultdict


# ---------------------------------------------------------------------------
# Geometry helpers (no external deps)
# ---------------------------------------------------------------------------

def _parse_gpkg_point(blob):
    """Extract (lon, lat) from a GeoPackage geometry blob (point geometries)."""
    if blob is None:
        return None
    if blob[:2] != b"GP":
        return None
    flags = blob[3]
    env = (flags >> 1) & 0x07
    env_bytes = {0: 0, 1: 32, 2: 48, 3: 48, 4: 64}.get(env, 0)
    wkb = blob[8 + env_bytes:]
    endian = "<" if wkb[0] == 1 else ">"
    gtype = struct.unpack(endian + "I", wkb[1:5])[0]
    if gtype % 1000 != 1:          # not a Point
        return None
    x, y = struct.unpack(endian + "dd", wkb[5:21])
    return (x, y)


def haversine_m(a, b):
    """Great-circle distance in metres between two (lon, lat) tuples."""
    (lo1, la1), (lo2, la2) = a, b
    r = 6371000.0
    p = math.pi / 180.0
    dla = (la2 - la1) * p
    dlo = (lo2 - lo1) * p
    h = (math.sin(dla / 2) ** 2
         + math.cos(la1 * p) * math.cos(la2 * p) * math.sin(dlo / 2) ** 2)
    return 2 * r * math.asin(min(1.0, math.sqrt(h)))


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------

# Candidate column names we understand, in priority order.
_NAME_KEYS = ("stop_name", "name", "raw_name")
_DIR_KEYS = ("direction_id", "double", "dir", "dir_bin")


_SAFE_IDENT = re.compile(r"^[A-Za-z_][A-Za-z0-9_]*$")


def _safe_ident(name):
    """Reject anything that is not a bare SQL identifier.

    SQLite cannot bind table/column names, so they must be interpolated into the
    query text. The layer name comes from a user-supplied CLI argument
    (path.gpkg:layer), so it is validated here before it can reach the SQL.
    Column names come from PRAGMA table_info but are validated too, so a
    maliciously crafted GeoPackage cannot inject either.
    """
    if not isinstance(name, str) or not _SAFE_IDENT.match(name):
        raise SystemExit(f"Unsafe SQL identifier rejected: {name!r}")
    return name


def _load_gpkg(path, layer):
    con = sqlite3.connect(path)
    try:
        if layer is None:
            layer = "transit_stops"
        layer = _safe_ident(layer)
        cols = [r[1] for r in con.execute(f'PRAGMA table_info("{layer}")')]  # nosec B608 — identifier validated by _safe_ident
        if not cols:
            raise SystemExit(f"Layer '{layer}' not found in {path}")
        gcol = "geom" if "geom" in cols else ("geometry" if "geometry" in cols else None)
        if gcol is None:
            raise SystemExit(f"No geometry column in layer '{layer}'")
        attrs = [_safe_ident(c) for c in cols if c != gcol]
        gcol = _safe_ident(gcol)
        rows = []
        # nosec B608 — every identifier below passed _safe_ident (bare
        # [A-Za-z_][A-Za-z0-9_]* only); no user value is interpolated.
        q = f'SELECT {",".join(chr(34)+c+chr(34) for c in attrs)},"{gcol}" FROM "{layer}"'  # nosec B608
        for r in con.execute(q):
            d = dict(zip(attrs, r[:-1]))
            d["_xy"] = _parse_gpkg_point(r[-1])
            rows.append(d)
        return rows
    finally:
        con.close()


def _load_geojson(path):
    with open(path, "r", encoding="utf-8") as f:
        gj = json.load(f)
    rows = []
    for feat in gj.get("features", []):
        props = dict(feat.get("properties") or {})
        geom = feat.get("geometry")
        xy = None
        if geom and geom.get("type") == "Point":
            c = geom.get("coordinates") or []
            if len(c) >= 2:
                xy = (c[0], c[1])
        props["_xy"] = xy
        rows.append(props)
    return rows


def load_layer(spec):
    """spec = 'path', 'path.gpkg' or 'path.gpkg:layer'. Returns rows with _xy.

    A trailing ':layer' is only honoured for GeoPackage inputs, and only when
    the part before the last ':' ends in '.gpkg' — this avoids mistaking a
    Windows drive letter (C:\\...) for a layer separator.
    """
    low = spec.lower()
    if ".gpkg" in low:
        base, sep, lyr = spec.rpartition(":")
        if sep and base.lower().endswith(".gpkg"):
            return _load_gpkg(base, lyr or None)
        return _load_gpkg(spec, None)
    if low.endswith((".geojson", ".json")):
        return _load_geojson(spec)
    raise SystemExit(f"Unrecognised input (need .gpkg or .geojson): {spec}")


# ---------------------------------------------------------------------------
# Metric helpers
# ---------------------------------------------------------------------------

def _first_key(row, keys):
    for k in keys:
        if k in row:
            return k
    return None


def _name_of(row, name_key):
    return (row.get(name_key) or "").strip() if name_key else ""


def _norm_name(s):
    s = (s or "")
    s = unicodedata.normalize("NFC", s)
    s = " ".join(s.split())            # collapse + trim whitespace
    return s.lower()


def percentiles(values, ps):
    if not values:
        return {p: None for p in ps}
    v = sorted(values)
    n = len(v)
    return {p: round(v[min(n - 1, int(p / 100.0 * n))], 1) for p in ps}


def nearest_within(rows):
    """List of each point's nearest-neighbour distance (m) within the layer."""
    pts = [r["_xy"] for r in rows if r["_xy"]]
    out = []
    for i, a in enumerate(pts):
        best = math.inf
        for j, b in enumerate(pts):
            if i == j:
                continue
            d = haversine_m(a, b)
            if d < best:
                best = d
        if best < math.inf:
            out.append(best)
    return out


def cross_nearest(a_rows, b_rows):
    """For each A point, distance to nearest B point (m)."""
    bpts = [b["_xy"] for b in b_rows if b["_xy"]]
    out = []
    for a in a_rows:
        if not a["_xy"] or not bpts:
            continue
        out.append(min(haversine_m(a["_xy"], bp) for bp in bpts))
    return out


def collapse_places(rows, dist_m):
    """Union-find collapse of points within dist_m into physical places."""
    pts = [r for r in rows if r["_xy"]]
    parent = list(range(len(pts)))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i in range(len(pts)):
        for j in range(i + 1, len(pts)):
            if haversine_m(pts[i]["_xy"], pts[j]["_xy"]) <= dist_m:
                parent[find(i)] = find(j)
    groups = defaultdict(list)
    for i in range(len(pts)):
        groups[find(i)].append(pts[i])
    return list(groups.values())


def pairing_stats(rows, pair_max_m, name_key):
    """% of stops that have a neighbour within pair_max_m and share its name."""
    pts = [(r["_xy"], _norm_name(_name_of(r, name_key))) for r in rows if r["_xy"]]
    have = same = diff = 0
    seps = []
    for i, (xy, nm) in enumerate(pts):
        best = math.inf
        bnm = None
        for j, (xy2, nm2) in enumerate(pts):
            if i == j:
                continue
            d = haversine_m(xy, xy2)
            if d < best:
                best, bnm = d, nm2
        if best <= pair_max_m:
            have += 1
            seps.append(best)
            if bnm == nm:
                same += 1
            else:
                diff += 1
    n = len(pts) or 1
    return {
        "n_stops": len(pts),
        "with_neighbor": have,
        "with_neighbor_pct": round(100.0 * have / n, 1),
        "neighbor_same_name": same,
        "neighbor_diff_name": diff,
        "same_name_frac_of_paired": round(same / have, 3) if have else None,
        "median_sep_m": percentiles(seps, [50])[50],
    }


def naming_stats(rows, name_key):
    names = [(r.get(name_key) if name_key else None) or "" for r in rows]
    blank = sum(1 for x in names if not x.strip() or x.strip().lower() == "unnamed")
    whitespace = sum(1 for x in names if x != x.strip() or "  " in x)
    nonblank = [x.strip() for x in names if x.strip() and x.strip().lower() != "unnamed"]
    norm_counts = Counter(_norm_name(x) for x in nonblank)
    dup_groups = sum(1 for v in norm_counts.values() if v > 1)
    return {
        "name_key": name_key,
        "total": len(names),
        "blank_or_unnamed": blank,
        "with_whitespace_issues": whitespace,
        "distinct_names": len(norm_counts),
        "duplicate_name_groups": dup_groups,
    }


def direction_stats(rows):
    dk = _first_key(rows[0], _DIR_KEYS) if rows else None
    if not dk:
        return {"direction_key": None}
    return {"direction_key": dk, "distribution": dict(Counter(str(r.get(dk)) for r in rows))}


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------

def build_report(candidate, reference, args):
    cand_name_key = _first_key(candidate[0], _NAME_KEYS) if candidate else None
    ref_name_key = _first_key(reference[0], _NAME_KEYS) if reference else None

    cand_places = collapse_places(candidate, args.place_m)
    ref_places = collapse_places(reference, args.place_m)

    r2c = cross_nearest(reference, candidate)
    c2r = cross_nearest(candidate, reference)

    rep = {
        "params": {
            "place_m": args.place_m,
            "pair_max_m": args.pair_m,
            "coverage_tight_m": args.coverage_m,
        },
        "counts": {
            "candidate_stops": len([r for r in candidate if r["_xy"]]),
            "reference_stops": len([r for r in reference if r["_xy"]]),
            "candidate_places": len(cand_places),
            "reference_places": len(ref_places),
        },
        "spacing_nn_m": {
            "candidate": percentiles(nearest_within(candidate), [5, 25, 50, 75, 95]),
            "reference": percentiles(nearest_within(reference), [5, 25, 50, 75, 95]),
        },
        "pairing": {
            "candidate": pairing_stats(candidate, args.pair_m, cand_name_key),
            "reference": pairing_stats(reference, args.pair_m, ref_name_key),
        },
        "coverage": {
            "reference_to_candidate": {
                "pctiles_m": percentiles(r2c, [50, 75, 90, 95]),
                "within_tight": sum(1 for d in r2c if d <= args.coverage_m),
                "of": len(r2c),
                "pct": round(100.0 * sum(1 for d in r2c if d <= args.coverage_m) / (len(r2c) or 1), 1),
            },
            "candidate_to_reference": {
                "pctiles_m": percentiles(c2r, [50, 75, 90, 95]),
                "within_tight": sum(1 for d in c2r if d <= args.coverage_m),
                "of": len(c2r),
                "pct": round(100.0 * sum(1 for d in c2r if d <= args.coverage_m) / (len(c2r) or 1), 1),
            },
        },
        "naming": {
            "candidate": naming_stats(candidate, cand_name_key),
            "reference": naming_stats(reference, ref_name_key),
        },
        "direction": {
            "candidate": direction_stats(candidate),
            "reference": direction_stats(reference),
        },
    }
    return rep


def _fmt_pct(d):
    return "{" + ", ".join(f"p{k}={v}" for k, v in d.items()) + "}"


def print_scorecard(rep, args):
    c = rep["counts"]
    role = "PARITY: layer B" if args.parity else "REFERENCE (human)"
    print("=" * 72)
    print("  TfC stop-creation validation scorecard")
    print(f"  place_m={rep['params']['place_m']}  pair_max_m={rep['params']['pair_max_m']}"
          f"  coverage_m={rep['params']['coverage_tight_m']}")
    print("=" * 72)

    print("\nCOUNTS                          candidate      " + role)
    print(f"  stops                         {c['candidate_stops']:>9}   {c['reference_stops']:>15}")
    print(f"  physical places (≤{rep['params']['place_m']:.0f}m)        {c['candidate_places']:>9}   {c['reference_places']:>15}")

    print("\nSPACING  nearest-neighbour (m)")
    print(f"  candidate  {_fmt_pct(rep['spacing_nn_m']['candidate'])}")
    print(f"  reference  {_fmt_pct(rep['spacing_nn_m']['reference'])}")

    pc, pr = rep["pairing"]["candidate"], rep["pairing"]["reference"]
    print(f"\nPAIRING (neighbour within {rep['params']['pair_max_m']:.0f}m)")
    print(f"  candidate : {pc['with_neighbor']}/{pc['n_stops']} = {pc['with_neighbor_pct']}%"
          f"  same-name {pc['same_name_frac_of_paired']}  median sep {pc['median_sep_m']}m")
    print(f"  reference : {pr['with_neighbor']}/{pr['n_stops']} = {pr['with_neighbor_pct']}%"
          f"  same-name {pr['same_name_frac_of_paired']}  median sep {pr['median_sep_m']}m")

    cov = rep["coverage"]
    cm = rep["params"]["coverage_tight_m"]
    print(f"\nCOVERAGE (within {cm:.0f}m)")
    r2c = cov["reference_to_candidate"]
    c2r = cov["candidate_to_reference"]
    print(f"  reference → candidate : {r2c['pct']}% ({r2c['within_tight']}/{r2c['of']})"
          f"   dist {_fmt_pct(r2c['pctiles_m'])}")
    print(f"  candidate → reference : {c2r['pct']}% ({c2r['within_tight']}/{c2r['of']})"
          f"   dist {_fmt_pct(c2r['pctiles_m'])}")

    nc, nr = rep["naming"]["candidate"], rep["naming"]["reference"]
    print("\nNAMING                          candidate      " + role)
    print(f"  name column                   {str(nc['name_key']):>9}   {str(nr['name_key']):>15}")
    print(f"  blank / 'Unnamed'             {nc['blank_or_unnamed']:>9}   {nr['blank_or_unnamed']:>15}")
    print(f"  whitespace issues             {nc['with_whitespace_issues']:>9}   {nr['with_whitespace_issues']:>15}")
    print(f"  distinct names                {nc['distinct_names']:>9}   {nr['distinct_names']:>15}")

    dc, dr = rep["direction"]["candidate"], rep["direction"]["reference"]
    print("\nDIRECTION")
    print(f"  candidate key={dc.get('direction_key')}  {dc.get('distribution','')}")
    print(f"  reference key={dr.get('direction_key')}  {dr.get('distribution','')}")
    print("=" * 72)


def main(argv=None):
    ap = argparse.ArgumentParser(description="Compare candidate vs reference stops layers.")
    ap.add_argument("-c", "--candidate", required=True,
                    help="candidate layer: path.gpkg[:layer] or path.geojson")
    ap.add_argument("-r", "--reference", required=True,
                    help="reference layer: path.gpkg[:layer] or path.geojson")
    ap.add_argument("--place-m", dest="place_m", type=float, default=40.0,
                    help="collapse radius for counting physical places (default 40)")
    ap.add_argument("--pair-m", dest="pair_m", type=float, default=40.0,
                    help="max separation to count a directional pair (default 40)")
    ap.add_argument("--coverage-m", dest="coverage_m", type=float, default=30.0,
                    help="tight coverage radius (default 30)")
    ap.add_argument("--parity", action="store_true",
                    help="both inputs are candidate implementations (SQL vs geopandas)")
    ap.add_argument("--json", dest="json_out", default=None,
                    help="also write machine-readable JSON to this path")
    args = ap.parse_args(argv)

    candidate = load_layer(args.candidate)
    reference = load_layer(args.reference)
    if not candidate:
        raise SystemExit("Candidate layer has no rows.")
    if not reference:
        raise SystemExit("Reference layer has no rows.")

    rep = build_report(candidate, reference, args)
    print_scorecard(rep, args)

    if args.json_out:
        with open(args.json_out, "w", encoding="utf-8") as f:
            json.dump(rep, f, indent=2, ensure_ascii=False)
        print(f"\nJSON written to {args.json_out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
