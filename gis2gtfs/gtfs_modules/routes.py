'''
Main Steps:
1. Load trips data (trips.geojson):
Contains route-level attributes (e.g., route_id, agency_id, route_type, etc.)
2. Deduplicate:
Retain only one record per route_id (remove duplicates caused by origin/destination trip direction)
3. Remap agency_id: some SDI exports store transit_trips.agency_id as the agency
   surrogate gid (e.g., 1) rather than the text code (e.g., 'Daladala'). GTFS
   requires routes.agency_id to match agency.agency_id, so we look up the code
   from agency.csv (gid -> agency_id) and overwrite the column when the raw
   value is numeric.
4. Add route_short_name:
Rename route_short to route_short_name
5. Include conditional logic:
If continuous_dropoff_pickup == TRUE, add continuous_pickup and continuous_drop_off columns.
Otherwise, exclude them.
6. Export routes.txt in GTFS format.
'''
import os
import pandas as pd
import geopandas as gpd


def _build_agency_gid_to_code(data_raw_dir):
    """Return a dict mapping agency gid (as str) -> agency_id code.

    Returns an empty dict if agency.csv is missing or lacks the required cols.
    """
    agency_csv = os.path.join(data_raw_dir, "agency.csv")
    if not os.path.exists(agency_csv):
        return {}
    ag = pd.read_csv(agency_csv, encoding="utf-8")
    if "agency_id" not in ag.columns or "gid" not in ag.columns:
        return {}
    return {str(g): str(a) for g, a in zip(ag["gid"], ag["agency_id"]) if pd.notna(g) and pd.notna(a)}


def generate(data_dir, data_raw_dir, continuous_dropoff_pickup=True):
    """
    Generate GTFS routes.txt file.

    Parameters:
        data_dir (str): Output GTFS folder
        data_raw_dir (str): Folder containing raw trips.geojson
        continuous_dropoff_pickup (bool): Whether to include continuous pickup/dropoff columns
    """
    print("Loading trip data from trips.geojson...")

    trips_path = os.path.join(data_raw_dir, "trips.geojson")
    trips_gdf = gpd.read_file(trips_path)

    # Remove duplicate route_ids (keep first occurrence)
    routes_df = trips_gdf.drop_duplicates(subset=["route_id"]).copy()

    # Drop geometry
    routes_df = routes_df.drop(columns="geometry", errors="ignore")

    # Remap agency_id (gid) -> agency code if needed. Some exports store the
    # surrogate gid; agency.txt always writes the text code, so unless we
    # translate we get a foreign_key_violation in the validator.
    gid_to_code = _build_agency_gid_to_code(data_raw_dir)
    if gid_to_code and "agency_id" in routes_df.columns:
        raw = routes_df["agency_id"].astype(str)
        routes_df["agency_id"] = raw.map(lambda v: gid_to_code.get(v, v))

    # Rename route_short to route_short_name
    routes_df.rename(columns={"route_short": "route_short_name"}, inplace=True)

    # Select columns
    if continuous_dropoff_pickup:
        required_columns = [
            "route_id", "agency_id", "route_long", "route_short_name",
            "route_type", "continuous_pickup", "continuous_drop_off"
        ]
    else:
        required_columns = [
            "route_id", "agency_id", "route_long", "route_short_name", "route_type"
        ]

    routes_df = routes_df[required_columns].rename(columns={"route_long": "route_long_name"})

    output_file = os.path.join(data_dir, "routes.txt")
    routes_df.to_csv(output_file, index=False, encoding='utf-8')
    print("routes.txt written.")
