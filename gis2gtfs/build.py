import os
from qgis.core import QgsProviderRegistry

# Import GTFS module functions
from .gtfs_modules.download_db_data import download_db_data
from .gtfs_modules.estimate_missing_times import run_estimate_missing_times
from .gtfs_modules.feed_info import generate as generate_feed_info
from .gtfs_modules.agency import generate as generate_agency
from .gtfs_modules.calendar import generate as generate_calendar
from .gtfs_modules.frequencies import generate as generate_frequencies
from .gtfs_modules.routes import generate as generate_routes
from .gtfs_modules.shapes import generate as generate_shapes
from .gtfs_modules.stops import generate as generate_stops
from .gtfs_modules.trips import generate as generate_trips
from .gtfs_modules.stop_times import generate as generate_stop_times


def run_gtfs_pipeline(
    sdi_mode,
    conn_name,
    gpkg_path,
    data_raw_dir,
    data_dir,
    feed_name,
    feed_start_date,
    feed_end_date,
    feed_version,
    feed_lang,
    start_date,
    end_date,
    service_id,
    continuous_dropoff_pickup=True,
    dwell_time_sec=15,
    db_host="sdi.data.transportforcairo.com",
    db_port=5432,
):
    """
    Main GTFS build pipeline using connection name from QGIS UI.
    Converts SDI GIS data into GTFS-compliant text files.

    Args:
        data_raw_dir (str): Folder for raw extracted data
        data_dir (str): Output folder for GTFS .txt files
        feed_name (str): GTFS feed version name (also used as feed_version if not specified)
        feed_start_date (int): YYYYMMDD
        feed_end_date (int): YYYYMMDD
        feed_version (str): GTFS feed version (shown to user apps)
        feed_lang (str): Feed language (e.g., 'en')
        start_date (int): Calendar start date (YYYYMMDD)
        end_date (int): Calendar end date (YYYYMMDD)
        service_id (str): Base GTFS service ID (e.g., 'Ground')
        continuous_dropoff_pickup (bool): If True, includes informal boarding flags
        dwell_time_sec (int): Dwell time per stop in seconds
        db_host (str): PostgreSQL host
        db_port (int): PostgreSQL port (default 5432)
    """
    print("🚀 Starting GTFS generation pipeline...")

    # --- Step 0: Download DB data
    print("📥 Downloading SDI data...", conn_name if sdi_mode == 'postgres' else gpkg_path)
    download_db_data(
        sdi_mode=sdi_mode,
        conn_name=conn_name,
        gpkg_path=gpkg_path,
        output_dir=data_raw_dir
    )
    print("✅ Data downloaded")

    # --- Step 1: Estimate missing times
    print("⏱ Estimating missing travel times...")
    run_estimate_missing_times(
        input_csv=os.path.join(data_raw_dir, "travel_times_trackpoints.csv"),
        output_csv=os.path.join(data_raw_dir, "travel_times_trackpoints_filled_na.csv"),
        vehicle=1  # assuming vehicle-type specific data
    )
    print("✅ Travel times calculated")

    # --- GTFS Generation Scripts
    print("🧱 Generating GTFS components...")

    # commented part seems clearer but I'm not sure about generate()
    generate_feed_info(data_dir, feed_name, feed_start_date, feed_end_date, feed_version, feed_lang)
    print("✅ feed_info.txt DONE")

    generate_agency(data_dir, data_raw_dir)
    print("✅ agency.txt DONE")

    # calendar.py now derives services dynamically from intervals.mon..sun
    # (see the helper in gtfs_modules/_interval_expansion.py). The
    # service_id parameter acts as a legacy fallback -- it's only used
    # when the raw intervals.csv has no day columns (pre-refresh
    # exporter). In the new pipeline the generated calendar.txt has one
    # row per distinct day pattern (svc_daily, svc_weekday, etc.).
    generate_calendar(data_dir, data_raw_dir, start_date, end_date,
                      fallback_service_id=f"{service_id}_Daily")
    print("✅ calendar.txt DONE")

    generate_frequencies(data_dir, data_raw_dir)
    print("✅ frequencies.txt DONE")

    generate_routes(data_dir, data_raw_dir, continuous_dropoff_pickup=continuous_dropoff_pickup)
    print("✅ routes.txt DONE")

    generate_shapes(data_dir, data_raw_dir)
    print("✅ shapes.txt DONE")

    generate_stops(data_dir, data_raw_dir)
    print("✅ stops.txt DONE")

    generate_trips(data_dir, data_raw_dir)
    print("✅ trips.txt DONE")

    generate_stop_times(data_dir, data_raw_dir, dwell_time_sec=dwell_time_sec)
    print("✅ stop_times.txt DONE")

    print("🎉 GTFS pipeline completed successfully!")