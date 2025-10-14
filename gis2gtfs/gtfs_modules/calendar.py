'''
This script:
1. Creates a static calendar table manually in memory using tibble()
2. Specifies two service types:
"Ground_Daily": operates every day
"Ground_Weekdays": operates only Monday to Thursday and Sunday (not Friday/Saturday)
3. Includes GTFS-required fields:
monday through sunday (binary flags)
start_date, end_date (from global variables)
service_id
4. Exports it to calendar.txt in GTFS format (CSV with no index)
'''
import os
import pandas as pd


def generate(data_dir, start_date, end_date, service_id):
    """
    Generate GTFS calendar.txt file.

    Parameters:
        data_dir (str): Directory to save the output GTFS file
        start_date (int): Start date in YYYYMMDD format
        end_date (int): End date in YYYYMMDD format
        service_id (str): Base service ID to be used (e.g. "Ground")
    """
    output_file = os.path.join(data_dir, "calendar.txt")

    print("🗓️ Generating calendar.txt...")

    calendar_df = pd.DataFrame([
        {
            "service_id": "Ground_Daily",
            "monday": 1, "tuesday": 1, "wednesday": 1, "thursday": 1,
            "friday": 1, "saturday": 1, "sunday": 1,
            "start_date": start_date,
            "end_date": end_date
        },
        {
            "service_id": "Ground_Weekdays",
            "monday": 1, "tuesday": 1, "wednesday": 1, "thursday": 1,
            "friday": 0, "saturday": 0, "sunday": 1,
            "start_date": start_date,
            "end_date": end_date
        }
    ])

    calendar_df.to_csv(output_file, index=False)
    print("✅ calendar.txt written.")
