'''
This is a simple script that:
1. Loads the previously extracted agency.csv (from the database)
2. Selects only the required GTFS fields:
agency_id
agency_name
agency_url
agency_timezone
3. Exports the result as a comma-delimited GTFS file: agency.txt
'''

import os
import pandas as pd


def generate(data_dir, data_raw_dir):
    """
    Create GTFS agency.txt file from raw agency.csv.

    Parameters:
        data_dir (str): Output directory for GTFS files
        data_raw_dir (str): Directory where raw CSV data is stored
    """
    input_file = os.path.join(data_raw_dir, "agency.csv")
    output_file = os.path.join(data_dir, "agency.txt")

    print("📥 Reading agency.csv...")
    agency_df = pd.read_csv(input_file)

    # # 🚨 TEMP PATCH: Fill missing agency_name
    # # THIS PART IS FOR TESTING PURPOSES 
    # # This block should be removed once SDI data is manually corrected
    # if "agency_url" in agency_df.columns:
    #     agency_df["agency_url"] = agency_df["agency_url"].fillna("https://transportforcairo.com/")
    # if "agency_timezone" in agency_df.columns:
    #     agency_df["agency_timezone"] = agency_df["agency_timezone"].fillna("Africa/Cairo")
    # # 🚨 END OF TEMP PATCH    

    print("📤 Writing agency.txt...")
    agency_df = agency_df[["agency_id", "agency_name", "agency_url", "agency_timezone"]]
    agency_df.to_csv(output_file, index=False)

    print("✅ agency.txt written.")
