'''
This script:
1. Creates a static GTFS feed_info table containing metadata about the feed.
2. Pulls in the following variables from the global config:
feed_start_date, feed_end_date — date range of the feed
feed_name — used as feed_version
3. Hardcodes:
feed_publisher_name = "Transport for Cairo"
feed_publisher_url = "http://transportforcairo.com"
feed_lang = "en"
4. Writes the output as a comma-delimited file feed_info.txt, per GTFS specification
'''
import os
import pandas as pd


def generate(data_dir, feed_name, feed_start_date, feed_end_date, feed_version=None, feed_lang="en"):
    """
    Generate GTFS feed_info.txt file.

    Parameters:
        data_dir (str): Directory to save the GTFS output file
        feed_name (str): Used as default feed_version if not separately specified
        feed_start_date (int): Feed validity start date in YYYYMMDD
        feed_end_date (int): Feed validity end date in YYYYMMDD
        feed_version (str, optional): GTFS feed version (defaults to feed_name)
        feed_lang (str): Language of the feed (default is 'en')
    """
    output_file = os.path.join(data_dir, "feed_info.txt")
    version = feed_version or feed_name

    print("🧾 Generating feed_info.txt...")

    feed_info_df = pd.DataFrame([{
        "feed_publisher_name": "Transport for Cairo",
        "feed_publisher_url": "http://transportforcairo.com",
        "feed_start_date": feed_start_date,
        "feed_end_date": feed_end_date,
        "feed_version": version,
        "feed_lang": feed_lang
    }])

    feed_info_df.to_csv(output_file, index=False)
    print("✅ feed_info.txt written.")
