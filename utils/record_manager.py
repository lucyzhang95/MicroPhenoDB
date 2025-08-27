"""Record Manager for MicroPhenoDB Parser
Manages caching and processing of the final parsed records.
"""
import os
import sys

from MicroPhenoDB_parser import MicroPhenoDBParser
from utils.cache_helper import CacheHelper


class RecordCacheManager:
    """Manages caching of the final processed records."""

    def __init__(self, cache_helper: CacheHelper):
        self.cache_helper = cache_helper

    def cache_records(self, records: list, filename_base="microphenodb_parsed_records"):
        """
        Caches a given list of records to both pickle and JSON formats.
        This method is responsible for SAVING, not generating.
        """
        if not records:
            print("‼️ Warning: No records provided to cache.")
            return

        print(f"\n▶️ Caching {len(records)} records...")
        pickle_filename = f"{filename_base}.pkl"
        json_filename = f"{filename_base}.json"

        self.cache_helper.save_pickle(records, pickle_filename)
        self.cache_helper.save_json(records, json_filename)
        print("🎉 Full MicroPhenoDB parsed records cached successfully.")

    def load_cached_records(self, filename_base="microphenodb_parsed_records"):
        """
        Loads cached records from pickle format.
        """
        pickle_filename = f"{filename_base}.pkl"
        try:
            records = self.cache_helper.load_pickle(pickle_filename)
            print(f"✅ Loaded {len(records)} cached records from {pickle_filename}")
            return records
        except FileNotFoundError:
            print(f"❌ No cached records found at {pickle_filename}")
            return None

    @staticmethod
    def run_data_pipeline():
        """
        Orchestrates the entire process of parsing and caching data.
        """

        current_dir = os.path.dirname(os.path.abspath(__file__))
        parent_dir = os.path.dirname(current_dir)
        if parent_dir not in sys.path:
            sys.path.insert(0, parent_dir)

        cache_helper = CacheHelper(cache_dir="cache")
        parser = MicroPhenoDBParser(data_dir="downloads")
        record_cacher = RecordCacheManager(cache_helper)

        records_list = list(parser.load_microphenodb_data())
        final_records = sorted(records_list, key=lambda record: record["_id"])
        record_cacher.cache_records(final_records)

        return final_records
