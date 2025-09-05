"""Record Manager for MicroPhenoDB Parser
Manages caching and processing of the final parsed records.
"""
import glob
import os

from MicroPhenoDB_parser import MicroPhenoDBParser
from utils.cache_helper import CacheHelper, RecordHelper


class RecordCacheManager:
    """Manages caching of the final processed records."""

    def __init__(self):
        self.cache_helper = CacheHelper()
        self.rec_helper = RecordHelper()

    def cache_records(self, records: list, filename_base="microphenodb_parsed_records"):
        """Cache the final processed records to a JSONL file."""
        if not records:
            print("***Warning: No records provided to cache.")
            return

        print(f"\n-> Processing {len(records)} records...")

        jsonl_filename = f"{filename_base}.jsonl"
        jsonl_path = self.rec_helper.save_jsonl(records, jsonl_filename)

        print("Records cached and exported successfully.")
        return jsonl_path

    def load_cached_records(self, filename_base="microphenodb_parsed_records"):
        """Load cached records from the most recent record JSONL file."""
        jsonl_pattern = os.path.join("records", f"{filename_base}_*.jsonl")
        jsonl_files = glob.glob(jsonl_pattern)

        if not jsonl_files:
            print(f"!!! No JSONL files found matching pattern: {filename_base}_*.jsonl")
            return None

        most_recent_jsonl = max(jsonl_files, key=os.path.getmtime)
        jsonl_filename = os.path.basename(most_recent_jsonl)

        try:
            records = self.rec_helper.load_jsonl(jsonl_filename)
            print(f"-> Loaded {len(records)} records from JSONL: {jsonl_filename}")
            return records
        except Exception as e:
            print(f"!!! Failed to load JSONL file {jsonl_filename}: {e}")
            return None

    @staticmethod
    def run_data_pipeline():
        """Run the full data pipeline: parse, process, cache, and export records."""
        parser = MicroPhenoDBParser()
        records_list = list(parser.load_microphenodb_data())
        final_records = sorted(records_list, key=lambda record: record["_id"])

        record_manager = RecordCacheManager()

        jsonl_path = record_manager.cache_records(final_records)

        return final_records, jsonl_path
