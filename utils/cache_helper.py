import json
import os
import pickle
from datetime import datetime
from typing import Any, Dict, List


class CacheHelper:
    """Handles low-level saving and loading of data to/from the filesystem.
    pkl is only used for intermediate data storage."""

    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    CACHE_DIR = os.path.join(SCRIPT_DIR, "cache")
    DATA_DIR = os.path.join(SCRIPT_DIR, "parsed_records")

    def __init__(self, data_dir=DATA_DIR, cache_dir=CACHE_DIR):
        self.cache_dir = os.path.join(os.getcwd(), cache_dir)
        os.makedirs(self.cache_dir, exist_ok=True)
        self.data_dir = os.path.join(os.getcwd(), data_dir)
        os.makedirs(self.data_dir, exist_ok=True)

    def save_pickle(self, obj, f_name):
        """Saves an object to a pickle file."""
        with open(os.path.join(self.cache_dir, f_name), "wb") as out_f:
            pickle.dump(obj, out_f)

    def load_pickle(self, f_name):
        """Loads an object from a pickle file."""
        path = os.path.join(self.cache_dir, f_name)
        if os.path.exists(path):
            with open(path, "rb") as in_f:
                return pickle.load(in_f)
        return None

    def save_json(self, obj, f_name):
        """Saves an object to a JSON file."""
        with open(os.path.join(self.data_dir, f_name), "w") as out_f:
            json.dump(obj, out_f, indent=4)

    def load_json(self, f_name):
        """Loads an object from a JSON file."""
        path = os.path.join(self.data_dir, f_name)
        if os.path.exists(path):
            with open(path, "r") as in_f:
                return json.load(in_f)
        return None

    def save_jsonl(self, records: List[Dict[str, Any]], f_name: str):
        """Saves records to a JSONL file with standardized formatting."""
        if not records:
            print("‼️ Warning: No records provided for JSONL export.")
            return None

        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        base_name = f_name.replace(".jsonl", "") if f_name.endswith(".jsonl") else f_name
        jsonl_filename = f"{base_name}_{timestamp}.jsonl"
        jsonl_path = os.path.join(self.data_dir, jsonl_filename)

        exported_count = 0
        with open(jsonl_path, "w", encoding="utf-8") as f:
            for record in records:
                clean_record = self._standardize_record(record)
                json.dump(clean_record, f, ensure_ascii=False, separators=(",", ":"))
                f.write("\n")
                exported_count += 1

        print(f"✅ Exported {exported_count} records to {jsonl_path}")
        return jsonl_path

    def load_jsonl(self, f_name: str) -> List[Dict[str, Any]]:
        """Load records from JSONL file."""
        path = os.path.join(self.data_dir, f_name)

        if not os.path.exists(path):
            return []

        records = []
        with open(path, "r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if line:
                    try:
                        records.append(json.loads(line))
                    except json.JSONDecodeError:
                        continue

        print(f"✅ Loaded {len(records)} records from {f_name}")
        return records

    def _standardize_record(self, record: Dict[str, Any]) -> Dict[str, Any]:
        """
        SIMPLIFIED: Only essential transformations.
        Main goal: Convert xrefs from dict to list (your key requirement).
        """
        clean_record = record.copy()

        for entity_key in ["subject", "object"]:
            if entity_key in clean_record and isinstance(clean_record[entity_key], dict):
                entity = clean_record[entity_key].copy()
                if "xrefs" in entity and isinstance(entity["xrefs"], dict):
                    entity["xrefs"] = [f"{k}:{v}" for k, v in entity["xrefs"].items() if v]
                clean_record[entity_key] = entity

        return clean_record
