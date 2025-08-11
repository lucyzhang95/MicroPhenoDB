import json
import os
import pickle


class CacheHelper:
    """Handles low-level saving and loading of data to/from the filesystem."""

    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
    CACHE_DIR = os.path.join(SCRIPT_DIR, "cache")

    def __init__(self, cache_dir=CACHE_DIR):
        self.cache_dir = os.path.join(os.getcwd(), cache_dir)
        os.makedirs(self.cache_dir, exist_ok=True)

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
        with open(os.path.join(self.cache_dir, f_name), "w") as out_f:
            json.dump(obj, out_f, indent=4)

    def load_json(self, f_name):
        """Loads an object from a JSON file."""
        path = os.path.join(self.cache_dir, f_name)
        if os.path.exists(path):
            with open(path, "r") as in_f:
                return json.load(in_f)
        return None