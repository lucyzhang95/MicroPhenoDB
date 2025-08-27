"""Data Cache Pipeline for MicroPhenoDB Parser
Orchestrates the caching of all data types.
"""

import os

from utils.cache_manager import CacheManager


class DataCachePipeline:
    """Pipeline for caching data for the parser."""

    def __init__(self, cache_dir="cache", downloads_dir="downloads"):
        self.downloads_dir = downloads_dir
        self.cache_dir = cache_dir
        self.ncit_path = os.path.join(self.downloads_dir, "NCIT.txt")
        self.cache_manager = CacheManager(cache_dir)

    def run(self):
        """Runs the full data caching pipeline."""
        print("▶️ Caching taxon data...")
        self._cache_taxon_data()

        print("\n▶️ Caching disease data...")
        self._cache_disease_data()

        print("\n▶️ Caching publication data...")
        self._cache_pubmed_data()

        print("\n▶️ Caching anatomical entity to Uberon mapping...")
        self._cache_anatomical_data()

        print("\n🎉 All pre-data successfully cached.\n")

    def _cache_taxon_data(self):
        """Caches all taxon-related mappings and information."""
        self.cache_manager.get_or_cache_ncits2taxids_mapping()
        self.cache_manager.get_or_cache_ete3_taxon_name2taxid()
        self.cache_manager.get_or_cache_entrez_taxon_name2taxid()
        self.cache_manager.get_or_cache_rapidfuzz_taxon_name2taxid()
        self.cache_manager.get_or_cache_bt_taxon_name2taxid()
        self.cache_manager.get_or_cache_manual_taxon_name2taxid()
        self.cache_manager.get_or_cache_taxon_info()

    def _cache_disease_data(self):
        """Caches all disease-related mappings and information."""
        self.cache_manager.get_or_cache_disease_name2efo()
        self.cache_manager.get_or_cache_disease_info()

    def _cache_pubmed_data(self):
        """Caches publication metadata from PubMed."""
        self.cache_manager.get_or_cache_pubmed_metadata()

    def _cache_anatomical_data(self):
        """Caches anatomical entity to Uberon ID mapping."""
        self.cache_manager.get_or_cache_anatomical_entity2uberon()
