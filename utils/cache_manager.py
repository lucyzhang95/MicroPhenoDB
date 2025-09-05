"""Cache Manager for MicroPhenoDB Parser
Handles caching and retrieval of various data mappings.
"""

import functools
import logging
import os
import re
from collections import defaultdict

from dotenv import load_dotenv

from manual_annotations.anatomy2uberon import MANUAL_ANATOMICAL_ENTITY_MAPPING_OVERRIDE
from manual_annotations.disease_name2id import MANUAL_MAP_DISEASE_INFO
from manual_annotations.ncit2taxid import (
    MANUAL_TAXID_MAPPING_OVERRIDES,
    MANUAL_TAXID_MAPPING_PATCHES,
)
from manual_annotations.taxon_name2taxid import (
    MANUAL_MAP_UNMAPPED_TAXON_NAMES,
    OVERRIDE_BT_mapped_TAXID,
)
from utils.cache_helper import CacheHelper
from utils.ontology_mapper import RapidFuzzUtils, Text2TermUtils
from utils.ontology_services import (
    BioThingsService,
    EbiTaxonomyService,
    EntrezTaxonomyService,
    ETE3TaxonomyService,
    PubMedService,
    UberonService,
)
from utils.reader import FileReader
from utils.text_preprocessors import TextSemanticPreprocessor, TextStructurePreprocessor

log = logging.getLogger(__name__)


class OntologyNameProcessor(TextStructurePreprocessor, TextSemanticPreprocessor):
    """Applies a series of rules to clean and standardize ontology names."""

    def preprocess_taxon_name(self, name: str) -> str:
        """Preprocesses a taxon name by applying a series of cleaning rules."""
        name = name.lower().strip()
        name = self.remove_special_char(name)
        name = self.remove_colon4name(name)
        name = self.remove_dot4name_except_in_sp(name)
        name = self.remove_hyphen4name(name)
        name = self.remove_parentheses(name)
        name = self.split_name_by_slash(name)
        name = self.process_comma(name)
        name = self.remove_and_in_name(name)
        name = self.remove_strain_in_taxon_name(name)
        name = self.remove_type_in_taxon_name(name)
        name = self.remove_group_in_taxon_name(name)
        name = self.remove_serovar_in_taxon_name(name)
        name = self.remove_pre_postfix_in_taxon_name(name)
        name = self.remove_spp_in_name(name)
        name = self.expand_taxon_name_abbrev(name)
        name = self.replace_taxon_name(name)
        return name

    def preprocess_disease_name(self, name: str) -> str:
        """Preprocesses a disease name by applying a series of cleaning rules."""
        name = name.lower().strip()
        name = self.remove_non_english_chars_in_name(name)
        name = self.remove_parentheses(name)
        name = self.remove_in_and_one_word_after_in_name(name)
        name = self.split_on_conjunction_in_name(name, "with", "right")
        name = self.split_on_conjunction_in_name(name, "and", "right")
        name = self.split_on_conjunction_in_name(name, "among", "left")
        name = self.remove_word_related_in_name(name)
        name = self.remove_leading_non_in_name(name)
        return name

    def convert_preprocessed_name2dict(self, names: list, preprocessor_func) -> dict:
        """
        :param names:
        :param preprocessor_func: a function to preprocess the names,
        e.g., preprocess_taxon_name or preprocess_disease_name
        :return:
        {"original_name": "preprocessed_name"}
        """
        return {original_name: preprocessor_func(original_name) for original_name in names}


class NCIt2TaxidMapper:
    """Handles the mapping of names to taxonomic identifiers using various services."""

    def __init__(self):
        self.ncit_service = EbiTaxonomyService()

    def _get_ncit_code(self, ncit_file_path) -> list:
        """Extracts NCIT codes from NCIT.txt file."""
        reader = FileReader()
        return sorted(
            [
                line[2].split("_")[1]
                for line in reader.read_file(ncit_file_path, has_header=False)
                if "NCIT" in line[2]
            ]
        )

    def ncits2taxids(self, ncit_path=None):
        """Maps NCIT codes to NCBI Taxonomy IDs."""
        if ncit_path is None:
            ncit_path = os.path.join("downloads", "NCIT.txt")
        ncit_codes = self._get_ncit_code(ncit_path)
        ncits2taxids, notfound_ncit = self.ncit_service.async_run_ncit_codes_to_taxids(ncit_codes)

        for name, patch_info in MANUAL_TAXID_MAPPING_PATCHES.items():
            if name in notfound_ncit:
                notfound_ncit[name].update(
                    {
                        "id": f"NCBITaxon:{patch_info['taxid']}",
                        "taxid": patch_info["taxid"],
                        "description": notfound_ncit[name]["description"]
                        if notfound_ncit[name]["description"]
                        else patch_info.get("description", ""),
                    }
                )

        ncits2taxids.update(notfound_ncit)

        for name, override_taxid in MANUAL_TAXID_MAPPING_OVERRIDES.items():
            if name in ncits2taxids:
                ncits2taxids[name].update(
                    {
                        "id": f"NCBITaxon:{override_taxid['taxid']}",
                        "taxid": override_taxid["taxid"],
                        "description": ncits2taxids[name]["description"],
                        "mapping_tool": "ebi_ols"
                        if ncits2taxids[name]["description"]
                        else override_taxid.get("description", ""),
                    }
                )

        return ncits2taxids


class CacheManager(CacheHelper):
    load_dotenv(dotenv_path=".env")

    def __init__(self, cache_dir="cache", data_dir="downloads"):
        """Initializes the CacheManager"""
        super().__init__(cache_dir)
        self.data_dir = data_dir
        self.cache_dir = cache_dir
        self.file_reader = FileReader()
        self.id_mapper = NCIt2TaxidMapper()
        self.name_processor = OntologyNameProcessor()
        self.ete3_service = ETE3TaxonomyService()
        self.entrez_service = EntrezTaxonomyService()
        self.bt_service = BioThingsService()
        self.t2t_utils = Text2TermUtils()
        self.pubmed_service = PubMedService(email=os.getenv("EMAIL_ADDRESS"))
        self.uberon_service = UberonService()
        self.taxon_automated_mapping_sources = [
            self.get_or_cache_ete3_taxon_name2taxid,
            self.get_or_cache_entrez_taxon_name2taxid,
            self.get_or_cache_rapidfuzz_taxon_name2taxid,
            self.get_or_cache_bt_taxon_name2taxid,
        ]
        self.all_taxon_mapping_sources = self.taxon_automated_mapping_sources + [
            self.get_or_cache_manual_taxon_name2taxid
        ]

    @staticmethod
    def _cache_decorator(cache_f_name):
        """Decorator to cache the result of a function call."""

        def decorator(func):
            """Decorator to cache the result of a function call."""

            @functools.wraps(func)
            def wrapper(self, *args, **kwargs):
                """Wrapper function to handle caching logic."""
                cache_f_path = os.path.join(self.cache_dir, cache_f_name)
                if os.path.exists(cache_f_path):
                    print(f"[DONE] {cache_f_name} already cached. Loading...")
                    return self.load_pickle(cache_f_name)
                else:
                    print(f">>> Caching {cache_f_name}...")
                    result = func(self, *args, **kwargs)
                    self.save_pickle(result, cache_f_name)
                    print(f"[DONE] {cache_f_name} successfully cached.")
                    return result

            return wrapper

        return decorator

    def _get_unmapped_taxon_names(self, initial_names_set, list_of_mapping_funcs):
        """Filters an initial set of names against sources of mapped names."""
        mapped_names = set()
        for func in list_of_mapping_funcs:
            mapped_names.update(func().keys())
        return sorted(list(initial_names_set - mapped_names))

    def _get_all_taxon_names(self, file_path=None):
        if file_path is None:
            file_path = os.path.join(self.data_dir, "core_table.txt")
        core_data = self.file_reader.read_file(file_path)
        return sorted(list(set(line[1].lower().strip() for line in core_data if line)))

    @_cache_decorator("ncit2taxid.pkl")
    def get_or_cache_ncits2taxids_mapping(self):
        """Checks if the NCIT to NCBI Taxonomy ID mapping is cached."""
        return self.id_mapper.ncits2taxids()

    def _get_taxon_names_for_id_map(self):
        """Gets taxon names that need to be mapped to NCBI Taxonomy IDs."""
        core_taxon_names = self._get_all_taxon_names()
        cached_ncit2taxids = self.get_or_cache_ncits2taxids_mapping()
        return sorted(list(set(core_taxon_names) - set(cached_ncit2taxids.keys())))

    def _preprocessed_taxon_names_mapping(self) -> dict:
        """Preprocesses taxon names for mapping."""
        taxon_names = self._get_taxon_names_for_id_map()
        preprocessed_taxon_names = self.name_processor.convert_preprocessed_name2dict(
            taxon_names, self.name_processor.preprocess_taxon_name
        )
        return preprocessed_taxon_names

    def _get_taxon_names_for_ete3_mapping(self):
        """Gets taxon names that need to be mapped to NCBI Taxonomy IDs using ETE3."""
        processed_taxon_names_mapping = self._preprocessed_taxon_names_mapping()
        taxon_names4ete3 = list(set(processed_taxon_names_mapping.values()))
        return sorted(taxon_names4ete3)

    @_cache_decorator("ete3_taxon_name2taxid.pkl")
    def get_or_cache_ete3_taxon_name2taxid(self):
        """Caches the ETE3 taxon name to NCBI Taxonomy ID mapping."""
        taxon_names4ete3 = self._get_taxon_names_for_ete3_mapping()
        return self.ete3_service.ete3_taxon_name2taxid(taxon_names4ete3)

    def _get_taxon_names_for_entrez_mapping(self):
        """Gets taxon names that need to be mapped to NCBI Taxonomy IDs using Entrez."""
        initial_names = set(self._get_taxon_names_for_ete3_mapping())
        taxon_names_to_map = self._get_unmapped_taxon_names(
            initial_names, [self.get_or_cache_ete3_taxon_name2taxid]
        )
        return sorted(list(taxon_names_to_map))

    @_cache_decorator("entrez_taxon_name2taxid.pkl")
    def get_or_cache_entrez_taxon_name2taxid(self):
        """Caches the Entrez taxon name to NCBI Taxonomy ID mapping."""
        taxon_names_to_map = self._get_taxon_names_for_entrez_mapping()
        return self.entrez_service.async_run_entrez_taxon_names2taxids(taxon_names_to_map)

    def _get_taxon_names_for_rapidfuzz_mapping(self):
        """Gets taxon names that need to be mapped to NCBI Taxonomy IDs using Rapidfuzz."""
        initial_names = set(self._get_taxon_names_for_ete3_mapping())
        taxon_names_to_map = self._get_unmapped_taxon_names(
            initial_names,
            [self.get_or_cache_ete3_taxon_name2taxid, self.get_or_cache_entrez_taxon_name2taxid],
        )
        return sorted(list(taxon_names_to_map))

    @_cache_decorator("rapidfuzz_taxon_name2taxid.pkl")
    def get_or_cache_rapidfuzz_taxon_name2taxid(self):
        """Caches the BioThings taxon name to NCBI Taxonomy ID mapping."""
        taxon_names_to_map = self._get_taxon_names_for_rapidfuzz_mapping()
        return RapidFuzzUtils().fuzzy_matched_name2taxid(taxon_names_to_map)

    def _get_taxon_names_for_bt_mapping(self):
        """Gets taxon names that need to be mapped to NCBI Taxonomy IDs using BioThings."""
        initial_names = set(self._get_taxon_names_for_ete3_mapping())
        taxon_names_to_map = self._get_unmapped_taxon_names(
            initial_names,
            [
                self.get_or_cache_ete3_taxon_name2taxid,
                self.get_or_cache_entrez_taxon_name2taxid,
                self.get_or_cache_rapidfuzz_taxon_name2taxid,
            ],
        )
        return sorted(list(taxon_names_to_map))

    @_cache_decorator("bt_taxon_name2taxid.pkl")
    def get_or_cache_bt_taxon_name2taxid(self):
        """Caches the BioThings taxon name to NCBI Taxonomy ID mapping."""
        taxon_names_to_map = self._get_taxon_names_for_bt_mapping()
        bt_mapped = self.bt_service.query_bt_taxon_name2taxid(taxon_names_to_map)

        for name, override_info in OVERRIDE_BT_mapped_TAXID.items():
            if name in bt_mapped:
                bt_mapped[name].update({"taxid": override_info["taxid"]})

        if "microorganism" in bt_mapped:
            del bt_mapped["microorganism"]
        return bt_mapped

    def _get_last_unmapped_taxon_names(self):
        """Gets taxon names not mapped by any service."""
        initial_names = set(self._get_taxon_names_for_ete3_mapping())
        unmapped_taxon_names = self._get_unmapped_taxon_names(
            initial_names, self.taxon_automated_mapping_sources
        )
        return sorted(list(unmapped_taxon_names))

    @_cache_decorator("manual_taxon_name2taxid.pkl")
    def get_or_cache_manual_taxon_name2taxid(self):
        """Caches the manual mapping of taxon names to NCBI Taxonomy IDs."""
        manual_mapping = {
            name: {"taxid": taxid, "mapping_tool": "manual"}
            for name, taxid in MANUAL_MAP_UNMAPPED_TAXON_NAMES.items()
        }
        return manual_mapping

    def _convert_preprocessed_name2original_name(self) -> dict:
        """Converts preprocessed taxon names to their original names and updates the mapping with taxid."""
        processed_taxon_name_mapping = self._preprocessed_taxon_names_mapping()

        all_mapped_names = {}
        for source_func in self.all_taxon_mapping_sources:
            all_mapped_names.update(source_func())

        original_name_mapping = {}
        for original_name, preprocessed_name in processed_taxon_name_mapping.items():
            if preprocessed_name in all_mapped_names:
                info = all_mapped_names[preprocessed_name]
                original_name_mapping[original_name] = {
                    **info,
                    "name": preprocessed_name,
                }
        return original_name_mapping

    def _get_taxids_for_taxon_info(self) -> list:
        """Extracts taxids from the taxon information."""
        mapped_taxon_names = self._convert_preprocessed_name2original_name()
        ncit2taxid = self.get_or_cache_ncits2taxids_mapping()
        all_mapped_taxon_names = {**ncit2taxid, **mapped_taxon_names}
        taxids = [info["taxid"] for info in all_mapped_taxon_names.values() if "taxid" in info]
        return sorted(list(set(taxids)))

    def query_bt_taxon_info(self) -> dict:
        """Queries BioThings for taxon information based on taxids."""
        taxids = self._get_taxids_for_taxon_info()
        taxon_info = self.bt_service.query_bt_taxon_info(taxids)
        return taxon_info

    def map_original_name2taxon_info(self) -> dict:
        """Maps original taxon names to enriched taxon information."""
        taxon_info = self.query_bt_taxon_info()

        mapped_taxon_names = self._convert_preprocessed_name2original_name()
        ncit2taxid = self.get_or_cache_ncits2taxids_mapping()
        all_name_mappings = {**mapped_taxon_names, **ncit2taxid}

        final_mapping = {}
        for original_name, source_info in all_name_mappings.items():
            taxid = source_info.get("taxid")
            if taxid and str(taxid) in taxon_info:
                enriched_info = {
                    **source_info,
                    **taxon_info[str(taxid)],
                    "original_name": original_name,
                }
                final_mapping[original_name] = enriched_info

        return final_mapping

    @_cache_decorator("original_taxon_name_and_info.pkl")
    def get_or_cache_taxon_info(self):
        """Caches and retrieves taxon information."""
        return self.map_original_name2taxon_info()

    def _get_disease_efo_mapping(self, efo_f_path=None):
        """Reads the EFO file and returns a mapping of scientific names to disease information."""
        if efo_f_path is None:
            efo_f_path = os.path.join(self.data_dir, "EFO.txt")
        reader = FileReader()
        efo_map = {}
        for line in reader.read_file(efo_f_path):
            d_name, sci_d_name, _id, desc = (
                line[0].lower().strip(),
                line[1].lower().strip(),
                line[2],
                line[3].strip(),
            )
            _id = _id.replace("_", ":").lower()
            if ":" in _id:
                prefix = _id.split(":")[0].strip().lower()
                efo_map[sci_d_name] = {
                    "id": _id.upper() if "orphanet" not in _id else _id,
                    "name": sci_d_name,
                    "original_name": d_name,
                    "description": f"{desc}[{prefix.upper()}]",
                    "type": "biolink:Disease",
                    "xrefs": {prefix: _id.upper() if "orphanet" not in _id else _id},
                }
        return efo_map

    @_cache_decorator("efo_disease_name2id.pkl")
    def get_or_cache_disease_name2efo(self):
        """Caches and retrieves the disease name to EFO mapping."""
        return self._get_disease_efo_mapping()

    @_cache_decorator("text2term_disease_name2id.pkl")
    def get_or_cache_text2term_disease_name2id(self, disease_names_to_map):
        """Caches the Text2Term disease name to ID mapping."""
        return self.t2t_utils.text2term_name2id(disease_names_to_map)

    def _merge_xrefs(self, target_record: dict, source_record: dict) -> None:
        """Merges xrefs from a source record into a target record without overwriting."""
        source_xrefs = source_record.get("xrefs")
        if not source_xrefs:
            return

        target_xrefs = target_record.setdefault("xrefs", {})
        for key, value in source_xrefs.items():
            target_xrefs.setdefault(key, value)

    def _ensure_xref_from_id(self, record: dict) -> None:
        """Adds a self-referential xref from the record's ID if no xrefs exist."""
        if record.get("xrefs"):
            return

        record_id = record.get("id")
        if record_id and ":" in record_id:
            prefix = record_id.split(":", 1)[0].upper()
            record["xrefs"] = {prefix: record_id.upper()}

    def _resolve_new_disease_names(self, names_to_resolve: set[str]) -> tuple[dict, dict]:
        """Identify and preprocess unresolved disease names."""
        if not names_to_resolve:
            return {}, {}

        preprocess = self.name_processor.preprocess_disease_name
        original_to_preprocessed = {name: preprocess(name) for name in names_to_resolve}
        unique_preprocessed_names = sorted(set(original_to_preprocessed.values()))

        preprocessed_to_id_map = self.get_or_cache_text2term_disease_name2id(
            unique_preprocessed_names
        )
        disease_ids = sorted(
            list({v["id"] for v in preprocessed_to_id_map.values() if v and v.get("id")})
        )
        id_to_info_map = self.bt_service.query_bt_disease_info(disease_ids)

        return preprocessed_to_id_map, id_to_info_map

    def get_text2term_disease_id_and_bt_info(self) -> dict:
        """Maps disease names to their corresponding ontology IDs and metadata."""
        core_path = os.path.join(self.data_dir, "core_table.txt")
        core_data = self.file_reader.read_file(core_path)
        original_names = {
            line[2].lower().strip()
            for line in core_data
            if line[2].lower() not in ("null", "not foundthogenic")
        }

        efo_mapped = self.get_or_cache_disease_name2efo()
        efo_lookup = {**efo_mapped}
        for record in efo_mapped.values():
            if original_name := record.get("original_name"):
                efo_lookup.setdefault(original_name, record)

        final_mapping = {}
        names_to_resolve = set()
        for name in original_names:
            if name in efo_lookup:
                record = efo_lookup[name].copy()
                record["original_name"] = name
                final_mapping[name] = record
            else:
                names_to_resolve.add(name)

        preprocessed_to_id_map, id_to_info_map = self._resolve_new_disease_names(names_to_resolve)

        preprocess = self.name_processor.preprocess_disease_name
        missing_reasons = defaultdict(dict)

        for original_name in names_to_resolve:
            preprocessed_name = preprocess(original_name)
            id_info = preprocessed_to_id_map.get(preprocessed_name) or {}
            disease_id = id_info.get("id")
            if not disease_id:
                missing_reasons[original_name] = {
                    "reason": "no_text2term_id",
                    "preprocessed_name": preprocessed_name,
                }
                continue

            info = id_to_info_map.get(disease_id)
            if info:
                final_record = info.copy()
                final_record["id"] = disease_id
                final_record["original_name"] = original_name
                self._merge_xrefs(final_record, id_info)
                self._ensure_xref_from_id(final_record)
                final_mapping[original_name] = final_record
            else:
                final_mapping[original_name] = {
                    "id": disease_id,
                    "original_name": original_name,
                }
                missing_reasons[original_name] = {
                    "reason": "no_bt_info_for_t2t_id",
                    "id": disease_id,
                    "preprocessed_name": preprocessed_name,
                }

        return final_mapping

    @_cache_decorator("original_disease_name_and_info.pkl")
    def get_or_cache_disease_info(self):
        """Caches total mapped disease with descriptions."""
        t2t_mapped = self.get_text2term_disease_id_and_bt_info()
        efo_mapped = self.get_or_cache_disease_name2efo()
        final_disease_info = {**t2t_mapped, **efo_mapped, **MANUAL_MAP_DISEASE_INFO}
        return final_disease_info

    def _get_pubmed_metadata(self, file_path=None):
        """Reads the PubMed metadata file and returns a mapping of PMIDs to metadata."""
        if file_path is None:
            file_path = os.path.join(self.data_dir, "core_table.txt")
        core_data = self.file_reader.read_file(file_path)
        pmids = [line[4] for line in core_data if re.match(r"\b\d+\b", line[4])]
        pub_metadata = self.pubmed_service.query_pubmed_metadata(pmids)
        return pub_metadata

    @_cache_decorator("pubmed_metadata.pkl")
    def get_or_cache_pubmed_metadata(self):
        """Caches and retrieves PubMed metadata."""
        return self._get_pubmed_metadata()

    def _map_anatomical_entity2uberon(self, file_path=None) -> list:
        if file_path is None:
            file_path = os.path.join(self.data_dir, "core_table.txt")
        core_data = self.file_reader.read_file(file_path)
        anatomical_entities = sorted(
            list(
                set(
                    line[6].lower().strip()
                    for line in core_data
                    if line[6] and line[6].lower().strip() not in ("other", "unknown")
                )
            )
        )
        uberon_mapped = self.uberon_service.async_run_anatomical_entities_to_uberon_ids(
            anatomical_entities
        )

        for entity, _ in uberon_mapped.items():
            if entity in MANUAL_ANATOMICAL_ENTITY_MAPPING_OVERRIDE:
                manual_info = MANUAL_ANATOMICAL_ENTITY_MAPPING_OVERRIDE[entity]
                uberon_mapped[entity].update(
                    {
                        "id": manual_info["id"],
                        "name": manual_info["name"],
                        "original_name": manual_info["original_name"],
                        "category": manual_info["category"],
                    }
                )
        return uberon_mapped

    @_cache_decorator("uberon_anatomical_entity2uberon.pkl")
    def get_or_cache_anatomical_entity2uberon(self):
        """Caches and retrieves the anatomical entity to Uberon ID mapping."""
        return self._map_anatomical_entity2uberon()
