import functools
import os
import re
import uuid

from dotenv import load_dotenv

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
from utils.file_reader import FileReader
from utils.ontology_mapper import RapidFuzzUtils, Text2TermUtils
from utils.ontology_services import (
    BioThingsService,
    EbiTaxonomyService,
    EntrezTaxonomyService,
    ETE3TaxonomyService,
    PubMedService,
)
from utils.text_preprocessors import TextSemanticPreprocessor, TextStructurePreprocessor


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
        """Maps NCIT codes to NCBI Taxonomy IDs.
        Final mapping derived 582 unique NCIT codes to 568 NCBI Taxonomy IDs.
        :param ncit_path: Path to the NCIT.txt file.
        :return: A dictionary mapping NCIT codes to NCBI Taxonomy IDs.
        e.g.,
        'clostridiales xiii': {'description': 'A bacterial family of uncertain placement in the phylum...[NCIT]',
        'xrefs': {'ncit': 'C85925'},
        'id': 'NCBITaxon:189325',
        'taxid': 189325},
        'mapping_tool': 'ebi_ols'}
        """
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

        ncits2taxids.update(notfound_ncit)  # add notfound NCIT dicts to the main mapping

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
                    print(f"✅ {cache_f_name} already cached. Loading...")
                    return self.load_pickle(cache_f_name)
                else:
                    print(f"▶️ Caching {cache_f_name}...")
                    result = func(self, *args, **kwargs)
                    self.save_pickle(result, cache_f_name)
                    print(f"🎉 {cache_f_name} successfully cached.")
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
        """Gets taxon names that need to be mapped to NCBI Taxonomy IDs.
        1767 names in core table
        515 names in NCIT.txt
        1252 names left
        """
        core_taxon_names = self._get_all_taxon_names()
        cached_ncit2taxids = self.get_or_cache_ncits2taxids_mapping()
        return sorted(list(set(core_taxon_names) - set(cached_ncit2taxids.keys())))

    def _preprocessed_taxon_names_mapping(self) -> dict:
        """Preprocesses taxon names for mapping.
        1201 names after preprocessing, 51 names merged

        :return: A dictionary mapping preprocessed taxon names to their original names.
        """
        taxon_names = self._get_taxon_names_for_id_map()  # already extracted ncit mapped named
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

        # override some taxids with manual mapping
        for name, override_info in OVERRIDE_BT_mapped_TAXID.items():
            if name in bt_mapped:
                bt_mapped[name].update({"taxid": override_info["taxid"]})

        if "microorganism" in bt_mapped:
            del bt_mapped["microorganism"]  # remove "microorganism" as it is too generic
        return bt_mapped

    def _get_last_unmapped_taxon_names(self):
        """Gets taxon names that were not mapped by any service."""
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
        """

        taxon_info:
        '42862': {'id': 'NCBITaxon:42862',
        'taxid': 42862,
        'name': 'rickettsia felis',
        'parent_taxid': 114277,
        'lineage': [42862,
                    114277,
                    780,
                    33988,
                    775,
                    766,
                    28211,
                    1224,
                    3379134,
                    2,
                    131567,
                    1],
        'rank': 'species'},

        ncit2taxid:
        'clostridiales xiii': {'description': 'A bacterial family of uncertain placement in the phylum Firmicutes and the order Clostridiales that is used to classify the genus Mogibacterium.[NCIT]',
        'xrefs': {'ncit': 'C85925'},
        'id': 'NCBITaxon:189325',
        'taxid': 189325},
        'mapping_tool': 'ebi_ols'}

        """
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
            _id = _id.replace("_", ":")
            if ":" in _id:
                prefix = _id.split(":")[0].strip().lower()
                efo_map[sci_d_name] = {
                    "id": _id.upper(),
                    "name": sci_d_name,
                    "original_name": d_name,
                    "description": f"{desc}[{prefix.upper()}]",
                    "type": "biolink:Disease",
                    "xrefs": {prefix: _id.upper()},
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

    def get_text2term_disease_id_and_bt_info(self) -> dict:
        """Map disease names to disease ontology and info."""
        core_data = self.file_reader.read_file(os.path.join(self.data_dir, "core_table.txt"))
        core_disease_names = [
            line[2].lower().strip()
            for line in core_data
            if line[2].lower() != "null" and line[2].lower() != "not foundthogenic"
        ]
        efo_mapped = self.get_or_cache_disease_name2efo()
        initial_disease_names = sorted(list(set(core_disease_names) - set(efo_mapped.keys())))

        original_to_preprocessed_map = self.name_processor.convert_preprocessed_name2dict(
            initial_disease_names, self.name_processor.preprocess_disease_name
        )
        preprocessed_names_to_map = sorted(list(set(original_to_preprocessed_map.values())))

        preprocessed_name_to_id_map = self.get_or_cache_text2term_disease_name2id(
            preprocessed_names_to_map
        )
        unique_disease_ids = sorted(
            list(set(v["id"] for v in preprocessed_name_to_id_map.values() if v.get("id")))
        )
        disease_info = self.bt_service.query_bt_disease_info(unique_disease_ids)

        final_mapping = {}
        for original_name, preprocessed_name in original_to_preprocessed_map.items():
            id_info = preprocessed_name_to_id_map.get(preprocessed_name)
            disease_id = id_info.get("id") if id_info else None
            if disease_id and disease_id in disease_info:
                info = disease_info[disease_id].copy()
                info["original_name"] = original_name
                final_mapping[original_name] = info

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


class DataCachePipeline:
    """Pipeline for caching data for the parser."""

    def __init__(self, cache_dir="cache", downloads_dir="downloads"):
        self.downloads_dir = downloads_dir
        self.cache_dir = cache_dir
        self.ncit_path = os.path.join(self.downloads_dir, "NCIT.txt")
        self.cache_manager = CacheManager(cache_dir)

    def run_cache_taxon_names2taxids(self):
        """Caches the NCIT to NCBI Taxonomy ID mappings.
        1700 taxids in total combined from all sources. 1615 unique taxids.
        2 has no taxids
        """
        # taxon preprocessing and mapping
        ncit2taxid_mapping = self.cache_manager.get_or_cache_ncits2taxids_mapping()
        ete3_taxon_name2taxid = self.cache_manager.get_or_cache_ete3_taxon_name2taxid()
        entrez_taxon_name2taxid = self.cache_manager.get_or_cache_entrez_taxon_name2taxid()
        rapidfuzz_taxon_name2taxid = self.cache_manager.get_or_cache_rapidfuzz_taxon_name2taxid()
        bt_taxon_name2taxid = self.cache_manager.get_or_cache_bt_taxon_name2taxid()
        manual_taxon_name2taxid = self.cache_manager.get_or_cache_manual_taxon_name2taxid()
        taxon_info = self.cache_manager.get_or_cache_taxon_info()

        # disease mapping
        efo_disease_name2id = self.cache_manager.get_or_cache_disease_name2efo()
        disease_info = self.cache_manager.get_or_cache_disease_info()

        # pubmed metadata
        pubmed_metadata = self.cache_manager.get_or_cache_pubmed_metadata()

        return (
            ncit2taxid_mapping,
            ete3_taxon_name2taxid,
            entrez_taxon_name2taxid,
            rapidfuzz_taxon_name2taxid,
            bt_taxon_name2taxid,
            manual_taxon_name2taxid,
            taxon_info,
            efo_disease_name2id,
            disease_info,
            pubmed_metadata,
        )


class MicroPhenoDBParser:
    """Orchestrates the entire data processing pipeline for MicroPhenoDB."""

    def __init__(self, data_dir="downloads"):
        self.data_dir = data_dir
        self.file_reader = FileReader()
        self.name_processor = OntologyNameProcessor()
        self.cache_manager = CacheManager()
        self.cache_pipeline = DataCachePipeline()

    def _get_file_path(self, filename):
        return os.path.join(self.data_dir, filename)

    def _get_organism_type(self, node) -> str:
        type_map = {
            2: "biolink:Bacterium",
            2157: "Archaea",
            2759: "Eukaryote",
            10239: "biolink:Virus",
        }
        for taxid, biolink_type in type_map.items():
            if taxid in node.get("lineage", []):
                return biolink_type
        return "Other"

    def _get_microbe_node(self, name, taxon_map):
        node = taxon_map.get(name.lower())
        if node:
            node = node.copy()
            node["type"] = "biolink:OrganismTaxon"
            node["organism_type"] = self._get_organism_type(node)
            node.pop("mapping_tool", None)
            return node
        return None

    def _get_disease_node(self, name, disease_map):
        return (
            disease_map.get(name.lower())
            if name.lower() != "null" and name.lower() != "not foundthogenic"
            else None
        )

    def _get_association_node(self, score, position, qualifier, pmid, pub_map):
        assoc = {
            "predicate": "biolink:OrganismalEntityAsAModelOfDiseaseAssociation",
            "type": "biolink:associated_with",
            "score": float(score),
            "anatomical_entity": position.lower(),
            "aggregator_knowledge_source": "infores:MicroPhenoDB",
            "evidence_type": "ECO:0000305",  # manual assertion
            "publication": pub_map.get(pmid),
        }
        if qualifier and qualifier.lower() != "tendency":
            assoc["qualifier"] = qualifier.lower()
        return assoc

    def _remove_empty_values(self, obj):
        """
        Recursively removes keys/items that have None or "empty" values
        (e.g., '', {}, []).
        """
        if isinstance(obj, dict):
            return {
                k: self._remove_empty_values(v)
                for k, v in obj.items()
                if v not in (None, {}, [], "")
            }
        if isinstance(obj, list):
            return [self._remove_empty_values(v) for v in obj if v not in (None, {}, [], "")]
        return obj

    def load_microphenodb_data(self):
        """Loads and yields the final processed data records."""
        (
            _,
            _,
            _,
            _,
            _,
            _,
            taxon_map,
            _,
            disease_map,
            pub_map,
        ) = self.cache_pipeline.run_cache_taxon_names2taxids()

        core_f_path = self._get_file_path("core_table.txt")
        for line in self.file_reader.read_file(core_f_path):
            _, organism_name, disease_name, score, pmid, _, position, qualifier = (
                part.strip() for part in line
            )

            subject_node = self._get_microbe_node(organism_name, taxon_map)
            subject_node = self._remove_empty_values(subject_node)

            object_node = self._get_disease_node(disease_name, disease_map)
            object_node = self._remove_empty_values(object_node)
            if object_node:
                description = object_node.get("description")
                if description and "null" in str(description).lower():
                    del object_node["description"]

            if not subject_node or not object_node:
                continue

            association_node = self._get_association_node(score, position, qualifier, pmid, pub_map)
            association_node = self._remove_empty_values(association_node)

            yield {
                "_id": str(uuid.uuid4()),
                "association": association_node,
                "object": object_node,
                "subject": subject_node,
            }


class RecordCacheManager:
    """Manages caching of the final processed records."""

    def __init__(self, cache_helper: CacheHelper):
        self.cache_helper = cache_helper
        self.parser = MicroPhenoDBParser()

    def cache_records(self, records: list, filename_base="microphenodb_parsed_records"):
        """
        Caches a given list of records to both pickle and JSON formats.
        This method is now only responsible for SAVING, not generating.
        """
        if not records:
            print("❗️Warning: No records provided to cache.")
            return

        print(f"Caching {len(records)} records...")
        pickle_filename = f"{filename_base}.pkl"
        json_filename = f"{filename_base}.json"

        self.cache_helper.save_pickle(records, pickle_filename)
        self.cache_helper.save_json(records, json_filename)
        print("✅ Final records cached successfully.")

    @staticmethod
    def run_data_pipeline():
        """
        Orchestrates the entire process of parsing and caching data.
        """
        cache_helper = CacheHelper(cache_dir="cache")
        parser = MicroPhenoDBParser(data_dir="downloads")
        record_cacher = RecordCacheManager(cache_helper)
        records_list = list(parser.load_microphenodb_data())
        final_records = sorted(records_list, key=lambda record: record["_id"])
        record_cacher.cache_records(final_records)


if __name__ == "__main__":
    RecordCacheManager.run_data_pipeline()
    print("Full MicroPhenoDB parsed records cached successfully.")
