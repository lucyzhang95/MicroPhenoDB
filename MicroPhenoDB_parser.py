"""MicroPhenoDB Parser Class
Contains the main parser class for processing MicroPhenoDB data.
"""

import os
import uuid

from utils.cache_manager import CacheManager
from utils.cache_pipeline import DataCachePipeline
from utils.reader import FileReader


class MicroPhenoDBParser:
    """Orchestrates the entire data processing pipeline for MicroPhenoDB."""

    def __init__(self, data_dir="downloads"):
        self.data_dir = data_dir
        self.file_reader = FileReader()
        self.cache_manager = CacheManager()
        self.cache_pipeline = DataCachePipeline()

    def _get_file_path(self, filename):
        return os.path.join(self.data_dir, filename)

    @staticmethod
    def _get_organism_type(node) -> str:
        """Determine an organism type based on taxonomic lineage."""
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
        """Create a microbe node from taxon mapping."""
        node = taxon_map.get(name.lower())
        if node:
            node = node.copy()
            node["type"] = "biolink:OrganismTaxon"
            node["organism_type"] = self._get_organism_type(node)
            node.pop("mapping_tool", None)
            return node
        return None

    @staticmethod
    def _get_disease_node(name, disease_map):
        """Get the disease node from disease mapping."""
        return (
            disease_map.get(name.lower())
            if name.lower() != "null" and name.lower() != "not foundthogenic"
            else None
        )

    @staticmethod
    def _get_association_node(score, position, qualifier, pmid, pub_map, anatomical_map):
        """Create an association node with all relevant metadata."""
        assoc = {
            "category": "biolink:OrganismalEntityAsAModelOfDiseaseAssociation",
            "predicate": "biolink:affects",
            "primary_knowledge_source": ["infores:Disbiome", "infores:HMDAD"],
            "aggregator_knowledge_source": "infores:MicroPhenoDB",
            "has_evidence": "ECO:0000305",  # manual assertion
            "agent_type": "biolink:manual_agent",
            "score": float(score),
            "publication": pub_map.get(pmid),
            "anatomical_entity": anatomical_map.get(
                position.lower(),
                {"original_name": position.lower(), "type": "biolink:AnatomicalEntity"},
            ),
        }
        if qualifier and qualifier.lower() != "tendency":
            assoc["qualifier"] = qualifier.lower()
        return assoc

    def _remove_empty_values(self, obj):
        """Recursively removes keys/items that have None or empty values."""
        if isinstance(obj, dict):
            return {
                k: self._remove_empty_values(v)
                for k, v in obj.items()
                if v not in (None, {}, [], "")
            }
        if isinstance(obj, list):
            return [self._remove_empty_values(v) for v in obj if v not in (None, {}, [], "")]
        return obj

    def _remove_trailing_spaces(self, data):
        """Recursively removes leading/trailing whitespace from all string values."""
        if isinstance(data, dict):
            return {key: self._remove_trailing_spaces(value) for key, value in data.items()}
        elif isinstance(data, list):
            return [self._remove_trailing_spaces(item) for item in data]
        elif isinstance(data, str):
            return data.strip()
        else:
            return data

    @staticmethod
    def _remove_duplicate_records(records):
        """
        Removes duplicate records based on the '_id' field.
        If a record with the same '_id' exists with the same pmid but different scores,
        take the one with the highest score.
        """
        unique_records = {}
        for record in records:
            _id = record.get("_id")
            if not _id:
                continue

            if _id not in unique_records:
                unique_records[_id] = record
            else:
                new_assoc = record.get("association", {})
                new_pub = new_assoc.get("publication", {})
                new_pmid = new_pub.get("pmid")
                new_score = new_assoc.get("score")

                existing_record = unique_records[_id]
                existing_assoc = existing_record.get("association", {})
                existing_pub = existing_assoc.get("publication", {})
                existing_pmid = existing_pub.get("pmid")
                existing_score = existing_assoc.get("score")

                if (
                        new_pmid is not None
                        and new_score is not None
                        and new_pmid == existing_pmid
                        and abs(new_score) > abs(existing_score)
                ):
                    unique_records[_id] = record

        return list(unique_records.values())

    def load_microphenodb_data(self):
        """Loads and yields the final processed data records."""
        self.cache_pipeline.run()

        taxon_data = self.cache_manager.get_or_cache_taxon_info()
        disease_data = self.cache_manager.get_or_cache_disease_info()
        pubmed_data = self.cache_manager.get_or_cache_pubmed_metadata()
        anatomical_data = self.cache_manager.get_or_cache_anatomical_entity2uberon()

        core_f_path = self._get_file_path("core_table.txt")

        records = []
        for line in self.file_reader.read_file(core_f_path):
            _, organism_name, disease_name, score, pmid, _, position, qualifier = (
                part.strip() for part in line
            )

            # subject node (microbe/taxon)
            subject_node = self._get_microbe_node(organism_name, taxon_data)
            subject_node = self._remove_empty_values(subject_node)

            # object node (disease)
            object_node = self._get_disease_node(disease_name, disease_data)
            object_node = self._remove_empty_values(object_node)
            if object_node:
                description = object_node.get("description")
                if description and "null" in str(description).lower():
                    del object_node["description"]

            # skip records without a valid subject or object
            if not subject_node or not object_node:
                continue

            # association node
            association_node = self._get_association_node(
                score, position, qualifier, pmid, pubmed_data, anatomical_data
            )
            association_node = self._remove_empty_values(association_node)

            _id = (
                f"{subject_node['id'].split(':')[1]}"
                f"_{association_node['predicate'].split(':')[1]}"
                f"_{object_node['id'].split(':')[1]}"
                if subject_node and object_node
                else str(uuid.uuid4())
            )

            records.append(
                {
                    "_id": _id,
                    "association": association_node,
                    "object": object_node,
                    "subject": subject_node,
                }
            )

        print(f"Loaded a total of {len(records)} records with duplicated _id from {core_f_path}.")

        # remove duplicates and clean data
        unique_records = self._remove_duplicate_records(records)
        for record in unique_records:
            cleaned_record = self._remove_trailing_spaces(record)
            yield cleaned_record
