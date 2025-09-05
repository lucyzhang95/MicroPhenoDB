"""Record Statistics Report Generator
Generates comprehensive statistics for MicroPhenoDB parsed records.
"""

import json
import os
from collections import Counter, defaultdict
from datetime import datetime
from typing import Any, Dict, List

from utils.record_manager import RecordCacheManager


class RecordStatsReporter:
    """Generates comprehensive statistics for parsed records."""

    def __init__(self, records: List[Dict[str, Any]]):
        self.records = records
        self.stats = {}

    def _extract_curie_prefix(self, curie_id: str) -> str:
        """Extract prefix from CURIE format (e.g., 'NCBITaxon:123' -> 'NCBITaxon')"""
        if isinstance(curie_id, str) and ":" in curie_id:
            return curie_id.split(":", 1)[0]
        return "Unknown"

    def _safe_get_nested(self, record: dict, keys: list, default=None):
        """Safely get nested dictionary values."""
        current = record
        for key in keys:
            if isinstance(current, dict) and key in current:
                current = current[key]
            else:
                return default
        return current

    def generate_basic_stats(self):
        """Generate basic record count statistics."""
        self.stats["basic_stats"] = {
            "total_record_count": len(self.records),
            "total_unique_ids": len(set(record.get("_id", "") for record in self.records)),
        }

    def generate_id_duplication_analysis(self):
        """Analyze _id duplication patterns."""
        id_counts = Counter(record.get("_id", "") for record in self.records)
        duplicated_ids = {id_val: count for id_val, count in id_counts.items() if count > 1}

        duplication_analysis = {
            "total_duplicate_id_groups": len(duplicated_ids),
            "duplicate_records_count": sum(duplicated_ids.values()) - len(duplicated_ids),
            "duplicated_ids_detail": dict(
                sorted(duplicated_ids.items(), key=lambda x: x[1], reverse=True)
            ),
        }

        # sample duplicate records for investigation
        if duplicated_ids:
            sample_duplicate_id = list(duplicated_ids.keys())[0]
            sample_records = [r for r in self.records if r.get("_id") == sample_duplicate_id]
            duplication_analysis["sample_duplicate_records"] = sample_records[:3]  # only first 3

        self.stats["id_duplication_analysis"] = duplication_analysis

    def generate_subject_stats(self):
        """Generate subject-related statistics."""
        subject_ids = [self._safe_get_nested(r, ["subject", "id"]) for r in self.records]
        subject_ids = [sid for sid in subject_ids if sid]

        # CURIE analysis
        curie_prefixes = [self._extract_curie_prefix(sid) for sid in subject_ids]
        curie_counts = Counter(curie_prefixes)

        # organism type analysis
        organism_types = [
            self._safe_get_nested(r, ["subject", "organism_type"]) for r in self.records
        ]
        organism_types = [ot for ot in organism_types if ot]
        organism_type_counts = Counter(organism_types)

        # rank analysis
        ranks = [self._safe_get_nested(r, ["subject", "rank"]) for r in self.records]
        ranks = [rank for rank in ranks if rank]
        rank_counts = Counter(ranks)

        # description analysis
        descriptions = [self._safe_get_nested(r, ["subject", "description"]) for r in self.records]
        descriptions = [desc for desc in descriptions if desc]

        # xrefs analysis
        xrefs_analysis = defaultdict(int)
        for record in self.records:
            subject_xrefs = self._safe_get_nested(record, ["subject", "xrefs"], [])
            if isinstance(subject_xrefs, list):
                for curie in subject_xrefs:
                    prefix = self._extract_curie_prefix(curie)
                    xrefs_analysis[prefix] += 1

        # sample subject records
        sample_subjects = []
        seen_curie_prefixes = set()
        for record in self.records:
            subject = record.get("subject", {})
            subject_id = subject.get("id", "")
            prefix = self._extract_curie_prefix(subject_id)
            if prefix not in seen_curie_prefixes and len(sample_subjects) < 5:
                sample_subjects.append(subject)
                seen_curie_prefixes.add(prefix)

        self.stats["subject_stats"] = {
            "total_subjects_with_id": len(subject_ids),
            "unique_subject_ids": len(set(subject_ids)),
            "curie_prefixes": {
                "counts": dict(curie_counts),
                "total_unique_prefixes": len(curie_counts),
            },
            "organism_types": {
                "counts": dict(organism_type_counts),
                "total_unique_types": len(organism_type_counts),
            },
            "taxonomic_ranks": {
                "counts": dict(rank_counts),
                "total_unique_ranks": len(rank_counts),
            },
            "descriptions": {
                "total_with_description": len(descriptions),
                "description_coverage_percentage": (len(descriptions) / len(self.records)) * 100
                if self.records
                else 0,
            },
            "xrefs": {
                "xref_sources": dict(xrefs_analysis),
                "total_unique_xref_sources": len(xrefs_analysis),
            },
            "sample_subjects": sample_subjects,
        }

    def generate_object_stats(self):
        """Generate object (disease) related statistics."""
        object_ids = [self._safe_get_nested(r, ["object", "id"]) for r in self.records]
        object_ids = [oid for oid in object_ids if oid]

        # CURIE analysis
        curie_prefixes = [self._extract_curie_prefix(oid) for oid in object_ids]
        curie_counts = Counter(curie_prefixes)

        # description analysis
        descriptions = [self._safe_get_nested(r, ["object", "description"]) for r in self.records]
        descriptions = [desc for desc in descriptions if desc]

        # xrefs analysis
        xrefs_analysis = defaultdict(int)
        for record in self.records:
            subject_xrefs = self._safe_get_nested(record, ["object", "xrefs"], [])
            if isinstance(subject_xrefs, list):
                for curie in subject_xrefs:
                    prefix = self._extract_curie_prefix(curie)
                    xrefs_analysis[prefix] += 1

        # sample object records
        sample_objects = []
        seen_curie_prefixes = set()
        for record in self.records:
            obj = record.get("object", {})
            object_id = obj.get("id", "")
            prefix = self._extract_curie_prefix(object_id)
            if prefix not in seen_curie_prefixes and len(sample_objects) < 5:
                sample_objects.append(obj)
                seen_curie_prefixes.add(prefix)

        self.stats["object_stats"] = {
            "total_objects_with_id": len(object_ids),
            "unique_object_ids": len(set(object_ids)),
            "curie_prefixes": {
                "counts": dict(curie_counts),
                "total_unique_prefixes": len(curie_counts),
            },
            "descriptions": {
                "total_with_description": len(descriptions),
                "description_coverage_percentage": (len(descriptions) / len(self.records)) * 100
                if self.records
                else 0,
            },
            "xrefs": {
                "xref_sources": dict(xrefs_analysis),
                "total_unique_xref_sources": len(xrefs_analysis),
            },
            "sample_objects": sample_objects,
        }

    def generate_association_stats(self):
        """Generate association-related statistics."""
        # anatomical entity analysis
        anatomical_entities = []
        for record in self.records:
            anat_entity = self._safe_get_nested(record, ["association", "anatomical_entity"])
            if anat_entity:
                anatomical_entities.append(anat_entity)

        anatomical_original_names = [
            ae.get("original_name", "") for ae in anatomical_entities if isinstance(ae, dict)
        ]
        anatomical_counts = Counter(anatomical_original_names)

        # publication (PMID) analysis
        pmids = []
        publications = []
        for record in self.records:
            pub = self._safe_get_nested(record, ["association", "publication"])
            if pub and isinstance(pub, dict):
                publications.append(pub)
                pmid = pub.get("pmid")
                if pmid:
                    pmids.append(pmid)

        pmid_counts = Counter(pmids)

        # score analysis
        scores = [self._safe_get_nested(r, ["association", "score"]) for r in self.records]
        scores = [s for s in scores if s is not None]

        score_ranges = {
            "negative": len([s for s in scores if s < 0]),
            "zero": len([s for s in scores if s == 0]),
            "positive": len([s for s in scores if s > 0]),
        }

        # qualifier analysis
        qualifiers = [self._safe_get_nested(r, ["association", "qualifier"]) for r in self.records]
        qualifiers = [q for q in qualifiers if q]
        qualifier_counts = Counter(qualifiers)

        self.stats["association_stats"] = {
            "anatomical_entities": {
                "total_with_anatomical_entity": len(anatomical_entities),
                "unique_anatomical_entities": len(set(anatomical_original_names)),
                "anatomical_entity_counts": dict(anatomical_counts.most_common(20)),  # Top 20
                "sample_anatomical_entities": anatomical_entities[:5],
            },
            "publications": {
                "total_with_publication": len(publications),
                "total_unique_pmids": len(set(pmids)),
                "pmid_frequency_distribution": {
                    "top_10_most_cited_pmids": dict(pmid_counts.most_common(10))
                },
                "sample_publications": publications[:3],
            },
            "scores": {
                "total_with_score": len(scores),
                "score_ranges": score_ranges,
                "score_statistics": {
                    "min": min(scores) if scores else None,
                    "max": max(scores) if scores else None,
                    "mean": sum(scores) / len(scores) if scores else None,
                },
            },
            "qualifiers": {
                "total_with_qualifier": len(qualifiers),
                "unique_qualifiers": len(set(qualifiers)),
                "qualifier_counts": dict(qualifier_counts),
            },
        }

    def generate_comprehensive_report(self) -> Dict[str, Any]:
        """Generate comprehensive statistics report."""
        print("\nGenerating comprehensive record statistics...")

        # add metadata
        self.stats["metadata"] = {
            "report_version": "1.0",
            "analysis_date": datetime.now().isoformat(),
            "total_records_analyzed": len(self.records),
        }

        self.generate_basic_stats()
        self.generate_id_duplication_analysis()
        self.generate_subject_stats()
        self.generate_object_stats()
        self.generate_association_stats()

        return self.stats

    def save_report(self, output_path: str = "reports/record_stats.json"):
        """Save the statistics report to a JSON file."""
        os.makedirs(os.path.dirname(output_path), exist_ok=True)

        with open(output_path, "w", encoding="utf-8") as f:
            json.dump(self.stats, f, indent=2, ensure_ascii=False)

        print(f"-> Statistics report saved to: {output_path}")


def generate_record_statistics_report():
    """Main function to generate and save a record statistics report."""
    record_manager = RecordCacheManager()

    records = record_manager.load_cached_records()
    if not records:
        print("!!! No cached records found.")

    reporter = RecordStatsReporter(records)
    stats = reporter.generate_comprehensive_report()

    reporter.save_report()

    print("\nRECORD STATISTICS SUMMARY:")
    print(f"Total Records: {stats['basic_stats']['total_record_count']}")
    print(f"Unique IDs: {stats['basic_stats']['total_unique_ids']}")
    print(f"Duplicate ID Groups: {stats['id_duplication_analysis']['total_duplicate_id_groups']}")
    print(
        f"Subject CURIE Prefixes: {stats['subject_stats']['curie_prefixes']['total_unique_prefixes']}"
    )
    print(
        f"Object CURIE Prefixes: {stats['object_stats']['curie_prefixes']['total_unique_prefixes']}"
    )
    print(f"Unique PMIDs: {stats['association_stats']['publications']['total_unique_pmids']}")

    return stats


if __name__ == "__main__":
    generate_record_statistics_report()
