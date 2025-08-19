import tarfile

import text2term
from rapidfuzz import fuzz, process


class RapidFuzzUtils:
    """Utilities for fuzzy string matching using RapidFuzz."""

    def __init__(self, tar_path="taxdump.tar.gz"):
        """Loads and parses the NCBI taxonomy data upon initialization."""
        print(f"Initializing RapidFuzzUtils and loading data from '{tar_path}'...")

        self.ref_name_to_taxid = self._parse_names_dmp_from_taxdump(tar_path)
        self.ref_names = list(self.ref_name_to_taxid.keys())
        print(f"Ready. Loaded {len(self.ref_names)} reference names.")

    def _parse_names_dmp_from_taxdump(self, tar_path, f_name="names.dmp") -> dict:
        """Parses the names.dmp file from the NCBI Taxonomy database dump."""
        keep_classes = {
            "scientific name",
            "synonym",
            "equivalent name",
            "genbank synonym",
            "genbank anamorph",
        }
        name2taxid = {}
        with tarfile.open(tar_path, "r:gz") as tar_f:
            member = tar_f.getmember(f_name)
            with tar_f.extractfile(member) as in_f:
                for line in in_f:
                    parts = [part.strip().decode("utf-8") for part in line.strip().split(b"|")]
                    if (
                        len(parts) >= 4 and parts[3] in keep_classes
                    ):  # parts[0] == taxid, parts[1] == name
                        name2taxid[parts[1].lower()] = int(parts[0])
        return name2taxid

    def fuzzy_match(self, query_names: list, scorer=fuzz.token_sort_ratio, score_cutoff=90) -> dict:
        """
        Performs fuzzy matching of query names against the pre-loaded reference names.
        """
        matches = {}
        for name in query_names:
            match = process.extractOne(
                name, self.ref_names, scorer=scorer, score_cutoff=score_cutoff
            )
            if match:
                matched_name, score, _ = match
                matches[name] = {
                    "matched_name": matched_name,
                    "score": score,
                    "mapping_tool": "rapidfuzz",
                }
        return matches

    def fuzzy_matched_name2taxid(self, query_names: list) -> dict:
        """
        Finds the best fuzzy match for each query name and returns its corresponding taxid.
        """
        fuzzy_matches = self.fuzzy_match(query_names)

        for _original_name, match_info in fuzzy_matches.items():
            matched_name = match_info["matched_name"]
            taxid = self.ref_name_to_taxid.get(matched_name)
            if taxid:
                match_info["taxid"] = taxid

        return fuzzy_matches


class Text2TermUtils:
    def text2term_name2id(
        self, disease_names, ontology="MONDO", url="http://purl.obolibrary.org/obo/mondo.owl"
    ):
        """Maps disease names to ontology identifiers using text2term."""
        if not text2term.cache_exists(ontology):
            text2term.cache_ontology(ontology_url=url, ontology_acronym=ontology)
        df = text2term.map_terms(
            list(set(disease_names)), ontology, use_cache=True, min_score=0.8, max_mappings=1
        )
        if df is None or df.empty:
            return {}
        df = df[~df["Mapped Term CURIE"].astype(str).str.contains("NCBITAXON", na=False)]

        text2term_mapped = {}
        for _, row in df.iterrows():
            source_term = row["Source Term"]
            mapped_curie = row["Mapped Term CURIE"]

            text2term_mapped[source_term] = {
                "id": mapped_curie,
                "mapping_tool": "text2term",
                "xrefs": {mapped_curie.split(':')[0].lower(): mapped_curie},
            }
        return text2term_mapped
