import asyncio
import csv
import json
import os
import pickle
import re
import ssl
import tarfile
import uuid

import aiohttp
import biothings_client as bt
import chardet
import text2term
from Bio import Entrez
from dotenv import load_dotenv
from ete3 import NCBITaxa
from rapidfuzz import fuzz, process
from tqdm.asyncio import tqdm

load_dotenv()


class CacheHelper:
    """Handles low-level saving and loading of data to/from the filesystem."""

    def __init__(self, cache_dir="cache"):
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


class FileReader:
    """Handles reading and decoding of input files."""

    def _detect_encoding(self, in_file_path):
        """Detects the character encoding of a file."""
        with open(in_file_path, "rb") as f:
            raw = f.read(5000)
        return chardet.detect(raw)["encoding"]

    def read_file(self, in_file_path, has_header=True):
        """Reads a file and yields its lines, handling encoding issues."""
        encoding = self._detect_encoding(in_file_path)
        if encoding == "ascii":
            encoding = "utf-8"
        try:
            with open(in_file_path, "r", encoding=encoding, errors="ignore") as in_f:
                reader = csv.reader(in_f, delimiter="\t")
                if has_header:
                    try:
                        next(reader)
                    except StopIteration:
                        return
                for line in reader:
                    yield line
        except UnicodeDecodeError as e:
            print(f"Unicode error with {encoding} on file {in_file_path}: {e}")
        except FileNotFoundError:
            print(f"Error: File not found at {in_file_path}")


class TextStructurePreprocessor:
    """Rules for cleaning the structure of text strings."""

    def remove_special_char(self, name: str) -> str:
        """Removes special characters from the name."""
        return re.sub(r"[?!#*&+]", "", name).strip()

    def remove_colon4name(self, name: str) -> str:
        """Replaces colons with spaces in the name."""
        return re.sub(r":", " ", name).strip()

    def remove_dot4name_except_in_sp(self, name: str) -> str:
        """Removes dots from the name except in 'sp.' or 'spp.'."""
        name = re.sub(r"\b(sp|spp)\.", r"\1__dot__", name)
        numeric_matches = re.findall(r"\d+(?:\.\d+)+", name)
        for match in numeric_matches:
            protected = match.replace(".", "__dot__")
            name = name.replace(match, protected)
        name = name.replace(".", "")
        return name.replace("__dot__", ".").strip()

    def remove_hyphen4name(self, name: str) -> str:
        """Replaces hyphens with spaces in the name, handling specific cases."""
        name = name.replace("butyrate-producing", "BUTYRATEPRODUCING")
        name = re.sub(r"(?<=[a-z])-(?=\d)", " ", name)
        name = re.sub(r"(?<=[a-z])-(?=(like|associated|related|positive|negative|)\b)", " ", name)
        return name.replace("BUTYRATEPRODUCING", "butyrate-producing").strip()

    def split_name_by_slash(self, name: str) -> str:
        """Splits the name by slashes, keeping the first part."""
        return re.split(r"(?<=[a-zA-Z])/\s*(?=[a-zA-Z])", name)[0].strip()

    def split_on_conjunction_in_name(self, name, keyword, prefer):
        """Splits the name on a conjunction, returning the preferred part."""
        pattern = r"\b" + re.escape(keyword) + r"\b"
        parts = [p.strip() for p in re.split(pattern, name, maxsplit=1) if p.strip()]
        if not parts:
            return name.strip()
        if prefer == "left":
            return parts[0]
        if prefer == "right":
            return parts[1] if len(parts) > 1 else parts[0]
        return name

    def remove_parentheses(self, name: str) -> str:
        """Removes parentheses and their contents from the name."""
        return re.sub(r"\s*\(.+\)", "", name).strip()

    def process_comma(self, name: str) -> str:
        """Processes commas in the name, removing trailing words."""
        name = re.sub(r",\s*(and\s+)?[a-z]{1,3}\b", "", name)
        return name.split(",")[0].strip()


class TextSemanticPreprocessor:
    """Rules for cleaning the semantics of text strings."""

    def remove_non_english_chars_in_name(self, name: str) -> str:
        """Removes non-ASCII characters from the name."""
        return re.sub(r"[^\x00-\x7F]+", "", name)

    def remove_and_in_name(self, name: str) -> str:
        """Removes 'and' from the name, keeping the first part."""
        return name.split(" and ")[0].strip()

    def remove_in_and_one_word_after_in_name(self, name: str) -> str:
        """Removes 'in' followed by a single word from the name."""
        name = re.sub(r"\bin\b\s+\w+\b", "", name)
        return re.sub(r"\s{2,}", " ", name).strip()

    def remove_word_related_in_name(self, name: str) -> str:
        """Removes 'related' and similar terms from the name."""
        return re.sub(r"\b\w+-related\s*", "", name).strip()

    def remove_leading_non_in_name(self, name: str) -> str:
        """Removes leading non-taxonomic terms from the name."""
        return re.sub(r"^(?:non\w+\s+)+", "", name).strip()

    def remove_spp_in_name(self, name: str) -> str:
        """Removes 'sp.' or 'spp.' from the name."""
        return re.sub(r"\bspps?\b.*", "", name).strip()

    def remove_strain_in_taxon_name(self, name: str) -> str:
        """Removes 'strain' or 'strains' from the name."""
        return re.sub(r"\s{2,}", " ", re.sub(r"\bstrains?\b", "", name).strip())

    def remove_type_in_taxon_name(self, name: str) -> str:
        """Removes 'type' or 'types' from the name."""
        return re.sub(r"\btypes?\b", "", name).strip()

    def remove_group_in_taxon_name(self, name: str) -> str:
        """Removes 'group' or 'groups' from the name."""
        name = re.sub(r"\bgroups?\s+[a-z]\b", "", name)
        name = re.sub(r"\b(groups?|subgroup)\b$", "", name)
        return re.sub(r"\b(groups?|subgroup)\b(?=\s+[a-z])", "", name).strip()

    def remove_serovar_in_taxon_name(self, name: str) -> str:
        """Removes 'serovar' or 'serovars' from the name."""
        return re.sub(r"\bserovars?.*\b", "", name).strip()

    def remove_pre_postfix_in_taxon_name(self, name: str) -> str:
        """Removes common prefixes and suffixes from the name."""
        name = re.sub(r"^\s*b\s+", "", name)
        name = re.sub(r"\bstains?\b", "", name)
        name = re.sub(
            r"(coagulase negative|(?:non)?hemolytic|sensu lato|complex(?:es)?|incertae sedis|rapid growers)",
            "",
            name,
        )
        return re.sub(r"\s{2,}", " ", name).strip()

    def expand_taxon_name_abbrev(self, name: str) -> str:
        """Expands common abbreviations in taxon names."""
        expansions = {
            r"\be histolytica\b": "entamoeba histolytica",
            r"\bp multocida\b": "pasteurella multocida",
            r"\bhsv\b": "herpes simplex virus",
            r"\bhpv\b": "human papillomavirus",
            r"\bhiv\b": "human immunodeficiency virus",
            r"\blcm\b": "lymphocytic choriomeningitis",
            r"\bcluster\b": "family",
            r"\bspirochaeta\b": "borrelia",
            r"\bpiv\b": "parainfluenza virus",
            r"\btm7\b": "candidatus saccharimonadota",
            r"\brubella\b": "rubella virus",
            r"\bmumps\b": "mumps virus",
            r"\bntm\b": "mycobacteriales",  # changed from "non-tuberculous mycobacteria" to "mycobacteriales"
            r"^\bsr1\b$": "candidatus absconditibacteriota",
            r"zygomycetes": "mucoromycota",  # zygomycetes is an obsolete term for mucoromycota and zoopagomycota.
        }
        for pattern, replacement in expansions.items():
            name = re.sub(pattern, replacement, name)
        return re.sub(r"\s{2,}", " ", name).strip()

    def replace_taxon_name(self, name: str) -> str:
        """Standardizes plural or variant taxon names to their singular form and correct spelling."""
        replacements = {
            r"streptococci": "streptococcus",
            r"lactobacilli": "lactobacillus",
            r"enterococci": "enterococcus",
            r"staphylococci": "staphylococcus",
            r"coxsackie": "coxsackievirus",
            r"gemellales": "gemella",
            r"\bparainfluenza virus\b": "orthorubulavirus",
            r"\bparainfluenza\b": "orthorubulavirus",  # make rank broader to include all parainfluenza viruses
            r"\bpapovavirus\b": "papillomavirus",  # papovavirus obsolete term for both papillomavirus and polyomavirus
        }
        for old, new in replacements.items():
            name = re.sub(rf"\b{old}\b", new, name, flags=re.IGNORECASE)
        return name.strip()


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


class EbiTaxonomyService:
    """Handles NCIt API services and related tasks.
    This service fetches NCIT codes and their corresponding NCBI Taxonomy IDs from the EBI OLS API.
    """

    def get_ncit_code(self, ncit_file_path) -> list:
        """Extracts NCIT codes from NCIT.txt file."""
        reader = FileReader()
        return sorted(
            [
                line[2].split("_")[1]
                for line in reader.read_file(ncit_file_path, has_header=False)
                if "NCIT" in line[2]
            ]
        )

    async def async_query_ebi_ncit_code_to_taxid(self, session, ncit_code: str):
        """Map NCIT identifier to NCBI Taxid using EBI API.

        :param session: aiohttp session for making requests
        :param ncit_code: a NCIT code e.g., "C125969"
        :return:
        """
        url = f"https://www.ebi.ac.uk/ols4/api/ontologies/ncit/terms?iri=http://purl.obolibrary.org/obo/NCIT_{ncit_code}"
        try:
            async with session.get(url, timeout=10) as resp:
                if resp.status != 200:
                    print(f"Failed to connect ebi API: status {resp.status}")
                    return None
                data = await resp.json()
                terms = data.get("_embedded", {}).get("terms", [{}])[0]
                name = terms.get("label", "").lower()
                description = next(iter(terms.get("description", [])), "")
                annot = terms.get("annotation", {})
                if "NCBI_Taxon_ID" in annot:
                    taxid = annot["NCBI_Taxon_ID"][0]
                    return name, {
                        "id": f"NCBITaxon:{int(taxid)}",
                        "taxid": int(taxid),
                        "description": f"{description}[NCIT]" if description else "",
                        "xrefs": {"ncit": ncit_code},
                    }
                return name, {
                    "description": f"{description}[NCIT]" if description else "",
                    "xrefs": {"ncit": ncit_code},
                }
        except aiohttp.ClientError as e:
            print(f"EBI API request failed for NCIT_{ncit_code}: {e}")
            return None

    async def async_query_ebi_ncit_codes_to_taxids(self, ncit_codes: list):
        """Map NCIT identifiers to NCBI Taxids using EBI API

        :param ncit_codes: a list of NCIT codes e.g., ["C85924", "C83526", ...]
        :return ncit2taxid: a dictionary mapping NCIT codes to taxids
        {'Trichostrongylus colubriformis': {'taxid': 6319, 'ncit': 'C125969', 'description': 'A species of...'}}
        :return notfound_ncit: a dictionary with NCIT codes failed to map taxid
        {'Trypanosoma brucei gambiense': {'ncit': 'C125975', 'description': 'A species of parasitic protozoa...'}}
        """
        ncit2taxids, notfound_ncit = {}, {}
        async with aiohttp.ClientSession() as session:
            tasks = [
                self.async_query_ebi_ncit_code_to_taxid(session, ncit_code)
                for ncit_code in ncit_codes
            ]
            results = await tqdm.gather(*tasks, desc="Querying taxids from NCIT codes...")
        for result in results:
            if result:
                name, info = result
                if "taxid" in info:
                    ncit2taxids[name] = info
                else:
                    notfound_ncit[name] = info
        return ncit2taxids, notfound_ncit

    def async_run_ncit_codes_to_taxids(self, ncit_codes: list):
        """Maps NCIT codes to NCBI Taxonomy IDs."""
        ncit2taxids, notfound_ncit = asyncio.run(
            self.async_query_ebi_ncit_codes_to_taxids(ncit_codes)
        )
        print("✅ EBI NCIT Codes to TaxID mapping completed!")
        return ncit2taxids, notfound_ncit


class ETE3TaxonomyService:
    """Handles interactions with local NCBI Taxonomy database file0s via ETE3."""

    def __init__(self):
        if not os.path.exists("taxdump.tar.gz"):
            print("Taxonomy database not found. Downloading and updating...")
            ssl._create_default_https_context = ssl._create_unverified_context
            self.ncbi_taxa = NCBITaxa()
            self.ncbi_taxa.update_taxonomy_database()
            print("NCBI Taxonomy Database update complete.")
        else:
            print("NCBI Taxonomy database found. Loading into memory...")
            self.ncbi_taxa = NCBITaxa()
        print("ETE3TaxonomyService is initialized and ready.")

    def ete3_taxon_name2taxid(self, taxon_names: list) -> dict:
        """Maps taxon names to NCBI Taxonomy IDs using ETE3."""
        name2taxid = self.ncbi_taxa.get_name_translator(sorted(list(set(taxon_names))))
        print("✅ ETE3 Taxonomy Name to TaxID mapping completed!")
        return {
            name: {"taxid": int(taxid_list[0]), "mapping_tool": "ete3"}
            for name, taxid_list in name2taxid.items()
            if taxid_list
        }


class EntrezTaxonomyService:
    """Handles interactions with NCBI Entrez API."""

    EMAIL = os.getenv("EMAIL_ADDRESS")

    def __init__(self):
        Entrez.email = self.EMAIL
        self.semaphore = asyncio.Semaphore(3)

    async def async_query_entrez_taxon_name2taxid(
        self, taxon_name: str, max_retries=3, initial_backoff=1
    ):
        """Asynchronously queries NCBI Entrez for a taxon name."""
        for attempt in range(max_retries):
            try:
                async with self.semaphore:

                    def blocking_entrez_call():
                        """Performs a blocking call to the Entrez API."""
                        handle = Entrez.esearch(
                            db="taxonomy", term=taxon_name, retmode="xml", retmax=1
                        )
                        rec = Entrez.read(handle)
                        handle.close()
                        return rec

                    rec = await asyncio.to_thread(blocking_entrez_call)
                    await asyncio.sleep(1)

                    if rec and rec.get("IdList"):
                        return taxon_name, {
                            "taxid": int(rec["IdList"][0]),
                            "mapping_tool": "entrez",
                        }
                    else:
                        return None

            except Exception as e:
                print(f"Attempt {attempt + 1}/{max_retries} failed for '{taxon_name}': {e}")
                if attempt + 1 == max_retries:
                    break

                backoff_time = initial_backoff * (2**attempt)
                print(f"Retrying in {backoff_time} seconds...")
                await asyncio.sleep(backoff_time)

        print(f"All retries failed for '{taxon_name}'.")
        return None

    async def async_query_entrez_taxon_names2taxids(self, taxon_names: list) -> list:
        """Runs Entrez queries for a list of names concurrently."""
        tasks = [self.async_query_entrez_taxon_name2taxid(name) for name in taxon_names]
        results = await tqdm.gather(*tasks, desc="Querying Entrez Taxonomy Names...")
        return dict(res for res in results if res)

    def async_run_entrez_taxon_names2taxids(self, taxon_names: list) -> dict:
        results = asyncio.run(self.async_query_entrez_taxon_names2taxids(taxon_names))
        print("✅ Entrez Taxonomy Name to TaxID mapping completed!")
        return results


class BioThingsService:
    """Handles interactions with BioThings API for taxonomic information."""

    def query_bt_taxon_name2taxid(self, taxon_names: list) -> dict:
        """Maps taxon names to NCBI Taxonomy IDs using BioThings API."""
        client = bt.get_client("taxon")
        results = client.querymany(
            list(set(taxon_names)), scopes=["scientific_name", "other_names"], fields=["taxid"]
        )
        return {
            item["query"]: {"taxid": int(item["taxid"]), "mapping_tool": "biothings"}
            for item in results
            if "notfound" not in item
        }

    def query_bt_taxon_info(self, taxids: list) -> dict:
        """Fetches taxon information from BioThings API."""
        if not taxids:
            return {}
        client = bt.get_client("taxon")
        infos = client.gettaxa(
            list(set(taxids)), fields=["scientific_name", "parent_taxid", "lineage", "rank"]
        )
        return {
            info["_id"]: {
                "id": f"NCBITaxon:{int(info['_id'])}",
                "taxid": int(info["_id"]),
                "name": info.get("scientific_name"),
                "parent_taxid": info.get("parent_taxid"),
                "lineage": info.get("lineage"),
                "rank": info.get("rank"),
            }
            for info in infos
            if "notfound" not in info
        }

    def query_bt_disease_info(self, ids: list) -> dict:
        """Fetches disease information from BioThings API using Mondo ontology."""
        if not ids:
            return {}
        client = bt.get_client("disease")
        results = client.querymany(
            list(set(ids)),
            scopes=["mondo.xrefs.hp", "mondo.xrefs.efo", "mondo.mondo"],
            fields=["mondo.definition", "mondo.label"],
        )
        d_info = {}
        for info in results:
            if "notfound" not in info:
                mondo = info.get("mondo", {})
                d_info[info["query"]] = {
                    "id": info["_id"],
                    "name": mondo.get("label"),
                    "description": mondo.get("definition"),
                    "type": "biolink:Disease",
                }
        return d_info


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
            }
        return text2term_mapped


class PubMedService:
    """Handles interactions with NCBI PubMed service via Entrez."""

    def __init__(self, email):
        Entrez.email = email

    def query_pubmed_metadata(self, pmids):
        """Fetches PubMed metadata for a list of PMIDs."""
        if not pmids:
            return {}
        handle = Entrez.efetch(db="pubmed", id=",".join(map(str, set(pmids))), retmode="xml")
        records = Entrez.read(handle)
        handle.close()
        result = {}
        for article in records.get("PubmedArticle", []):
            try:
                pmid = str(article["MedlineCitation"]["PMID"])
                article_data = article["MedlineCitation"]["Article"]
                title = article_data.get("ArticleTitle", "")
                abstract = " ".join(
                    str(p) for p in article_data.get("Abstract", {}).get("AbstractText", [])
                )
                doi = next(
                    (
                        str(el)
                        for el in article.get("PubmedData", {}).get("ArticleIdList", [])
                        if el.attributes.get("IdType") == "doi"
                    ),
                    "",
                )
                result[pmid] = {
                    "pmid": int(pmid),
                    "name": title,
                    "summary": f"{abstract} [abstract]" if abstract else "",
                    "doi": doi,
                    "type": "biolink:Publication",
                }
            except (KeyError, IndexError) as e:
                print(f"Failed to parse PubMed article: {e}")
        return result


class NCIt2TaxidMapper:
    """Handles the mapping of names to taxonomic identifiers using various services."""

    def __init__(self):
        self.ncit_service = EbiTaxonomyService()

    # names in notfound_ncit
    _MANUAL_TAXID_MAPPING_PATCHES = {
        "human parainfluenza virus": {"taxid": 2905673},
        "trichoderma": {"taxid": 5543},
        "trichomonas vaginalis": {"taxid": 5722},
        "trypanosoma brucei gambiense": {"taxid": 31285},
        "malassezia furfur": {"taxid": 55194},
        "clostridium cluster iv": {"taxid": 1689151},
        "clostridium cluster xvi": {"taxid": 543347},
        "powassan virus": {"taxid": 11083},
        "mycobacterium xenopi": {"taxid": 1789},
        "clostridiales xi": {"taxid": 186804},
        "clostridiales xiii": {"taxid": 189325},
        "alpha-amylase (aspergillus oryzae)": {
            "taxid": 5062,
            "description": "A fungus used in East Asia to saccharify rice, sweet potato, and barley in the making of alcoholic beverages such as sake and shochu, and also to ferment soybeans for making soy sauce and miso. It is one of the different koji molds used for food fermentation.[Wikipedia]",
        },
        "japanese encephalitis virus strain nakayama-nih antigen (formaldehyde inactivated)": {
            "taxid": 11076,
            "description": "A virus from the family Flaviviridae, part of the Japanese encephalitis serocomplex of nine genetically and antigenically related viruses, some of which are particularly severe in horses, and four of which, including West Nile virus, are known to infect humans. The enveloped virus is closely related to the West Nile virus and the St. Louis encephalitis virus. The positive sense single-stranded RNA genome is packaged in the capsid which is formed by the capsid protein.[Wikipedia]",
        },
    }

    # existing taxids in ncits2taxids
    _MANUAL_TAXID_MAPPING_OVERRIDES = {
        "bacteroides dorei": {"taxid": 357276},  # NCIT mapped 357256 is wrong
        "ruminococcaceae": {"taxid": 216572},  # NCIT mapped 541000 is obsolete
        "peptococcus": {"taxid": 2740},  # NCIT mapped 2840 is wrong
        "mumps virus": {"taxid": 2560602},  # NCIT mapped 11161 is obsolete
        "lymphocytic choriomeningitis virus": {"taxid": 3052303},  # NCIT mapped 11623 is obsolete
        "trichosporon": {"taxid": 5552},  # NCIT mapped 599816 is obsolete
        "bacillus cereus": {"taxid": 1396},  # NCIT mapped 13968 is wrong
        "bacillaceae": {"taxid": 186817},  # NCIT mapped 18681 is wrong
        "actinobaculum": {
            "taxid": 76833,
            "description": "Actinobaculum is a bacterial genus in the family Actinomycetaceae.[Wikipedia]",
        },
        "rikenellaceae": {
            "taxid": 171550,
            "description": "Rikenellaceae is a family of Gram-negative bacteria described by Noel R. Krieg in 2015. It contains nine genera, five of which are validly published by the International Code of Nomenclature of Prokaryotes.[2] Bacteria with 16S ribosomal RNA highly similar to the Rikenella genus, as compared to the larger taxonomic order Bacteroidales, are classified in this family. This family consists of non-motile, rod-shaped bacteria that are tolerant of bile. Most Rikenellaceae species have been identified in the gastrointestinal tract microbiomes of various animals. Bacteria of this taxonomic family are elevated in the gut microbiomes of mice that are leptin-resistant obese and diabetic.[4] However, Rikenellaceae bacteria are depleted in the gut microbiomes of obese American adults, leading to reduced synthesis of butyrate and disrupted metabolism. Gut microbiomes with elevated levels of Rikenellaceae bacteria are associated with lupus and Alzheimer's disease in mice and colorectal cancer in humans.[Wikipedia]",
        },
        "bacillus coagulans": {
            "taxid": 1398,
            "description": "H. coagulans is a Gram-positive, catalase-positive, spore-forming, motile, facultative anaerobe rod that measures approximately 0.9 μm by 3.0 μm to 5.0 μm. It may appear Gram negative when entering the stationary phase of growth. The optimum temperature for growth is 50 °C (122 °F); the range of temperatures tolerated is 30–55 °C (86–131 °F). IMViC tests VP and MR (methyl red) are positive.[Wikipedia]",
        },
        "aspergillus clavatus": {
            "taxid": 5057,
            "description": "Aspergillus clavatus is a species of fungus in the genus Aspergillus with conidia dimensions 3–4.5 x 2.5–4.5 μm. It is found in soil and animal manure. The fungus was first described scientifically in 1834 by the French mycologist John Baptiste Henri Joseph Desmazières. The fungus can produce the toxin patulin, which may be associated with disease in humans and animals. This species is only occasionally pathogenic. Other sources have identified many species of Aspergillus as producing dry, hydrophobic spores that are easily inhaled by humans and animals. Due to the small size of the spores, about 70% of spores of A. fumigatus are able to penetrate into the trachea and primary bronchi and close to 1% into alveoli. Inhalation of spores of Aspergillus is a health risk. A. clavatus is allergenic, causing the occupational hypersensitivity pneumonitis known as malt-worker's lung.[Wikipedia]",
        },
        "alternaria alternata": {
            "taxid": 5599,
            "description": "Alternaria alternata is a fungus causing leaf spots, rots, and blights on many plant parts, and other diseases. It is an opportunistic pathogen on over 380 host species of plant. It can also cause upper respiratory tract infections and asthma in humans with compromised immunity.[Wikipedia]",
        },
        "sarcina": {
            "taxid": 1266,
            "description": "Sarcina is a genus of gram-positive cocci bacteria in the family Clostridiaceae. A synthesizer of microbial cellulose, various members of the genus are human flora and may be found in the skin and large intestine. The genus takes its name from the Latin word 'sarcina,' meaning pack or bundle, after the cuboidal (2x2x2) cellular associations they form during division along three planes. The genus's type species is Sarcina ventriculi, a variety found on the surface of cereal seeds, in soil, mud, and in the stomachs of humans, rabbits, and guinea pigs.[Wikipedia]",
        },
    }

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
        ncit_codes = self.ncit_service.get_ncit_code(ncit_path)
        ncits2taxids, notfound_ncit = self.ncit_service.async_run_ncit_codes_to_taxids(ncit_codes)

        for name, patch_info in self._MANUAL_TAXID_MAPPING_PATCHES.items():
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

        for name, override_taxid in self._MANUAL_TAXID_MAPPING_OVERRIDES.items():
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

    def get_or_cache_ncits2taxids_mapping(self):
        """Checks if the NCIT to NCBI Taxonomy ID mapping is cached."""
        cache_f_name = "ncit2taxid.pkl"
        cache_f_path = os.path.join(self.cache_dir, cache_f_name)

        if os.path.exists(cache_f_path):
            print("NCIT to NCBI Taxonomy mapping already cached. Loading...")
            return self.load_pickle(cache_f_name)
        else:
            print("Caching NCIT to NCBI Taxonomy ID mapping...")
            ncits2taxids = self.id_mapper.ncits2taxids()
            self.save_pickle(ncits2taxids, cache_f_name)
            print("✅ NCIT to NCBI Taxonomy ID mapping successfully cached.")
            return ncits2taxids

    def _get_all_taxon_names(self, file_path=None):
        if file_path is None:
            file_path = os.path.join(self.data_dir, "core_table.txt")
        core_data = self.file_reader.read_file(file_path)
        return sorted(list(set(line[1].lower().strip() for line in core_data if line)))

    def _get_taxon_names_for_id_map(self):
        """Gets taxon names that need to be mapped to NCBI Taxonomy IDs.
        1767 names in core table
        515 names in NCIT.txt
        1252 names left
        """
        core_taxon_names = self._get_all_taxon_names()
        cached_ncit2taxids = self.get_or_cache_ncits2taxids_mapping()
        ncit_taxon_names = set(cached_ncit2taxids.keys())
        taxon_names_to_map = set(core_taxon_names) - ncit_taxon_names
        return sorted(list(taxon_names_to_map))

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

    def get_or_cache_ete3_taxon_name2taxid(self):
        """Caches the ETE3 taxon name to NCBI Taxonomy ID mapping."""
        cache_f_name = "ete3_taxon_name2taxid.pkl"
        cache_f_path = os.path.join(self.cache_dir, cache_f_name)

        if os.path.exists(cache_f_path):
            print("ETE3 Taxon Name to TaxID mapping already cached. Loading...")
            return self.load_pickle(cache_f_name)
        else:
            print("Caching ETE3 Taxon Name to TaxID mapping...")
            taxon_names4ete3 = self._get_taxon_names_for_ete3_mapping()
            ete3_mapped = self.ete3_service.ete3_taxon_name2taxid(taxon_names4ete3)
            self.save_pickle(ete3_mapped, cache_f_name)
            print("✅ ETE3 Taxon Name to TaxID mapping successfully cached.")
            return ete3_mapped

    def _get_taxon_names_for_entrez_mapping(self):
        """Gets taxon names that need to be mapped to NCBI Taxonomy IDs using Entrez."""
        taxon_names = self._get_taxon_names_for_ete3_mapping()
        cached_ete3_taxon_names = self.get_or_cache_ete3_taxon_name2taxid()
        ete3_taxon_names = set(cached_ete3_taxon_names.keys())
        taxon_names_to_map = set(taxon_names) - ete3_taxon_names
        return sorted(list(taxon_names_to_map))

    def get_or_cache_entrez_taxon_name2taxid(self):
        """Caches the Entrez taxon name to NCBI Taxonomy ID mapping."""
        cache_f_name = "entrez_taxon_name2taxid.pkl"
        cache_f_path = os.path.join(self.cache_dir, cache_f_name)

        if os.path.exists(cache_f_path):
            print("Entrez Taxon Name to TaxID mapping already cached. Loading...")
            return self.load_pickle(cache_f_name)
        else:
            print("Caching Entrez Taxon Name to TaxID mapping...")
            taxon_names_to_map = self._get_taxon_names_for_entrez_mapping()
            entrez_mapped = self.entrez_service.async_run_entrez_taxon_names2taxids(
                taxon_names_to_map
            )
            self.save_pickle(entrez_mapped, cache_f_name)
            print("✅ Entrez Taxon Name to TaxID mapping successfully cached.")
            return entrez_mapped

    def _get_taxon_names_for_rapidfuzz_mapping(self):
        """Gets taxon names that need to be mapped to NCBI Taxonomy IDs using Rapidfuzz."""
        taxon_names = self._get_taxon_names_for_ete3_mapping()
        cached_ete3_taxon_names = self.get_or_cache_ete3_taxon_name2taxid()
        ete3_taxon_names = set(cached_ete3_taxon_names.keys())
        cached_entrez_taxon_names = self.get_or_cache_entrez_taxon_name2taxid()
        entrez_taxon_names = set(cached_entrez_taxon_names.keys())

        taxon_names_to_map = set(taxon_names) - ete3_taxon_names - entrez_taxon_names
        return sorted(list(taxon_names_to_map))

    def get_or_cache_rapidfuzz_taxon_name2taxid(self):
        """Caches the BioThings taxon name to NCBI Taxonomy ID mapping."""
        cache_f_name = "rapidfuzz_taxon_name2taxid.pkl"
        cache_f_path = os.path.join(self.cache_dir, cache_f_name)

        if os.path.exists(cache_f_path):
            print("RapidFuzz Taxon Name to TaxID mapping already cached. Loading...")
            return self.load_pickle(cache_f_name)
        else:
            print("Caching RapidFuzz Taxon Name to TaxID mapping...")
            taxon_names_to_map = self._get_taxon_names_for_rapidfuzz_mapping()
            rapidfuzz_mapped = RapidFuzzUtils().fuzzy_matched_name2taxid(taxon_names_to_map)
            self.save_pickle(rapidfuzz_mapped, cache_f_name)
            print("✅ RapidFuzz Taxon Name to TaxID mapping successfully cached.")
            return rapidfuzz_mapped

    def _get_taxon_names_for_bt_mapping(self):
        """Gets taxon names that need to be mapped to NCBI Taxonomy IDs using BioThings."""
        taxon_names = self._get_taxon_names_for_ete3_mapping()
        cached_ete3_taxon_names = self.get_or_cache_ete3_taxon_name2taxid()
        ete3_taxon_names = set(cached_ete3_taxon_names.keys())
        cached_entrez_taxon_names = self.get_or_cache_entrez_taxon_name2taxid()
        entrez_taxon_names = set(cached_entrez_taxon_names.keys())
        cached_rapid_fuzz_taxon_names = self.get_or_cache_rapidfuzz_taxon_name2taxid()
        rapidfuzz_taxon_names = set(cached_rapid_fuzz_taxon_names.keys())

        taxon_names_to_map = (
            set(taxon_names) - ete3_taxon_names - entrez_taxon_names - rapidfuzz_taxon_names
        )
        return sorted(list(taxon_names_to_map))

    OVERRIDE_BT_mapped_TAXID = {
        "influenza": {"taxid": 11320},  # influenza a virus
        "bifidobacterium infantis": {"taxid": 1682},
        "powassan": {"taxid": 11083},  # powassan virus
        "rubella virus virus": {"taxid": 11041},  # rubella virus
        "st louis encephalitis": {"taxid": 11080},  # st. louis encephalitis virus
        "yeasts": {"taxid": 5206},  # cryptococcus
    }

    def get_or_cache_bt_taxon_name2taxid(self):
        """Caches the BioThings taxon name to NCBI Taxonomy ID mapping."""
        cache_f_name = "bt_taxon_name2taxid.pkl"
        cache_f_path = os.path.join(self.cache_dir, cache_f_name)

        if os.path.exists(cache_f_path):
            print("BioThings Taxon Name to TaxID mapping already cached. Loading...")
            return self.load_pickle(cache_f_name)
        else:
            print("Caching BioThings Taxon Name to TaxID mapping...")
            taxon_names_to_map = self._get_taxon_names_for_bt_mapping()
            bt_mapped = self.bt_service.query_bt_taxon_name2taxid(taxon_names_to_map)

            # override some taxids with manual mapping
            for name, override_info in self.OVERRIDE_BT_mapped_TAXID.items():
                if name in bt_mapped:
                    bt_mapped[name].update({"taxid": override_info["taxid"]})

            if "microorganism" in bt_mapped:
                del bt_mapped["microorganism"]  # remove "microorganism" as it is too generic
            self.save_pickle(bt_mapped, cache_f_name)
            print("✅ BioThings Taxon Name to TaxID mapping successfully cached.")
            return bt_mapped

    def _get_unmapped_taxon_names(self):
        """Gets taxon names that were not mapped by any service."""
        taxon_names = self._get_taxon_names_for_ete3_mapping()
        cached_ete3_taxon_names = self.get_or_cache_ete3_taxon_name2taxid()
        ete3_taxon_names = set(cached_ete3_taxon_names.keys())
        cached_entrez_taxon_names = self.get_or_cache_entrez_taxon_name2taxid()
        entrez_taxon_names = set(cached_entrez_taxon_names.keys())
        cached_rapid_fuzz_taxon_names = self.get_or_cache_rapidfuzz_taxon_name2taxid()
        rapidfuzz_taxon_names = set(cached_rapid_fuzz_taxon_names.keys())
        cached_bt_taxon_names = self.get_or_cache_bt_taxon_name2taxid()
        bt_taxon_names = set(cached_bt_taxon_names.keys())

        unmapped_taxon_names = (
            set(taxon_names)
            - set(ete3_taxon_names)
            - set(entrez_taxon_names)
            - set(rapidfuzz_taxon_names)
            - set(bt_taxon_names)
        )
        return sorted(list(unmapped_taxon_names))

    _MANUAL_MAP_UNMAPPED_TAXON_NAMES = {
        "actinomyces neuii subspecies neuii": 144053,
        "actinomyces radingae": 131110,
        "actinomyces turicensis": 131111,
        "candidate division candidatus saccharimonadota genomosp": 239137,
        "candidate division candidatus saccharimonadota single cell isolate tm7b": 447455,
        "candidate division candidatus saccharimonadota single cell isolate tm7c": 447456,
        "clostridia family i": 31979,
        "clostridium family xiva": 543317,
        "clostridium family xviii": 189325,
        "clostridium lituseburense": 1537,
        "clostridium xi": 186804,
        "clostridium xivb": 543317,
        "coxsackievirus a virus": 12066,
        "creutzfeldt jakob disease": 36469,
        "cryptococcus albidus": 100951,
        "cysticercosis": 6204,
        "escherichia vulneris": 566,
        "eubacterium tortuosum": 39494,
        "hookworms cestodes": 6157,
        "neisseriagonorrhoeae": 485,
        "orthorubulavirus 1c4": 2560526,
        "ovine jaagziekte virus": 11746,
        "prevotella multisaccharivorax": 310514,
        "prevotella nanceiensis": 425941,
        "saccharomyces castellii": 27288,
        "syphilis": 160,
        "trophyrema": 2039,
        "uncultured clostridiales ii": 186801,
        "vulvovaginal candidiasis": 5476,
    }

    def get_or_cache_manual_taxon_name2taxid(self):
        """Caches the manual mapping of taxon names to NCBI Taxonomy IDs."""
        cache_f_name = "manual_taxon_name2taxid.pkl"
        cache_f_path = os.path.join(self.cache_dir, cache_f_name)

        if os.path.exists(cache_f_path):
            print("Manual Taxon Name to TaxID mapping already cached. Loading...")
            return self.load_pickle(cache_f_name)
        else:
            print("Caching Manual Taxon Name to TaxID mapping...")
            manual_mapping = {
                name: {"taxid": taxid, "mapping_tool": "manual"}
                for name, taxid in self._MANUAL_MAP_UNMAPPED_TAXON_NAMES.items()
            }
            self.save_pickle(manual_mapping, cache_f_name)
            return manual_mapping

    # TODO: count the output record of this function
    def _convert_preprocessed_name2original_name(self) -> dict:
        """Converts preprocessed taxon names to their original names and updates the mapping with taxid."""
        processed_taxon_name_mapping = self._preprocessed_taxon_names_mapping()
        cached_ete3_taxon_names = self.get_or_cache_ete3_taxon_name2taxid()
        cached_entrez_taxon_names = self.get_or_cache_entrez_taxon_name2taxid()
        cached_rapidfuzz_taxon_names = self.get_or_cache_rapidfuzz_taxon_name2taxid()
        cached_bt_taxon_names = self.get_or_cache_bt_taxon_name2taxid()
        cached_manual_taxon_names = self.get_or_cache_manual_taxon_name2taxid()

        mapped_names = {
            **cached_ete3_taxon_names,
            **cached_entrez_taxon_names,
            **cached_rapidfuzz_taxon_names,
            **cached_bt_taxon_names,
            **cached_manual_taxon_names,
        }

        original_name_mapping = {}
        for original_name, preprocessed_name in processed_taxon_name_mapping.items():
            if preprocessed_name in mapped_names:
                info = mapped_names[preprocessed_name]
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

    def get_or_cache_taxon_info(self):
        """Caches and retrieves taxon information."""
        cache_f_name = "original_taxon_name_and_info.pkl"
        cache_f_path = os.path.join(self.cache_dir, cache_f_name)

        if os.path.exists(cache_f_path):
            print("Original taxon name and info already cached. Loading...")
            return self.load_pickle(cache_f_name)
        else:
            print("Caching original taxon name and info...")
            taxon_info = self.map_original_name2taxon_info()
            self.save_pickle(taxon_info, cache_f_name)
            print("✅ Original taxon name and info successfully cached.")
            return taxon_info

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

    def get_or_cache_disease_name2efo(self):
        """Caches and retrieves the disease name to EFO mapping."""
        cache_f_name = "efo_disease_name2id.pkl"
        cache_f_path = os.path.join(self.cache_dir, cache_f_name)

        if os.path.exists(cache_f_path):
            print("Disease name to EFO mapping already cached. Loading...")
            return self.load_pickle(cache_f_name)
        else:
            print("Caching disease name to EFO mapping...")
            efo_map = self._get_disease_efo_mapping()
            self.save_pickle(efo_map, cache_f_name)
            print("✅ Disease name to EFO mapping successfully cached.")
            return efo_map

    def get_or_cache_text2term_disease_name2id(self, disease_names_to_map):
        """Caches the Text2Term disease name to ID mapping."""
        cache_f_name = "text2term_disease_name2id.pkl"
        cache_f_path = os.path.join(self.cache_dir, cache_f_name)

        if os.path.exists(cache_f_path):
            print("Text2Term Disease Name to ID mapping already cached. Loading...")
            return self.load_pickle(cache_f_name)
        else:
            print("Caching Text2Term Disease Name to ID mapping...")
            text2term_mapped = self.t2t_utils.text2term_name2id(disease_names_to_map)
            self.save_pickle(text2term_mapped, cache_f_name)
            print("✅ Text2Term Disease Name to ID mapping successfully cached.")
            return text2term_mapped

    MANUAL_MAP_DISEASE_INFO = {
        "arthritis": {
            "id": "HP:0001369",
            "name": "arthritis",
            "original_name": "arthritis",
            "description": "Inflammation of a joint.[HP]",
            "type": "biolink:Disease",
            "xrefs": {"hp": "HP:0001369"},
        },
        "virus respiratory tract infection": {
            "id": "HP:0011947",
            "name": "respiratory tract infection",
            "original_name": "virus respiratory tract infection",
            "description": "An infection of the upper or lower respiratory tract.[HP]",
            "type": "biolink:Disease",
            "xrefs": {"hp": "HP:0011947"},
        },
        "urethritis": {
            "id": "HP:0500006",
            "name": "urethritis",
            "original_name": "urethritis",
            "description": "Urethritis is characterized by discharge, dysuria and/or urethral discomfort but may be asymptomatic. Common etiologies include gonococcal urethritis as well as infection by chlamydia trachomatis, mycoplasma genitalium, ureaplasma urealyticum, trichomonas vaginalis, anaerobes, herpes simplex virus, and adenovirus.[HP]",
            "type": "biolink:Disease",
            "xrefs": {"hp": "HP:0500006"},
        },
        "lungsaspergillosis": {
            "id": "MONDO:0000266",
            "name": "pulmonary aspergilloma",
            "original_name": "lungsaspergillosis",
            "description": "Aspergillosis is an infection, growth, or allergic response caused by the Aspergillus fungus. There are several different kinds of aspergillosis. One kind is allergic bronchopulmonary aspergillosis (also called ABPA), a condition where the fungus causes allergic respiratory symptoms similar to asthma, such as wheezing and coughing, but does not actually invade and destroy tissue. Another kind of aspergillosis is invasive aspergillosis. This infection usually affects people with weakened immune systems due to cancer, AIDS, leukemia, organ transplantation, chemotherapy, or other conditions or events that reduce the number of normal white blood cells. In this condition, the fungus invades and damages tissues in the body. Invasive aspergillosis most commonly affects the lungs, but can also cause infection in many other organs and can spread throughout the body (commonly affecting the kidneys and brain). Aspergilloma, a growth (fungus ball) that develops in an area of previous lung disease such as tuberculosis or lung abscess, is a third kind of aspergillosis. This type of aspergillosis is composed of a tangled mass of fungus fibers, blood clots, and white blood cells. The fungus ball gradually enlarges, destroying lung tissue in the process, but usually does not spread to other areas.[MONDO]",
            "type": "biolink:Disease",
            "xrefs": {"mondo": "MONDO:0000266"},
        },
        "gastric cancer": {
            "id": "MONDO:0001056",
            "name": "gastric cancer",
            "original_name": "gastric cancer",
            "description": "A primary or metastatic malignant neoplasm involving the stomach.[MONDO]",
            "type": "biolink:Disease",
            "xrefs": {"mondo": "MONDO:0001056"},
        },
        "vincent angina bacteria": {
            "id": "MONDO:0006865",
            "name": "necrotizing ulcerative gingivitis",
            "original_name": "vincent angina bacteria",
            "description": "A bacterial infectious process affecting the gums. It is characterized by the development of necrotic, ulcerated, and painful lesions with creation of pseudomembranes extending along the gingival margins.[MONDO]",
            "type": "biolink:Disease",
            "xrefs": {"mondo": "MONDO:0006865"},
        },
        "juvenile idiopathic arthritis": {
            "id": "MONDO:0011429",
            "name": "juvenile idiopathic arthritis",
            "original_name": "juvenile idiopathic arthritis",
            "description": "Juvenile idiopathic arthritis (JIA) is the term used to describe a group of inflammatory articular disorders of unknown cause that begin before the age of 16 and last over 6 weeks. The term juvenile idiopathic arthritis was chosen to signify the absence of any known mechanism underlying the disorder and to highlight the necessity of excluding other types of arthritis occurring in well defined diseases (in particular arthritis occurring in association with infectious, inflammatory and haematooncologic diseases).[MONDO]",
            "type": "biolink:Disease",
            "xrefs": {"mondo": "MONDO:0011429"},
        },
        "cryptosporidiosis": {
            "id": "MONDO:0015474",
            "name": "cryptosporidiosis",
            "original_name": "cryptosporidiosis",
            "description": "Intestinal infection with organisms of the genus Cryptosporidium. It occurs in both animals and humans. Symptoms include severe diarrhea.[MONDO]",
            "type": "biolink:Disease",
            "xrefs": {"mondo": "MONDO:0015474"},
        },
        "hemoglobinopathies": {
            "id": "MONDO:0019050",
            "name": "inherited hemoglobinopathy",
            "original_name": "hemoglobinopathies",
            "description": "An inherited disorder characterized by structural alterations of a globin chain within the hemoglobin molecule.[MONDO]",
            "type": "biolink:Disease",
            "xrefs": {"mondo": "MONDO:0019050"},
        },
        "acute hepatitis a virus infection": {
            "id": "MONDO:0005790",
            "name": "hepatitis a virus infection",
            "original_name": "acute hepatitis a virus infection",
            "description": "Acute inflammation of the liver caused by the hepatitis A virus. It is highly contagious and usually contracted through close contact with an infected individual or their feces, contaminated food or water.[MONDO]",
            "type": "biolink:Disease",
            "xrefs": {"mondo": "MONDO:0005790"},
        },
    }

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

    def get_or_cache_disease_info(self):
        """Caches total mapped disease with descriptions."""
        cache_f_name = "original_disease_name_and_info.pkl"
        cache_f_path = os.path.join(self.cache_dir, cache_f_name)

        if os.path.exists(cache_f_path):
            print("Disease Info already cached. Loading...")
            return self.load_pickle(cache_f_name)
        else:
            print("Caching Disease Info...")
            t2t_mapped = self.get_text2term_disease_id_and_bt_info()
            efo_mapped = self.get_or_cache_disease_name2efo()
            final_disease_info = {**t2t_mapped, **efo_mapped, **self.MANUAL_MAP_DISEASE_INFO}
            self.save_pickle(final_disease_info, cache_f_name)
            print("✅ Disease Info successfully cached.")
            return final_disease_info

    def _get_pubmed_metadata(self, file_path=None):
        """Reads the PubMed metadata file and returns a mapping of PMIDs to metadata."""
        if file_path is None:
            file_path = os.path.join(self.data_dir, "core_table.txt")
        core_data = self.file_reader.read_file(file_path)
        pmids = [line[4] for line in core_data if re.match(r"\b\d+\b", line[4])]
        pub_metadata = self.pubmed_service.query_pubmed_metadata(pmids)
        return pub_metadata

    def get_or_cache_pubmed_metadata(self):
        """Caches and retrieves PubMed metadata."""
        cache_f_name = "pubmed_metadata.pkl"
        cache_f_path = os.path.join(self.cache_dir, cache_f_name)

        if os.path.exists(cache_f_path):
            print("PubMed metadata already cached. Loading...")
            return self.load_pickle(cache_f_name)
        else:
            print("Caching PubMed metadata...")
            pub_metadata = self._get_pubmed_metadata()
            self.save_pickle(pub_metadata, cache_f_name)
            print("✅ PubMed metadata successfully cached.")
            return pub_metadata


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
            return node
        return None

    def _get_disease_node(self, name, disease_map):
        return (
            disease_map.get(name.lower())
            if name.lower() != "null" and name.lower() != "not foundthogenic"
            else {"original_name": name, "type": "biolink:Disease"}
        )

    def _get_association_node(self, score, position, qualifier, pmid, pub_map):
        assoc = {
            "predicate": "biolink:OrganismalEntityAsAModelOfDiseaseAssociation",
            "type": "biolink:associated_with",
            "score": float(score),
            "anatomical_entity": position.lower(),
            "aggregator_knowledge_source": "infores:MicroPhenoDB",
            "evidence_type": "ECO:0000305",  # manual assertion
            "publication": pub_map.get(pmid)
            if pmid in pub_map
            else {
                "pmid": int(pmid) if "NCIT" not in pmid and pmid != "#N/A" else None,
                "type": "biolink:Publication",
            },
        }
        if qualifier and qualifier.lower() != "tendency":
            assoc["qualifier"] = qualifier.lower()
        return assoc

    def _remove_empty_none_values(self, obj):
        if isinstance(obj, dict):
            cleaned = {}
            for k, v in obj.items():
                v_clean = self._remove_empty_none_values(v)
                if v_clean not in (None, {}, []):
                    cleaned[k] = v_clean
            return cleaned

        if isinstance(obj, list):
            cleaned_list = []
            for v in obj:
                v_clean = self._remove_empty_none_values(v)
                if v_clean not in (None, {}, []):
                    cleaned_list.append(v_clean)
            return cleaned_list
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
            self._remove_empty_none_values(subject_node)
            object_node = self._get_disease_node(disease_name, disease_map)
            self._remove_empty_none_values(object_node)
            if not subject_node or not object_node:
                continue
            association_node = self._get_association_node(score, position, qualifier, pmid, pub_map)
            self._remove_empty_none_values(association_node)
            yield {
                "_id": str(uuid.uuid4()),
                "subject": subject_node,
                "object": object_node,
                "association": association_node,
            }


class RecordCacheManager:
    """Manages caching of the final processed records."""

    def __init__(self, cache_helper: CacheHelper):
        self.cache_helper = cache_helper
        self.parser = MicroPhenoDBParser()

    def cache_records(self, records: list, filename_base="microphenodb"):
        """
        Caches a given list of records to both pickle and JSON formats.
        This method is now only responsible for SAVING, not generating.
        """
        if not records:
            print("Warning: No records provided to cache.")
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

        final_records = list(parser.load_microphenodb_data())
        record_cacher.cache_records(final_records)


if __name__ == "__main__":
    RecordCacheManager.run_data_pipeline()
    print("Full MicroPhenoDB parsed records cached successfully.")
