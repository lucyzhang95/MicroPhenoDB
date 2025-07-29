import asyncio
import csv
import json
import os
import pickle
import re
import ssl
import tarfile
import time
import uuid
from collections import Counter

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
        """ "Splits the name by slashes, keeping the first part."""
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
        }
        for pattern, replacement in expansions.items():
            name = re.sub(pattern, replacement, name)
        return re.sub(r"\s{2,}", " ", name).strip()


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

    def convert_preprocessed_name2dict(self, names: list) -> dict:
        """Converts a list of names to a dictionary with preprocessed names as keys."""
        return {original_name: self.preprocess_taxon_name(original_name) for original_name in names}


class NCBITaxonomyService:
    """Handles interactions with local NCBI Taxonomy database files."""

    def __init__(self):
        self.ncbi_taxa = None

    def initialize_local_ncbi_taxa(self):
        """Initializes the NCBITaxa from ETE3.
        This method ensures that the local NCBITaxa instance is created and updated via ETE3.
        """
        if self.ncbi_taxa is None:
            ssl._create_default_https_context = ssl._create_unverified_context
            self.ncbi_taxa = NCBITaxa()
            self.ncbi_taxa.update_taxonomy_database()

    def parse_names_dmp_from_taxdump(self, tar_path, f_name="names.dmp", keep_classes=None):
        """Parses the names.dmp file from the NCBI Taxonomy database dump."""
        if keep_classes is None:
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
                    if len(parts) >= 4 and parts[3] in keep_classes:
                        # parts[0] is taxid, part[1] is name, parts[3] is class
                        name2taxid[parts[1].lower()] = int(parts[0])
        return name2taxid

    def get_taxon_info_from_bt(self, taxids: list) -> dict:
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


class EntrezService:
    """Handles interactions with NCBI Entrez API."""

    def __init__(self, email):
        Entrez.email = email

    def entrez_taxon_name2taxid(self, taxon_name: str) -> dict:
        """Queries NCBI Entrez for a taxon name and returns its taxid."""
        try:
            handle = Entrez.esearch(db="taxonomy", term=taxon_name, retmode="xml", retmax=1)
            record = Entrez.read(handle)
            handle.close()
            if record["IdList"]:
                return {taxon_name: {"taxid": int(record["IdList"][0])}}
        except Exception as e:
            print(f"Entrez query failed for '{taxon_name}': {e}")
        return {}


class NCItService:
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

    async def query_taxid_from_ncit_code(self, session, ncit_code: str):
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

    async def query_taxids_from_ncit_codes(self, ncit_codes: list):
        """Map NCIT identifiers to NCBI Taxids using EBI API

        :param ncit_codes: a list of NCIT codes e.g., ["C85924", "C83526", ...]
        :return ncit2taxid: a dictionary mapping NCIT codes to taxids
        {'Trichostrongylus colubriformis': {'taxid': 6319, 'ncit': 'C125969', 'description': 'A species of parasitic...'}}
        :return notfound_ncit: a dictionary with NCIT codes failed to map taxid
        {'Trypanosoma brucei gambiense': {'ncit': 'C125975', 'description': 'A species of parasitic flagellate protozoa...'}}
        """
        ncit2taxids, notfound_ncit = {}, {}
        async with aiohttp.ClientSession() as session:
            tasks = [
                self.query_taxid_from_ncit_code(session, ncit_code) for ncit_code in ncit_codes
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

    def map_ncits_to_taxids(self, ncit_codes: list):
        """Maps NCIT codes to NCBI Taxonomy IDs."""
        ncit2taxids, notfound_ncit = asyncio.run(self.query_taxids_from_ncit_codes(ncit_codes))
        return ncit2taxids, notfound_ncit


class PubMedService:
    """Handles interactions with NCBI PubMed service via Entrez."""

    def __init__(self, email):
        Entrez.email = email

    def get_pubmed_metadata(self, pmids):
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
                    "id": f"PMID:{pmid}",
                    "pmid": int(pmid),
                    "name": title,
                    "summary": f"{abstract} [abstract]" if abstract else "",
                    "doi": doi,
                    "type": "biolink:Publication",
                }
            except (KeyError, IndexError) as e:
                print(f"Failed to parse PubMed article: {e}")
        return result


class DiseaseUtils:
    """Utilities for fetching and processing disease information."""

    def get_efo_disease_info(self, efo_file_path):
        """Reads the EFO file and returns a mapping of scientific names to disease information."""
        reader = FileReader()
        efo_map = {}
        for line in reader.read_file(efo_file_path):
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

    def bt_get_disease_info(self, ids: list) -> dict:
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


class IDMapper:
    """Handles the mapping of names to taxonomic identifiers using various services."""

    def __init__(self):
        self.ncit_service = NCItService()

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
        'clostridiales xiii': {'description': 'A bacterial family of uncertain placement in the phylum Firmicutes and the order Clostridiales that is used to classify the genus Mogibacterium.[NCIT]',
        'xrefs': {'ncit': 'C85925'},
        'id': 'NCBITaxon:189325',
        'taxid': 189325},
        'mapping_tool': 'ebi_ols'}
        """
        if ncit_path is None:
            ncit_path = os.path.join("downloads", "NCIT.txt")
        ncit_codes = self.ncit_service.get_ncit_code(ncit_path)
        ncits2taxids, notfound_ncit = self.ncit_service.map_ncits_to_taxids(ncit_codes)

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
                        "description": ncits2taxids[name]["description"]
                        if ncits2taxids[name]["description"]
                        else override_taxid.get("description", ""),
                    }
                )

        ncits2taxids["mapping_tool"] = "ebi_ols"

        return ncits2taxids


class NameMapper:
    def __init__(self, email):
        self.entrez_service = EntrezService(email)
        self.ncit_service = NCItService()
        self.ncbi_tax_service = NCBITaxonomyService()

    def ete3_taxon_name2taxid(self, taxon_names: list) -> dict:
        """Maps taxon names to NCBI Taxonomy IDs using ETE3."""
        self.ncbi_tax_service.initialize_local_ncbi_taxa()
        name2taxid = self.ncbi_tax_service.ncbi_taxa.get_name_translator(list(set(taxon_names)))
        return {
            name: {"taxid": int(taxid_list[0])}
            for name, taxid_list in name2taxid.items()
            if taxid_list
        }

    def entrez_batch_name2taxid(self, taxon_names: list, sleep=0.34) -> dict:
        """Maps taxon names to NCBI Taxonomy IDs using Entrez in batches."""
        mapped = {}
        for name in set(taxon_names):
            mapped.update(self.entrez_service.entrez_taxon_name2taxid(name))
            time.sleep(sleep)
        return mapped

    def bt_name2taxid(self, taxon_names: list) -> dict:
        """Maps taxon names to NCBI Taxonomy IDs using BioThings API."""
        client = bt.get_client("taxon")
        results = client.querymany(
            list(set(taxon_names)), scopes=["scientific_name", "other_names"], fields=["taxid"]
        )
        return {
            item["query"]: {"taxid": int(item["taxid"])}
            for item in results
            if "notfound" not in item
        }

    def fuzzy_match(
        self, query_names: list, ref_names: list, scorer=fuzz.token_sort_ratio, score_cutoff=90
    ) -> dict:
        """Performs fuzzy matching of query names against reference names."""
        matches = {}
        for name in query_names:
            match = process.extractOne(name, ref_names, scorer=scorer, score_cutoff=score_cutoff)
            if match:
                matched_name, score, _ = match
                matches[name] = {"matched_name": matched_name, "score": score}
        return matches

    def fuzzy_matched_name2taxid(self, fuzzy_matches: dict, ref_name_to_taxid: dict) -> dict:
        """Maps fuzzy matched names to NCBI Taxonomy IDs."""
        for _, match_info in fuzzy_matches.items():
            matched_name = match_info["matched_name"]
            if matched_name in ref_name_to_taxid:
                match_info["taxid"] = int(ref_name_to_taxid[matched_name])
        return fuzzy_matches


class InfoMapper:
    """Handles the mapping of preprocessed names and NCIT codes to taxon information."""

    def __init__(self, email):
        self.entrez_service = EntrezService(email)
        self.ncit_service = NCItService()
        self.ncbi_tax_service = NCBITaxonomyService()

    def get_taxid_from_cache(self, cache_data: dict) -> dict:
        """Extracts taxid from cached data."""
        return {
            original_name: taxon_map_d["taxid"]
            for original_name, taxon_map_d in cache_data.items()
            if "taxid" in taxon_map_d
        }

    def map_preprocessed_name2taxon_info(self, taxid_dict: dict, bt_taxon_info: dict) -> dict:
        """Maps preprocessed names to taxon information."""
        return {
            name: bt_taxon_info[str(taxid)]
            for name, taxid in taxid_dict.items()
            if str(taxid) in bt_taxon_info
        }

    def map_original_name2taxon_info(
        self, name_map: dict, preprocessed_name2taxon_info: dict
    ) -> dict:
        """Maps original names to taxon information based on preprocessed names."""
        original_name2taxon_info = {}
        for ori_name, pre_name in name_map.items():
            key = pre_name if pre_name in preprocessed_name2taxon_info else ori_name
            if key in preprocessed_name2taxon_info:
                original_name2taxon_info[ori_name] = {
                    **preprocessed_name2taxon_info[key],
                    "original_name": ori_name,
                }
        return original_name2taxon_info

    def map_ncit2taxon_info(self, bt_taxon_info: dict, cached_ncit2taxid: dict) -> dict:
        """Maps NCIT codes to taxon information using BioThings API."""
        ncit2taxon_info = {}
        for name, ncit_entry in cached_ncit2taxid.items():
            taxid = ncit_entry.get("taxid")
            taxon_info = bt_taxon_info.get(str(taxid)) if taxid else None
            if taxon_info:
                combined = {**ncit_entry, **taxon_info}
                if "ncit" in combined:
                    combined["xrefs"] = {"ncit": combined.pop("ncit")}
                combined.pop("mapping_tool", None)
                ncit2taxon_info[name] = combined
        return ncit2taxon_info

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
        return dict(zip(df["Source Term"], df["Mapped Term CURIE"]))

    def map_bt_disease_info(
        self, disease_name2id: dict, disease_name_map: dict, disease_info: dict
    ) -> dict:
        """Maps disease names to their information using BioThings API."""
        final_d_mapping = {}
        for old_name, new_name in disease_name_map.items():
            _id = disease_name2id.get(new_name)
            if _id and _id in disease_info:
                d_info = disease_info[_id].copy()
                d_info["original_name"] = old_name
                final_d_mapping[old_name] = d_info
        return final_d_mapping


class CacheManager(CacheHelper):
    def __init__(self, cache_dir="cache"):
        super().__init__(cache_dir)
        self.cache_dir = cache_dir
        self.id_mapper = IDMapper()

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
            return ncits2taxids


class DataCachePipeline:
    """Pipeline for caching data for the parser."""

    def __init__(self, cache_dir="cache", downloads_dir="downloads"):
        self.downloads_dir = downloads_dir
        self.cache_dir = cache_dir
        self.ncit_path = os.path.join(self.downloads_dir, "NCIT.txt")
        self.cache_manager = CacheManager(cache_dir)

    def run(self):
        """Caches the NCIT to NCBI Taxonomy ID mappings."""
        ncit2taxid_mapping = self.cache_manager.get_or_cache_ncits2taxids_mapping()
        return ncit2taxid_mapping


class MicroPhenoDBParser:
    """Orchestrates the entire data processing pipeline for MicroPhenoDB."""

    def __init__(self, data_dir="downloads"):
        self.data_dir = data_dir
        self.file_reader = FileReader()
        self.name_processor = OntologyNameProcessor()
        self.email = os.getenv("EMAIL_ADDRESS", "user@example.com")
        self.name_mapper = NameMapper(self.email)
        self.info_mapper = InfoMapper()
        self.ncbi_tax_service = NCBITaxonomyService()
        self.pubmed_service = PubMedService(self.email)
        self.disease_utils = DiseaseUtils()
        self.cache_helper = CacheHelper()

    def _get_file_path(self, filename):
        return os.path.join(self.data_dir, filename)

    def _get_organism_type(self, node) -> str:
        type_map = {
            2: "biolink:Bacterium",
            2157: "biolink:Archaea",
            2759: "biolink:Eukaryote",
            10239: "biolink:Virus",
        }
        for taxid, biolink_type in type_map.items():
            if taxid in node.get("lineage", []):
                return biolink_type
        return "biolink:OrganismalEntity"

    def _get_all_taxon_names(self, file_path):
        return {line[1].lower().strip() for line in self.file_reader.read_file(file_path) if line}

    def _get_taxon_names_to_map(self, all_names, mapped_names):
        return list(all_names - set(mapped_names.keys()))

    def _cache_and_load_data(self):
        """Manages the caching and loading of all intermediate and final data."""
        # Taxon Mapping
        if not (taxon_map := self.cache_helper.load_pickle("original_taxon_name2taxid.pkl")):
            print("Taxon map not found in cache. Generating...")
            # 1. NCIT mapping
            ncit_codes = NCItService().get_ncit_code(self._get_file_path("NCIT.txt"))
            ncit_mapped = self.name_mapper.ncit2taxid(ncit_codes)
            # 2. Other mappings
            all_taxon_names = self._get_all_taxon_names(self._get_file_path("core_table.txt"))
            names_to_map = self._get_taxon_names_to_map(all_taxon_names, ncit_mapped)
            preprocessed_map = self.name_processor.convert_preprocessed_name2dict(names_to_map)
            unique_preprocessed = list(set(preprocessed_map.values()))
            # 3. Pipeline
            pipeline = {
                "ete3": self.name_mapper.ete3_taxon_name2taxid,
                "entrez": self.name_mapper.entrez_batch_name2taxid,
                "bt": self.name_mapper.bt_name2taxid,
            }
            all_mappings = {}
            remaining = unique_preprocessed
            for _, func in pipeline.items():
                results = func(remaining)
                all_mappings.update(results)
                remaining = [name for name in remaining if name not in results]
            # 4. Fuzzy
            taxdump_path = os.path.join(os.getcwd(), "taxdump.tar.gz")
            if not os.path.exists(taxdump_path):
                self.ncbi_tax_service.download_ncbi_taxdump()
            ref_map = self.ncbi_tax_service.parse_names_dmp_from_taxdump(taxdump_path)
            fuzzy_matches = self.name_mapper.fuzzy_match(remaining, list(ref_map.keys()))
            all_mappings.update(self.name_mapper.fuzzy_matched_name2taxid(fuzzy_matches, ref_map))
            # 5. Finalize
            final_taxid_map = {
                name: info["taxid"] for name, info in all_mappings.items() if "taxid" in info
            }
            final_taxid_map.update(
                {name: info["taxid"] for name, info in ncit_mapped.items() if "taxid" in info}
            )
            taxon_info_dump = self.ncbi_tax_service.get_taxon_info_from_bt(
                list(final_taxid_map.values())
            )
            taxon_map = self.info_mapper.map_ncit2taxon_info(taxon_info_dump, ncit_mapped)
            preprocessed_info = self.info_mapper.map_preprocessed_name2taxon_info(
                final_taxid_map, taxon_info_dump
            )
            taxon_map.update(
                self.info_mapper.map_original_name2taxon_info(preprocessed_map, preprocessed_info)
            )
            self.cache_helper.save_pickle(taxon_map, "original_taxon_name2taxid.pkl")

        # Disease Mapping
        if not (disease_map := self.cache_helper.load_pickle("original_disease_name2id.pkl")):
            print("Disease map not found in cache. Generating...")
            efo_map = self.disease_utils.get_efo_disease_info(self._get_file_path("EFO.txt"))
            core_disease_names = {
                line[2].lower().strip()
                for line in self.file_reader.read_file(self._get_file_path("core_table.txt"))
                if line[2].lower() not in ["null", "not foundthogenic"]
            }
            names_to_map = list(core_disease_names - set(efo_map.keys()))
            preprocessed_map = {
                name: self.name_processor.preprocess_disease_name(name) for name in names_to_map
            }
            name2id = self.info_mapper.text2term_name2id(list(set(preprocessed_map.values())))
            disease_info_dump = self.disease_utils.bt_get_disease_info(list(set(name2id.values())))
            disease_map = efo_map
            disease_map.update(
                self.info_mapper.map_bt_disease_info(name2id, preprocessed_map, disease_info_dump)
            )
            self.cache_helper.save_pickle(disease_map, "original_disease_name2id.pkl")

        # Publication Info
        if not (pub_map := self.cache_helper.load_pickle("publication_metadata.pkl")):
            print("Publication metadata not found in cache. Generating...")
            pmids = {
                line[4]
                for line in self.file_reader.read_file(self._get_file_path("core_table.txt"))
                if re.match(r"^\d+$", line[4])
            }
            pub_map = self.pubmed_service.get_pubmed_metadata(list(pmids))
            self.cache_helper.save_pickle(pub_map, "publication_metadata.pkl")

        return taxon_map, disease_map, pub_map

    def load_microphenodb_data(self):
        """Loads and yields the final processed data records."""
        taxon_map, disease_map, pub_map = self._cache_and_load_data()
        core_f_path = self._get_file_path("core_table.txt")
        for line in self.file_reader.read_file(core_f_path):
            organism_name, disease_name, score, pmid, _, _, position, qualifier = (
                part.strip() for part in line
            )
            subject_node = self._get_microbe_node(organism_name, taxon_map)
            object_node = self._get_disease_node(disease_name, disease_map)
            if not subject_node or not object_node:
                continue
            association_node = self._get_association_node(score, position, qualifier, pmid, pub_map)
            yield {
                "_id": str(uuid.uuid4()),
                "subject": subject_node,
                "object": object_node,
                "association": association_node,
            }

    def _get_microbe_node(self, name, taxon_map):
        node = taxon_map.get(name.lower())
        if node:
            node = node.copy()
            node["type"] = self._get_organism_type(node)
            return node
        return None

    def _get_disease_node(self, name, disease_map):
        return disease_map.get(name.lower())

    def _get_association_node(self, score, position, qualifier, pmid, pub_map):
        assoc = {
            "predicate": "biolink:OrganismalEntityAsAModelOfDiseaseAssociation",
            "type": "biolink:associated_with",
            "score": float(score),
            "anatomical_entity": position.lower(),
            "aggregator_knowledge_source": "infores:MicroPhenoDB",
            "evidence_type": "ECO:0000305",
            "publication": pub_map.get(pmid, {"id": f"PMID:{pmid}", "type": "biolink:Publication"}),
        }
        if qualifier and qualifier.lower() != "tendency":
            assoc["qualifier"] = qualifier.lower()
        return assoc


class RecordCacheManager:
    """Manages caching of the final processed records."""

    def __init__(self, cache_helper: CacheHelper):
        self.cache_helper = cache_helper

    def cache_microphenodb_entire_records(
        self, records, filename_base="microphenodb_parsed_records.pkl"
    ):
        """Caches the list of records to both pickle and JSON formats."""
        print(f"Caching {len(records)} final records...")
        self.cache_helper.save_pickle(records, f"{filename_base}.pkl")
        self.cache_helper.save_json(records, f"{filename_base}.json")
        print("Final records cached successfully.")


if __name__ == "__main__":
    pipeline = DataCachePipeline()
    data = pipeline.run()
