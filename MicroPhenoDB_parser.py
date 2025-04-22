import asyncio
import csv
import os
import pickle
import re
import ssl
import tarfile
import time

import aiohttp
import biothings_client as bt
import chardet
import requests
from Bio import Entrez
from ete3 import NCBITaxa
from rapidfuzz import fuzz, process

CACHE_DIR = os.path.join(os.getcwd(), "cache")
os.makedirs(CACHE_DIR, exist_ok=True)


def save_pickle(obj, f_name):
    """
    :param obj:
    :param f_name: files should only be existing in the cache directory
    :return:
    """
    f_path = os.path.join(CACHE_DIR, f_name)
    with open(f_path, "wb") as in_f:
        pickle.dump(obj, in_f)


def load_pickle(f_name):
    f_path = os.path.join(CACHE_DIR, f_name)
    if os.path.exists(f_path):
        with open(f_path, "rb") as in_f:
            return pickle.load(in_f)
    return None


def detect_encoding(in_file):
    with open(in_file, "rb") as f:
        raw = f.read(5000)
    return chardet.detect(raw)["encoding"]


def read_file(in_file, has_header=True):
    encoding = detect_encoding(in_file)
    if encoding == "ascii":
        encoding = "utf-8"
    try:
        with open(in_file, "r", encoding=encoding) as in_f:
            reader = csv.reader(in_f, delimiter="\t")
            if has_header:
                next(reader)
            for line in reader:
                yield line
    except UnicodeDecodeError as e:
        print(f"Unicode error with {encoding} on file {in_file}: {e}")


def download_ncbi_taxdump():
    ssl._create_default_https_context = ssl._create_unverified_context
    ncbi = NCBITaxa()
    ncbi.update_taxonomy_database()


def parse_names_dmp_from_taxdump(tar_path, f_name="names.dmp", keep_classes=None):
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
                name_parts = [part.strip().decode("utf-8") for part in line.strip().split(b"|")]
                if len(name_parts) < 4:
                    continue

                taxid, name_txt, _, name_class = name_parts[:4]
                if name_class in keep_classes:
                    name2taxid[name_txt.lower()] = int(taxid)
    return name2taxid


def get_ncit_code(in_file) -> list:
    """

    :param in_file:
    :return:
    ['C111133', 'C111204', 'C112407', 'C112408', 'C112409',...]
    """
    # 80 entries do not have ncit code, 582 unique ncit codes (755 with redundancy)
    # 80 + 755 = 835 which matches the source
    ncit_data = read_file(in_file, has_header=False)
    ncit_codes = [line[2].split("_")[1] for line in ncit_data if "NCIT" in line[2]]
    return ncit_codes


async def fetch_ncit_taxid(session, ncit_code: list, notfound_ncit: dict):
    iri = f"http://purl.obolibrary.org/obo/NCIT_{ncit_code}"
    url = "https://www.ebi.ac.uk/ols4/api/ontologies/ncit/terms"
    params = {"iri": iri}

    try:
        async with session.get(url, params=params, timeout=10) as resp:
            if resp.status != 200:
                print(f"Failed to connect ebi API: status {resp.status}")

            data = await resp.json()
            terms_list = data.get("_embedded").get("terms")
            if not terms_list:
                print(f"Failed to find terms for NCIT:{ncit_code}")
            terms = terms_list[0]
            annot = terms.get("annotation")
            name = terms.get("label").lower()
            description = next(iter(terms.get("description")), "")
            if "NCBI_Taxon_ID" in annot:
                taxid_list = annot.get("NCBI_Taxon_ID")
                taxid = taxid_list[0]
                return name, {
                    "taxid": int(taxid),
                    "ncit": ncit_code,
                    "description": f"{description}[NCIT]",  # add description src
                }
            else:
                notfound_ncit[name] = {"ncit": ncit_code, "description": f"{description}[NCIT]"}
                # print(f"NO NCBI Taxonomy ID found for NCIT:{ncit_code}")

    except requests.exceptions.RequestException as e:
        print(f"Failed to connect ebi API: {e}")
        return None


async def ncit2taxid(ncit_codes: list) -> [dict | dict]:
    """Map NCIT identifiers to NCBI Taxids using EBI API

    :param ncit_codes: a list of NCIT codes e.g., ["C85924", "C83526", ...]
    :return ncit2taxid: a dictionary mapping NCIT codes to taxids
    {'Trichostrongylus colubriformis': {'taxid': 6319, 'ncit': 'C125969', 'description': 'A species of parasitic...'}}
    :return notfound_ncit: a dictionary with NCIT codes failed to map taxid
    {'Trypanosoma brucei gambiense': {'ncit': 'C125975', 'description': 'A species of parasitic flagellate protozoa...'}}
    """
    ncit2taxids = {}
    notfound_ncit = {}

    async with aiohttp.ClientSession() as session:
        tasks = [fetch_ncit_taxid(session, ncit_code, notfound_ncit) for ncit_code in ncit_codes]
        results = await asyncio.gather(*tasks)

    for result in results:
        if result:
            name, taxid_info = result
            ncit2taxids[name] = taxid_info
    return ncit2taxids, notfound_ncit


def hard_code_ncit2taxid(ncit_codes: list) -> dict:
    """Manual map leftover NCIT identifiers to taxids
    There are a total of 15 NCIT identifiers need to be manually mapped
    2 key-values are removed due to non-microorganism property

    :param ncit_codes: a list of NCIT codes derived from NCIT.txt file
    :return: a dictionary mapping NCIT codes to taxids
    {'trypanosoma brucei gambiense': {
    'ncit': 'C125975',
    'description': 'A species of parasitic flagellate protozoa in the order Kinetoplastida. T. b. gambiense is transmitted by the tsetse fly and causes West African trypanosomiasis. Humans are the main reservoir for T. b. gambiense, but this species can also be found in animals.',
    'taxid': 31285}}
    """
    ncit2taxids, notfound_ncit = asyncio.run(ncit2taxid(ncit_codes))
    notfound_ncit["alpha-amylase (aspergillus oryzae)"][
        "description"
    ] = "A fungus used in East Asia to saccharify rice, sweet potato, and barley in the making of alcoholic beverages such as sake and shochu, and also to ferment soybeans for making soy sauce and miso. It is one of the different koji molds used for food fermentation.[Wikipedia]"
    notfound_ncit[
        "japanese encephalitis virus strain nakayama-nih antigen (formaldehyde inactivated)"
    ][
        "description"
    ] = "A virus from the family Flaviviridae, part of the Japanese encephalitis serocomplex of nine genetically and antigenically related viruses, some of which are particularly severe in horses, and four of which, including West Nile virus, are known to infect humans. The enveloped virus is closely related to the West Nile virus and the St. Louis encephalitis virus. The positive sense single-stranded RNA genome is packaged in the capsid which is formed by the capsid protein.[Wikipedia]"
    # Only 2 keys do not have taxids nor descriptions: Growth Hormone-Releasing Hormone Analogue and Metastatic Breast Carcinoma
    manual_taxid_map = {
        "powassan virus": 11083,
        "alpha-amylase (aspergillus oryzae)": 5062,
        "clostridiales xi": 189325,
        "trichoderma": 5543,
        "malassezia furfur": 55194,
        "clostridium cluster xvi": 543347,
        "trypanosoma brucei gambiense": 31285,
        "clostridium cluster iv": 1689151,
        "japanese encephalitis virus strain nakayama-nih antigen (formaldehyde inactivated)": 11076,
        "trichomonas vaginalis": 5722,
        "clostridiales xiii": 189325,
        "human parainfluenza Virus": 336871,
        "mycobacterium xenopi": 1789,
    }

    for name, taxid in manual_taxid_map.items():
        if name in notfound_ncit:
            notfound_ncit[name].update({"taxid": taxid})

    ncit2taxids.update(notfound_ncit)
    # manual change the taxid of bacteroides dorei, since the src mapping is wrong
    ncit2taxids["bacteroides dorei"]["taxid"] = 357276
    ncit2taxids["clostridiales xi"]["taxid"] = 884684
    return ncit2taxids


def cache_hard_code_ncit2taxid(ncit_codes: list, cache_file="ncit2taxid.pkl"):
    """
    can be used only once to cache ncit2taxid, or it will be overwritten every time
    :param ncit_codes:
    :param cache_file:
    :return: a dict with taxon names from the NCIT.txt file as keys and values as follows:
    {'clostridium cluster iv': {'ncit': 'C129412',
                               'description': 'A group of at least 8 bacterial species...
                               'taxid': 1689151,
                               'mapping_tool': 'obo'}, ...}
    """
    cached = load_pickle(cache_file)
    if cached:
        return cached
    result = hard_code_ncit2taxid(ncit_codes)
    for name in result:
        result[name]["mapping_tool"] = "obo"
    save_pickle(result, cache_file)
    return result


def get_all_taxon_names(in_file) -> list:
    """

    :param in_file:
    :return:
    ['clostridiaceae bacterium sk061', 'mycobacterium spp-rapid growers', 'alicyclobacillus',...]
    """
    core_data = read_file(in_file)
    core_taxon_names = [line[1].lower().strip() for line in core_data if line]
    return core_taxon_names


def get_taxon_names2map(in_file1, in_file2) -> dict:
    core_taxon_names = get_all_taxon_names(in_file1)
    ncit_mapped_names = load_pickle(in_file2)
    name2map = [name.strip() for name in core_taxon_names if name not in ncit_mapped_names]
    return name2map


def remove_special_char(name: str) -> str:
    name = re.sub(r"[?!#*&+]", "", name).strip()
    return name


def remove_colon4name(name: str) -> str:
    name = re.sub(r":", " ", name).strip()
    return name


def remove_dot4name_except_in_sp(name: str) -> str:
    name = re.sub(r"\b(sp|spp)\.", r"\1__dot__", name)
    numeric_matches = re.findall(r"\d+(?:\.\d+)+", name)
    for match in numeric_matches:
        protected = match.replace(".", "__dot__")
        name = name.replace(match, protected)
    name = name.replace(".", "")
    name = name.replace("__dot__", ".")
    return name.strip()


def remove_hyphen4name(name: str) -> str:
    name = name.replace("butyrate-producing", "BUTYRATEPRODUCING")
    name = re.sub(r"(?<=[a-z])-(?=\d)", " ", name)  # letter-digit
    name = re.sub(r"(?<=[a-z])-(?=(like|associated|related|positive|negative|)\b)", " ", name)
    name = name.replace("BUTYRATEPRODUCING", "butyrate-producing")
    return name.strip()


def split_name_by_slash(name: str) -> str:
    name = re.split(r"(?<=[a-zA-Z])/\s*(?=[a-zA-Z])", name)[0].strip()
    return name


def remove_parentheses(name: str) -> str:
    name = re.sub(r"\s*\(.+\)", "", name).strip()
    return name


def process_comma(name: str) -> str:
    name = re.sub(r",\s*(and\s+)?[a-z]{1,3}\b", "", name)
    name = name.split(",")[0].strip()
    return name.strip()


def remove_and_in_name(name: str) -> str:
    if " and " in name:
        name = name.split("and")[0].strip()
    return name


def remove_spp_in_name(name: str) -> str:
    name = re.sub(r"\bspps?\b.*", "", name).strip()
    return name


def remove_strain_in_taxon_name(name: str) -> str:
    name = re.sub(r"\bstrains?\b", "", name).strip()
    name = re.sub(r"\s{2,}", " ", name)
    return name


def remove_type_in_taxon_name(name: str) -> str:
    name = re.sub(r"\btypes?\b", "", name).strip()
    return name


def remove_group_in_taxon_name(name: str) -> str:
    name = re.sub(r"\bgroups?\s+[a-z]\b", "", name)  # group(s) [a-z]
    name = re.sub(r"\b(groups?|subgroup)\b$", "", name)  # sub/group(s) at the end
    name = re.sub(r"\b(groups?|subgroup)\b(?=\s+[a-z])", "", name)  # sub/group(s) in the middle
    return name.strip()


def remove_serovar_in_taxon_name(name: str) -> str:
    name = re.sub(r"\bserovars?.*\b", "", name).strip()
    return name


def remove_pre_postfix_in_taxon_name(name: str) -> str:
    name = re.sub(r"^\s*b\s+", "", name)
    name = re.sub(r"\bstains?\b", "", name)
    name = re.sub(
        r"(coagulase negative|(?:non)?hemolytic|sensu lato complexe|rapid growers|complex(?:es)?|incertae sedis)",
        "",
        name,
    )
    name = re.sub(r"\s{2,}", " ", name)
    return name.strip()


def expand_taxon_name_abbrev(name: str) -> str:
    name = re.sub(r"\be histolytica\b", "entamoeba histolytica", name)
    name = re.sub(r"\bp multocida\b", "pasteurella multocida", name)
    name = re.sub(r"\bhsv\b", "herpes simplex virus", name)
    name = re.sub(r"\bhpv\b", "human papillomavirus", name)
    name = re.sub(r"\bhiv\b", "human immunodeficiency virus", name)
    name = re.sub(r"\blcm\b", "lymphocytic choriomeningitis", name)
    name = re.sub(r"\bcluster\b", "family", name)
    name = re.sub(r"\bspirochaeta\b", "borrelia", name)
    name = re.sub(r"\s{2,}", " ", name)
    return name.strip()


def preprocess_taxon_name(name: str) -> str:
    name = name.strip()

    # character and structure cleaning
    name = remove_special_char(name)
    name = remove_colon4name(name)
    name = remove_dot4name_except_in_sp(name)
    name = remove_hyphen4name(name)
    name = remove_parentheses(name)
    name = split_name_by_slash(name)
    name = process_comma(name)

    # semantic cleaning
    name = remove_and_in_name(name)
    name = remove_strain_in_taxon_name(name)
    name = remove_type_in_taxon_name(name)
    name = remove_group_in_taxon_name(name)
    name = remove_serovar_in_taxon_name(name)
    name = remove_pre_postfix_in_taxon_name(name)
    name = remove_spp_in_name(name)
    name = expand_taxon_name_abbrev(name)

    return name.strip()


def convert_preprocessed_name2dict(names: list) -> dict:
    return {original_name: preprocess_taxon_name(original_name) for original_name in names}


def ete3_taxon_name2taxid(taxon_names: list) -> dict:
    """Use ete3 to map taxonomy names to NCBI taxonomy ids
    ete3 is good at mapping exact taxonomy names and fast without accessing API

    :param taxon_names: a list of taxon names (the values of the output of preprocess_taxon_name)
    :return: a dictionary mapping taxon names to taxid numbers
    {'human papillomavirus 11': 10580, 'veillonella sp.': 1926307, ...}
    """
    ete3_mapped = {}
    taxon_names = set(taxon_names)

    ncbi = NCBITaxa()
    ete3_name2taxid = ncbi.get_name_translator(taxon_names)
    for name, taxid in ete3_name2taxid.items():
        if taxid:
            ete3_mapped[name] = {"taxid": int(taxid[0])}
    return ete3_mapped


def cache_ete3_taxon_name2taxid(taxon_names: list, cache_file="ete3_name2taxid.pkl") -> dict:
    """

    :param taxon_names:
    :param cache_file:
    :return: a dict with first preprocessed names as keys and values as follows:
    {'varicella zoster virus': {'taxid': 10335, 'mapping_tool': 'ete3'},...}
    """
    cache = load_pickle(cache_file)
    if cache is None:
        cache = {}

    names2query = [name for name in taxon_names if name not in cache]
    if names2query:
        result = ete3_taxon_name2taxid(taxon_names)
        for name in result:
            result[name]["mapping_tool"] = "ete3"
        cache.update(result)
        save_pickle(cache, cache_file)
    return cache


def entrez_taxon_name2taxid(taxon_name: list) -> dict:
    """Map taxonomy names to NCBI taxonomy ids using entrez API
    Entrez is good at mapping reclassified taxonomy names that are outdated
    but very slow due to no batch query allowed, so expensive to query

    :param taxon_name:
    :return:
    """
    Entrez.email = "bazhang@scripps.edu"
    try:
        handle = Entrez.esearch(db="taxonomy", term=taxon_name, retmode="xml", retmax=1)
        record = Entrez.read(handle)
        handle.close()
        if record["IdList"]:
            taxid = int(record["IdList"][0])
            return {taxon_name: {"taxid": taxid}}
    except Exception as e:
        print(f"Entrez query failed for '{taxon_name}': {e}")


def entrez_batch_name2taxid(taxon_names: list, sleep=0.34) -> dict:
    entrez_mapped = {}
    for name in set(taxon_names):
        result = entrez_taxon_name2taxid(name)
        if result:
            entrez_mapped.update(result)
        time.sleep(sleep)
    return entrez_mapped


def cache_entrez_batch_name2taxid(taxon_names: list, cache_file="entrez_name2taxid.pkl") -> dict:
    """Manual function to batch query entrez API (1 at a time)

    :param taxon_names:
    :param cache_file:
    :return: a dict with first preprocessed names as keys and values are as follows:
    {'bacteriodes prevotella': {'taxid': 838, 'mapping_tool': 'entrez'},...}
    """
    cache = load_pickle(cache_file)
    if cache is None:
        cache = {}

    names2query = [name for name in taxon_names if name not in cache]
    if names2query:
        result = entrez_batch_name2taxid(names2query)
        for name in result:
            result[name]["mapping_tool"] = "entrez"
        cache.update(result)
        save_pickle(cache, cache_file)
    return cache


# TODO: Taxon name resolver (preprocess with special character only, ete3 first, entrez second, then detailed name preprocess, lastly using biothings...)
def bt_name2taxid(taxon_names: list) -> dict:
    """Map taxonomy names to NCBI taxonomy ids using BioThings API
    BioThings is okay at mapping taxonomy names with abbreviations in the "other_names" scope field
    However, some of the abbreviation mappings are wrong, e.g., "piv":{"taxid":204463,"mapping_source":"bt"}
    {204463: Comamonadaceae bacterium PIV-3D}...

    :param taxon_names:
    :return:
    {{'herpes simplex virus': 126283, 'mapping_source': 'bt'},...}
    """
    taxon_names = set(taxon_names)
    get_taxon = bt.get_client("taxon")
    taxon_info = get_taxon.querymany(
        taxon_names,
        scopes=["scientific_name", "other_names"],
        fields=["taxid"],
    )

    bte_mapped = {}
    for d in taxon_info:
        if "notfound" not in d:
            if d["query"] not in bte_mapped:
                bte_mapped[d["query"]] = {"taxid": int(d["taxid"])}
    return bte_mapped


def cache_bt_name2taxid(taxon_names: list, cache_file="bt_name2taxid.pkl") -> dict:
    """

    :param taxon_names:
    :param cache_file:
    :return: a dict with first preprocessed names as keys and values are as follows:
    {'herpes simplex virus': {'taxid': 126283, 'mapping_source': 'bt'},...}
    """
    cache = load_pickle(cache_file)
    if cache is None:
        cache = {}

    names2query = [name for name in taxon_names if name not in cache]
    if names2query:
        result = bt_name2taxid(names2query)
        for name in result:
            result[name]["mapping_source"] = "bt"
        cache.update(result)
        save_pickle(cache, cache_file)
    return cache


def fuzzy_match(
    name_query: list, name_reference: list, scorer=fuzz.token_sort_ratio, score_cutoff=None
) -> dict:
    """Map messy taxonomy names to NCBI taxonomy names using RapidFuzz
    Customizable with scorer options, cutoff score of 90 is used for MicroPhenoDB

    :param name_query:
    :param name_reference:
    :param scorer:
    :param score_cutoff:
    :return:
    {'calicivirusm': {'calicivirusm': 'calici virus', 'score': 91.66666666666666, 'mapping_tool': 'rapidfuzz'},...}
    """
    fuzzy_matched = {}
    for name in name_query:
        match = process.extractOne(name, name_reference, scorer=scorer, score_cutoff=score_cutoff)
        if match:
            matched_name, score, idx = match
            fuzzy_matched[name] = {name: matched_name, "score": score, "mapping_tool": "rapidfuzz"}
    return fuzzy_matched


def fuzzy_matched_name2taxid(fuzzy_matched_names: dict, ncbi_name_dmp: dict) -> dict:
    """Map fuzzy matched taxonomy names to NCBI taxonomy ids using NCBI taxdump.tar.gz as reference

    :param fuzzy_matched_names:
    {'calicivirusm': {'calicivirusm': 'calici virus', 'score': 91.66666666666666, 'mapping_tool': 'rapidfuzz'},...}
    :param ncbi_name_dmp:
    :return:
    {'calicivirusm': {'calicivirusm': 'calici virus', 'score': 91.66666666666666, 'mapping_tool': 'rapidfuzz', 'taxid': 1916234},...}
    """
    fuzz_name2taxid = {}
    for name, match in fuzzy_matched_names.items():
        if match[name] in ncbi_name_dmp:
            match["taxid"] = int(ncbi_name_dmp[match[name]])
            fuzz_name2taxid[name] = match
    return fuzz_name2taxid


def get_taxid_from_cache(cache_data):
    taxid_map = {}
    for original_name, taxon_map_d in cache_data.items():
        if "taxid" in taxon_map_d:
            taxid_map[original_name] = taxon_map_d["taxid"]
    return taxid_map


def get_taxon_info_from_bt(taxids) -> dict:
    """

    :param taxids:
    :return:
    '36855': {'taxid': 36855,
    'scientific_name': 'brucella canis',
    'parent_taxid': 234,
    'lineage': [36855,
                234,
                2826938,
                118882,
                356,
                28211,
                1224,
                3379134,
                2,
                131567,
                1],
    'rank': 'species'}
    }
    """
    taxids = set(taxids)
    get_taxon = bt.get_client("taxon")
    taxon_info = get_taxon.gettaxa(
        taxids, fields=["scientific_name", "parent_taxid", "lineage", "rank"]
    )

    taxon = {}
    for info in taxon_info:
        if "notfound" not in info.keys():
            taxon[info["query"]] = {
                "taxid": int(info["_id"]),
                "scientific_name": info["scientific_name"],
                "parent_taxid": int(info["parent_taxid"]),
                "lineage": info["lineage"],
                "rank": info["rank"],
            }
    return taxon


def map_preprocessed_name2taxon_info(taxid_dict: dict, bt_taxon_info: dict) -> dict:
    """

    :param taxid_dict:
    {'gulbenkiania': 397456, 'gymnopilus': 86085, 'haemophilus quentini': 123834,...}
    :param bt_taxon_info:
    '36855': {'taxid': 36855,
    'scientific_name': 'brucella canis',
    'parent_taxid': 234,
    'lineage': [36855,
                234,
                2826938,
                118882,
                356,
                28211,
                1224,
                3379134,
                2,
                131567,
                1],
    'rank': 'species'}
    }
    :return:
    {'brucella canis': {'taxid': 36855, 'scientific_name': 'brucella canis', 'parent_taxid': 234, 'lineage': [36855, 234, 2826938, 118882, 356, 28211, 1224, 3379134, 2, 131567, 1], 'rank': 'species'}}
    """
    return {
        name: bt_taxon_info[str(taxid)]
        for name, taxid in taxid_dict.items()
        if str(taxid) in bt_taxon_info
    }


def map_original_name2taxon_info(name_map: dict, preprocessed_name2taxon_info: dict) -> dict:
    """

    :param name_map:
    {'enterovirus (nonpolio)': 'enterovirus', 'haemophilus influenzae (nontypeable)': 'haemophilus influenzae', 'absidia': 'absidia',...}
    :param preprocessed_name2taxon_info:
    {'enterovirus': {'taxid': 12059, 'scientific_name': 'enterovirus', 'parent_taxid': 2946630, 'lineage': [12059, 2946630, 12058, 464095, 2732506, 2732408, 2732396, 2559587, 10239, 1], 'rank': 'genus'}}
    :return:
    {'enterovirus (nonpolio)': {'taxid': 12059, 'scientific_name': 'enterovirus', 'parent_taxid': 2946630, 'lineage': [12059, 2946630, 12058, 464095, 2732506, 2732408, 2732396, 2559587, 10239, 1], 'rank': 'genus'}}
    """
    return {
        ori_name: preprocessed_name2taxon_info[pre_name]
        for ori_name, pre_name in name_map.items()
        if pre_name in preprocessed_name2taxon_info
    }


def map_ncit2taxon_info(bt_taxon_info: dict, cached_ncit2taxid: dict) -> dict:
    """

    :param bt_taxon_info:
    :param cached_ncit2taxid:
    :return:
    {'bacteroides dorei': {'taxid': 357276, 'description': 'A species of anaerobic, gram-negative, rod shaped bacteria assigned to the phylum Bacteroidetes. This specis is non-spore forming, non-motile, does not hydrolyze esculin, is indole negative and nitrate is not reduced.[NCIT]', 'scientific_name': 'phocaeicola dorei', 'parent_taxid': 909656, 'lineage': [357276, 909656, 815, 171549, 200643, 976, 68336, 1783270, 3379134, 2, 131567, 1], 'rank': 'species', 'xrefs': {'ncit': 'C111133'}}}
    """
    ncit2taxon_info = {}
    for name, ncit2taxid in cached_ncit2taxid.items():
        if "taxid" in ncit2taxid:
            if str(ncit2taxid["taxid"]) in bt_taxon_info:
                updated_ncit2taxid = ncit2taxid.copy()
                updated_ncit2taxid.update(bt_taxon_info[str(ncit2taxid["taxid"])])
                ncit2taxon_info[name] = updated_ncit2taxid

    for name, taxon_info in ncit2taxon_info.items():
        if "taxid" and "ncit" in taxon_info:
            taxon_info["xrefs"] = {"ncit": taxon_info["ncit"]}
            taxon_info.pop("ncit")
            taxon_info.pop("mapping_tool")
    return ncit2taxon_info


if __name__ == "__main__":
    """
    in_f_ncit = os.path.join("downloads", "NCIT.txt")
    ncit_codes = get_ncit_code(in_f_ncit)
    # print(ncit_codes)
    # print(len(ncit_codes))

    # ncit2taxids = cache_hard_code_ncit2taxid(ncit_codes)  # 582 records after manual mapping
    # print(ncit2taxids)
    # print(f"Mapped NCIT taxon: {len(ncit2taxids)}")
    """

    in_f_core = os.path.join("downloads", "core_table.txt")
    total_taxon_names = get_all_taxon_names(in_f_core)
    print(f"Total taxon names: {len(total_taxon_names)}")  # 5529 with redundancy
    print(f"Total unique taxon names: {len(set(total_taxon_names))}")  # 1767 unique

    taxon_names = get_taxon_names2map(in_f_core, "ncit2taxid.pkl")
    # print(taxon_names)
    print(
        f"Total taxon names exclude NCIT covered to map with redundancy: {len(taxon_names)}"
    )  # 2450 redundant names
    print(
        f"Unique taxon names exclude NCIT covered to map: {len(set(taxon_names))}"
    )  # 1252 unique names need to be mapped

    # save_pickle(taxon_names, "MicrophenoDB_original_taxon_names.pkl")

    preprocessed_names = convert_preprocessed_name2dict(taxon_names)

    """
    names4ete3 = [new_name for old_name, new_name in preprocessed_names.items()]
    print(f"Unique names for ete3 after preprocessing: {len(set(names4ete3))}")

    # Now ete3 has 1031/ hit, 221/ unique names need to map
    # cache_ete3_mapped = cache_ete3_taxon_name2taxid(names4ete3)
    # print(f"Cached ete3 mapped: {len(cache_ete3_mapped)}")
    ete3_mapped = load_pickle("ete3_name2taxid.pkl")
    # print(ete3_mapped)
    print(f"Cached ete3 mapped: {len(ete3_mapped)}")  # 1031 names mapped

    # Now entrez has 56/169 mapped, 113/169 no hit
    names4entrez = [
        new_name for old_name, new_name in preprocessed_names.items() if new_name not in ete3_mapped
    ]
    print(f"Names to map for entrez: {len(set(names4entrez))}")  # 169 names to map for entrez
    # cache_entrez_mapped = cache_entrez_batch_name2taxid(names4entrez)
    # print(cached_entrez_mapped)
    # print(f"Cached entrez mapped: {len(set(cache_entrez_mapped))}")
    entrez_mapped = load_pickle("entrez_name2taxid.pkl")
    # print(entrez_mapped)
    print(f"Cached entrez mapped: {len(set(entrez_mapped))}")  # 56 names mapped

    # biothings: 1/113 with 1 hit, 21/113 with dup hits, 91/113 no hit
    names4bt = [
        new_name
        for old_name, new_name in preprocessed_names.items()
        if new_name not in ete3_mapped and new_name not in entrez_mapped
    ]
    print(f"Names to map for bt: {len(set(names4bt))}")  # 113 names to map for bt
    # cache_bt_mapped = cache_bt_name2taxid(names4bt)
    # print(cache_bt_mapped)
    # print(f"Cached bt mapped: {len(cache_bt_mapped)}")
    bt_mapped = load_pickle("bt_name2taxid.pkl")
    # print(bt_mapped)
    print(f"Cached bt mapped: {len(bt_mapped)}")  # 22 names mapped

    # after bt mapping, 91 no hit
    names4fuzz = [
        new_name
        for old_name, new_name in preprocessed_names.items()
        if new_name not in ete3_mapped
        and new_name not in entrez_mapped
        and new_name not in bt_mapped
    ]
    print(f"Taxon names to fuzzy match: {len(set(names4fuzz))}")  # 91 names for fuzzy match

    if not os.path.exists("cache/ncbi_taxdump.pkl"):
        tar_path = os.path.join(os.getcwd(), "taxdump.tar.gz")
        taxdump = parse_names_dmp_from_taxdump(tar_path)
        save_pickle(taxdump, "ncbi_taxdump.pkl")

    ncbi_name_dmp = load_pickle("ncbi_taxdump.pkl")
    ncbi_taxon_names = [name for name, taxon in ncbi_name_dmp.items()]

    # fuzzy-matched names: 53 mapped > 90 score cutoff
    # fuzzy_matched = fuzzy_match(names4fuzz, ncbi_taxon_names, score_cutoff=90)
    # print(f"Fuzzy matched names with a cutoff score of 90: {len(fuzzy_matched)}")

    # fuzzy_matched_taxid = fuzzy_matched_name2taxid(fuzzy_matched, ncbi_name_dmp)
    # print(fuzzy_matched_taxid)
    # save_pickle(fuzzy_matched_taxid, "rapidfuzz_name2taxid.pkl")
    cache_fuzzy_matched_taxid = load_pickle("rapidfuzz_name2taxid.pkl")
    print(f"Fuzzy matched names with a cutoff score of 90: {len(cache_fuzzy_matched_taxid)}")

    # 38 no hit
    no_hit = [name for name in names4fuzz if name not in cache_fuzzy_matched_taxid]
    print(f"No hit after preprocess, ete3, entrez, bt, and fuzzy match: {len(set(no_hit))}")
    """

    # Start querying taxids using BT and get more taxon info
    ncit_cached = load_pickle("ncit2taxid.pkl")
    ete3_cached = load_pickle("ete3_name2taxid.pkl")
    entrez_cached = load_pickle("entrez_name2taxid.pkl")
    bt_cached = load_pickle("bt_name2taxid.pkl")
    fuzz_cached = load_pickle("rapidfuzz_name2taxid.pkl")

    ncit_taxid = get_taxid_from_cache(ncit_cached)
    ete3_taxid = get_taxid_from_cache(ete3_cached)
    entrez_taxid = get_taxid_from_cache(entrez_cached)
    bt_taxid = get_taxid_from_cache(bt_cached)
    fuzz_taxid = get_taxid_from_cache(fuzz_cached)

    taxid_dicts = [ncit_taxid, ete3_taxid, entrez_taxid, bt_taxid, fuzz_taxid]
    combined_taxids = {name: taxid for d in taxid_dicts for name, taxid in d.items()}
    print(f"Combined taxids from all cache: {len(combined_taxids)}")

    taxids = [taxid for name, taxid in combined_taxids.items()]
    print(f"Combined unique taxids: {len(set(taxids))}")
    taxon_info = get_taxon_info_from_bt(taxids)
    print(f"Taxon info: {len(taxon_info)}, {len(set(taxids)) - len(taxon_info)} less")

    preprocessed_name_map = map_preprocessed_name2taxon_info(combined_taxids, taxon_info)
    # print(preprocessed_name_map)
    original_name_map = map_original_name2taxon_info(preprocessed_names, preprocessed_name_map)
    # print(original_name_map)
    # save_pickle(original_name_map, "original_name2taxon_info.pkl")

    ncit_name_map = map_ncit2taxon_info(taxon_info, ncit_cached)
