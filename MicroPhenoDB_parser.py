import asyncio
import csv
import os
import pickle
import re
import time

import aiohttp
import biothings_client as bt
import chardet
import requests
from Bio import Entrez
from ete3 import NCBITaxa

CACHE_DIR = os.path.join(os.getcwd(), "cache")
os.makedirs(CACHE_DIR, exist_ok=True)


def save_pickle(obj, f_name):
    f_path = os.path.join(CACHE_DIR, f_name)
    with open(f_path, "wb") as f:
        pickle.dump(obj, f)


def load_pickle(f_name):
    f_path = os.path.join(CACHE_DIR, f_name)
    if os.path.exists(f_path):
        with open(f_path, "rb") as f:
            return pickle.load(f)


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


def get_ncit_code(in_file):
    # 80 entries do not have ncit code, 582 unique ncit codes (755 with redundancy)
    # 80 + 755 = 835 which matches the source
    ncit_data = read_file(in_file, has_header=False)
    ncit_codes = [line[2].split("_")[1] for line in ncit_data if "NCIT" in line[2]]
    return ncit_codes


async def fetch_ncit_taxid(session, ncit_code, notfound_ncit) -> [dict]:
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
                notfound_ncit[name] = {"ncit": ncit_code, "description": description}
                # print(f"NO NCBI Taxonomy ID found for NCIT:{ncit_code}")

    except requests.exceptions.RequestException as e:
        print(f"Failed to connect ebi API: {e}")


async def ncit2taxid(ncit_codes) -> [dict | dict]:
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


def hard_code_ncit2taxid(ncit_codes) -> dict:
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


def cached_hard_code_ncit2taxids(ncit_codes, cache_file="ncit2taxids.pkl"):
    cached = load_pickle(cache_file)
    if cached:
        return cached
    result = hard_code_ncit2taxid(ncit_codes)
    for name in result:
        result[name]["mapping_source"] = "obo"
    save_pickle(result, cache_file)
    return result


def get_taxon_names(in_file, ncit_codes) -> list:
    obj = read_file(in_file)
    mapped_taxids = hard_code_ncit2taxid(ncit_codes)
    names4map = []

    for line in obj:
        names = line[1].lower()
        if names not in mapped_taxids:
            names4map.append(names)
    return names4map


def remove_special_char(name):
    name = re.sub(r"[?!#*&+]", "", name).strip()
    return name


def remove_colon4name(name):
    name = re.sub(r":", " ", name).strip()
    return name


def remove_dot4name_except_in_sp(name):
    name = re.sub(r"\b(sp|spp)\.", r"\1__dot__", name)
    numeric_matches = re.findall(r"\d+(?:\.\d+)+", name)
    for match in numeric_matches:
        protected = match.replace(".", "__dot__")
        name = name.replace(match, protected)
    name = name.replace(".", "")
    name = name.replace("__dot__", ".")
    return name.strip()


def remove_hyphen4name(name):
    name = name.replace("butyrate-producing", "BUTYRATEPRODUCING")
    name = re.sub(r"(?<=[a-z])-(?=\d)", " ", name)  # letter-digit
    name = re.sub(r"(?<=[a-z])-(?=(like|associated|related|positive|negative|)\b)", " ", name)
    name = name.replace("BUTYRATEPRODUCING", "butyrate-producing")
    return name.strip()


def split_name_by_slash(name):
    name = re.split(r"(?<=[a-zA-Z])/\s*(?=[a-zA-Z])", name)[0].strip()
    return name


def remove_parentheses(name):
    name = re.sub(r"\s*\(.+\)", "", name).strip()
    return name


def process_comma(name):
    name = re.sub(r",\s*(and\s+)?[a-z]{1,3}\b", "", name)
    name = name.split(",")[0].strip()
    return name.strip()


def remove_and_in_name(name):
    name = name.split("and")[0].strip()
    return name


def preprocess_taxon_name(names):
    name_map = {}
    for old_name in names:
        new_name = old_name.strip()
        new_name = remove_special_char(new_name)
        new_name = remove_colon4name(new_name)
        new_name = remove_dot4name_except_in_sp(new_name)
        new_name = remove_hyphen4name(new_name)
        new_name = remove_parentheses(new_name)
        new_name = split_name_by_slash(new_name)
        new_name = process_comma(new_name)
        new_name = remove_and_in_name(new_name)
        name_map[old_name] = new_name
    return name_map


def ete_taxon_name2taxid(taxon_names):
    """Use ete3 to map taxonomy names to NCBI taxonomy ids
    ete3 is good at mapping exact taxonomy names and fast without accessing API

    :param taxon_names: a list of taxon names (the values of the output of preprocess_taxon_name)
    :return: a dictionary mapping taxon names to taxid numbers
    {'human papillomavirus 11': 10580, 'veillonella sp.': 1926307, ...}
    """
    taxon_names = set(taxon_names)

    ncbi = NCBITaxa()
    ete3_name2taxid = ncbi.get_name_translator(taxon_names)
    ete3_mapped = {name: taxid[0] for name, taxid in ete3_name2taxid.items() if taxid}
    return ete3_mapped


def entrez_taxon_name2taxid(taxon_name):
    try:
        handle = Entrez.esearch(db="taxonomy", term=taxon_name, retmode="xml", retmax=1)
        record = Entrez.read(handle)
        handle.close()
        if record["IdList"]:
            taxid = int(record["IdList"][0])
            return {taxon_name: taxid}
    except Exception as e:
        print(f"Entrez query failed for '{taxon_name}': {e}")


def entrez_batch_name2taxid(taxon_names, sleep=0.34):
    mapping = {}
    failed = []

    for name in taxon_names:
        result, failure = entrez_taxon_name2taxid(name)
        mapping.update(result)
        if failure:
            failed.append(failure)
        time.sleep(sleep)

    return mapping, failed


# TODO: Taxon name resolver (preprocess with special character only, ete3 first, entrez second, then detailed name preprocess, lastly using biothings...)
def name2taxid(names):
    names = set(names)
    get_taxon = bt.get_client("taxon")
    taxon_info = get_taxon.querymany(
        names,
        scopes="scientific_name",
        fields=["_id", "scientific_name", "lineage", "parent_taxid", "rank"],
    )
    notfound = []
    for d in taxon_info:
        if "notfound" in d:
            notfound.append(d["query"])
    return notfound


if __name__ == "__main__":
    """
    in_f_ncit = os.path.join("downloads", "NCIT.txt")
    NCIT = get_ncit_code(in_f_ncit)
    # print(len(NCIT))

    # 567 records in mapped
    # print(taxids)
    new_taxids = cached_hard_code_ncit2taxids(NCIT)  # 582 records after manual mapping
    # print(new_taxids)
    print(f"Mapped NCIT taxon: {len(new_taxids)}")
    """

"""
    in_f_core = os.path.join("downloads", "core_table.txt")
    names4map = get_taxon_names(in_f_core, NCIT)
    # print(names4map)
    print(len(set(names4map)))  # 1259 names need to use biothings to map, 515 names mapped
    # (1196 names need to be mapped if I want to get 95% retrieval rate)

    processed_name_map = preprocess_taxon_name(names4map)
    print(f"Unique processed name map: {len(processed_name_map)}")
    # print(processed_name_map)

    # the processed name can be the same for different pre-processed names
    # need to save the old name so that I can map these back
    name_query = [new_name for old_name, new_name in processed_name_map.items()]
    print(f"Unique Names for query: {len(set(name_query))}")  # 1182 new name after preprocessing

    ncbi = NCBITaxa()
    ete3_name2taxid = ncbi.get_name_translator(
        name_query
    )  # ete3 successfully mapped 1016 taxon, 166 no hit
    ete3_mapped = {name: taxid[0] for name, taxid in ete3_name2taxid.items() if taxid}
    print(f"ete3 mapped: {len(ete3_mapped)}")
    name4entrez = [name for name in name_query if name not in ete3_mapped]
    # print(set(name4entrez))
    print(len(set(name4entrez)))  # 166 no hit

    # 206 found unique hits, 799 has dup hits, and 177 has no hit
    name2taxid = name2taxid(name_query)
    print(name2taxid)

    # use Entrez to map the 166 notfound names, 57 mapped, 109 no hit
    # biothings: 40 with 1 hit, 34 found dup hits, and 92 no hit (out of 166 names)
    # Currently mapped 1184 names ~ 94% of the retrieval rate
"""
