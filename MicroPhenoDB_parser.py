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


def get_ncit_code(in_file):
    # 80 entries do not have ncit code, 582 unique ncit codes (755 with redundancy)
    # 80 + 755 = 835 which matches the source
    ncit_data = read_file(in_file, has_header=False)
    ncit_codes = [line[2].split("_")[1] for line in ncit_data if "NCIT" in line[2]]
    return ncit_codes


async def fetch_ncit_taxid(session, ncit_code, notfound_ncit):
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
        return None


async def ncit2taxid(ncit_codes):
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


def cache_hard_code_ncit2taxid(ncit_codes, cache_file="ncit2taxid.pkl"):
    """
    can be used only once to cache ncit2taxid, or it will be overwritten every time
    :param ncit_codes:
    :param cache_file:
    :return: a dict with taxon names from the NCIT.txt file as keys and values as follows:
    {'clostridium cluster iv': {'ncit': 'C129412',
                               'description': 'A group of at least 8 bacterial species...
                               'taxid': 1689151,
                               'mapping_source': 'obo'}, ...}
    """
    cached = load_pickle(cache_file)
    if cached:
        return cached
    result = hard_code_ncit2taxid(ncit_codes)
    for name in result:
        result[name]["mapping_source"] = "obo"
    save_pickle(result, cache_file)
    return result


def get_all_taxon_names(in_file) -> list:
    core_data = read_file(in_file)
    core_taxon_names = [line[1].lower().strip() for line in core_data if line]
    return core_taxon_names


def get_taxon_names2map(in_file1, in_file2):
    core_taxon_names = get_all_taxon_names(in_file1)
    ncit_mapped_names = load_pickle(in_file2)
    name2map = [name.strip() for name in core_taxon_names if name not in ncit_mapped_names]
    return name2map


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
    if " and " in name:
        name = name.split("and")[0].strip()
    return name


def remove_spp_in_name(name):
    name = re.sub(r"\bspps?\b.*", "", name).strip()
    return name


def remove_strain_in_taxon_name(name):
    name = re.sub(r"\bstrains?\b", "", name).strip()
    name = re.sub(r"\s{2,}", " ", name)
    return name


def remove_type_in_taxon_name(name):
    name = re.sub(r"\btypes?\b", "", name).strip()
    return name


def remove_group_in_taxon_name(name):
    name = re.sub(r"\bgroups?\s+[a-z]\b", "", name)  # group(s) [a-z]
    name = re.sub(r"\b(groups?|subgroup)\b$", "", name)  # sub/group(s) at the end
    name = re.sub(r"\b(groups?|subgroup)\b(?=\s+[a-z])", "", name)  # sub/group(s) in the middle
    return name.strip()


def remove_serovar_in_taxon_name(name):
    name = re.sub(r"\bserovars?.*\b", "", name).strip()
    return name


def remove_pre_postfix_in_taxon_name(name):
    name = re.sub(r"^\s*b\s+", "", name)
    name = re.sub(r"\bstains?\b", "", name)
    name = re.sub(
        r"(coagulase negative|(?:non)?hemolytic|sensu lato complexe|rapid growers|complex(?:es)?|incertae sedis)",
        "",
        name,
    )
    name = re.sub(r"\s{2,}", " ", name)
    return name.strip()


def expand_taxon_name_abbrev(name):
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


def preprocess_taxon_name(name):
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


def convert_preprocessed_name2dict(names):
    return {original_name: preprocess_taxon_name(original_name) for original_name in names}


def ete3_taxon_name2taxid(taxon_names):
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
            ete3_mapped[name] = {name: int(taxid[0])}
    return ete3_mapped


def cache_ete3_taxon_name2taxid(taxon_names, cache_file="ete3_name2taxid.pkl"):
    """

    :param taxon_names:
    :param cache_file:
    :return: a dict with first preprocessed names as keys and values as follows:
    {'varicella zoster virus': {'varicella zoster virus': 10335, 'mapping_source': 'ete3'},...}
    """
    cache = load_pickle(cache_file)
    if cache is None:
        cache = {}

    names2query = [name for name in taxon_names if name not in cache]
    if names2query:
        result = ete3_taxon_name2taxid(taxon_names)
        for name in result:
            result[name]["mapping_source"] = "ete3"
        cache.update(result)
        save_pickle(cache, cache_file)
    return cache


def entrez_taxon_name2taxid(taxon_name):
    Entrez.email = "bazhang@scripps.edu"
    try:
        handle = Entrez.esearch(db="taxonomy", term=taxon_name, retmode="xml", retmax=1)
        record = Entrez.read(handle)
        handle.close()
        if record["IdList"]:
            taxid = int(record["IdList"][0])
            return {taxon_name: {taxon_name: taxid}}
    except Exception as e:
        print(f"Entrez query failed for '{taxon_name}': {e}")


def entrez_batch_name2taxid(taxon_names, sleep=0.34):
    entrez_mapped = {}
    for name in set(taxon_names):
        result = entrez_taxon_name2taxid(name)
        if result:
            entrez_mapped.update(result)
        time.sleep(sleep)
    return entrez_mapped


def cache_entrez_batch_name2taxid(taxon_names, cache_file="entrez_name2taxid.pkl"):
    """

    :param taxon_names:
    :param cache_file:
    :return: a dict with first preprocessed names as keys and values are as follows:
    {'bacteriodes prevotella': {'bacteriodes prevotella': 838, 'mapping_source': 'entrez'},...}
    """
    cache = load_pickle(cache_file)
    if cache is None:
        cache = {}

    names2query = [name for name in taxon_names if name not in cache]
    if names2query:
        result = entrez_batch_name2taxid(names2query)
        for name in result:
            result[name]["mapping_source"] = "entrez"
        cache.update(result)
        save_pickle(cache, cache_file)
    return cache


# TODO: Taxon name resolver (preprocess with special character only, ete3 first, entrez second, then detailed name preprocess, lastly using biothings...)
def bte_name2taxid(taxon_names):
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
                bte_mapped[d["query"]] = {d["query"]: int(d["taxid"])}
    return bte_mapped


def cache_bte_name2taxid(taxon_names, cache_file="bte_name2taxid.pkl"):
    """

    :param taxon_names:
    :param cache_file:
    :return: a dict with first preprocessed names as keys and values are as follows:
    {'herpes simplex virus': {'herpes simplex virus': 126283, 'mapping_source': 'bte'},...}
    """
    cache = load_pickle(cache_file)
    if cache is None:
        cache = {}

    names2query = [name for name in taxon_names if name not in cache]
    if names2query:
        result = bte_name2taxid(names2query)
        for name in result:
            result[name]["mapping_source"] = "bte"
        cache.update(result)
        save_pickle(cache, cache_file)
    return cache


def fuzzy_match(
    name_query: list, name_reference: list, scorer=fuzz.token_sort_ratio, score_cutoff=None
) -> dict:
    fuzzy_matched = {}
    for name in name_query:
        match = process.extractOne(name, name_reference, scorer=scorer, score_cutoff=score_cutoff)
        if match:
            matched_name, score, idx = match
            fuzzy_matched[name] = {name: matched_name, "score": score, "mapping_tool": "rapidfuzz"}
    return fuzzy_matched


if __name__ == "__main__":
    in_f_ncit = os.path.join("downloads", "NCIT.txt")
    ncit_codes = get_ncit_code(in_f_ncit)
    # print(ncit_codes)
    # print(len(ncit_codes))

    # 567 records in mapped
    # print(taxids)
    # ncit2taxids = hard_code_ncit2taxid(ncit_codes)
    # print(ncit2taxids)
    # new_taxids = cached_hard_code_ncit2taxid(ncit_codes)  # 582 records after manual mapping
    # print(new_taxids)
    # print(f"Mapped NCIT taxon: {len(new_taxids)}")

    in_f_core = os.path.join("downloads", "core_table.txt")
    total_taxon_names = get_all_taxon_names(in_f_core)
    print(f"Total taxon names: {len(total_taxon_names)}")  # 5529 with redundancy
    print(f"Total unique taxon names: {len(set(total_taxon_names))}")  # 1774 unique

    taxon_names = get_taxon_names2map(in_f_core, "ncit2taxid.pkl")
    # print(taxon_names)
    print(
        f"Total taxon names exclude NCIT covered to map with redundancy: {len(taxon_names)}"
    )  # 2456 redundant names
    print(
        f"Unique taxon names exclude NCIT covered to map: {len(set(taxon_names))}"
    )  # 1259 unique names need to be mapped
    # # (1196 names need to be mapped if I want to get 95% retrieval rate)

    preprocessed_names = convert_preprocessed_name2dict(taxon_names)
    preprocessed_names4map = [new_name for old_name, new_name in preprocessed_names.items()]
    print(
        f"Unique names after preprocess1: {len(set(preprocessed_names4map))}"
    )  # 1244 unique names after preprocessing

    # TODO: need to rerun cache and everything below after combine preprocess into one fn;
    #  Need to add fuzz_match fn before rerun all mapping tools
    #  or use the mapping tools first, then the no-hit names can be further processed using text2term to get taxid directly with cutoff score
    #  or use the text2term directly to get taxid with cutoff score
    # Now ete3 has 969/1244 hit, 275/1259 unique names need to map
    # ete3_mapped = ete3_taxon_name2taxid(preprocessed_names2map)
    # cache_ete3_mapped = cache_ete3_taxon_name2taxid(preprocessed_names2map)
    # print(cache_ete3_mapped)
    # print(f"Cached ete3 mapped: {len(cache_ete3_mapped)}")
    ete3_mapped = load_pickle("ete3_name2taxid.pkl")
    print(f"Cached ete3 mapped: {len(ete3_mapped)}")  # 969 mapped

    # Now entrez has 56/275 mapped, 219/275 no hit
    names4entrez = [
        new_name for old_name, new_name in preprocessed_names.items() if new_name not in ete3_mapped
    ]
    print(f"Names to map for entrez: {len(set(names4entrez))}")  # 275 names to map for entrez
    # entrez_mapped = entrez_batch_name2taxid(names4entrez)
    # cache_entrez_mapped = cache_entrez_batch_name2taxid(names4entrez)
    # print(cached_entrez_mapped)
    # print(f"Cached entrez mapped: {len(set(cache_entrez_mapped))}")
    entrez_mapped = load_pickle("entrez_name2taxid.pkl")
    print(f"Cached entrez mapped: {len(set(entrez_mapped))}")  # 56 mapped

    # biothings: 2/219 with 1 hit, 22/219 found dup hits, and 195/219 no hit (24/219 mapped)
    names4bte = [
        new_name
        for old_name, new_name in preprocessed_names.items()
        if new_name not in ete3_mapped and new_name not in entrez_mapped
    ]
    print(f"Names to map for bte: {len(set(names4bte))}")  # 219 names to map for bte
    # bte_mapped = bte_name2taxid(names4bte)
    # cache_bte_mapped = cache_bte_name2taxid(names4bte)
    # print(cache_bte_mapped)
    # print(f"Cached bte mapped: {len(cache_bte_mapped)}")
    bte_mapped = load_pickle("bte_name2taxid.pkl")
    print(f"Cached bte mapped: {len(bte_mapped)}")  # 24 mapped

    # after bte mapping, 195 no hit
    names4preprocess = [
        new_name
        for old_name, new_name in preprocessed_names.items()
        if new_name not in ete3_mapped
        and new_name not in entrez_mapped
        and new_name not in bte_mapped
    ]
    print(
        f"Taxon names to preprocess2: {len(set(names4preprocess))}"
    )  # 195 names to map after preprocess2

    # Currently mapped 1049/1259 names ~ 84% of the retrieval rate

    # tar_path = os.path.join(os.getcwd(), "taxdump.tar.gz")
    # taxdump = parse_names_dmp_from_taxdump(tar_path)
    # cache_taxdump = save_pickle(taxdump, "ncbi_taxdump.pkl")
    cached_taxdump = load_pickle("ncbi_taxdump.pkl")
