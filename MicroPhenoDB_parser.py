import asyncio
import csv
import os
import re

import aiohttp
import biothings_client as bt
import chardet
import requests
from ete3 import NCBITaxa


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
    obj = read_file(in_file, has_header=False)
    seen = set()
    # 80 entries do not have ncit code, 582 unique ncit codes (755 with redundancy)
    # 80 + 755 = 835 which matches the source
    for line in obj:
        if len(line) > 2 and "NCIT" in line[2]:
            ncit_code = line[2].split("_")[1]
            if ncit_code not in seen:
                seen.add(ncit_code)
                yield ncit_code


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
                    "description": f"{description} [NCIT]",  # add description src
                }
            else:
                notfound_ncit[name] = {"ncit": ncit_code, "description": description}
                # print(f"NO NCBI Taxonomy ID found for NCIT:{ncit_code}")

    except requests.exceptions.RequestException as e:
        print(f"Failed to connect ebi API: {e}")


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


def hard_code_ncit2taxid(ncit_codes):
    """Manual map leftover NCIT identifiers to taxids
    There are a total of 15 NCIT identifiers need to be manually mapped
    2 key-values are removed due to non-microorganism property

    :param ncit_codes: a list of NCIT codes
    :return: a dictionary mapping NCIT codes to taxids
    {'Trypanosoma brucei gambiense': {'ncit': 'C125975', 'description': 'A species of parasitic flagellate protozoa...'}}
    """
    ncit2taxids, notfound_ncit = asyncio.run(ncit2taxid(ncit_codes))
    notfound_ncit["alpha-amylase (aspergillus oryzae)"][
        "description"
    ] = "A fungus used in East Asia to saccharify rice, sweet potato, and barley in the making of alcoholic beverages such as sake and shōchū, and also to ferment soybeans for making soy sauce and miso. It is one of the different koji molds used for food fermentation. [Wikipedia]"
    notfound_ncit[
        "japanese encephalitis virus strain nakayama-nih antigen (formaldehyde inactivated)"
    ][
        "description"
    ] = "A virus from the family Flaviviridae, part of the Japanese encephalitis serocomplex of nine genetically and antigenically related viruses, some of which are particularly severe in horses, and four of which, including West Nile virus, are known to infect humans.[13] The enveloped virus is closely related to the West Nile virus and the St. Louis encephalitis virus. The positive sense single-stranded RNA genome is packaged in the capsid which is formed by the capsid protein. [Wikipedia]"
    # 2 keys do not have taxids nor descriptions: Growth Hormone-Releasing Hormone Analogue and Metastatic Breast Carcinoma
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
    # manual change the taxid of Bacteroides dorei, since the src mapping is wrong
    ncit2taxids["bacteroides dorei"]["taxid"] = 357276
    ncit2taxids["clostridiales xi"]["taxid"] = 884684
    return ncit2taxids


def get_taxon_names(in_file, ncit_codes):
    obj = read_file(in_file)
    mapped_taxids = hard_code_ncit2taxid(ncit_codes)
    names4map = []

    for line in obj:
        names = line[1].lower()
        if names not in mapped_taxids:
            names4map.append(names)
    return names4map


def preprocess_taxon_name(names):
    name_map = {}
    manual = {
        "butyrate-producing bacterium ssc/2": "anaerostipes hadrus",
        "neisseria weaverii": "neisseria weaveri",
        "lcm virus": "mammarenavirus choriomeningitidis",
        "erysipelatoclostridium ramosum": "thomasclavelia ramosa",
    }

    for old_name in set(names):
        if old_name in manual:
            name_map[old_name] = manual[old_name]
            continue

        name = old_name.lower().strip()
        name = re.split(r"[(/,]| and | namely | such as ", name)[0].strip()
        name = re.sub(r"\bsensu lato\b", "", name)
        name = re.sub(r"\bcomplex\b", "", name)
        name = re.sub(r"\bcluster\b.*", "", name)
        name = re.sub(r"\bgroup\b.*", "", name)
        name = re.sub(r"\bsubgroup\b.*", "", name).strip()
        name = re.sub(r"\bclade\b", "", name)
        name = re.sub(r"\bincertae sedis\b", "", name)
        name = re.sub(r"\btypes? \d+[a-zA-Z]*\b", "", name)
        name = re.sub(r"\bserovars? \w+\b", "", name)
        name = re.sub(r"\bsubsp(ecies)? \w+\b", "", name)
        name = re.sub(r"\bserogroup? \w+\b", "", name)
        name = re.sub(r"\bstrain\b", "", name)
        name = re.sub(r"\bspp\b.*", "", name).strip()
        name = name.replace("??", "").replace("?", "").replace(":", "")
        name = re.sub(r"\bb\.\s*", "", name)
        name = re.sub(r"\be\.\s*", "entamoeba ", name)
        name = re.sub(r"\bhsv[-\s]*\d*", "herpes simplex virus", name)
        name = re.sub(r"\bhpv[-\s]*\d*", "human papillomavirus", name)
        name = re.sub(r"\bhmpv\b", "human metapneumovirus", name)
        name = re.sub(r"\bebv\b", "epstein-barr virus", name)
        name = re.sub(r"\bp\.\s*", "pasteurella multocida", name)
        name = re.sub(r"\s+", " ", name).strip()

        name_map[old_name] = name

    return name_map


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
    in_f_ncit = os.path.join("downloads", "NCIT.txt")
    NCIT = [ncit for ncit in get_ncit_code(in_f_ncit)]
    # print(len(NCIT))
    #
    # # NCIT = ["C85924", "C83526"]
    # taxids, notfound = asyncio.run(ncit2taxid(NCIT))  # 567 records in mapped
    # print(taxids)
    # print(len(taxids))
    # new_taxids = hard_code_ncit2taxid(NCIT)  # 582 records after manual mapping
    # print(len(new_taxids))
    # print(len(notfound))

    in_f_core = os.path.join("downloads", "core_table.txt")
    names4map = get_taxon_names(in_f_core, NCIT)
    # print(names4map)
    # print(len(set(names4map)))  # 1259 names need to use biothings to map, 515 names mapped

    processed_name_map = preprocess_taxon_name(names4map)
    # print(processed_name_map)

    name_query = [new_name for old_name, new_name in processed_name_map.items()]
    ncbi = NCBITaxa()
    ete3_name2taxid = ncbi.get_name_translator(
        name_query
    )  # ete3 successfully mapped 1003 taxon, 256 no hit
    # print(ete3_name2taxid)
    # print(len(ete3_name2taxid))

    # 264 found unique hits, 792 has dup hits, and 203 has no hit
    name2taxid = name2taxid(name_query)
    print(name2taxid)

    # use entrez to map the 264 notfound names, 127 mapped, 137 no hit
