import asyncio
import csv
import os

import aiohttp
import chardet
import requests


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
    notfound_ncit["Alpha-Amylase (Aspergillus oryzae)"][
        "description"
    ] = "A fungus used in East Asia to saccharify rice, sweet potato, and barley in the making of alcoholic beverages such as sake and shōchū, and also to ferment soybeans for making soy sauce and miso. It is one of the different koji molds used for food fermentation. [Wikipedia]"
    notfound_ncit[
        "Japanese encephalitis Virus Strain Nakayama-NIH Antigen (Formaldehyde Inactivated)"
    ][
        "description"
    ] = "A virus from the family Flaviviridae, part of the Japanese encephalitis serocomplex of nine genetically and antigenically related viruses, some of which are particularly severe in horses, and four of which, including West Nile virus, are known to infect humans.[13] The enveloped virus is closely related to the West Nile virus and the St. Louis encephalitis virus. The positive sense single-stranded RNA genome is packaged in the capsid which is formed by the capsid protein. [Wikipedia]"
    # 2 keys do not have taxids nor descriptions: Growth Hormone-Releasing Hormone Analogue and Metastatic Breast Carcinoma
    manual_taxid_map = {
        "Powassan Virus": 11083,
        "Alpha-Amylase (Aspergillus oryzae)": 5062,
        "Clostridiales XI": 189325,
        "Trichoderma": 5543,
        "Malassezia furfur": 55194,
        "Clostridium Cluster XVI": 543347,
        "Trypanosoma brucei gambiense": 31285,
        "Clostridium Cluster IV": 1689151,
        "Japanese encephalitis Virus Strain Nakayama-NIH Antigen (Formaldehyde Inactivated)": 11076,
        "Trichomonas vaginalis": 5722,
        "Clostridiales XIII": 189325,
        "Human Parainfluenza Virus": 336871,
        "Mycobacterium xenopi": 1789,
    }

    for name, taxid in manual_taxid_map.items():
        if name in notfound_ncit:
            notfound_ncit[name].update({"taxid": taxid})

    ncit2taxids.update(notfound_ncit)
    # manual change the taxid of Bacteroides dorei, since the src mapping is wrong
    ncit2taxids["Bacteroides dorei"]["taxid"] = 357276
    # TODO: need to doublecheck if the name is "Clostridiales XI" or "Clostridiales xi"
    ncit2taxids["Clostridiales XI"]["taxid"] = 884684
    return ncit2taxids


if __name__ == "__main__":
    filename = os.path.join("downloads", "NCIT.txt")
    NCIT = [ncit for ncit in get_ncit_code(filename)]
    print(len(NCIT))

    # NCIT = ["C85924", "C83526"]
    # taxids, notfound = asyncio.run(ncit2taxid(NCIT))  # 567 records in mapped
    # print(len(taxids))
    # new_taxids = hard_code_ncit2taxid(NCIT)  # 582 records after manual mapping
    # print(len(new_taxids))
    # print(len(notfound))
    # mismatch for C111133 from the original database, so need to manually change the taxid to 357276
    # Need to hard-code for {"C111133": 357276, "C85924": 884684}
    # Only 512 taxon are shared between NCIT.txt (580) and core_table.txt (1774)
    # There are 450 shared diseases between core_table.txt (500) and EFO.txt (515)
    # Unique disease prefix: {'Orphanet', 'DOID', 'EFO', 'HP'}
    # Diseases do not have ids: 49
