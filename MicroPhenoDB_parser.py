import asyncio
import csv
import os

import aiohttp
import requests


def read_txt_file(in_file):
    with open(in_file) as in_f:
        reader = csv.reader(in_f, delimiter="\t")
        for line in reader:
            yield line


def get_ncit_code(in_file) -> [list, list]:
    obj = read_txt_file(in_file)
    ncit_codes = []
    names = []
    # 80 entries do not have ncit code, 582 unique ncit codes (755 with redundancy)
    # 80 + 755 = 835 which matches the source
    for line in obj:
        if len(line) > 2 and "NCIT" in line[2]:
            ncit_codes.append(line[2].split("_")[1])
        else:
            names.append(line[0])
    return ncit_codes, names


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
            name = terms.get("label")
            description = next(iter(terms.get("description")), "")
            if "NCBI_Taxon_ID" in annot:
                taxid_list = annot.get("NCBI_Taxon_ID")
                taxid = taxid_list[0]
                return name, {
                    "taxid": int(taxid),
                    "ncit": ncit_code,
                    "description": f"{description} [NCIT]", # add description src
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
    ncit_codes = set(ncit_codes)
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
    NCIT, _ = get_ncit_code(filename)
    # NCIT = ["C85924", "C83526"]
    taxids, notfound = asyncio.run(ncit2taxid(NCIT)) # 567 records in mapped
    print(len(taxids))
    new_taxids = hard_code_ncit2taxid(NCIT) # 582 records after manual mapping
    print(len(new_taxids))
    # print(len(notfound))
    # mismatch for C111133 from the original database, so need to manually change the taxid to 357276
    # Need to hard-coded for {"C111133": 357276, "C85924": 884684}
    # Only 512 taxon are shared between NCIT.txt (580) and core_table.txt (1774)
    # There are 450 shared diseases between core_table.txt (500) and EFO.txt (480)
    # In EFO file, there are
    # TODO: find shared disease name in core_table.txt and EFO.txt to make sure there is >80% overlap
