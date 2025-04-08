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
                    "description": description,
                }
            else:
                notfound_ncit[name] = {"ncit": ncit_code, "description": description}
                # print(f"NO NCBI Taxonomy ID found for NCIT:{ncit_code}")

    except requests.exceptions.RequestException as e:
        print(f"Failed to connect ebi API: {e}")


async def ncit2taxid(ncit_codes):
    """

    :param ncit_codes:
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


def hard_code_ncit2taxid():
    ncit2taxids, notfound_ncit = asyncio.run(ncit2taxid(NCIT))


if __name__ == "__main__":
    filename = os.path.join("downloads", "NCIT.txt")
    NCIT, _ = get_ncit_code(filename)
    # NCIT = ["C85924", "C83526"]
    taxids, notfound = asyncio.run(ncit2taxid(NCIT))
    print(notfound)
    # print(len(notfound))
    # 15 organisms NCIT cannot map to taxid
    # mismatch for C111133 from the original database, so need to manually change the taxid to 357276
    # Need to hard-coded for {"C111133": 357276, "C85924": 884684}
