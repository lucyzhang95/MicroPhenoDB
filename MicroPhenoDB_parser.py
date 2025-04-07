import csv
import os
import requests
import aiohttp
import asyncio


def read_txt_file(in_file):
    with open(in_file) as in_f:
        reader = csv.reader(in_f, delimiter="\t")
        for line in reader:
            yield line


def get_ncit_code(in_file) -> list:
    obj = read_txt_file(in_file)
    ncit_codes = []
    # 80 entries do not have ncit code, 582 unique ncit codes (755 with redundancy)
    for line in obj:
        if len(line) > 2 and "NCIT" in line[2]:
            ncit_codes.append(line[2].split("_")[1])
    return ncit_codes


def ncit2taxid(ncit_code: list) -> dict:
    """Convert NCIT code to NCBI taxid

    :param ncit_code: a list of NCIT codes
    :return: a dictionary of ncit_code to taxid
    """
    ncit_code = set(ncit_code)
    ncit2taxid = {}
    for code in ncit_code:
        iri = f"http://purl.obolibrary.org/obo/NCIT_{code}"
        url = "https://www.ebi.ac.uk/ols4/api/ontologies/ncit/terms"
        params = {"iri": iri}

        try:
            response = requests.get(url, params=params, timeout=3)
            response.raise_for_status()
        except requests.RequestException as e:
            print(f"Failed to fetch NCIT code {code}: {e}")
            continue

        terms_list = response.json().get("_embedded", {}).get("terms", [])
        if not terms_list:
            print(f"Failed to find terms for NCIT: {code}")
            continue

        terms = terms_list[0]
        annot = terms.get("annotation", {})
        taxid_list = annot.get("NCBI_Taxon_ID", [])

        if not taxid_list:
            print(f"No NCBI Taxonomy ID found for NCIT: {code}")
            continue
        taxid = taxid_list[0]
        name = terms.get("label", "Unknown")
        description = next(iter(terms.get("description", [])), "")

        ncit2taxid[name] = {
            "taxid": int(taxid),
            "ncit": code,
            "description": description,
        }
    return ncit2taxid


if __name__ == "__main__":
    filename = os.path.join("downloads", "NCIT.txt")
    NCIT = get_ncit_code(filename)
    # NCIT = ["C85924", "C83526"]
    taxids = ncit2taxid(NCIT)
    print(taxids)
    # mismatch for C111133 from the original database, so need to manually change the taxid to 357276
    # taxids["C111133"] = 357276
    # print(taxids)
