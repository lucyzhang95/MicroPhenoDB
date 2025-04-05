import requests


def ncit2taxid(ncit_code: list) -> dict[str, int]:
    """Convert NCIT code to NCBI taxid

    :param ncit_code: a list of NCIT codes
    :return: a dictionary of ncit_code to taxid
    """
    taxids_mapped = {}
    for code in ncit_code:
        iri = f"http://purl.obolibrary.org/obo/NCIT_{code}"
        url = "https://www.ebi.ac.uk/ols4/api/ontologies/ncit/terms"
        params = {"iri": iri}

        response = requests.get(url, params=params)
        if not response.ok:
            print(f"Failed to fetch NCIT:{ncit_code}")

        terms = response.json().get("_embedded", {}).get("terms", [])
        if not terms:
            print(f"No term found for NCIT:{ncit_code}")

        annotations = terms[0].get("annotation", {})
        taxid_list = annotations.get("NCBI_Taxon_ID", [])

        if taxid_list:
            taxids_mapped[code] = int(taxid_list[0])
        else:
            print(f"No NCBI Taxonomy ID found for NCIT:{ncit_code}")
    return taxids_mapped

NCTI = ["C111133", "C111204", "C112407"]
taxid = ncit2taxid(NCTI)
print(taxid)
# mismatch for C111133 from the original database, so need to manually change the taxid to 357276
taxid["C111133"] = 357276
print(taxid)