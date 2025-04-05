import requests


def ncit2taxid(ncit_code: list) -> dict[str, int]:
    """Convert NCIT code to NCBI taxid

    :param ncit_code: a list of NCIT codes
    :return: a dictionary of ncit_code to taxid
    """
    ncit2taxid = {}
    for code in ncit_code:
        iri = f"http://purl.obolibrary.org/obo/NCIT_{code}"
        url = "https://www.ebi.ac.uk/ols4/api/ontologies/ncit/terms"
        params = {"iri": iri}

        response = requests.get(url, params=params)
        if response.status_code == requests.codes.ok:
            annot = response.json().get("_embedded").get("terms")[0].get("annotation")
            taxid = annot.get("NCBI_Taxon_ID")[0]
            if taxid:
                ncit2taxid[code] = int(taxid)
            else:
                print(f"No NCBI Taxonomy ID found for NCIT:{ncit_code}")
        else:
            print(f"Failed to fetch NCIT:{ncit_code}")
    return ncit2taxid


NCTI = ["C111133", "C111204", "C112407"]
taxids = ncit2taxid(NCTI)
print(taxids)
# mismatch for C111133 from the original database, so need to manually change the taxid to 357276
taxids["C111133"] = 357276
print(taxids)
