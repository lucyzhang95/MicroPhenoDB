import csv
import requests


def read_txt_file(in_file):
    with open(in_file) as in_f:
        reader = csv.reader(in_f, delimiter='\t')
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

        response = requests.get(url, params=params)
        if response.status_code == requests.codes.ok:
            terms =  response.json().get("_embedded").get("terms")[0]
            taxid = terms.get("annotation").get("NCBI_Taxon_ID")[0]
            name = terms.get("label")
            for desc in terms.get("description"):
                if taxid:
                    ncit2taxid[name] = {
                        "taxid": int(taxid),
                        "ncit": code,
                        "description": desc
                    }
                else:
                    print(f"No NCBI Taxonomy ID found for NCIT:{ncit_code}")
        else:
            print(f"Failed to fetch NCIT:{ncit_code}")
    return ncit2taxid


if __name__ == '__main__':
    NCTI = ["C111133", "C111204", "C112407"]
    taxids = ncit2taxid(NCTI)
    print(taxids)
    # mismatch for C111133 from the original database, so need to manually change the taxid to 357276
    # taxids["C111133"] = 357276
    # print(taxids)
