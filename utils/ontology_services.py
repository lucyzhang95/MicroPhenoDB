import asyncio
import os
import ssl

import aiohttp
import biothings_client as bt
from Bio import Entrez
from dotenv import load_dotenv
from ete3 import NCBITaxa
from tqdm.asyncio import tqdm

load_dotenv()


class EbiTaxonomyService:
    """Handles NCIt API services and related tasks.
    This service fetches NCIT codes and their corresponding NCBI Taxonomy IDs from the EBI OLS API.
    """

    async def async_query_ebi_ncit_code_to_taxid(self, session, ncit_code: str):
        """Map NCIT identifier to NCBI Taxid using EBI API.

        :param session: aiohttp session for making requests
        :param ncit_code: a NCIT code e.g., "C125969"
        :return:
        """
        url = f"https://www.ebi.ac.uk/ols4/api/ontologies/ncit/terms?iri=http://purl.obolibrary.org/obo/NCIT_{ncit_code}"
        try:
            async with session.get(url, timeout=10) as resp:
                if resp.status != 200:
                    print(f"!!! Failed to connect ebi API: status {resp.status}")
                    return None
                data = await resp.json()
                terms = data.get("_embedded", {}).get("terms", [{}])[0]
                name = terms.get("label", "").lower()
                description = next(iter(terms.get("description", [])), "")
                annot = terms.get("annotation", {})
                if "NCBI_Taxon_ID" in annot:
                    taxid = annot["NCBI_Taxon_ID"][0]
                    return name, {
                        "id": f"NCBITaxon:{int(taxid)}",
                        "taxid": int(taxid),
                        "description": f"{description}[NCIT]" if description else "",
                        "xrefs": {"ncit": f"NCIT:{ncit_code}"},
                    }
                return name, {
                    "description": f"{description}[NCIT]" if description else "",
                    "xrefs": {"ncit": f"NCIT:{ncit_code}"},
                }
        except aiohttp.ClientError as e:
            print(f"!!! EBI API request failed for NCIT_{ncit_code}: {e}")
            return None

    async def async_query_ebi_ncit_codes_to_taxids(self, ncit_codes: list):
        """Map NCIT identifiers to NCBI Taxids using EBI API

        :param ncit_codes: a list of NCIT codes e.g., ["C85924", "C83526", ...]
        :return ncit2taxid: a dictionary mapping NCIT codes to taxids
        {'Trichostrongylus colubriformis': {'taxid': 6319, 'ncit': 'C125969', 'description': 'A species of...'}}
        :return notfound_ncit: a dictionary with NCIT codes failed to map taxid
        {'Trypanosoma brucei gambiense': {'ncit': 'C125975', 'description': 'A species of parasitic protozoa...'}}
        """
        ncit2taxids, notfound_ncit = {}, {}
        async with aiohttp.ClientSession() as session:
            tasks = [
                self.async_query_ebi_ncit_code_to_taxid(session, ncit_code)
                for ncit_code in ncit_codes
            ]
            results = await tqdm.gather(*tasks, desc="Querying taxids from NCIT codes...")
        for result in results:
            if result:
                name, info = result
                if "taxid" in info:
                    ncit2taxids[name] = info
                else:
                    notfound_ncit[name] = info
        return ncit2taxids, notfound_ncit

    def async_run_ncit_codes_to_taxids(self, ncit_codes: list):
        """Maps NCIT codes to NCBI Taxonomy IDs."""
        ncit2taxids, notfound_ncit = asyncio.run(
            self.async_query_ebi_ncit_codes_to_taxids(ncit_codes)
        )
        print("[DONE] EBI NCIT Codes to TaxID mapping completed!")
        return ncit2taxids, notfound_ncit


class ETE3TaxonomyService:
    """Handles interactions with local NCBI Taxonomy database file0s via ETE3."""

    _instance = None

    def __new__(cls, *args, **kwargs):
        if cls._instance is None:
            cls._instance = super(ETE3TaxonomyService, cls).__new__(cls)
        return cls._instance

    def __init__(self):
        if hasattr(self, "_initialized") and self._initialized:
            return

        print("-> Initializing ETE3TaxonomyService for the first time...")
        self.ncbi_taxa = NCBITaxa()

        if not os.path.exists("taxdump.tar.gz"):
            print("Taxonomy database not found. Downloading and updating...")
            ssl._create_default_https_context = ssl._create_unverified_context
            self.ncbi_taxa.update_taxonomy_database()
            print("[DONE] NCBI Taxonomy Database update complete.")
        else:
            print(">>> NCBI Taxonomy database found and loaded.")

        print("---------- ETE3TaxonomyService is now ready ----------\n")

        self._initialized = True

    def ete3_taxon_name2taxid(self, taxon_names: list) -> dict:
        """Maps taxon names to NCBI Taxonomy IDs using ETE3."""
        name2taxid = self.ncbi_taxa.get_name_translator(sorted(list(set(taxon_names))))
        print("[DONE] ETE3 Taxonomy Name to TaxID mapping completed!")
        return {
            name: {"taxid": int(taxid_list[0]), "mapping_tool": "ete3"}
            for name, taxid_list in name2taxid.items()
            if taxid_list
        }


class EntrezTaxonomyService:
    """Handles interactions with NCBI Entrez API."""

    EMAIL = os.getenv("EMAIL_ADDRESS")

    def __init__(self):
        Entrez.email = self.EMAIL
        self.semaphore = asyncio.Semaphore(3)

    async def async_query_entrez_taxon_name2taxid(
            self, taxon_name: str, max_retries=3, initial_backoff=1
    ):
        """Asynchronously queries NCBI Entrez for a taxon name."""
        for attempt in range(max_retries):
            try:
                async with self.semaphore:

                    def blocking_entrez_call():
                        """Performs a blocking call to the Entrez API."""
                        handle = Entrez.esearch(
                            db="taxonomy", term=taxon_name, retmode="xml", retmax=1
                        )
                        rec = Entrez.read(handle)
                        handle.close()
                        return rec

                    rec = await asyncio.to_thread(blocking_entrez_call)
                    await asyncio.sleep(1)

                    if rec and rec.get("IdList"):
                        return taxon_name, {
                            "taxid": int(rec["IdList"][0]),
                            "mapping_tool": "entrez",
                        }
                    else:
                        return None

            except Exception as e:
                print(f"Attempt {attempt + 1}/{max_retries} failed for '{taxon_name}': {e}")
                if attempt + 1 == max_retries:
                    break

                backoff_time = initial_backoff * (2 ** attempt)
                print(f"Retrying in {backoff_time} seconds...")
                await asyncio.sleep(backoff_time)

        print(f"!!! All retries failed for '{taxon_name}'.")
        return None

    async def async_query_entrez_taxon_names2taxids(self, taxon_names: list) -> list:
        """Runs Entrez queries for a list of names concurrently."""
        tasks = [self.async_query_entrez_taxon_name2taxid(name) for name in taxon_names]
        results = await tqdm.gather(*tasks, desc="Querying Entrez Taxonomy Names...")
        return dict(res for res in results if res)

    def async_run_entrez_taxon_names2taxids(self, taxon_names: list) -> dict:
        results = asyncio.run(self.async_query_entrez_taxon_names2taxids(taxon_names))
        print("[DONE] Entrez Taxonomy Name to TaxID mapping completed!")
        return results


class BioThingsService:
    """Handles interactions with BioThings API for taxonomic information."""

    def query_bt_taxon_name2taxid(self, taxon_names: list) -> dict:
        """Maps taxon names to NCBI Taxonomy IDs using BioThings API."""
        client = bt.get_client("taxon")
        results = client.querymany(
            list(set(taxon_names)), scopes=["scientific_name", "other_names"], fields=["taxid"]
        )
        return {
            item["query"]: {"taxid": int(item["taxid"]), "mapping_tool": "biothings"}
            for item in results
            if "notfound" not in item
        }

    def query_bt_taxon_info(self, taxids: list) -> dict:
        """Fetches taxon information from BioThings API."""
        if not taxids:
            return {}
        client = bt.get_client("taxon")
        infos = client.gettaxa(
            list(set(taxids)), fields=["scientific_name", "parent_taxid", "lineage", "rank"]
        )
        return {
            info["_id"]: {
                "id": f"NCBITaxon:{int(info['_id'])}",
                "taxid": int(info["_id"]),
                "name": info.get("scientific_name"),
                "parent_taxid": info.get("parent_taxid"),
                "lineage": info.get("lineage"),
                "rank": info.get("rank"),
            }
            for info in infos
            if "notfound" not in info
        }

    def query_bt_disease_info(self, ids: list) -> dict:
        """Fetches disease information from BioThings API using Mondo ontology."""
        if not ids:
            return {}
        client = bt.get_client("disease")
        results = client.querymany(
            list(set(ids)),
            scopes=["mondo.xrefs.hp", "mondo.xrefs.efo", "mondo.mondo"],
            fields=["mondo.definition", "mondo.label"],
        )
        d_info = {}
        for info in results:
            if "notfound" not in info:
                mondo = info.get("mondo", {})
                d_info[info["query"]] = {
                    "id": info["_id"],
                    "name": mondo.get("label"),
                    "description": mondo.get("definition"),
                }
        return d_info


class PubMedService:
    """Handles interactions with NCBI PubMed service via Entrez."""

    def __init__(self, email):
        Entrez.email = email

    def query_pubmed_metadata(self, pmids):
        """Fetches PubMed metadata for a list of PMIDs."""
        if not pmids:
            return {}
        handle = Entrez.efetch(db="pubmed", id=",".join(map(str, set(pmids))), retmode="xml")
        records = Entrez.read(handle)
        handle.close()
        result = {}
        for article in records.get("PubmedArticle", []):
            try:
                pmid = str(article["MedlineCitation"]["PMID"])
                article_data = article["MedlineCitation"]["Article"]
                title = article_data.get("ArticleTitle", "")
                abstract = " ".join(
                    str(p) for p in article_data.get("Abstract", {}).get("AbstractText", [])
                )
                doi = next(
                    (
                        str(el)
                        for el in article.get("PubmedData", {}).get("ArticleIdList", [])
                        if el.attributes.get("IdType") == "doi"
                    ),
                    "",
                )
                result[pmid] = {
                    "pmid": f"PMID:{int(pmid)}",
                    "name": title,
                    "summary": f"{abstract} [abstract]" if abstract else "",
                    "doi": f"doi:{doi}" if doi else "",
                    "type": "biolink:Publication",
                }
            except (KeyError, IndexError) as e:
                print(f"!! Failed to parse PubMed article: {e}")
        return result


class UberonService:
    """Handles anatomy EBI OLS service."""

    async def async_query_anatomical_entity_to_uberon_id(
            self, session, term, url, match_type="exact", ontology="uberon", rows=1
    ):
        """Async helper function that now accepts a URL."""
        params = {"q": term, "ontology": ontology, "rows": rows, "start": 0}

        if match_type == "exact":
            params["exact"] = "true"
        elif match_type == "partial":
            params["exact"] = "false"
            params["queryFields"] = "label"
        elif match_type == "fuzzy":
            params["q"] = term + "~"
            params["exact"] = "false"
        else:
            raise ValueError("match_type must be 'exact', 'partial', or 'fuzzy'.")

        try:
            async with session.get(url, params=params) as resp:
                resp.raise_for_status()
                data = await resp.json()

                for doc in data.get("response", {}).get("docs", []):
                    iri = doc.get("iri")
                    if iri and "UBERON_" in iri:
                        uberon_id = iri.split("/")[-1].replace("_", ":")

                        return term, {
                            "id": uberon_id,
                            "name": doc.get("label"),
                            "original_name": term,
                            "type": "biolink:AnatomicalEntity",
                        }
        except aiohttp.ClientError as e:
            print(f"!! Error fetching term '{term}': {e}")
        return term, None

    async def async_anatomical_entities_to_uberon_ids(
            self, terms, match_type="exact", base_url="https://www.ebi.ac.uk/ols4/api/search"
    ):
        """
        Query OLS for a list of terms, using the provided base_url.
        """
        async with aiohttp.ClientSession() as session:
            tasks = [
                self.async_query_anatomical_entity_to_uberon_id(
                    session, term, url=base_url, match_type=match_type
                )
                for term in terms
            ]

            all_results = await tqdm.gather(*tasks, desc="Querying UBERON IDs...")

            return {term: info for term, info in all_results if info is not None}

    def async_run_anatomical_entities_to_uberon_ids(self, terms, match_type="exact"):
        base_url = "https://www.ebi.ac.uk/ols4/api/search"
        results = asyncio.run(
            self.async_anatomical_entities_to_uberon_ids(
                terms, match_type=match_type, base_url=base_url
            )
        )
        print("[DONE] UBERON ID mapping completed!")
        return results
