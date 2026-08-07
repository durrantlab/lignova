r"""Implementation of the UniChem API parser class."""

import gzip
import io
import urllib.request
from typing import Any

import pandas as pd
from loguru import logger
from wqm.api import BaseAPI, ClientConfig
from wqm.api.errors import PermanentAPIError, TransientAPIError


class UniChemAPI(BaseAPI):
    """
    class for UniChem calls via the REST API and FTP bulk mappings.

    Tasks
        - search :  search for compounds given an ID in one source, find its equivalents in other _SOURCES.
        - connectivity : connectivity search using InChI / InChIKey / sourceID.
        - mapping" : download whole-source mapping files from the FTP, save them as Parquet
    """

    _TASKS = ["search", "connectivity", "mapping"]
    _FTP_BASE = "https://ftp.ebi.ac.uk/pub/databases/chembl/UniChem"
    _BASE_URL = "https://www.ebi.ac.uk/unichem/api/v1"
    _SOURCES = {
        "chembl": 1,
        "drugbank": 2,
        "rcsb_pdb": 3,
        "gtopdb": 4,
        "pdbe": 5,
        "chebi": 7,
        "fdasrs": 14,
        "surechembl": 15,
        "hmdb": 18,
        "pubchem": 22,
        "nmrshiftdb2": 24,
        "molport": 28,
        "bindingdb": 31,
        "comptox": 32,
        "lipidmaps": 33,
        "drugcentral": 34,
        "brenda": 37,
        "rhea": 38,
        "swisslipids": 41,
        "probes_and_drugs": 49,
        "ccdc": 50,
    }

    def __init__(
        self,
        base_url: str | None = None,
        task: str = "search",
        output_dir: str | None = None,
        **kwargs: Any,
    ):
        r"""Initialize the UniChemAPI object.

        Args:
            base_url: Override for the base URL. Falls back to _BASE_URL
            task: Task type: search, connectivity, or mapping. Defaults to search.
            output_dir: Directory to save output files. Defaults to None.
        """
        if task not in self._TASKS:
            raise ValueError(f"Invalid task. Choose {self._TASKS}")
        self.task = task
        self.output_dir = output_dir
        config = ClientConfig(headers={"Content-Type": "application/json"})
        super().__init__(self._default_url(), config=config, **kwargs)

    def _default_url(self) -> str:
        r"""Ensures the base URL is set according to the task type, defaulting to the appropriate API endpoint or FTP base URL."""
        if self.task == "search":
            return f"{self._BASE_URL}/compounds"
        elif self.task == "connectivity":
            return f"{self._BASE_URL}/connectivity"
        return self._FTP_BASE

    def _update_sources(self) -> dict[str, int]:
        """Fetch the latest _SOURCES from the API and update the _SOURCES dictionary. if the FTP fetch fails, fallback to the default _SOURCES.
        Returns:
            dict: Updated _SOURCES dictionary with source names as keys and their corresponding IDs as values.
        """
        url_base = f"{self._FTP_BASE}/data/table_dumps"
        file_path = f"{url_base}/source.tsv.gz"
        try:
            with urllib.request.urlopen(file_path) as response:
                buf = io.BytesIO(response.read())

            with gzip.GzipFile(fileobj=buf) as gz:
                df = pd.read_csv(gz, sep="\t", dtype=str)
            new_src = {}
            for _, row in df.iterrows():
                name = row["NAME"].strip().lower()
                src_id = int(row["SRC_ID"])
                new_src[name] = src_id

            UniChemAPI._SOURCES = new_src
            self._SOURCES = UniChemAPI._SOURCES
            logger.debug(f"Updated sources from FTP: {len(new_src)} sources found.")
        except Exception as e:
            logger.warning(
                f"Failed to update sources from FTP: {e},using default sources."
            )
        return self._SOURCES

    @staticmethod
    def _mapping_url(src_a: int, src_b: int) -> str:
        lo, hi = sorted((src_a, src_b))
        return (
            f"{UniChemAPI._FTP_BASE}/data/wholeSourceMapping/src_id{lo}/"
            f"src{lo}src{hi}.txt.gz"
        )

    async def _get_remote_etag(self, url: str) -> str:
        try:
            resp = await self._request("HEAD", url)
        except (PermanentAPIError, TransientAPIError):
            return ""
        return resp.headers.get("ETag", "") or resp.headers.get("Last-Modified", "")
