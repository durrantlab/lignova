# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 University of Pittsburgh — Of the Commonwealth System of Higher Education
# Source: https://github.com/durrantlab/lignova

"""Implementation of the PubChem API parser class.
https://pubchemdocs.ncbi.nlm.nih.gov/pug-rest
"""

import asyncio
from collections.abc import Awaitable, Callable
from typing import Any

import httpx
from aiolimiter import AsyncLimiter
from loguru import logger
from reqadence.api import APIResponseType, BaseAPI, ClientConfig
from reqadence.api.errors import PermanentAPIError

from .model import DEFAULT_PROPERTIES, AssayInfo, CompoundProperties

PUBCHEM_BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"


class PubChemAPI(BaseAPI):
    """Asynchronous Client for PubChem REST API."""

    def __init__(
        self,
        base_url: str = PUBCHEM_BASE_URL,
        response_format: APIResponseType = APIResponseType.JSON,
        client_factory: Callable[..., httpx.AsyncClient] = httpx.AsyncClient,
        rate_limiter: AsyncLimiter | None = None,
        sleep: Callable[[float], Awaitable[None]] = asyncio.sleep,
    ) -> None:
        """Initialize the PubChemAPI client.

        Args:
            base_url : Base URL for the PubChem API. Defaults to the PubChem
                REST root [PUBCHEM_BASE_URL].
            response_format : The response format for the API, JSON or XML.
                Defaults to JSON.
            client_factory : httpx client factory. Pass a cache policy's factory
                to enable caching. Defaults to a plain httpx.AsyncClient.
            rate_limiter : AsyncLimiter instance for rate limiting. Defaults to
                5 requests per second when None.
            sleep: A callable for sleeping between retries. Defaults to
                asyncio.sleep.
        """
        config = ClientConfig(response_format=response_format)
        super().__init__(
            base_url,
            config=config,
            client_factory=client_factory,
            rate_limiter=rate_limiter,
            sleep=sleep,
        )
        self.response_format: str = str(config.response_format)

    async def get_assay(self, aid: int) -> AssayInfo | None:
        """Fetch the concise assay information for a given PubChem Assay ID (AID).

        Args:
            aid : PubChem Assay ID.

        Returns:
            An AssayInfo object containing the assay information, or None if not
            found.
        """
        url = f"assay/aid/{aid}/concise/{self.response_format}"
        try:
            data = await self._get_json(url)
        except PermanentAPIError:
            logger.warning("HTTP failed to retrieve data for Assay ID: {aid}", aid=aid)
            return None

        if data is None:
            logger.warning("No data available to retrieve for assay {aid}.", aid=aid)
            return None
        if not isinstance(data, dict):
            logger.warning(
                "Unexpected response format for assay {aid}: {t}",
                aid=aid,
                t=type(data).__name__,
            )
            return None
        if "Fault" in data:
            fault_msg = data["Fault"].get("Message", "Unknown fault")
            logger.info(
                "AID {aid} not available on PubChem: {fault_msg}",
                aid=aid,
                fault_msg=fault_msg,
            )
            return None
        return AssayInfo.from_concise(aid, data)

    async def _get_cids_info(
        self, cid: list[int], properties: list[str]
    ) -> dict[int, dict[str, Any]]:
        r"""Get compound property information from PubChem API of one or more CIDs.

        Args:
            cid : List of PubChem Compound IDs.
            properties : list of properties to retrieve.

        Returns:
            A dictionary mapping each CID to its compound information, or an
            empty dict.
        """
        if len(properties) == 0:
            logger.error("No properties provided.")
            raise ValueError("No properties provided.")
        if not cid:
            return {}
        cid_str = ",".join(str(c) for c in cid)
        properties_str = ",".join(properties)
        url = f"compound/cid/{cid_str}/property/{properties_str}/{self.response_format}"
        try:
            data = await self._get_json(url)
        except PermanentAPIError:
            logger.warning(
                "HTTP failed to retrieve information for CID {cid}.", cid=cid
            )
            return {}
        if data is None:
            logger.warning(
                "No data available to retrieve information for CID {cid}.", cid=cid
            )
            return {}
        results: dict[int, dict[str, Any]] = {}
        if "PropertyTable" in data and "Properties" in data["PropertyTable"]:
            for properties_data in data["PropertyTable"]["Properties"]:
                entry_cid = properties_data.get("CID")
                if entry_cid is None:
                    continue
                compound_info: dict[str, Any] = {}
                for prop in properties:
                    if str(prop) in properties_data:
                        compound_info[str(prop)] = properties_data[str(prop)]
                    else:
                        logger.warning(
                            "Failed to find {prop} information for CID {cid}.",
                            prop=prop,
                            cid=entry_cid,
                        )
                        compound_info[str(prop)] = None
                results[int(entry_cid)] = compound_info
        return results

    async def enrich_cid_properties(
        self, assay_data: AssayInfo, extra_properties: list[str] | None = None
    ) -> AssayInfo:
        """Fetch properties for each unique CID and enrich the assay in place.

        Args:
            assay_data : AssayInfo object containing the assay data.
            extra_properties : Additional properties to retrieve beyond the default set of properties declared in
                DEFAULT_PROPERTIES.

        Returns:
            The same AssayInfo, with each record's properties populated.
        """
        properties = list(dict.fromkeys(DEFAULT_PROPERTIES + (extra_properties or [])))
        results = await self._get_cids_info(assay_data.unique_cids, properties)
        by_cid: dict[int, CompoundProperties | None] = {
            cid: CompoundProperties.model_validate(info) if info else None
            for cid, info in results.items()
        }
        for r in assay_data.records:
            if r.cid is not None:
                r.properties = by_cid.get(r.cid)
        return assay_data

    async def get_enriched_assay(
        self, aid: int, extra_properties: list[str] | None = None
    ) -> AssayInfo | None:
        """Fetch and enrich an assay in one step.

        Args:
            aid : PubChem Assay ID.
            extra_properties : Additional properties to retrieve on top of
                DEFAULT_PROPERTIES.

        Returns:
            An enriched AssayInfo, or None if the assay could not be fetched.
        """
        assay_data = await self.get_assay(aid)
        if assay_data is None:
            return None
        return await self.enrich_cid_properties(assay_data, extra_properties)
