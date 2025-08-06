r""" Implementation of the PubChem API parser class.
https://pubchemdocs.ncbi.nlm.nih.gov/pug-rest
"""

import time
from collections.abc import Iterable

import requests
from loguru import logger


class PubChemAPI:
    r"""Class for parsing PubChem API.

    Args:
        api_key (str): PubChem API key.

    Attributes:
        api_key (str): PubChem API key.
        format (str): Data format (JSON).

    """

    def __init__(
        self,
        api_key: str = "https://pubchem.ncbi.nlm.nih.gov/rest/pug",
        retrieve_format: str = "JSON",
    ):
        r"""Initialize the PubChemAPI class.
        Args:
            api_key : PubChem API key.
                Defaults to "https://pubchem.ncbi.nlm.nih.gov/rest/pug".
            retrieve_format : Data format to retrieve.
                Defaults to "JSON".
        """
        self.api_key = api_key
        self.retrieve_format = retrieve_format

    def get_cids(self, aid: int, active: bool = True) -> Iterable[int]:
        r"""Get compound IDs from PubChem API.

        Args:
            aid : PubChem Assay ID.
            active : Boolean to filter active compounds.
                Defaults to True.
        Returns:
            A list of compound IDs.
        """
        url = f"{self.api_key}/assay/aid/{str(aid)}/cids/{self.retrieve_format}"
        if active:
            url += "?cids_type=active"
        else:
            url += "?cids_type=inactive"
        response = requests.get(url, timeout=30)
        if response.status_code == 200:
            data = response.json()
            logger.debug(f"Response data: {data}")
            cids = []
            if "InformationList" in data and "Information" in data["InformationList"]:
                for entry in data["InformationList"]["Information"]:
                    if "CID" in entry:
                        cids.extend(entry["CID"])
                return list(set(cids))
        elif response.status_code == 204:
            status = "active" if active else "inactive"
            logger.warning(f"No {status} CIDs found for Assay ID: {aid}.")
        elif response.status_code == 503:
            logger.warning(
                f"Service unavailable for Assay ID: {aid}. Retrying in 5 seconds."
            )
            time.sleep(5)
            return self.get_cids(aid, active)
        else:
            print(f"Failed to fetch CIDs. Status code: {response.status_code}")

    def get_cids_info(self, cid: int, properties: list[str]) -> dict[str, str]:
        r"""Get compound information from PubChem API.

        Args:
            cid :PubChem Compound ID.
            properties : list of properties to retrieve.

        Returns:
            A dictionary with the Compound information.
        """
        if len(properties) == 0:
            logger.error("No properties provided.")
            raise ValueError("No properties provided.")
        cids_str = str(cid)
        properties_str = ",".join(properties)
        url = f"{self.api_key}/compound/cid/{cids_str}/property/{properties_str}/{self.retrieve_format}"
        response = requests.get(url, timeout=30)
        if response.status_code == 200:
            response = response.json()
            logger.debug(f"Retrieved compound information for CID {cid}.")
            logger.debug(f"Properties: {response['PropertyTable']['Properties']}")
            compound_info = {}
            if "Properties" in response["PropertyTable"]:
                properties_data = response["PropertyTable"]["Properties"][0]
                for prop in properties:
                    if str(prop) in properties_data:
                        logger.debug(f"{prop}: {properties_data[str(prop)]}")
                        compound_info[str(prop)] = properties_data[str(prop)]
                    else:
                        logger.warning(
                            f"Failed to find {prop} information for CID {cid}."
                        )
                        compound_info[str(prop)] = ""
                print(f"Compound information for CID {cid}: {compound_info}")
                return compound_info
        elif response.status_code == 204:
            logger.warning(
                f"No information found for CID {cid}. Status code: {response.status_code}"
            )
            return {}
        elif response.status_code == 503:
            # wait for 5 seconds and retry again
            logger.warning(f"Service unavailable for CID {cid}. Retrying in 5 seconds.")
            time.sleep(5)
            return self.get_cids_info(cid, properties)
        logger.warning(f"Failed to retrieve compound information for CID {cid}.")
        return {}

    def get_binding_affinity(
        self, aid: int, cid: list[int]
    ) -> dict[int, dict[str, str]]:
        r"""Get binding affinity information from PubChem API.

        Args:
            aid : PubChem Assay ID.
            cid : list of PubChem Compound IDs.

        Returns:
            A dictionary with the Binding affinity information.
        """
        url = f"{self.api_key}/assay/aid/{str(aid)}/concise/{self.retrieve_format}"
        response = requests.get(url, timeout=30)
        data = response.json()
        if response.status_code == 200:
            if "Table" in data and "Row" in data["Table"]:
                columns = data["Table"]["Columns"]["Column"]
                rows = data["Table"]["Row"]
                # Extract columns and rows with "Activity" in the column name
                activity_data = {}
                for row in rows:
                    cid_value = None
                    cid_data = {}
                    for column, cell in zip(columns, row["Cell"]):
                        if column == "CID":
                            logger.debug(f"Extracting CID value: {cell}")
                            cid_value = int(cell)
                        if "Activity" in column:
                            cid_data[column] = cell
                            logger.debug(cid_data)
                    if cid_value in cid:
                        activity_data[cid_value] = cid_data
                return activity_data
            logger.warning(
                f"Failed to retrieve binding affinity information for CID {cid}."
            )
            return {}
        elif response.status_code == 204:
            logger.warning(
                f"No binding affinity information found for CID {cid}. Status code: {response.status_code}"
            )
            return {}
        elif response.status_code == 503:
            # wait for 5 seconds and retry again
            logger.warning(f"Service unavailable for CID {cid}. Retrying in 5 seconds.")
            time.sleep(5)
            return self.get_binding_affinity(aid, cid)
        else:
            logger.warning(
                f"Failed to retrieve binding affinity information for CID {cid}."
            )
            return {}
