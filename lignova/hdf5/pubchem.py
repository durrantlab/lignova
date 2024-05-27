import requests
from loguru import logger

r""" Implementation of the PubChem API parser class.
https://pubchemdocs.ncbi.nlm.nih.gov/pug-rest
"""


class PubChemAPI:
    r"""Class for parsing PubChem API.

    Parameters:
    ----------
        api_key (str): PubChem API key.

    Attributes:
    ----------
        api_key (str): PubChem API key.
        format (str): Data format (JSON).

    """

    def __init__(
        self,
        api_key: str = "https://pubchem.ncbi.nlm.nih.gov/rest/pug",
        format: str = "JSON",
    ):
        self.api_key = api_key
        self.format = format

    def get_cids_info(self, cid: int, properties: list) -> dict:
        r"""Get compound information from PubChem API.

        Parameters:
        ----------
            cid : int
                PubChem Compound ID.
            properties : list
                List of properties to retrieve.

        Returns:
        ----------
            dict: Compound information.
        """
        if len(properties) == 0:
            logger.error("No properties provided.")
            raise ValueError("No properties provided.")
        cids_str = str(cid)
        properties_str = ",".join(properties)
        url = f"{self.api_key}/compound/cid/{str(cid)}/property/{properties_str}/{self.format}"
        response = requests.get(url, timeout=10)
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
                return compound_info
            else:
                logger.warning(f"Failed to find properties information for CID {cid}.")
                return {}
        else:
            logger.warning(f"Failed to retrieve compound information for CID {cid}.")
            return {}

    def get_binding_affinity(self, aid: int, cid: list[int]) -> dict:
        r"""Get binding affinity information from PubChem API.

        Parameters:
        ----------
            aid : int
                PubChem Assay ID.
            cid : List[int]
                List of PubChem Compound IDs.
        Returns:
        ----------
            dict: Binding affinity information.
        """
        url = f"{self.api_key}/assay/aid/{str(aid)}/concise/{self.format}"
        response = requests.get(url)
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
            else:
                logger.warning(
                    f"Failed to retrieve binding affinity information for CID {cid}."
                )
                return {}
        else:
            logger.warning(
                f"Failed to retrieve binding affinity information for CID {cid}."
            )
            return {}
