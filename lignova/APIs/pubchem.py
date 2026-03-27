r"""Implementation of the PubChem API parser class.
https://pubchemdocs.ncbi.nlm.nih.gov/pug-rest
"""

from loguru import logger

from .base import BaseAPI


class PubChemAPI(BaseAPI):
    r"""Class for parsing PubChem API.

    Args:
        api_key (str): PubChem API key.

    Attributes:
        api_key (str): PubChem API key.
        format (str): Data format (JSON).

    """

    _BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

    def get_cids(self, aid: int, active: bool = True) -> list[str]:
        r"""Get compound IDs from PubChem API.

        Args:
            aid : PubChem Assay ID.
            active : Boolean to filter active compounds.
                Defaults to True.
        Returns:
            A list of compound IDs.
        """
        url = f"{self.base_url}/assay/aid/{str(aid)}/cids/{self.response_format}"
        if active:
            url += "?cids_type=active"
        else:
            url += "?cids_type=inactive"
        data = self._get_json(url)
        cids: list[int] = []
        if (
            data
            and "InformationList" in data
            and "Information" in data["InformationList"]
        ):
            for entry in data["InformationList"]["Information"]:
                if "CID" in entry:
                    cids.extend(entry["CID"])
            return list(set(cids))

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
        url = f"{self.base_url}/compound/cid/{cids_str}/property/{properties_str}/{self.response_format}"
        data = self._get_json(url)
        if data is None:
            logger.warning(f"Failed to retrieve information for CID {cid}.")
            return {}
        compound_info: dict[str, str] = {}
        if "PropertyTable" in data and "Properties" in data["PropertyTable"]:
            properties_data = data["PropertyTable"]["Properties"][0]
            for prop in properties:
                if str(prop) in properties_data:
                    compound_info[str(prop)] = properties_data[str(prop)]
                else:
                    logger.warning(f"Failed to find {prop} information for CID {cid}.")
                    compound_info[str(prop)] = ""
            return compound_info

    def get_binding_affinity(
        self, aid: int, cid: list[str]
    ) -> tuple[dict[int, None | float], str]:
        r"""Get binding affinity information from PubChem API.

        Args:
            aid : PubChem Assay ID.
            cid : list of PubChem Compound IDs.

        Returns:
            A dictionary with the Binding affinity information.
            str: The type of binding affinity (e.g., IC50, Ki).
        """
        url = f"{self.base_url}/assay/aid/{str(aid)}/concise/{self.response_format}"
        data = self._get_json(url)
        if data is None:
            logger.warning(
                f"Failed to retrieve binding affinity information for CID {cid}."
            )
            return {}, None
        if "Table" in data and "Row" in data["Table"]:
            columns = data["Table"]["Columns"]["Column"]
            # from columns get the index of "CID" and "Activity Value [uM]"
            if (
                "CID" not in columns
                or "Activity Value [uM]" not in columns
                or "Activity Name" not in columns
            ):
                logger.warning(f"No CID or Activity column found for AID {aid}.")
                return {}, None
            cid_index = columns.index("CID")
            activity_index = columns.index("Activity Value [uM]")
            activity_name_index = columns.index("Activity Name")
            # Extract columns and rows with "Activity" in the column name
            activity_data = {}
            for row in data["Table"]["Row"]:
                cid_value = row["Cell"][cid_index]
                activity = row["Cell"][activity_index]
                activity_name = row["Cell"][activity_name_index]
                if cid_value in cid:
                    try:
                        activity_data[int(cid_value)] = float(activity)
                    except ValueError:
                        activity_data[int(cid_value)] = None
            return activity_data, activity_name

    def get_pubmed_id(self, aid: int) -> int | None:
        r"""Get binding affinity information from PubChem API.

        Args:
            aid : PubChem Assay ID.

        Returns:
            A list with the PubMed IDs.
        """
        url = f"{self.base_url}/assay/aid/{str(aid)}/concise/{self.response_format}"
        data = self._get_json(url)
        if data is None:
            logger.warning(f"Failed to retrieve PubMed ID information for AID {aid}.")
            return None

        if "Table" in data and "Row" in data["Table"]:
            columns = data["Table"]["Columns"]["Column"]
            if "PubMed ID" not in columns:
                logger.warning(f"No PubMed ID column found for AID {aid}.")
                return None
            else:
                pubmed_ids = columns.index("PubMed ID")
                rows = data["Table"]["Row"]
                pubmed_ids = rows[0]["Cell"][pubmed_ids]
                if pubmed_ids == "":
                    return None
                return int(pubmed_ids)
