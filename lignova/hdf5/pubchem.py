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

    """

    def __init__(self, api_key: str, format: str = "json"):
        self.api_key = api_key
        self.format = format

    def get_smiles(self, cid: int, properties: list) -> dict:
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
        pass

    def get_binding_affinity(self, aid: int) -> dict:
        r"""Get binding affinity information from PubChem API.

        Parameters:
        ----------
            aid : int
                PubChem Assay ID.

        Returns:
        ----------
            dict: Binding affinity information.
        """
        pass
