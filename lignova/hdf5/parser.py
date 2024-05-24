r"""Class for parsing HDF5 files."""
from typing import List, Tuple

import h5py
import numpy as np


class HDF5Parser:
    r"""Class for parsing HDF5 files.

    Parameters:
    ----------
        file_path (str): Path to the HDF5 file.

    Attributes:
    ----------
        file_path (str): Path to the HDF5 file.
    """

    def __init__(self, file_path: str):
        self.file_path = file_path

    def read(self, dataset: str) -> np.ndarray | List[np.ndarray]:
        r"""Read a dataset from the HDF5 file.

        Parameters:
        ----------
            dataset:str
                Name of the dataset to read.

        Returns:
        ----------
            np.ndarray | List[np.ndarray]: Data from the dataset.
        """
        pass

    def write(self, dataset: str, data: np.ndarray | List[np.ndarray]) -> None:
        r"""Write data to the HDF5 file.

        Parameters:
        ----------
            dataset : str
                Name of the dataset to write.
            data : np.ndarray | List[np.ndarray]
                Data to write.
        """
        pass

    def read_attributes(self, dataset: str) -> Tuple[dict, dict]:
        r"""Read attributes from the dataset.

        Parameters:
        ----------
            dataset : str
                Name of the dataset to read attributes from.

        Returns:
        ----------
            Tuple[dict, dict]: Tuple containing attributes and their values.
        """
        pass

    def write_attributes(self, dataset: str, attributes: dict) -> None:
        r"""Write attributes to the dataset.

        Parameters:
        ----------
            dataset : str
                Name of the dataset to write attributes to i.e. the group.
            attributes : dict
                Attributes to write.
        """
        pass
