r"""Class for parsing HDF5 files."""
from typing import List, Tuple

import os

import h5py
import numpy as np
from loguru import logger


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

    def create(self) -> None:
        r"""Create an HDF5 file."""
        h5py.File(self.file_path, "w")

    def read(self, dataset: str) -> np.ndarray | List[np.ndarray]:
        r"""Read a dataset from the HDF5 file.

        Parameters:
        ----------
            dataset:str
                Path to the dataset to read.

        Returns:
        ----------
            np.ndarray | List[np.ndarray]: Data from the dataset.
        """
        # check if the file exists
        if not os.path.exists(self.file_path):
            raise FileNotFoundError(f"File {self.file_path} not found.")
        try:
            with h5py.File(self.file_path, "r") as hdf5_file:
                logger.debug(f"File keys {list(hdf5_file.keys())}.")
                if dataset not in list(hdf5_file.keys()):
                    raise KeyError(f"Dataset {dataset} not found.")
                data = hdf5_file[dataset][()]
        except Exception as e:
            raise e
        return data

    def write(self, dataset: str, data: np.ndarray | List[np.ndarray]) -> None:
        r"""Write data to the HDF5 file.

        Parameters:
        ----------
            dataset : str
                Name of the dataset to write. If the dataset does not exist, it will be created.
            data : np.ndarray | List[np.ndarray]
                Data to write. If the dataset already exists, the data will be appended to the dataset.
        """
        try:
            with h5py.File(self.file_path, "r+") as hdf5_file:
                if dataset not in hdf5_file.keys():
                    logger.warning(f"Dataset {dataset} not found. Creating dataset.")
                    hdf5_file.create_dataset(dataset, data=data, shape=data.shape)
                else:
                    hdf5_file[dataset].resize(
                        (hdf5_file[dataset].shape[0] + data.shape[0]), axis=0
                    )
                    hdf5_file[dataset][-data.shape[0] :] = data

        except Exception as e:
            raise e

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
        # check if the file exists
        if not os.path.exists(self.file_path):
            raise FileNotFoundError(f"File {self.file_path} not found.")
        try:
            with h5py.File(self.file_path, "r") as hdf5_file:
                if dataset not in hdf5_file.keys():
                    raise KeyError(f"Dataset {dataset} not found.")
                attributes = dict(hdf5_file[dataset].attrs)
        except Exception as e:
            raise e
        return attributes

    def write_attributes(self, dataset: str, attributes: dict) -> None:
        r"""Write attributes to the dataset.

        Parameters:
        ----------
            dataset : str
                Name of the dataset to write attributes to i.e. the group.
                If the dataset does not exist, it will be created.
            attributes : dict
                Attributes to write.
        """
        # check if the file exists
        if not os.path.exists(self.file_path):
            raise FileNotFoundError(f"File {self.file_path} not found.")
        try:
            with h5py.File(self.file_path, "r+") as hdf5_file:
                if dataset not in hdf5_file.keys():
                    logger.warning(f"Dataset {dataset} not found. Creating dataset.")
                    hdf5_file.create_dataset(dataset, shape=len(attributes))

                for key, value in attributes.items():
                    hdf5_file[dataset].attrs[key] = value
        except Exception as e:
            raise e

    def find_file_stats(self):
        r"""save the file stats to a file in the same directory as the hdf5 file"""
        # check if the file exists
        if not os.path.exists(self.file_path):
            raise FileNotFoundError(f"File {self.file_path} not found.")
        try:
            # save the file stats to a file in the same directory as the hdf5 file
            with open(self.file_path.replace(".hdf5", "_stats.txt"), "w") as f:
                with h5py.File(self.file_path, "r") as hdf5_file:
                    f.write(f"File Path: {self.file_path}\n")
                    f.write(f"Dataset Count: {len(list(hdf5_file.keys()))}\n")
                    for dataset in list(hdf5_file.keys()):
                        f.write(f"Dataset: {dataset}\n")
                        f.write(f"Attributes: {dict(hdf5_file[dataset].attrs)}\n")
                        f.write(f"Shape: {hdf5_file[dataset].shape}\n")
        except Exception as e:
            raise e
