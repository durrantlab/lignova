r"""Class for parsing HDF5 files."""

from typing import List

import os

import h5py
import numpy as np
from loguru import logger

from .files import FormatManager


class HDF5Parser(FormatManager):
    r"""Class for parsing HDF5 files.

    Args:
        file_path (str): Path to the HDF5 file.

    Attributes:
    ----------
        file_path (str): Path to the HDF5 file.
    """

    def create(self) -> None:
        r"""Create an HDF5 file."""
        h5py.File(self.file_path, "w")

    def read(self, path: str) -> np.ndarray | List[np.ndarray]:
        r"""Read a dataset or a group from the HDF5 file.

        Args:
            path: Path to the dataset or group to read.
        Returns:
            A numpy array or a list of numpy arrays with the dataset.
        """

        # check if the file exists
        if not os.path.exists(self.file_path):
            raise FileNotFoundError(f"File {self.file_path} not found.")
        try:
            with h5py.File(self.file_path, "r") as hdf5_file:
                data = hdf5_file[path]
                logger.debug(f"data: {data}")
                if isinstance(data, h5py.Dataset):
                    return data[()]
                if isinstance(data, h5py.Group):
                    return list(data.keys())
        except Exception as e:
            raise e

    def write(self, data: np.ndarray | List[np.ndarray], data_scheme: str) -> None:
        r"""Write data to the HDF5 file.

        Args:
            data : Data to write.
                If the dataset already exists, the data will be appended to the dataset.
            data_scheme : Name of the dataset to write.
                If the dataset does not exist, it will be created.
        Returns:
            None
        """
        try:
            with h5py.File(self.file_path, "r+") as hdf5_file:
                if data_scheme not in hdf5_file.keys():
                    logger.warning(
                        f"Dataset {data_scheme} not found. Creating dataset."
                    )
                    hdf5_file.create_dataset(data_scheme, data=data, shape=data.shape)
                else:
                    hdf5_file[data_scheme].resize(
                        (hdf5_file[data_scheme].shape[0] + data.shape[0]), axis=0
                    )
                    hdf5_file[data_scheme][-data.shape[0] :] = data

        except Exception as e:
            raise e

    def read_attributes(self, dataset: str) -> dict:
        r"""Read attributes from the dataset.

        Args:
            dataset : Name of the dataset to read attributes from.

        Returns:
            Dictionary containing attributes and their values.
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

        Args:
            dataset : Name of the dataset to write attributes to i.e. the group.
                If the dataset does not exist, it will be created.
            attributes : Attributes to write.
        Returns:
            None
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

    def find_file_stats(self) -> None:
        r"""save the file stats to a file in the same directory as the hdf5 file"""
        # check if the file exists
        if not os.path.exists(self.file_path):
            raise FileNotFoundError(f"File {self.file_path} not found.")
        try:
            # save the file stats to a file in the same directory as the hdf5 file
            with open(
                self.file_path.replace(".hdf5", "_stats.txt"), "w", encoding="utf-8"
            ) as txt_file:
                with h5py.File(self.file_path, "r") as hdf5_file:
                    txt_file.write(f"File Path: {self.file_path}\n")
                    txt_file.write(f"Dataset Count: {len(list(hdf5_file.keys()))}\n")
                    for dataset in list(hdf5_file.keys()):
                        txt_file.write(f"Dataset: {dataset}\n")
                        txt_file.write(
                            f"Attributes: {dict(hdf5_file[dataset].attrs)}\n"
                        )
                        txt_file.write(f"Shape: {hdf5_file[dataset].shape}\n")
        except Exception as e:
            raise e

    def is_path_valid(self, path: str) -> bool:
        r"""Check if the path exists in the HDF5 file.
        Note that this is a recursive function. Thus takes a long time to run on large files.

        Args:
            path: Path to check in the HDF5 file.

        Returns:
            True if the path exists in the HDF5 file, False otherwise.
        """
        valid_path = False

        def valid_visitor(name, node):
            nonlocal valid_path
            if name == path:
                valid_path = True

        with h5py.File(self.file_path, "r", encoding="utf-8") as hdf5_file:
            hdf5_file.visititems(valid_visitor)
            return valid_path
