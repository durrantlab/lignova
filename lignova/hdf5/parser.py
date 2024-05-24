r"""Class for parsing HDF5 files."""
from typing import List, Tuple, Union

import h5py
import numpy as np


class HDF5Parser:
    r"""Class for parsing HDF5 files.

    Args:
        file_path (str): Path to the HDF5 file.

    Attributes:
        file_path (str): Path to the HDF5 file.
    """

    def __init__(self, file_path: str):
        self.file_path = file_path
