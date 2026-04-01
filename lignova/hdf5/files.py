r"""Implementation for a base class for file parsers."""

import os
from abc import ABC, abstractmethod
from collections.abc import Iterable

import pyarrow as pa
from loguru import logger


class FormatManager(ABC):
    r"""Abstract class for FormatManager ."""

    def __init__(self, file_path: str, schema: pa.Schema | None = None) -> None:
        r"""Initialize the FormatManager class.
        Args:
            file_path: the path to the file.
            schema: the schema of the file. Default is None.
        Returns:
            None
        """
        # check if file exists
        if not os.path.exists(file_path):
            # raise warning if file does not exis
            logger.warning(f"File {file_path} not found. Creating file.")

        self.file_path = file_path
        self.file_extension = os.path.splitext(file_path)[1]
        self.file_name = os.path.basename(file_path)
        self.file_dir = os.path.dirname(file_path)
        self.schema = schema

    @abstractmethod
    def create(self) -> None:
        r"""Create a file."""
        raise NotImplementedError()

    @abstractmethod
    def read(
        self,
        columns: str | list | None = None,
    ) -> Iterable:
        r"""Read data from a file."""
        raise NotImplementedError()

    @abstractmethod
    def write(self, data: Iterable) -> None:
        r"""Write data to a file."""
        raise NotImplementedError()
