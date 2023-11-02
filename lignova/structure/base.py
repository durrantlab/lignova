r"""Implementation of structure class."""
from typing import Union

from abc import ABC, abstractmethod
import os

class Structure(ABC):
    r"""Base class for any physical structure."""

    def __init__(self, file_path: Union[str, None] = None, file_id: Union[str, None] = None):
        # TODO:DONE add in other file info as optional parameters.
        self.file_path = file_path
        if file_id is None and file_path is not None:
            self.file_name = os.path.basename(file_path)
            self.file_id, self.file_ext = os.path.splitext(self.file_name)
            self.file_ext = self.file_ext.lstrip(".")
        elif file_id is not None and self.file_path is not None:
            self.file_id = file_id
            self.file_ext = os.path.splitext(self.file_path)[1].lstrip(".")


    @abstractmethod
    def load(
        self,
        file_path: Union[str, None] = None,
        write: bool = False,
        write_path: Union[None, str] = None,
        pdb_id: Union[str, None] = None,
    ) -> None:
        r"""Load structural data from a variety of sources."""
        raise NotImplementedError()


class Prepared:
    def get_info(self):
        """This function gives information about the prepared structure."""
        # TODO: Write function that gets all attributes arbitrarily and creates
        # dictionary to return.
