r"""Implements methods for reading and writing files."""

from typing import Union

import os
from tempfile import NamedTemporaryFile

import MDAnalysis as mda
import pandas as pd
from loguru import logger


def get_file_ext(file_path: str) -> str:
    r"""Get the file extension from a file path.
    Parameters
    ----------
    file_path
        Path to file.
        Returns
        -------
        File extension.
    """
    return os.path.splitext(file_path)[-1]


def write_text(
    text: str | pd.DataFrame | mda.Universe,
    write_path: Union[None, str] = None,
    file_ext: Union[None, str] = None,
) -> str:
    r"""General method to write text files.

    Parameters
    ----------
    text
        Text to write to file.
    write_path
        Path to write to file. If ``None``, then a ``NamedTemporaryFile`` will
        be created instead.
    file_ext
        Specify the file extension if ``write_path`` is ``None``.

    Returns
    -------

        Path to file that was just written.
    """
    if write_path is None:
        # Use delete=False to keep the file after closing
        with NamedTemporaryFile(
            mode="w+", encoding="utf-8", suffix=file_ext, delete=False
        ) as temp_file:
            write_path = temp_file.name
            logger.info(f"Writing to temporary file: {write_path}")
            if isinstance(text, pd.DataFrame):
                text.to_csv(temp_file, index=False, header=True)
            elif isinstance(text, mda.core.groups.AtomGroup):
                text.write(temp_file.name)
            else:
                temp_file.write(text)
    else:
        logger.info(f"Writing to file: {write_path}")
        with open(write_path, "w", encoding="utf-8") as file:
            file.write(text)

    return write_path
