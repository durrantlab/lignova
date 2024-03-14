"""Implementation of utility functions for the analysis module."""
from typing import TextIO, Union

import os
import subprocess

from loguru import logger

from ..docking.contexts import GlideContext


def interconvert_mae_sdf(
    test_file: Union[str, TextIO],
    output_filename: str,
    ntruct: Union[int, str, None] = None,
    context: GlideContext = GlideContext.get_current(),
):
    r"""Convert ligand(s) to SDF format.

    Parameters
    ----------
    test_file : Union[str, TextIO]
        Test file name or file object.
    output_filename : str
        Output file name.
    ntruct : Union[int, str, None]
        Number of structures to convert. Default 1:5 i.e the first 5 structures.
    context : GlideContext
        Docking context to run the command.
    """
    # GET THE path of the file extension using the os.path.splitext() function
    # if the file extension is .sdf, then the file is in SDF format
    # if the file extension is .mae, then the file is in MAE format
    # if the file extension is neither .sdf nor .mae, then the file format is not supported
    if not os.path.exists(test_file):
        logger.error(f"File {test_file} does not exist.")
        return
    logger.info(os.path.splitext(test_file)[1])
    if os.path.splitext(test_file)[1] == ".sdf":
        logger.info(f"Input file is in SDF format.Converting to MAE format.")
        fileformat = "-isdf"
        outformat = "-omae"
    elif (
        os.path.splitext(test_file)[1] == ".maegz"
        or os.path.splitext(test_file)[1] == ".mae"
    ):
        logger.info(f"Input file is in MAE format.Converting to SDF format.")
        fileformat = "-imae"
        outformat = "-osdf"
    else:
        logger.error(
            "Input file format not supported. Please provide a file in MAE or SDF format."
        )
        return

    command = [
        context.command + "/utilities/sdconvert",
        fileformat,
        test_file,
        outformat,
        output_filename,
    ]
    if ntruct is not None:
        command.extend(["-n", str(ntruct)])
    else:
        command.extend(["-all"])
    logger.info(f"Running command: {' '.join(command)}")
    try:
        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
        )
        stdout, stderr = process.communicate()
        if process.returncode == 0:
            logger.info("File format conversion completed")
            logger.info(f"Output:\n{stdout}")
        else:
            logger.error("File format conversion failed ")
            logger.error(f"Error Output:\n{stderr}")
            raise subprocess.CalledProcessError(process.returncode, " ".join(command))
    except Exception as e:
        logger.error(f"An error occurred during rmsd calculation: {str(e)}")
        raise e
