r"""Implementation for different Schrodinger's suplementary functions."""
from typing import Iterable

import glob
import os
import shutil
import subprocess

from loguru import logger

from ..docking.contexts import GlideContext


def get_complexes(
    input_file: str,
    context: str = GlideContext.get_current(),
):
    r"""Convert docking file to complexes format and saving in write_dir.
    Parameters
    ----------
    input_file : str
        Path to Input file name.
    context : str
        Glide context. Default is GlideContext.get_current().
    """
    filename = os.path.basename(input_file)
    command = [
        context.command + "/run",
        "pv_convert.py",
        "-mode",
        "merge",
        "-o",
        os.path.splitext(filename)[0],
        input_file,
    ]
    try:
        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
        )
        stdout, stderr = process.communicate()
        if process.returncode == 0:
            logger.info(f"Converted {os.path.basename(input_file)} to complexes format")
            # the complexes file is written to the same directory as the input file find it and return it
            logger.debug(os.path.splitext(input_file))
            # check if the input file name has pv in it
            if "pv" in os.path.splitext(input_file)[0]:
                # then the complexes file name won't have pv in it
                logger.debug(os.path.splitext(input_file)[0].replace("pv", ""))
                complexes_file = glob.glob(
                    os.path.splitext(input_file)[0].replace("_pv", "")
                    + "*_complex.maegz"
                )[0]
            else:
                complexes_file = glob.glob(
                    os.path.splitext(input_file)[0] + "*_complex.maegz"
                )[0]
            # rename it to the same name as the input file
            os.rename(
                complexes_file, os.path.splitext(input_file)[0] + "_complexes.maegz"
            )
            # move the renamed file to the context write_dir
            shutil.move(
                os.path.splitext(input_file)[0] + "_complexes.maegz",
                os.path.join(
                    context.write_dir,
                    os.path.splitext(os.path.basename(input_file))[0]
                    + "_complexes.maegz",
                ),
            )
        else:
            logger.error(
                f"Complex generation failed for {os.path.basename(input_file)}"
            )
            logger.error(f"Error Output:\n{stderr}")
            raise subprocess.CalledProcessError(process.returncode, " ".join(command))
    except Exception as e:
        logger.error(f"An error occurred during Complex generation: {str(e)}")
        raise e


def convert_to_pdb(
    input_file: str,
    context: str = GlideContext.get_current(),
    n_structures: list = [1, 2],
):
    r"""Convert docking file to pdb format.
    Parameters
    ----------
    input_file : str
        Path to Input file name.

    context : str
        Glide context. Default is GlideContext.get_current().
    n_structures : int
        Number of structures to convert. Default is 1 ,2.
    """
    command = [
        context.command + "/utilities/structconvert",
        "-use_component_dict",
        "-n",
        ",".join(map(str, n_structures)),
        input_file,
        os.path.join(
            context.write_dir,
            f"{os.path.splitext(os.path.basename(input_file))[0]}.pdb",
        ),
    ]
    try:
        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
        )
        stdout, stderr = process.communicate()
        if process.returncode == 0:
            logger.info(f"Converted {input_file} to pdb format")
            logger.info(f"Output:\n{stdout}")
        else:
            logger.error(f"Conversion failed for {input_file}")
            logger.error(f"Error Output:\n{stderr}")
            raise subprocess.CalledProcessError(process.returncode, " ".join(command))
    except Exception as e:
        logger.error(f"An error occurred during conversion: {str(e)}")
        raise e
