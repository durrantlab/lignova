r"""Implementation for different Schrodinger's suplementary functions."""
from typing import Iterable, Union

import os
import subprocess

from loguru import logger

from ..docking.contexts import GlideContext


def manipulate_complexes(
    input_file: str,
    context: str = GlideContext.get_current(),
    outfile_name: Union[None, str] = "manipulate_outp.maegz",
    mode: Union[None, str] = "merge",
):
    r"""Convert docking file to complexes format and saving in write_dir.
    Parameters
    ----------
    input_file : str
        Path to Input file name.
    context : str
        Glide context. Default is GlideContext.get_current().
    outfile_name : str
        Output file name. Default is input file name with "manipulate_outp.maegz" as suffix.
    mode : str
        Mode to convert file. Default is "merge". Options are ["merge" - combine PV/EPV structures into complexes,
        ,"split_pv" - extract receptor and ligands from complexes and save as PV
        ,"split_epv" - extract receptor and ligands from complexes and save as EPV
        ,"split_ligand" - extract ligands from complexes
        ,"split_receptor" - extract receptors from complexes
        ,"pv_to_epv" - join multiple PV files into single EPV file
        ,"epv_to_pv" - split EPV file into multiple PV files ]
    """
    # check if "out" in the outfile_name and if so throw an error
    if "-out" in outfile_name:
        logger.error("Output file name cannot contain '-out'")
        raise ValueError("Output file name cannot contain '-out'")
    filename = os.path.basename(input_file)
    if "_pv" in os.path.splitext(filename)[0]:
        filename = os.path.splitext(filename)[0].replace("_pv", "")
    elif "_epv" in os.path.splitext(filename)[0]:
        filename = os.path.splitext(filename)[0].replace("_epv", "")
    else:
        filename = os.path.splitext(filename)[0]
    # check if mode is not from the list of options
    if mode not in [
        "merge",
        "split_pv",
        "split_epv",
        "split_ligand",
        "split_receptor",
        "pv_to_epv",
        "epv_to_pv",
    ]:
        logger.error(f"Mode {mode} not in the list of options")
        raise ValueError(f"Mode {mode} not in the list of options")
    if outfile_name != "manipulate_outp.maegz":
        new_filename = os.path.join(context.write_dir, outfile_name)
    else:
        new_filename = os.path.join(context.write_dir, filename + "_" + outfile_name)
    # check if new_filename exists and if so append a number to the filename
    if os.path.exists(new_filename):
        i = 1
        while os.path.exists(new_filename):
            new_filename = os.path.join(
                context.write_dir, filename + f"_{i}_" + outfile_name
            )
            i += 1
    command = [
        context.command + "/run",
        "pv_convert.py",
        "-mode",
        mode,
        "-o",
        new_filename,
        input_file,
    ]
    logger.debug(f"Running command: {' '.join(command)}")
    try:
        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
        )
        stdout, stderr = process.communicate()
        logger.debug(f"Output:\n{stdout}")
        print((all("skipping" in line.lower() for line in stdout.split("\n")[:-1])))
        # check if the filename has _pv or _epv and if so remove it
        if process.returncode == 0 and (
            (all("skipping" in line.lower() for line in stdout.split("\n")[:-1]))
            is False
            or len(stdout.split("\n")) > 0
        ):
            logger.info(f"Converted {os.path.basename(input_file)}")
            # Find the generated complexes file
            if mode not in ["split_epv", "pv_to_epv", "epv_to_pv"]:
                complexes_file = None
                for file in os.listdir(os.path.dirname(input_file)):
                    logger.debug(f"Checking file: {file}")
                    if (
                        file.startswith(filename)
                        and file.endswith(".maegz")
                        and file != os.path.basename(input_file)
                        and "-out" in file
                    ):
                        complexes_file = os.path.join(os.path.dirname(input_file), file)
                        print(complexes_file)
                        break
                if complexes_file:
                    # Rename the complexes file to match the input file name
                    os.rename(complexes_file, new_filename)
                    logger.info(f"Converted file saved at: {new_filename}")
                else:
                    logger.error(
                        f"Failed to find the generated complexes file for {os.path.basename(input_file)}"
                    )
        else:
            logger.error(f"file manipulation failed for {os.path.basename(input_file)}")
            logger.error(f"Error Output:\n{stderr}")
            raise subprocess.CalledProcessError(process.returncode, " ".join(command))
    except Exception as e:
        logger.error(
            f"An error occurred during file manipulation: {stderr if stderr else stdout}"
        )
        raise e


# TODO:MAKE THIS UNIVERSAL FOR ALL FORMATS LATER
def convert_to_pdb(
    input_file: str,
    context: str = GlideContext.get_current(),
    n_structures: list = list(range(1, 4)),
):
    r"""Convert docking file to pdb format.
    Parameters
    ----------
    input_file : str
        Path to Input file name.

    context : str
        Glide context. Default is GlideContext.get_current().
    n_structures : int
        Number of structures to convert. Default is the first 3 structures.
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
