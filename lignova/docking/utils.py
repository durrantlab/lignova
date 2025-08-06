r"""Implementation for different Schrodinger's suplementary functions."""

import glob
import os
import subprocess

from loguru import logger

from ..docking.contexts import GlideContext


def manipulate_complexes(
    input_file: str,
    context: GlideContext = GlideContext.get_current(),
    outfile_name: None | str = "manipulate_outp.maegz",
    mode: None | str = "merge",
) -> None:
    r"""Convert docking file to complexes format and saving in write_dir.

    Args:
        input_file : Path to Input file name.
        context : the Glide context. Default is GlideContext.get_current().
        outfile_name : Output file name. Default is input file name with "manipulate_outp.maegz" as suffix.
        mode : Mode to convert file. Default is "merge".
            Options are ["merge" - combine PV/EPV structures into complexes,
            ,"split_pv" - extract receptor and ligands from complexes and save as PV
            ,"split_epv" - extract receptor and ligands from complexes and save as EPV
            ,"split_ligand" - extract ligands from complexes
            ,"split_receptor" - extract receptors from complexes
            ,"pv_to_epv" - join multiple PV files into single EPV file
            ,"epv_to_pv" - split EPV file into multiple PV files ]
    Returns:
        None
    """
    # check if "out" in the outfile_name and if so throw an error
    if "-out" in outfile_name:
        logger.error("Output file name cannot contain '-out'")
        raise ValueError("Output file name cannot contain '-out'")
    filename = os.path.basename(input_file).split(".")[0]
    logger.debug(f"Filename: {filename}")
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
    if outfile_name == "manipulate_outp.maegz":
        new_filename = filename + "_" + outfile_name
    else:
        new_filename = outfile_name

    # check if new_filename exists and if so append a number to the filename
    if os.path.exists(os.path.join(context.write_dir, new_filename)):
        i = 1
        while os.path.exists(os.path.join(context.write_dir, new_filename)):
            new_filename = os.path.join(
                context.write_dir, f"{filename}_{i}_{outfile_name}"
            )
            i += 1

    command = [context.command + "/run", "pv_convert.py", "-mode", mode]
    if mode in ["split_pv", "split_epv", "split_ligand", "split_receptor"]:
        command.extend(["-lig_last_mol"])
    command.extend(["-o", os.path.join(context.write_dir, new_filename), input_file])
    try:
        with subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
        ) as process:
            stdout, stderr = process.communicate()

        logger.debug(f"Output:\n{stdout}")
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
                # find the file with the same name as the input file
                # but with "-out" in the name and .maegz extension
                directory = os.path.dirname(input_file)
                logger.debug(f"looking for the file in: {directory}")
                for file in glob.glob(os.path.join(directory, "*")):
                    filename = os.path.basename(file)
                    # logger.debug(f"Checking file: {file}")
                    if (
                        filename.startswith(filename)
                        and (filename.endswith(".maegz") or filename.endswith(".mae"))
                        and filename != os.path.basename(input_file)
                        and "-out" in file
                    ):
                        logger.debug(f"found it: {file}")
                        complexes_file = file
                        break
                if complexes_file:
                    # Rename the complexes file to match the input file name
                    os.rename(
                        complexes_file, os.path.join(context.write_dir, new_filename)
                    )
                    logger.info(
                        f"Converted file saved at: {os.path.join(context.write_dir,new_filename)}"
                    )
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


# NOTE:IN THE FUTURE MAKE THIS UNIVERSAL FOR ALL FORMATS LATER
def convert_to_pdb(
    input_file: str,
    context: str = GlideContext.get_current(),
    n_structures: list | None = None,
    reorder: bool = False,
):
    r"""Convert docking file to pdb format.

    Args:
        input_file : Path to Input file name.
        context : Glide context. Default is GlideContext.get_current().
        n_structures : List of structure indices to convert. If None, converts all structures.
        reorder : if False, the atom indexing is preserved. Default is False.
    Returns:
        None
    """
    command = [context.command + "/utilities/structconvert", "-no_component_dict"]
    if not reorder:
        command.append("-no_renum")
    if n_structures:
        command.extend(["-n", ",".join(map(str, n_structures))])
    # Break down the long line
    output_path = os.path.join(
        context.write_dir, f"{os.path.splitext(os.path.basename(input_file))[0]}.pdb"
    )
    command.extend([input_file, output_path])

    try:
        with subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
        ) as process:
            stdout, stderr = process.communicate()

        logger.debug(f"Output:\n{stdout}")
        logger.debug(f"Error Output:\n{stderr}")
        if process.returncode == 0:
            logger.info(f"Converted {input_file} to pdb format")
            logger.info(f"Output:\n{stdout}")
        # Break down the long line
        elif (
            "Each Structure is converted independently and written to separate files"
            in stderr
        ):
            output_file_new = os.path.join(
                context.write_dir,
                f"{os.path.splitext(os.path.basename(input_file))[0]}-1.pdb",
            )
            logger.info(f"Converted files one by one: {output_file_new}")
        else:
            logger.error(f"Conversion failed for {input_file}")
            logger.error(f"Error Output:\n{stderr}")
            raise subprocess.CalledProcessError(process.returncode, " ".join(command))
    except Exception as e:
        logger.error(f"An error occurred during conversion: {str(e)}")
        raise e
