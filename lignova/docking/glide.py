"""Implements the Docking class."""
import glob
import os
import subprocess

from loguru import logger

from ..structure.ligand import Ligand, PreparedLigand
from ..structure.protein import PreparedProtein
from .contexts import GlideContext
from .docking import Docking


class Glide(Docking):
    r"""Perform docking with `glide_sif.py`.

    Command documentation can be found
    [here](https://www.schrodinger.com/
    sites/default/files/s3/release/2023-3/
    Documentation/html/utilities/program_utility_usage/glide_sif.html).
    """

    def __init__(self) -> None:
        # TODO:DONE? move this check with context in run
        # self.context= GlideContext.get_current()
        pass

    def run(self, target, ligand, context):
        r"""Dock ligand into protein grid."""
        # TODO:DONE Break apart into functions and update with new arguments
        # ensure that prepped_ligand and grid_file are defined and if not raise an error and exit
        if not prepped_ligand or not self.grid_file:
            logger.error("Prepared Ligand or Prepared Protein not defined.")
            raise SystemExit(1)
        jobname = prepped_ligand.lig_name + "_docking"
        command = [
            os.environ["SCHRODINGER"] + "/run",
            "glide_sif.py",
            "-gridfile",
            self.grid_file.file_path,
            "-ligandfile",
            prepped_ligand.file_path,
            "-calc_input_rms",
            "yes",
            "-forcefield",
            context.forcefield,
            "-precision",
            context.docking_protocol,
            "-nenhanced_sampling",
            context.n_enhanced_sampling,
            "-postdock_npose",
            context.postdock_nposes,
            os.path.join(context.write_dir, f"{jobname}.in"),
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
                logger.info(
                    f"Docking completed for {self.grid_file.protein_name} and {prepped_ligand.lig_name}"
                )
            else:
                logger.error(
                    f"Docking failed for {self.grid_file.protein_name} and {prepped_ligand.lig_name}"
                )
                logger.error(f"Error Output:\n{stderr}")
        except Exception as e:
            logger.error(f"An error occurred during docking: {str(e)}")
            raise e
        command = [os.environ["SCHRODINGER"] + "/glide", "-WAIT", os.path.join(context.write_dir,  f"{jobname}.in")]
        try:
            process = subprocess.Popen(
                command,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True,
            )
            stdout, stderr = process.communicate()
            if process.returncode == 0:
                logger.info(f"Glide docking completed for {jobname}")
            else:
                logger.error(f"Glide docking failed for {jobname}")
                logger.error(f"Error Output:\n{stderr}")
                raise subprocess.CalledProcessError(
                    process.returncode, " ".join(command)
                )
        except Exception as e:
            logger.error(f"An error occurred during Glide docking: {str(e)}")
            raise e

    @staticmethod
    def convert_to_mae(input_object, context):
        r"""Convert from the pdb format to mae format needed for Schrodinger"""
        command = [
            context.command + "/utilities/structconvert",
            input_object.file_path,
            os.path.join(context.write_dir,f"{input_object.file_id}.mae"),
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
                logger.info("Conversion completed successfully.")
            else:
                logger.error(f"Conversion failed with error:\n{stderr}")
        except Exception as e:
            logger.error(f"An error occurred during conversion: {str(e)}")
            raise e

    @staticmethod
    def PrepLigand(ligand, context):
        r"""Check the extension of the ligand file using the split function and
        Prepare ligands for docking using Schrödinger's LigPrep"""
        if ligand.file_ext == "pdb":
            Glide().convert_to_mae(ligand,context)
            ligand = Ligand(file_path=os.path.join(context.write_dir,f"{ligand.file_id}.mae"))
        logger.debug(ligand.file_path)
        if ligand.file_ext in ["sd", "mae", "smi", "csv"]:
            command = [
                context.command + "/ligprep",
                "-i" + ligand.file_ext,
                ligand.file_path,
                "-omae",
                os.path.join(context.write_dir, f"{ligand.file_id}_prepared.mae"),
                "-WAIT",
                "-epik",
                "-ph",
                context.lig_ph,
                "-pht",
                context.lig_pht,
                "-bff",
                context.lig_forcefield,
                "-ac",
                "-s",
                context.lig_stereoisomers,
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
                    logger.info(f"Ligand preparation completed for {ligand.file_id}")
                    prep = PreparedLigand(os.path.join(context.write_dir, f"{ligand.file_id}_prepared.mae"))
                else:
                    logger.error(f"Ligand preparation failed for {ligand.file_id}")
                    logger.error(f"Error Output:\n{stderr}")
                    raise subprocess.CalledProcessError(
                        process.returncode, " ".join(command)
                    )
            except Exception as e:
                logger.error(f"An error occurred during ligand preparation: {str(e)}")
                raise e
        else:
            logger.error(
                "The ligand file is not in the correct format. Please check the extension of the ligand file"
            )
            raise ValueError("Invalid ligand file format")
        return prep

    @staticmethod
    def PrepProtein(protein, context):
        r"""Prepare protein structures using Schrödinger's Protein Wizard"""
        command = [
            context.command + "/utilities/prepwizard",
            f"{protein.file_id}.mae",
            os.path.join(context.write_dir, f"{protein.file_id}_protein_prepared.mae"),
            f"{protein.file_id}_protein_prepared.mae",
            "-WAIT",
            "-epik_pH",
            context.epik_pH,
            "-epik_pHt",
            context.epik_pHt,
            "-propka_pH",
            context.propka_pH,
            "-r",
            context.rmsd,
            "-f".join(context.forcefield),
        ]
        logger.info(f"Preparing protein for PDB ID {protein.file_id}")
        try:
            process = subprocess.Popen(
                command,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True,
            )
            stdout, stderr = process.communicate()
            if process.returncode == 0:
                logger.info(f"Protein preparation completed for {protein.file_id}")
            else:
                logger.error(f"Protein preparation failed for {protein.file_id}")
                logger.error(f"Error Output:\n{stderr}")
                raise subprocess.CalledProcessError(
                    process.returncode, " ".join(command)
                )
        except Exception as e:
            logger.error(f"An error occurred during protein preparation: {str(e)}")
            raise e

        # Generate a grid around the ligand in the protein structure
        grid_command = [
            context.command + "/utilities/generate_glide_grids",
            "-rec_file",
            f"{protein.file_id}_protein_prepared.mae",
            "-lig_asl",
            "ligand",
            "-inner_box",
            context.inner_box,
            "-verbose",
            "-forcefield",
            context.forcefield,
            "-WAIT",
        ]
        logger.info(f"Generating grid for PDB ID {protein.file_id}")
        try:
            process = subprocess.Popen(
                grid_command,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True,
            )
            stdout, stderr = process.communicate()
            if process.returncode == 0:
                logger.info(f"Grid generation completed for PDB ID {protein.file_id}")
                # find the grid file and rename it to the protein name
                grid_file = glob.glob("generate-grids-gridgen.zip")[0]
                os.rename(grid_file, f"{protein.file_id}_grid.zip")
                os.rename("generate_glide_grids_run.log", f"{protein.file_id}_grid.log")
                prep = PreparedProtein(file_path=f"{protein.file_id}_grid.zip")
            else:
                logger.error(f"Grid generation failed for PDB ID {protein.file_id}")
                logger.error(f"Error Output:\n{stderr}")
                raise subprocess.CalledProcessError(
                    process.returncode, " ".join(grid_command)
                )
        except Exception as e:
            logger.error(f"Error while processing PDB ID {protein.file_id}: {e.stderr}")
        except Exception as e:
            logger.error(
                f"An unexpected error occurred for PDB ID {protein.file_id}: {str(e)}"
            )
            raise e
        return prep
