"""Implements the Docking class."""

import glob
import os
import shutil
import subprocess

from loguru import logger

from ..structure.ligand import Ligand, PreparedLigand
from ..structure.protein import PreparedProtein, Protein
from .contexts import GlideContext
from .docking import Docking


class Glide(Docking):
    r"""Perform docking with `glide_sif.py`.

    Command documentation can be found
    [here](https://www.schrodinger.com/sites/default/files/s3/release/2023-3/Documentation/html/utilities/program_utility_usage/glide_sif.html).
    """

    def __init__(self) -> None:
        # self.context= GlideContext.get_current()
        pass

    def run(
        self, target: PreparedProtein, ligand: PreparedLigand, context: GlideContext
    ) -> None:
        r"""Dock ligand into protein grid.

        Args:
            target : The protein structure to dock the ligand into
            ligand : The ligand structure to dock into the protein
            context : The context for the glide docking
        Returns:
            None
        """
        # ensure that prepped_ligand and grid_file are defined and if not raise an error and exit
        logger.info(ligand.file_path, ligand.file_id)
        jobname = str(ligand.file_id.split("_prepared")[0]) + "_docking"
        logger.info(jobname)
        command = [
            context.command + "/run",
            "glide_sif.py",
            "-gridfile",
            target.file_path,
            "-ligandfile",
            ligand.file_path,
            "-forcefield",
            context.forcefield,
            "-precision",
            context.docking_protocol,
            "-docking_method",
            context.docking_method,
            "-nenhanced_sampling",
            context.n_enhanced_sampling,
            "-postdock_npose",
            context.postdock_nposes,
            "-poses_per_lig",
            context.poses_per_lig,
            os.path.join(context.write_dir, jobname),
        ]
        try:
            with subprocess.Popen(
                command,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True,
            ) as process:
                _, stderr = process.communicate()
            if process.returncode == 0:
                logger.info(
                    f"Docking jobname file created for {target.file_id} and {ligand.file_id}"
                )
            else:
                logger.error(
                    f"Docking jobname file failed for {target.file_id} and {ligand.file_id}"
                )
                logger.error(f"Error Output:\n{stderr}")
        except Exception as e:
            logger.error(f"An error occurred during docking: {str(e)}")
            raise e
        command = [
            context.command + "/glide",
            "-WAIT",
            os.path.join(context.write_dir, f"{jobname}.in"),
        ]
        try:
            with subprocess.Popen(
                command,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True,
            ) as process:
                _, stderr = process.communicate()
            if process.returncode == 0:
                logger.info(f"Glide docking completed for {jobname}")
                shutil.move(
                    f"{jobname}.csv", os.path.join(context.write_dir, f"{jobname}.csv")
                )
                shutil.move(
                    f"{jobname}_pv.maegz",
                    os.path.join(context.write_dir, f"{jobname}_pv.maegz"),
                )
                shutil.move(
                    f"{jobname}.log", os.path.join(context.write_dir, f"{jobname}.log")
                )
                shutil.move(
                    f"{jobname}_skip.csv",
                    os.path.join(context.write_dir, f"{jobname}_skip.csv"),
                )
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
    def convert_to_mae(input_object: Ligand | Protein, context: GlideContext) -> None:
        r"""Convert from the pdb format to mae format needed for Schrodinger

        Args:
            input_object : The object to be converted to mae format
            context : The context for the glide docking
        Returns:
            None
        """
        command = [
            context.command + "/utilities/structconvert",
            input_object.file_path,
            os.path.join(context.write_dir, f"{input_object.file_id}.mae"),
        ]
        try:
            with subprocess.Popen(
                command,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True,
            ) as process:
                process.communicate()
            if process.returncode == 0:
                logger.info("Conversion completed successfully.")

        except Exception as e:
            logger.error(f"An error occurred during conversion: {str(e)}")
            raise e

    @staticmethod
    def PrepLigand(ligand: Ligand, context: GlideContext) -> PreparedLigand:
        r"""Check the extension of the ligand file using the split function and
        Prepare ligands for docking using Schrödinger's LigPrep

        Args:
            ligand : The ligand structure to be prepared for docking as a Ligand object
            context : The context for the glide docking
        Returns:
            The prepared ligand structure needed for docking as a PreparedLigand
        """
        if ligand.file_ext == "pdb":
            Glide().convert_to_mae(ligand, context)
            ligand = Ligand(
                file_path=os.path.join(context.write_dir, f"{ligand.file_id}.mae")
            )
        logger.debug(ligand.file_path)
        if ligand.file_ext in ["sd", "mae", "smi", "csv"]:
            command = [
                context.command + "/ligprep",
                "-i" + ligand.file_ext,
                ligand.file_path,
                "-omae",
                os.path.join(context.write_dir, f"{ligand.file_id}_prepared.mae"),
                "-ma",
                context.lig_max_mw,
                "-WAIT",
                "-epik" if context.lig_epik else "",
                "-ph",
                context.lig_ph,
                "-pht",
                context.lig_pht,
                "-bff",
                context.lig_forcefield,
                "-s",
                context.lig_stereoisomers,
            ]
            try:
                with subprocess.Popen(
                    command,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    universal_newlines=True,
                ) as process:
                    _, stderr = process.communicate()
                if process.returncode == 0:
                    logger.info(f"Ligand preparation completed for {ligand.file_id}")
                    # check if the prepared ligand file exists
                    if not os.path.exists(
                        os.path.join(
                            context.write_dir, f"{ligand.file_id}_prepared.mae"
                        )
                    ):
                        logger.error(f"Ligand preparation failed for {ligand.file_id}")
                        logger.error(f"Error Output:\n{stderr}")
                        raise subprocess.CalledProcessError(
                            process.returncode, " ".join(command)
                        )
                    prep = PreparedLigand(
                        file_path=os.path.join(
                            context.write_dir, f"{ligand.file_id}_prepared.mae"
                        )
                    )
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
                "The ligand file is not in the correct format. Please check the extension"
            )
            raise ValueError("Invalid ligand file format")
        return prep

    @staticmethod
    def PrepProtein(
        protein: Protein, context: GlideContext, ligand_asl: str | None = None
    ) -> PreparedProtein:
        r"""Prepare protein structures using Schrödinger's Protein Wizard

        Args:
            protein The protein structure to be prepared for docking as a Protein object
            context : The context for the glide docking
            ligand_asl : The name of the ligand in the protein structure to generate the grid around it
                (i.e ligand residue name). Set to None if the ligand is not known
        Returns:
            The prepared protein structure needed for docking as a PreparedProtein
        """
        command = [
            context.command + "/utilities/prepwizard",
            protein.file_path,
            os.path.join(context.write_dir, f"{protein.file_id}_protein_prepared.mae"),
            "-WAIT",
            "-fillsidechains" if context.fillsidechains else "",
            "-disulfides" if context.disulfides else "",
            "-rehtreat" if context.rehtreat else "",
            "-minimize_adj_h" if context.minimize_adj_h else "",
            "-epik_pH",
            context.epik_ph,
            "-epik_pHt",
            context.epik_pht,
            "-propka_pH",
            context.propka_ph,
            "-r",
            context.prot_rmsd,
            "-f" + context.forcefield,
            "-watdist",
            context.prot_watdist,
        ]
        if context.samplewater:
            command.extend(["-samplewater"])
        logger.info(f"Preparing protein for PDB ID {protein.file_id}")
        try:
            with subprocess.Popen(
                command,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True,
            ) as process:
                _, stderr = process.communicate()
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
        if ligand_asl is None:
            ligand_asl = "ligand"
        else:
            ligand_asl = f"res {ligand_asl.lower()}"
        # Generate a grid around the ligand in the protein structure
        grid_command = [
            context.command + "/utilities/generate_glide_grids",
            "-rec_file",
            os.path.join(context.write_dir, f"{protein.file_id}_protein_prepared.mae"),
            "-lig_asl",
            ligand_asl,
            "-inner_box",
            context.grid_innerbox,
            "-verbose",
            "-forcefield",
            context.forcefield,
            "-WAIT",
        ]
        logger.info(f"Generating grid for PDB ID {protein.file_id}")
        try:
            with subprocess.Popen(
                grid_command,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True,
            ) as process:
                _, stderr = process.communicate()
            if process.returncode == 0:
                logger.info(f"Grid generation completed for PDB ID {protein.file_id}")
                # find the grid file and rename it to the protein name
                grid_file = glob.glob("generate-grids-gridgen.zip")[0]
                os.rename(grid_file, f"{protein.file_id}_grid.zip")
                os.rename("generate_glide_grids_run.log", f"{protein.file_id}_grid.log")
                shutil.move(
                    f"{protein.file_id}_grid.zip",
                    os.path.join(context.write_dir, f"{protein.file_id}_grid.zip"),
                )
                shutil.move(
                    f"{protein.file_id}_grid.log",
                    os.path.join(context.write_dir, f"{protein.file_id}_grid.log"),
                )
                prep = PreparedProtein(file_path=f"{protein.file_id}_grid.zip")
            else:
                logger.error(f"Grid generation failed for PDB ID {protein.file_id}")
                logger.error(f"Error Output:\n{stderr}")
                raise subprocess.CalledProcessError(
                    process.returncode, " ".join(grid_command)
                )
        except Exception as e:
            logger.error(
                f"An unexpected error occurred for PDB ID {protein.file_id}: {str(e)}"
            )
            raise e
        return prep

    @staticmethod
    def sort_docking_results(docking_results: str, context: GlideContext) -> None:
        r"""Sort the docking maegz output based on the glide score

        Args:
            docking_results : The path to the docking results file
            context : The context for the glide docking
        Returns:
            None
        """
        command = [
            context.command + "/utilities/glide_sort",
            "-use_dscore",
            "-o",
            docking_results.replace(".maegz", "_sorted.maegz"),
            docking_results,
        ]
        try:
            with subprocess.Popen(
                command,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True,
            ) as process:
                _, stderr = process.communicate()
            if process.returncode == 0:
                logger.info(f"Sorting completed for {docking_results}")
            else:
                logger.error(f"Sorting failed for {docking_results}")
                logger.error(f"Error Output:\n{stderr}")
                raise subprocess.CalledProcessError(
                    process.returncode, " ".join(command)
                )
        except Exception as e:
            logger.error(f"An error occurred during sorting: {str(e)}")
            raise e
