r"""Implementation of combind."""

from typing import List, TextIO

import glob
import os
import subprocess

import pandas as pd
from loguru import logger

from .contexts.combind import CombindContext


class Combind(CombindContext):
    r"""Combind class for pose prediction/selection."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        if not self.validate():
            raise ValueError("Validation failed. Combind intialization aborted")
        self.activate = f'source {self.schrodinger_env + "/bin/activate"} && '

    def featurize(
        self,
        docking_filepaths: str | List[str],
        file_name: str,
        max_poses: int = 8500000,
        screen: bool = False,
    ) -> None:
        r"""Featurize the docking files.

        Args:
            docking_filepaths : Path to the docking file from GLIDE.
                Can be a single file or a list of two files.
            file_name : Name of the output file.
            max_poses : the maximum number of poses ComBind can consider. Default is 7000000
            screen : whether to the docking is for virtual screening (ComBindVS) or not.
                Default is False.
        Returns:
            None
        """
        if isinstance(docking_filepaths, str):
            docking_filepaths = [docking_filepaths]

        if len(docking_filepaths) == 0 or len(docking_filepaths) > 2:
            logger.error(
                "Invalid number of docking files provided. Must provide one or two files."
            )
            raise ValueError(
                "Invalid number of docking files provided. Must provide one or two files."
            )
        elif len(docking_filepaths) < 2 and screen:
            logger.error(
                "Screening requires two docking files. Please provide a second file."
            )
            raise ValueError(
                "Screening requires two docking files. Please provide a second file."
            )
        # make sure all the files exist
        for file in docking_filepaths:
            if not os.path.exists(file):
                logger.error(f"{file} does not exist.")
                raise FileNotFoundError(f"{file} does not exist.")
        command1 = [
            self.activate,
            self.command + "/combind",
            "featurize",
            "--max-poses",
            str(max_poses),
        ]
        if screen:
            command1.append("--no-mcss")
            command1.append("--screen")
        command1.append(os.path.join(self.work_dir, f"{file_name}_features"))
        command1.extend(docking_filepaths)
        command = " ".join(command1)

        try:
            process = subprocess.run(
                command,
                shell=True,
                executable="/bin/bash",
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=True,
            )
            if process.returncode == 0:
                logger.info(
                    "Schrodinger virtual environment activated & Featurization completed."
                )
                with open(
                    os.path.join(self.work_dir, f"{file_name}_features.log"),
                    "w",
                    encoding="utf-8",
                ) as file:
                    file.write(process.stdout.decode())
                # check if it created files in the work_dir with the words rmsd name gscore ifp
                # and delete them
                for file in os.listdir(self.work_dir):
                    if any(
                        keyword in file for keyword in ["rmsd", "name", "gscore", "ifp"]
                    ):
                        file_path = os.path.join(self.work_dir, file)
                        if os.path.isfile(file_path):
                            os.remove(file_path)
                            logger.info(f"Deleted file: {file_path}")
            else:
                error_message = (
                    "Failed to activate Schrodinger virtual environment"
                    + f"\n{process.stderr.decode()}."
                )
                logger.critical(error_message)
                raise NotImplementedError(error_message)
        except Exception as e:
            logger.error(f"Featurization failed. {str(e)}")
            raise e

    def select_pose(self, file_name: str | TextIO, features_dir: str) -> None:
        r"""Select the best pose from the docking file.

        Args:
            file_name : Name of the output file.
            features_dir : Name of the features directory.
        Returns:
            None
        """
        command1 = [
            self.activate,
            self.command + "/combind",
            "pose-prediction",
            features_dir,
            os.path.join(self.work_dir, f"{file_name}.csv"),
        ]
        command = " ".join(command1)
        try:
            process = subprocess.run(
                command,
                shell=True,
                executable="/bin/bash",
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=True,
            )
            if process.returncode == 0:
                logger.info(
                    "Schrodinger virtual environment activated & Pose prediction completed."
                )
            else:
                error_message = (
                    "Failed to activate Schrodinger virtual environment"
                    + f"\n{process.stderr.decode()}."
                )
                logger.critical(error_message)
                raise NotImplementedError(error_message)
        except Exception as e:
            logger.error(f"Pose prediction failed. {str(e)}")
            raise e

    def get_3d_top_pose(
        self,
        docking_filepath: str | TextIO,
        combind_csv: str | TextIO,
        extract_filename: str,
    ) -> None:
        r"""Get the top ligand pose after combind prediction.

        Args:
            docking_filepath: Path to the docking file from GLIDE.
            combind_csv: Path to the combind csv file.
            extract_filename: Name of the output file.
        Returns:
            None
        """
        command1 = [
            self.activate,
            self.command + "/combind",
            "extract-top-poses",
            combind_csv,
            docking_filepath,
        ]
        command = " ".join(command1)
        try:
            subprocess.run(
                command,
                shell=True,
                executable="/bin/bash",
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=True,
            )
            logger.info("Pose extraction from the docking file completed.")
            output_file = glob.glob(
                os.path.join(
                    self.work_dir,
                    f"{os.path.splitext(os.path.basename(combind_csv))[0]}_pv.maegz",
                )
            )
            os.rename(
                output_file[0],
                os.path.join(self.work_dir, f"{extract_filename}_pv.maegz"),
            )
        except subprocess.CalledProcessError as e:
            logger.error(f"Pose extraction failed: {e}")
            raise
        except Exception as e:
            logger.error(f"Pose extraction failed.\n {str(e)}")
            raise

    def compute_combind_score(self, features_dir: str, filename: str) -> None:
        r"""Compute the combind score.

        Args:
            features_dir : Path of the features directory.
            filename : Name of the output npy file.
        Returns:
            None
        """
        command1 = [
            self.activate,
            self.command + "/combind",
            "screen",
            os.path.join(self.work_dir, f"{filename}_screen.npy"),
            features_dir,
        ]
        command = " ".join(command1)
        try:
            process = subprocess.run(
                command,
                shell=True,
                executable="/bin/bash",
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=True,
            )
            if process.returncode == 0:
                logger.info("Combind score computation completed.")
            else:
                error_message = (
                    "Failed to compute the combind score"
                    + f"\n{process.stderr.decode()}."
                )
                logger.critical(error_message)
                raise NotImplementedError(error_message)
        except Exception as e:
            logger.error(f"Combind score computation failed.\n {str(e)}")
            raise e

    def apply_combind_score(
        self,
        docking_filepath: str | TextIO,
        combind_score_file: str | TextIO,
        output_filename: str | TextIO,
        sort: bool = True,
    ):
        r"""Apply the combind score to the docking file.

        Args:
            docking_filepath : Path to the docking file from GLIDE.
            combind_score_file :Path to the combind score file.
            output_filename :  name of the output file.
            sort : whether to sort the output file based on the combind score or not.
                Default is True.
        Returns:
            None
        """
        command1 = [
            self.activate,
            self.command + "/combind",
            "apply-scores",
            docking_filepath,
            combind_score_file,
            os.path.join(self.work_dir, f"{output_filename}_combind.maegz"),
        ]
        command = " ".join(command1)
        try:
            subprocess.run(
                command,
                shell=True,
                executable="/bin/bash",
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=True,
            )
            logger.info("Combind score application completed.")
            if sort:
                command = [
                    self.schrodinger + "/utilities/glide_sort",
                    "-use_prop_d",
                    "r_i_combind_score",
                    "-o",
                    os.path.join(
                        self.work_dir, f"{output_filename}_combind_sorted.maegz"
                    ),
                    os.path.join(self.work_dir, f"{output_filename}_combind.maegz"),
                ]
                result = subprocess.run(
                    command,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    universal_newlines=True,
                    check=True,
                )
                logger.info("Combind score application sorting completed.")
                with open(
                    os.path.join(self.work_dir, f"{output_filename}_combind_sort.log"),
                    "w",
                    encoding="utf-8",
                ) as file:
                    file.write(result.stdout)
        except subprocess.CalledProcessError as e:
            error_message = f"Failed to apply or sort the combind score: {e}"
            logger.critical(error_message)
            raise NotImplementedError(error_message) from e
        except Exception as e:
            logger.error(f"Combind score application failed.\n {str(e)}")
            raise

    def extract_data_csv(
        self,
        docking_file: str | TextIO,
        filename: str,
        filter_data: bool = True,
        top_poses: bool = False,
    ) -> str:
        r"""To extract the scores from schrodinger docking file including
            combind if run after apply_combind_score function.

        Args:
            docking_file :  Path to the docking file from GLIDE.
            filename : name of the output file.
            filter_data : whether to filter_data the docking file to include only scores or not.
                Default is True.
            top_poses : If True, extract the data of the top poses per ligand only.
                Default is False.
        Returns:
            file path where the csv file is saved.
        """
        # check if the docking file exists
        if not os.path.exists(docking_file):
            raise FileNotFoundError(f"{docking_file} does not exist.")
        command = [
            self.schrodinger + "/utilities/proplister",
            "-a",
            "-c",
            docking_file,
            "-o",
            os.path.join(self.work_dir, f"{filename}.csv"),
        ]
        with subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
        ) as process:
            _, stderr = process.communicate()
        if process.returncode == 0:
            logger.info("Data extraction completed.")
            if top_poses:
                # find the best pose for each ligand in terms of the highest combind score
                data = pd.read_csv(os.path.join(self.work_dir, f"{filename}.csv"))
                # add a column to called POSE with the loacation
                # of the row within the group of ligands so it would recount when the ligand changes
                data["POSE"] = data.groupby("s_m_title").cumcount() + 1
                # save the data
                # group by ligand and find the max combind score and get that row
                data = data.loc[data.groupby("s_m_title")["r_i_combind_score"].idxmax()]
                logger.debug(data)
                logger.info("Top poses extracted.")
                # save it with the same name
                data.to_csv(os.path.join(self.work_dir, f"{filename}.csv"), index=False)
            if filter_data:
                # read the csv file and filter_data the scores
                data = pd.read_csv(os.path.join(self.work_dir, f"{filename}.csv"))
                # if the data has a column called POSE keep it
                if "POSE" in data.columns:
                    data = data[
                        [
                            "s_m_title",
                            "i_i_glide_posenum",
                            "r_i_docking_score",
                            "r_i_combind_score",
                            "r_i_glide_emodel",
                            "r_i_glide_energy",
                            "r_i_glide_gscore",
                            "r_i_glide_ligand_efficiency",
                            "POSE",
                        ]
                    ]
                else:
                    data = data[
                        [
                            "s_m_title",
                            "i_i_glide_posenum",
                            "r_i_docking_score",
                            "r_i_combind_score",
                            "r_i_glide_emodel",
                            "r_i_glide_energy",
                            "r_i_glide_gscore",
                            "r_i_glide_ligand_efficiency",
                        ]
                    ]
                # save it with the same name
                data.to_csv(os.path.join(self.work_dir, f"{filename}.csv"), index=False)
                logger.info("Data filter_dataation completed.")
        else:
            error_message = "Failed to extract the data" + f"\n{stderr.decode()}."
            logger.critical(error_message)
            raise NotImplementedError(error_message)
        return os.path.join(self.work_dir, f"{filename}.csv")
