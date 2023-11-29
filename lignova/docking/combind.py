r"""Implementation of combind."""
from typing import TextIO, Union

import glob
import os
import subprocess

from loguru import logger

from .contexts.combind import CombindContext


class Combind(CombindContext):
    r"""Combind class for pose prediction/selection."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.activate = f'source {self.schrodinger_env +"/bin/activate"} && '

    def featurize(
        self,
        docking_filepath: Union[str, TextIO],
        file_name: Union[str],
    ):
        r"""Featurize the docking poses file.
        Parameters
        ----------
        docking_filepath : str, file-like object
            Path to the docking file.
        file_name : str,
            Name of the output file.
        """
        if docking_filepath is None or not os.path.exists(
            os.path.join(self.work_dir, docking_filepath)
        ):
            logger.error("Docking file is not found.")
            # find files with _pv.maegz extension in the work_dir using glob
            files = glob.glob(os.path.join(self.work_dir, "*_pv.maegz"))
            if len(files) == 0 or len(files) > 1:
                logger.error(
                    "No files or Multiple docking files found in the work_dir."
                )
                raise OSError(
                    "No files or Multiple docking files found in the work_dir."
                )
            docking_file = os.path.join(self.work_dir, str(files[0]))
        else:
            docking_file = os.path.join(self.work_dir, docking_filepath)
        command1 = [
            self.activate,
            self.command + "/combind",
            "featurize",
            os.path.join(self.work_dir, f"{file_name}_features"),
            docking_file,
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
                    "Schrodinger virtual environment activated & Featurization completed."
                )
                with open(
                    os.path.join(self.work_dir, f"{file_name}_features.log"),
                    "w",
                    encoding="utf-8",
                ) as file:
                    file.write(process.stdout.decode())
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

    def select_pose(self, file_name: Union[str, TextIO], features_dir: Union[str]):
        r"""Select the best pose from the docking file.
        Parameters
        ----------
        file_name : str, file-like object
            Name of the output file.
        features_dir : str
            name of the features directory.
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
        self, docking_filepath: Union[str, TextIO], combind_csv: Union[str, TextIO]
    ):
        r"""Get the top ligand pose after combind prediction.
        Parameters
        ----------
        docking_filepath : str, file-like object
            Path to the docking file from GLIDE.
        combind_csv : str, file-like object
            Path to the combind csv file.
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
            process = subprocess.run(
                command,
                shell=True,
                executable="/bin/bash",
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=True,
            )
            if process.returncode == 0:
                logger.info("Pose extraction from the docking file completed.")
            else:
                error_message = (
                    "Failed to activate extract the top pose from the docking file"
                    + f"\n{process.stderr.decode()}."
                )
                logger.critical(error_message)
                raise NotImplementedError(error_message)
        except Exception as e:
            logger.error(f"Pose extraction failed.\n {str(e)}")
            raise e

    def compute_combind_score(self, features_dir: Union[str], filename: Union[str]):
        r"""Compute the combind score.
        Parameters
        ----------
        features_dir : str
            Path of the features directory.
        filename : str
            Name of the output npy file.
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
        docking_filepath: Union[str, TextIO],
        combind_score_file: Union[str, TextIO],
        output_filename: Union[str, TextIO],
        sort: bool = True,
    ):
        r"""Apply the combind score to the docking file.
        Parameters
        ----------
        docking_filepath : str, file-like object
            Path to the docking file from GLIDE.
        combind_score_file : str, file-like object
            Path to the combind score file.
        output_filename : str, file-like object
            name of the output file.
        sort : bool, optional, default=True
            Sort the output file by the combind score.
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
            process = subprocess.run(
                command,
                shell=True,
                executable="/bin/bash",
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                check=True,
            )
            if process.returncode == 0:
                logger.info("Combind score application completed.")
                if sort:
                    command = [
                        self.schrodinger + "/utilities/glide_sort",
                        "-best_by_title",
                        "-use_prop_d",
                        "r_i_combind_score",
                        "-o",
                        os.path.join(
                            self.work_dir, f"{output_filename}_combind_sorted.maegz"
                        ),
                        os.path.join(self.work_dir, f"{output_filename}_combind.maegz"),
                    ]
                    process = subprocess.Popen(
                        command, stdout=subprocess.PIPE, stderr=subprocess.PIPE
                    )
                    stdout, stderr = process.communicate()
                    print(command)
                    print(stdout.decode())
                    print(stderr.decode())
                    if process.returncode == 0:
                        logger.info("Combind score application sorting completed.")
                        with open(
                            os.path.join(
                                self.work_dir, f"{output_filename}_combind_sort.log"
                            ),
                            "w",
                            encoding="utf-8",
                        ) as file:
                            file.write(stdout.decode())
                    else:
                        error_message = (
                            "Failed to sort the combind score application"
                            + f"\n{stderr.decode()}."
                        )
                        logger.critical(error_message)
                        raise NotImplementedError(error_message)
            else:
                error_message = (
                    "Failed to apply the combind score"
                    + f"\n{process.stderr.decode()}."
                )
                logger.critical(error_message)
                raise NotImplementedError(error_message)
        except Exception as e:
            logger.error(f"Combind score application failed.\n {str(e)}")
            raise e
