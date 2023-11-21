r"""Implementation of combind."""
from typing import Union

import glob
import os
import subprocess

from loguru import logger

from .contexts.combind import CombindContext


class Combind(CombindContext):
    r"""Combind class for pose prediction/selection."""

    def featurize(
        self,
        docking_filepath: Union[str, None] = None,
        file_name: Union[str, None] = None,
    ):
        r"""Featurize the docking poses file.
        Parameters
        ----------
        docking_filepath : str, optional
            Path to the docking file.
        file_name : str, optional
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
        command = [
            self.command,
            "/combind",
            "featurize",
            f"{file_name}_features",
            docking_file,
        ]
        try:
            process = subprocess.Popen(
                command, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
            stdout, stderr = process.communicate()
            if process.returncode == 0:
                logger.info("Featurization completed.")
            else:
                logger.error(f"Featurization failed\n{stderr}.")
        except Exception as e:
            logger.error(f"Featurization failed. {str(e)}")
            raise e

    def select_pose(
        self, file_name: Union[str, None] = None, features_dir: Union[str, None] = None
    ):
        r"""Select the best pose from the docking file.
        Parameters
        ----------
        file_name : str, optional
            Name of the output file.
        features_dir : str, optional
            name of the features directory.
        """
        command = [
            self.command,
            "/combind",
            "pose-prediction",
            self.work_dir,
            features_dir,
            f"{file_name}.csv",
        ]
        try:
            process = subprocess.Popen(
                command, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
            stdout, stderr = process.communicate()
            if process.returncode == 0:
                logger.info("Pose selection completed.")
            else:
                logger.error(f"Pose selection failed\n{stderr}.")
        except Exception as e:
            logger.error(f"Pose selection failed. {str(e)}")
            raise e
