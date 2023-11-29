r""" Implementation of the CombindContext class containing the configuration for Glide docking."""
from typing import Optional, Union

import os
import subprocess

from loguru import logger

_default_combind_context: "Optional[CombindContext]" = None

DEFAULT_COMMAND = "/home/mma121/PubChem_small/combind"
DEFAULT_WORK_DIR = "./tmp/6oav"
DEFAULT_SCHRODINGER = os.environ.get("SCHRODINGER", None)
DEFAULT_SCHRODINGER_ENV = "schrodinger.ve"


class CombindContext:
    r"""Singleton for combind pose selection configuration."""

    def __init__(
        self,
        command: Union[str, None] = None,
        schrodinger: Union[str, None] = None,
        work_dir: Union[str, None] = None,
        schrodinger_env: Union[str, None] = None,
    ):
        """Initialize the combind context."""
        self.work_dir = work_dir
        # Check if Schrödinger is installed or the $SCHRODINGER environment variable is set
        if not os.environ.get("SCHRODINGER") and schrodinger is None:
            error_message = "Schrödinger is not installed"
            "or the $SCHRODINGER environment variable is not set."
            logger.critical(error_message)
            raise NotImplementedError(error_message)

        self.schrodinger = schrodinger
        current_directory = os.getcwd()
        logger.debug(f"Current working directory before: {os.getcwd()}")
        os.chdir(work_dir)
        if not os.path.isdir(schrodinger_env):
            logger.debug(
                "Schrodinger virtual environment is not found. Creating in work directory..."
            )
            logger.debug(f"Current working directory after: {os.getcwd()}")
            command_run = [
                self.schrodinger + "/run",
                "schrodinger_virtualenv.py",
                schrodinger_env,
            ]
            process = subprocess.Popen(
                command_run,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            stderr = process.communicate()
            if process.returncode == 0:
                logger.info("Schrodinger virtual environment created.")
                self.schrodinger_env = os.path.join(os.getcwd(), schrodinger_env)
                os.chdir(path=current_directory)
            else:
                error_message = (
                    f"Schrodinger virtual environment failed\n{stderr[1].decode()}."
                )
                logger.critical(error_message)
                raise NotImplementedError(error_message)
        else:
            logger.info("Schrodinger virtual environment found.")
            self.schrodinger_env = os.path.join(self.work_dir, schrodinger_env)
            os.chdir(path=current_directory)
        if command is None:
            error_message = (
                "Please provide the path to combind. If you have not installed combind,"
                "please clone this repository:https://github.com/drorlab/combind.git"
            )
            logger.critical(error_message)
            raise NotImplementedError(error_message)
        logger.info("Setting up combind specific variables.")
        combind_variables = {
            "COMBINDHOME": command,
            "PATH": os.pathsep.join([os.environ.get("PATH"), command]),
        }
        for key, value in combind_variables.items():
            os.environ[key] = value
            # Get the value of COMBINDHOME from the environment
            combind_home = os.environ.get("COMBINDHOME")
            if combind_home not in os.environ.get("PATH"):
                os.environ["PATH"] += os.pathsep + combind_home
        # Activate Combind environment
        if os.environ.get("COMBINDHOME"):
            logger.info("Combind environment activated.")
            os.chdir(current_directory)
            logger.debug(os.getcwd())
            self.command = os.environ.get("COMBINDHOME")
            logger.info(f"Combind command: {self.command}")
        else:
            error_message = f"Combind environment activation failed\n{stderr[1]}."
            logger.critical(error_message)
            raise NotImplementedError(error_message)

    @staticmethod
    def get_current() -> "CombindContext":
        r"""Get or create a singleton context."""

        global _default_combind_context

        if _default_combind_context is None:
            _default_combind_context = CombindContext(
                command=DEFAULT_COMMAND,
                schrodinger=DEFAULT_SCHRODINGER,
                work_dir=DEFAULT_WORK_DIR,
                schrodinger_env=DEFAULT_SCHRODINGER_ENV,
            )
        return _default_combind_context

    @staticmethod
    def set_current(context: "CombindContext") -> None:
        """Set the current context."""
        global _default_combind_context
        _default_combind_context = context
