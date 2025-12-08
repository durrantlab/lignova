r""" Implementation of the CombindContext class containing the configuration for Glide docking."""

import os
import subprocess

from loguru import logger

DEFAULT_COMMAND = "/home/mma121/PubChem_small/combind"
DEFAULT_WORK_DIR = "./tmp/6oav"
DEFAULT_SCHRODINGER = os.environ.get("SCHRODINGER", None)
DEFAULT_SCHRODINGER_ENV = "schrodinger.ve"


class CombindContext:
    r"""Singleton for combind pose selection configuration."""

    def __init__(
        self,
        command: str | None = None,
        schrodinger: str | None = None,
        work_dir: str | None = None,
        schrodinger_env: str | None = None,
    ):
        """Initialize the combind context."""
        self.work_dir = work_dir
        self.schrodinger = schrodinger
        self.command = command
        self.schrodinger_env = schrodinger_env

    def validate(self) -> bool:
        """Validate the combind context parameters."""
        # Check if Schrödinger is installed or the $SCHRODINGER environment variable is set
        if not os.environ.get("SCHRODINGER") and self.schrodinger is None:
            logger.critical(
                "Schrödinger is not installed or the $SCHRODINGER environment variable is not set."
            )
            return False

        current_directory = os.getcwd()
        os.chdir(self.work_dir)
        if not os.path.isdir(self.schrodinger_env):
            logger.debug(
                "Schrodinger virtual environment is not found. Creating in work directory..."
            )
            command_run = [
                self.schrodinger + "/run",
                "schrodinger_virtualenv.py",
                self.schrodinger_env,
            ]
            process = subprocess.Popen(
                command_run,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            stderr = process.communicate()
            if process.returncode == 0:
                logger.info("Schrodinger virtual environment created.")
                self.schrodinger_env = os.path.join(os.getcwd(), self.schrodinger_env)
                os.chdir(path=current_directory)
            else:
                logger.critical(
                    f"Schrodinger virtual environment failed\n{stderr[1].decode()}."
                )
                return False
        else:
            logger.info("Schrodinger virtual environment found.")
            self.schrodinger_env = os.path.join(self.work_dir, self.schrodinger_env)
            os.chdir(path=current_directory)

        if self.command is None:
            logger.critical(
                "Please provide the path to combind. If you have not installed combind, "
                "please clone this repository: https://github.com/drorlab/combind.git"
            )
            return False

        logger.info("Setting up combind specific variables.")
        combind_variables = {
            "COMBINDHOME": self.command,
            "PATH": os.pathsep.join([os.environ.get("PATH"), self.command]),
        }

        for key, value in combind_variables.items():
            os.environ[key] = value
            combind_home = os.environ.get("COMBINDHOME")
            if combind_home not in os.environ.get("PATH"):
                os.environ["PATH"] += os.pathsep + combind_home

        # Activate Combind environment
        if os.environ.get("COMBINDHOME"):
            logger.info("Combind environment activated.")
            os.chdir(current_directory)
            self.command = os.environ.get("COMBINDHOME")
            logger.info(f"Combind command: {self.command}")
        else:
            logger.critical(f"Combind environment activation failed\n{stderr[1]}.")
            return False

        return True

    @staticmethod
    def get_current() -> "CombindContext":
        r"""Get or create a singleton context."""
        # pylint: disable-next=global-statement
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


_default_combind_context: CombindContext | None = None
