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
        self.schrodinger_env = schrodinger_env

        # Check if Schrödinger virtual environment is not found and create it
        if schrodinger_env and schrodinger_env not in os.listdir(self.work_dir):
            logger.debug(
                "Schrödinger virtual environment is not found. Creating in work directory..."
            )
            command = [
                self.schrodinger + "/run",
                "schrodinger_virtualenv.py",
                os.path.join(work_dir, schrodinger_env),
            ]
            process = subprocess.Popen(
                command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True
            )
            stderr = process.communicate()
            if process.returncode == 0:
                logger.info("Schrodinger virtual environment created.")
            else:
                error_message = (
                    f"Schrodinger virtual environment creation failed\n{stderr[1]}."
                )
                logger.critical(error_message)
                raise NotImplementedError(error_message)

        # Check if command is provided
        if command is None:
            error_message = (
                "Please provide the path to combind. If you have not installed combind,"
                "please clone this repository:https://github.com/drorlab/combind.git"
            )
            logger.critical(error_message)
            raise NotImplementedError(error_message)

        logger.info("Setting up combind specific variables.")

        # Update shebang line in setup.sh
        activate_combind = os.path.join(command, "setup.sh")
        with open(activate_combind, "r", encoding="utf-8") as file:
            script_content = file.read()

        shebang_line = (
            "#!/bin/bash\n" if not script_content.startswith("#!/bin/bash") else ""
        )
        modified_script_content = shebang_line + script_content

        with open(activate_combind, "w", encoding="utf-8") as file:
            file.write(modified_script_content)

        # Run chmod +x setup.sh
        chmod_command = ["chmod", "u+x", activate_combind]
        process = subprocess.Popen(
            chmod_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE
        )
        stderr = process.communicate()

        if process.returncode == 0:
            logger.info("chmod DONE.")
            current_directory = os.getcwd()
            os.chdir(command)
            source_command = f"source setup.sh && echo COMBINDHOME=$COMBINDHOME"
            process = subprocess.Popen(
                source_command,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                shell=True,
                executable="/bin/bash",
            )
            stdout, stderr = process.communicate()
            output_lines = stdout.decode().splitlines()
            for environment_variable in output_lines:
                key, value = environment_variable.split("=")
                print(key, value)
                os.environ[key] = value

            # Get the value of COMBINDHOME from the environment
            combind_home = os.environ.get("COMBINDHOME")
            print("COMBINDHOME:", combind_home)
            print("PATH:", os.environ.get("PATH"))
        # Activate Combind environment
        if os.environ.get("COMBINDHOME"):
            logger.info("Combind environment activated.")
            os.chdir(current_directory)
            print(os.getcwd())
            self.command = os.environ.get("COMBINDHOME")
            activate = [
                "source",
                os.path.join(work_dir, self.schrodinger_env, "bin/activate"),
            ]
            process = subprocess.Popen(
                activate, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True
            )
            stderr = process.communicate()

            if process.returncode != 0:
                error_message = (
                    f"Failed to activate Schrodinger virtual environment\n{stderr[1]}."
                )
                logger.critical(error_message)
                raise NotImplementedError(error_message)
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
