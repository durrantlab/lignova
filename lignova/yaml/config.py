r"""Implementation for yaml class to write configuration files."""

from typing import Any, Iterator

import os
from venv import logger

import yaml


class YamlConfig:
    """Class to handle YAML configuration files."""

    def __init__(self, file_path: str, data_dict: dict[str, Any] | None = None) -> None:
        """Initialize with the path to the YAML file.
        Args:
            file_path (str): Path to the YAML configuration file.
            data_dict (dict, optional): Dictionary to create the YAML file if it doesn't exist.
        """
        self.file_path = file_path
        if not os.path.exists(file_path) and data_dict is None:
            raise ValueError("Either file_path or dictionary must be provided.")
        elif data_dict is not None:
            self.data_dict = data_dict
            self.write_config(data_dict)
        else:
            # Ensure data_dict is populated by reading existing file
            self.data_dict = self.read_config()

    def read_config(self) -> dict[str, Any]:
        """Read the YAML configuration file and return its contents as a dictionary."""
        with open(self.file_path, "r") as file:
            config = yaml.safe_load(file) or {}
        return config

    def write_config(self, config: dict[str, Any]) -> None:
        """Write the given dictionary to the YAML configuration file."""
        with open(self.file_path, "w") as file:
            yaml.dump(config, file)
        self.data_dict = config

    def update_config(
        self,
        updates: dict[str, Any],
        nested: bool = False,
        parent_key: str | None = None,
    ) -> None:
        """Update the YAML configuration file with the given dictionary
            Allow both surface and deep updates.
        Args:
            updates (dict): Dictionary containing updates to apply.
            nested (bool): Whether to update a nested dictionary. Default is False.
            parent_key (str): The key of the parent dictionary to update if nested is True. Default is None.
        """
        config = self.data_dict
        if nested and parent_key is not None:
            if parent_key in config and isinstance(config[parent_key], dict):
                config[parent_key].update(updates)
            else:
                config[parent_key] = updates
        else:
            config.update(updates)
        self.write_config(config)

    def delete_key(self, key: str) -> None:
        """Delete a key from the YAML configuration file."""
        config = self.read_config()
        if key in config:
            del config[key]
            self.write_config(config)

    def _leaf_items(self, data: dict[str, Any]) -> Iterator[tuple[str, Any]]:
        """Yield (key, value) for all non-dict values in a nested dict."""
        for key, value in data.items():
            if isinstance(value, dict):
                # Go deeper: depth-first search
                yield from self._leaf_items(value)
            else:
                # Leaf: not a dict
                yield key, value

    def to_cli(self, data: dict[str, Any] | None = None) -> list[str]:
        """Convert the YAML configuration to a command-line argument string.
        Args:
            data (dict, optional): Dictionary to convert. If None, uses self.data_dict.
        """
        if data is None:
            data = self.data_dict
        arg: list[str] = []
        for key, value in self._leaf_items(data):
            if isinstance(value, bool):
                if value:
                    arg.append(f"--{key}")
            elif isinstance(value, list):
                # write the list as a list of comma separated values
                arg.append(f"--{key}")
                arg.extend(str(v) for v in value)
            elif value is None:
                continue

            else:
                arg.append(f"--{key} {value}")
        return arg
