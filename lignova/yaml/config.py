# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 University of Pittsburgh — Of the Commonwealth System of Higher Education
# Source: https://github.com/durrantlab/lignova

"""Implementation for yaml class to write configuration files."""

import os
from collections.abc import Iterator
from typing import Any

import yaml


class YamlConfig:
    """Class to handle YAML configuration files."""

    def __init__(self, file_path: str, data_dict: dict[str, Any] | None = None) -> None:
        """Initialize with the path to the YAML file.
        Args:
            file_path str: Path to the YAML configuration file.
            data_dict (dict | None) : Dictionary to create the YAML file if it doesn't exist.
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
        """Write the given dictionary to the YAML configuration file.
        Args:
            config (dict): Dictionary to write to the YAML file.
        """
        with open(self.file_path, "w") as file:
            yaml.dump(config, file)
        self.data_dict = config

    def validate(self) -> None:
        """Validate the YAML configuration file. This method can be overridden in subclasses to implement specific validation logic."""
        return

    def _deep_update(self, target: dict[str, Any], updates: dict[str, Any]) -> None:
        """Recursively merge updates into target without replacing whole sections.
        Args:
            target : Dictionary to update in place.
            updates: Dictionary of updates to merge.
        """
        for key, value in updates.items():
            if isinstance(value, dict) and isinstance(target.get(key), dict):
                self._deep_update(target[key], value)
            else:
                target[key] = value

    def update_config(
        self,
        updates: dict[str, Any],
        parent_key: str | tuple[str, ...] | None = None,
    ) -> None:
        """Update the YAML configuration file with the given dictionary.
            Allow both surface and deep updates.
        Args:
            updates : Dictionary containing updates to apply.
            parent_key : Key or key path to descend into. None updates at the top level.
        """
        config = self.data_dict
        if parent_key is not None:
            keys = (parent_key,) if isinstance(parent_key, str) else parent_key
            for key in keys:
                if not isinstance(config.get(key), dict):
                    config[key] = {}
                config = config[key]
        self._deep_update(config, updates)
        self.write_config(self.data_dict)
        self.validate()

    def delete_key(self, key: str) -> None:
        """Delete a key from the YAML configuration file.
        Args:
            key (str): The key to delete from the configuration.
        """
        config = self.read_config()
        if key in config:
            del config[key]
            self.write_config(config)

    def _leaf_items(self, data: dict[str, Any]) -> Iterator[tuple[str, Any]]:
        """Yield (key, value) for all non-dict values in a nested dict.
        Args:
            data (Dict): The dictionary to traverse.
        """
        for key, value in data.items():
            if isinstance(value, dict):
                # Go deeper: depth-first search
                yield from self._leaf_items(value)
            else:
                # Leaf: not a dict
                yield key, value

    def to_cli(
        self, data: dict[str, Any] | None = None, prefix: str = "--"
    ) -> list[str]:
        """Convert the YAML configuration to a command-line argument string.
        Args:
            data (Dict | None) : Dictionary to convert. If None, uses self.data_dict.
            prefix (str): Prefix for command-line arguments. Default is "--".
        """
        if data is None:
            data = self.data_dict
        arg: list[str] = []
        for key, value in self._leaf_items(data):
            if isinstance(value, bool):
                if value:
                    arg.append(f"{prefix}{key}")
            elif isinstance(value, list):
                # write the list as a list of comma separated values
                arg.append(f"{prefix}{key}")
                arg.extend(str(v) for v in value)
            elif value is None:
                continue

            else:
                arg.append(f"{prefix}{key} {value}")
        return arg
