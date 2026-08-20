#!/usr/bin/env python
# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 University of Pittsburgh — Of the Commonwealth System of Higher Education
# Source: https://github.com/durrantlab/lignova

"""
run PDB2PQR to protonate proteins
"""

import argparse
import os
import sys
from enum import StrEnum

import yaml
from loguru import logger

from lignova.preparation.meeko import Meeko
from lignova.preparation.mgltools import MglTools
from lignova.preparation.pdb2pqr import PDB2PQR
from lignova.structure.editing import (
    get_mda_universe,
    remove_hetatoms,
    write_mda_universe,
)
from lignova.yaml.meeko_config import MeekoConfig
from lignova.yaml.mgltools_config import MglToolsConfig
from lignova.yaml.protonation_config import ProtonationConfig


def strip_hetatm_pdb(pdb_in: str, pdb_out: str) -> str:
    """Write a copy of pdb_in with all HETATM records removed.
    Args:
        pdb_in : The path to the input PDB file.
        pdb_out : The path to the output PDB file with HETATM records removed.
    Returns:
        The path to the output PDB file with HETATM records removed."""
    u = get_mda_universe(pdb_in)
    protein = remove_hetatoms(u)
    if protein.n_atoms == 0:
        raise ValueError(f"{pdb_in} has no non-HETATM atoms to keep.")
    write_mda_universe(protein, pdb_out)
    return pdb_out


class PDBQT_Strategy(StrEnum):
    """An enum to specify how to generate the PDBQT file"""

    MEEKO = "meeko"
    MGLTOOLS = "mgltools"


_CAP_FLAGS = ("neutraln", "neutralc")
"""pdb2pqr options that add a hydrogen to a terminus, breaking Meeko's PQR path."""


def check_run_mgltools(method: PDBQT_Strategy) -> bool:
    """Check if MGLTools environment is available when needed (i.e., when the method is MGLTools).
    Args:
        method : The method to generate the PDBQT file.
    Returns:
        True if MGLTools environment is available or if the method is not MGLTools.
    """
    if method == PDBQT_Strategy.MGLTOOLS:
        try:
            from lignova.io import mgltools_env_exists
        except ImportError:
            return False
        try:
            return mgltools_env_exists()
        except Exception as e:
            logger.error(f"Failed to find MGLTools environment: {e}")
            return False
    return True


def has_caps(protonation_config: ProtonationConfig) -> bool:
    """Whether the protonation settings produce neutral (capped) termini.

    Args:
        protonation_config : The configuration used for the pdb2pqr run.

    Returns:
        True if either neutral-terminus flag is enabled.
    """
    general = protonation_config.data_dict.get("pdb2pqr", {}).get("general", {})
    return any(bool(general.get(flag)) for flag in _CAP_FLAGS)


def validate_methods(method: PDBQT_Strategy, methods_yaml: str) -> bool:
    config_class = MeekoConfig if method == PDBQT_Strategy.MEEKO else MglToolsConfig
    if os.path.exists(methods_yaml):
        with open(methods_yaml) as handle:
            raw_config = yaml.safe_load(handle) or {}
        if method.value not in raw_config:
            logger.error(
                "The configuration file {methods_yaml} does not contain settings "
                "for the specified method: {method_value}",
                methods_yaml=methods_yaml,
                method_value=method.value,
            )
            return False
        try:
            config_class(methods_yaml)
        except Exception as e:
            logger.error(f"Failed to load configuration from {methods_yaml}: {e}")
            return False
    return True


def working_w_defaults(pdbqt_conversion: str) -> bool:
    """Check if the provided pdbqt_conversion is the a path if so then we are working with the default configuration for the method.
    Args:
        pdbqt_conversion : The path to the configuration file for the method.
    Returns:
        True if the pdbqt_conversion is a path for non existing file, False otherwise."""
    if os.path.exists(pdbqt_conversion):
        return False
    return True


def resolve_methods(
    has_cap: bool,
    defaults: bool,
    config: MeekoConfig | MglToolsConfig,
    output_basepath: str,
    keep_hetatms: bool = False,
) -> str:
    """Resolve the method to generate the PDBQT file based on the provided configuration if possible.
    Args:
        has_cap : Whether the protein has capped termini.
        defaults : Whether the configuration is the default configuration for the method.
        config : The configuration object for the chosen method.
        output_basepath : The path to the output PQR and PDB files from PDB2PQR without the file extension.
        keep_hetatms : Whether to keep HETATM records in the output PDB file. Default is False to match MGLTools behavior under the default configuration
            when working with PDB files as input instead of PQR files
    Returns:
        The path to the generated PDBQT file.
    """
    pdb_path = f"{output_basepath}.pdb"
    pdbqt_base = output_basepath

    is_meeko = isinstance(config, MeekoConfig)

    if has_cap and is_meeko:
        input_path = pdb_path
        if not keep_hetatms:
            input_path = strip_hetatm_pdb(pdb_path, f"{output_basepath}_no_hetatm.pdb")
        else:
            logger.warning(
                "Keeping HETATM records; Meeko output may be incompatible downstream."
            )
        perc = config.data_dict["meeko"]["receptor_perception"]
        if perc.get("charge_model") == "read":
            logger.warning(
                "Capped termini force the PDB path, which has no charges to read; "
                "switching charge_model 'read' to 'gasteiger'."
            )
            perc["charge_model"] = "gasteiger"
        gen_obj = Meeko(input_path, pdbqt_base, config)
        return gen_obj.run()

    # Defaults (non-Meeko, or Meeko without caps) -> PDB path as before
    if defaults:
        input_path = pdb_path
        method_obj = Meeko if is_meeko else MglTools
        if not keep_hetatms:
            input_path = strip_hetatm_pdb(pdb_path, f"{output_basepath}_no_hetatm.pdb")
        gen_obj = method_obj(input_path, pdbqt_base, config)
        return gen_obj.run()

    # User config, no capped-Meeko conflict -> honor their input/output
    section = "meeko" if is_meeko else "mgltools"
    io = config.data_dict.get(section, {}).get("input_output", {})
    if is_meeko:
        user_input = (
            io.get("read_pqr") or io.get("read_pdb") or io.get("read_with_prody")
        )
        user_output = io.get("output_basename")
    else:
        user_input = io.get("receptor")
        user_output = io.get("outfile")

    if user_input is None:
        raise ValueError("Config has no input file set; cannot run.")

    gen_obj = (Meeko if is_meeko else MglTools)(user_input, user_output, config)
    return gen_obj.run()


def protonate(
    pdb_file: str,
    output_filepath: str,
    config_path: str,
    pdbqt_conversion: str,
    method: PDBQT_Strategy = PDBQT_Strategy.MEEKO,
    keep_hetatms: bool = False,
) -> None:
    """
    Protonate the protein using PDB2PQR and generate a PDBQT file needed for docking using the specified method.

    Args:
        pdb_file : The path to the pdb file
        output_filepath : The path to the output file
        config_path : The path to the PDB2PQR configuration file if not provided default config will be used
        pdbqt_conversion: The yaml file for the method to use for PDBQT conversion.
        method: the method to use to generate the PDBQT file
        keep_hetatms: Whether to keep HETATM records in the output PDB file. Default is False.
    """
    # check if the output directory exists if not create it
    output_dir = os.path.dirname(output_filepath)
    pdb_output_path = os.path.splitext(output_filepath)[0] + ".pdb"
    if not os.path.exists(output_dir):
        logger.warning(f"The directory {output_dir} does not exist.Creating it")
        os.makedirs(output_dir)
    protonation_config = ProtonationConfig(config_path, data_dict=None)
    protonation_config.update_config({"pdb-output": pdb_output_path})
    protonation_config.update_config({"whitespace": method == PDBQT_Strategy.MEEKO})
    pdb2pqr = PDB2PQR(
        pdb_file=pdb_file, outfile=output_filepath, config_obj=protonation_config
    )
    pdb2pqr.run()
    has_caps_flag = has_caps(protonation_config)
    if not validate_methods(method, pdbqt_conversion):
        logger.error(
            f"Invalid configuration for {method.value} in {pdbqt_conversion}. Please check the configuration file."
        )
        raise ValueError(
            f"Invalid configuration for {method.value} in {pdbqt_conversion}. Please check the configuration file."
        )
    if not check_run_mgltools(method):
        logger.error(
            "MGLTools environment not found. Please ensure MGLTools is installed or use Meeko instead."
        )
        raise EnvironmentError(
            "MGLTools environment not found. Please ensure MGLTools is installed or use Meeko instead ."
        )
    defaults = working_w_defaults(pdbqt_conversion)
    contig_obj = MeekoConfig if method == PDBQT_Strategy.MEEKO else MglToolsConfig
    m_config = contig_obj(pdbqt_conversion)
    output_basepath = os.path.splitext(output_filepath)[0]
    output_pdbqt_file = resolve_methods(
        has_caps_flag, defaults, m_config, output_basepath, keep_hetatms
    )
    return output_pdbqt_file


def run_cli():
    """Command line interface for protonation using PDB2PQR."""
    parser = argparse.ArgumentParser(description="Protonate proteins using PDB2PQR.")
    parser.add_argument(
        "-p",
        "--pdb_file",
        type=str,
        required=True,
        help="Path to the input pdb file.",
    )
    parser.add_argument(
        "-o",
        "--output_filepath",
        type=str,
        required=True,
        help="Path to the output protonated pqr file.",
    )
    parser.add_argument(
        "-c",
        "--config_path",
        type=str,
        required=False,
        default="pdb2pqr.yaml",
        help="Path to the pdb2pqr configuration file.",
    )
    parser.add_argument(
        "-m",
        "--method",
        type=PDBQT_Strategy,
        choices=list(PDBQT_Strategy),
        default=PDBQT_Strategy.MEEKO,
        help="Downstream PDBQT tool; controls PQR whitespace formatting.",
    )

    parser.add_argument(
        "-pq",
        "--pdbqt_conversion",
        type=str,
        required=True,
        help="The yaml file for the method to use for PDBQT conversion. Options are 'meeko' or 'mgltools'.",
    )

    parser.add_argument(
        "-k",
        "--keep_hetatms",
        action="store_true",
        default=False,
        help="Whether to keep HETATM records in the output PDB file. Default is False.",
    )
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    protonate(
        pdb_file=args.pdb_file,
        output_filepath=args.output_filepath,
        config_path=args.config_path,
        method=args.method,
        pdbqt_conversion=args.pdbqt_conversion,
        keep_hetatms=args.keep_hetatms,
    )


if __name__ == "__main__":
    run_cli()
