r""" Implementation of the GlideContext class containing the configuration for Glide docking."""

from typing import Optional, Union

import os

from loguru import logger

_default_glide_context: "Optional[GlideContext]" = None

DEFAULT_COMMAND = os.environ.get("SCHRODINGER", None)
DEFAULT_POSES_PER_LIG = "100"
DEFAULT_FORCEFIELD = "OPLS_2005"
DEFAULT_DOCKING_PROTOCOL = "SP"
DEFAULT_N_ENHANCED_SAMPLING = "4"
DEFAULT_LIG_PH = "7.0"
DEFAULT_LIG_PHT = "2.0"
DEFAULT_LIG_FORCEFIELD = "14"
DEFAULT_LIG_STERIOISOMERS = "32"
DEFAULT_EPIK_PH = "7.0"
DEFAULT_EPIK_PHT = "2.0"
DEFAULT_PROT_RMSD = "0.3"
DEFAULT_PROPKA_PH = "7.0"
DEFAULT_GRID_INNERBOX = "10"
DEFAULT_POSTDOCK_N_POSES = "100"
DEFAULT_WRITE_DIR = "./tmp/6oav"
DEFAULT_LIG_MAX_MW = "500"
DEFAULT_PROT_WATER_DIST = "5.0"
DEFAULT_FILLSIDECHAINS = "True"
DEFAULT_DISULFIDES = "True"
DEFAULT_REHTREAT = "True"
DEFAULT_SAMPLEWATER = "True"
DEFUALT_MINIMIZE_ADJ_H = "True"
DEFAULT_LIG_EPIK = "True"
DEFAULT_DOCKING_METHOD = "rigid"


# pylint: disable-next=too-many-instance-attributes
class GlideContext:
    r"""Singleton for Glide docking configuration."""

    # pylint: disable-next=too-many-arguments,too-many-locals
    def __init__(
        self,
        command: Union[str, None],
        forcefield: str,
        docking_protocol: Union[str, None],
        n_enhanced_sampling: Union[None, str],
        lig_ph: [str, None],
        lig_pht: Union[None, str],
        lig_forcefield: Union[str, None],
        lig_stereoisomers: Union[str, None],
        epik_ph: Union[None, str],
        epik_pht: Union[str, None],
        prot_rmsd: Union[str, None],
        propka_ph: Union[str, None],
        grid_innerbox: Union[None, str],
        postdock_nposes: Union[None, str],
        write_dir: Union[str, None],
        lig_max_mw: Union[str, None],
        prot_watdist: Union[str, None],
        poses_per_lig: Union[str, None],
        fillsidechains: Union[bool, None],
        disulfides: Union[bool, None],
        rehtreat: Union[bool, None],
        samplewater: Union[bool, None],
        minimize_adj_h: Union[bool, None],
        lig_epik: Union[bool, None],
        docking_method: Union[str, None],
    ):
        if not os.environ.get("SCHRODINGER") or command is None:
            logger.critical(
                "Schrödinger is not installed or the $SCHRODINGER environment variable is not set."
            )
            raise OSError(
                "Schrödinger is not installed or the $SCHRODINGER environment variable is not set."
            )
        self.command = command
        self.forcefield = forcefield
        self.docking_protocol = docking_protocol
        self.n_enhanced_sampling = n_enhanced_sampling
        self.lig_ph = lig_ph
        self.lig_pht = lig_pht
        self.lig_forcefield = lig_forcefield
        self.lig_stereoisomers = lig_stereoisomers
        self.epik_ph = epik_ph
        self.epik_pht = epik_pht
        self.prot_rmsd = prot_rmsd
        self.propka_ph = propka_ph
        self.grid_innerbox = grid_innerbox
        self.postdock_nposes = postdock_nposes
        self.write_dir = write_dir
        self.lig_max_mw = lig_max_mw
        self.prot_watdist = prot_watdist
        self.poses_per_lig = poses_per_lig
        self.fillsidechains = fillsidechains
        self.disulfides = disulfides
        self.rehtreat = rehtreat
        self.samplewater = samplewater
        self.minimize_adj_h = minimize_adj_h
        self.lig_epik = lig_epik
        self.docking_method = docking_method

    @staticmethod
    def get_current() -> "GlideContext":
        r"""Get or create a singleton context."""

        # pylint: disable-next=global-statement
        global _default_glide_context

        if _default_glide_context is None:
            _default_glide_context = GlideContext(
                command=DEFAULT_COMMAND,
                forcefield=DEFAULT_FORCEFIELD,
                docking_protocol=DEFAULT_DOCKING_PROTOCOL,
                n_enhanced_sampling=DEFAULT_N_ENHANCED_SAMPLING,
                lig_ph=DEFAULT_LIG_PH,
                lig_pht=DEFAULT_LIG_PHT,
                lig_forcefield=DEFAULT_LIG_FORCEFIELD,
                lig_stereoisomers=DEFAULT_LIG_STERIOISOMERS,
                epik_ph=DEFAULT_EPIK_PH,
                epik_pht=DEFAULT_EPIK_PHT,
                prot_rmsd=DEFAULT_PROT_RMSD,
                propka_ph=DEFAULT_PROPKA_PH,
                grid_innerbox=DEFAULT_GRID_INNERBOX,
                postdock_nposes=DEFAULT_POSTDOCK_N_POSES,
                write_dir=DEFAULT_WRITE_DIR,
                lig_max_mw=DEFAULT_LIG_MAX_MW,
                prot_watdist=DEFAULT_PROT_WATER_DIST,
                poses_per_lig=DEFAULT_POSES_PER_LIG,
                fillsidechains=DEFAULT_FILLSIDECHAINS,
                disulfides=DEFAULT_DISULFIDES,
                rehtreat=DEFAULT_REHTREAT,
                samplewater=DEFAULT_SAMPLEWATER,
                minimize_adj_h=DEFUALT_MINIMIZE_ADJ_H,
                lig_epik=DEFAULT_LIG_EPIK,
                docking_method=DEFAULT_DOCKING_METHOD,
            )
        return _default_glide_context

    @staticmethod
    def set_current(context: "GlideContext") -> None:
        """Set the current context."""
        global _default_glide_context
        _default_glide_context = context
