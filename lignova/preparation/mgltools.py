# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 University of Pittsburgh — Of the Commonwealth System of Higher Education
# Source: https://github.com/durrantlab/lignova

"""Implementation for receptor preparation using MGLTools pre-docking."""

import os
import re
import shutil
import string
from dataclasses import dataclass, replace

from loguru import logger

from lignova.io import run_mgltools_command
from lignova.yaml.mgltools_config import MglToolsConfig

_FLOAT_RE = re.compile(r"-?\d+\.\d+")
"""Matches any signed decimal, used to locate the numeric tail of a PQR line."""

_RECORD_RE = re.compile(r"^(HETATM|ATOM)\s*(\d+)")
"""Matches the record type and serial at the start of an atom line."""

_WATER_RESNAMES = frozenset({"HOH", "WAT", "TIP", "TIP3", "SOL"})
"""Residue names treated as water when stripping HETATM lines."""

_CHAIN_COL = 21
"""Zero-based index of the PDB chain identifier (spec column 22)."""

_RESSEQ_SLICE = slice(22, 27)
"""Zero-based slice covering resSeq plus insertion code (spec columns 23-27)."""

_WIDE_RESSEQ_DIGITS = 4
"""Residue-number width at which the chain column becomes adjacent."""


@dataclass(frozen=True, slots=True)
class PqrAtom:
    """The fields of a single PQR atom record, parsed out of a malformed line."""

    record: str
    """Record type, either ATOM or HETATM. PDB columns 1-6."""

    serial: int
    """Atom serial number, unique within the file. PDB columns 7-11."""

    atom_name: str
    """Atom name, e.g. CA, OXT, HD11. PDB columns 13-16."""

    resname: str
    """Residue name, e.g. LEU, HOH, M3A. PDB columns 18-20."""

    chain: str
    """Single-character chain identifier, or a space when absent. PDB column 22."""

    resseq: int
    """Residue sequence number. PDB columns 23-26."""

    icode: str
    """Insertion code for residues sharing a number, or a space. PDB column 27."""

    x: float
    """Orthogonal coordinate in angstroms. PDB columns 31-38, format 8.3f."""

    y: float
    """Orthogonal coordinate in angstroms. PDB columns 39-46, format 8.3f."""

    z: float
    """Orthogonal coordinate in angstroms. PDB columns 47-54, format 8.3f."""

    charge: float
    """Partial atomic charge in electrons, in the PDB occupancy column."""

    radius: float
    """Atomic radius in angstroms, in the PDB B-factor column."""


def _is_atom_record(line: str) -> bool:
    """Whether a line is an ATOM/HETATM record long enough to carry a resSeq.

    Args:
        line : A raw line from a PQR, PDB or PDBQT file.

    Returns:
        True if the line has positional fields through column 27.
    """
    return line.startswith(("ATOM", "HETATM")) and len(line) > _RESSEQ_SLICE.stop - 1


def _fmt_float(val: float, width: int, max_dec: int, need_space: bool = True) -> str:
    """Format a float into a fixed-width column, shedding decimals as needed.

    Args:
        val : Value to format.
        width : Exact column width the result must occupy.
        max_dec : Preferred number of decimal places; reduced until it fits.
        need_space : Require a leading space so the field cannot abut the
            previous one. Disable for the first field of a group.

    Returns:
        The formatted value, width-exact where possible.
    """
    for dec in range(max_dec, 0, -1):
        s = f"{val:{width}.{dec}f}"
        if len(s) == width and (not need_space or s[0] == " "):
            return s
    s = f"{val:{width}.0f}"
    if len(s) == width:
        return s
    return f"{val:{width}.{max_dec}f}"


def _is_water_line(line: str) -> bool:
    """Whether an atom line belongs to a water molecule.

    Args:
        line : A raw ATOM/HETATM line.

    Returns:
        True if any whitespace-separated token is a known water residue name.
    """
    return bool(set(line.split()) & _WATER_RESNAMES)


def _should_strip_hetatm(line: str, cleanup: str) -> bool:
    """Whether a HETATM line would be removed by prepare_receptor4.py's -U flag.

    Operates on the raw line via whitespace splitting, so it is format-agnostic
    and safe for PDB as well as PQR input.

    Args:
        line : A raw ATOM/HETATM line.
        cleanup : The cleanup mode passed to prepare_receptor4.py (-U flag).

    Returns:
        True if the line should be dropped.
    """
    if not line.startswith("HETATM"):
        return False
    if _is_water_line(line):
        return "waters" in cleanup
    return "nonstdres" in cleanup


def strip_hetatm_lines(lines: list[str], cleanup: str) -> tuple[list[str], int, int]:
    """Drop HETATM lines that the -U cleanup mode would remove anyway.

    Stripping these up front matters because MGLTools often cannot parse HETATM
    records at all, so it fails before the cleanup would have applied. No line
    is reformatted, so this is safe for any fixed-column format.

    Args:
        lines : The file's lines, without terminators.
        cleanup : The cleanup mode passed to prepare_receptor4.py (-U flag).

    Returns:
        Tuple of (kept lines, waters removed, other HETATMs removed).
    """
    kept: list[str] = []
    removed_waters = 0
    removed_hetatm = 0

    for line in lines:
        if _should_strip_hetatm(line, cleanup):
            if _is_water_line(line):
                removed_waters += 1
            else:
                removed_hetatm += 1
            continue
        kept.append(line)

    return kept, removed_waters, removed_hetatm


def has_fused_resseq(lines: list[str]) -> bool:
    """Whether a chain ID abuts a four-digit residue number on any atom line.

    Three-digit numbers right-pad into columns 23-26 and leave column 22
    visibly separate, so they parse correctly and need no intervention.

    Args:
        lines : The file's lines.

    Returns:
        True if at least one line would be corrupted by Mgltools PQR parser.
    """
    for line in lines:
        if not _is_atom_record(line) or line[_CHAIN_COL] == " ":
            continue
        resseq = line[22:26].strip()
        if len(resseq) >= _WIDE_RESSEQ_DIGITS and resseq.isdigit():
            return True
    return False


def residue_sequence(lines: list[str]) -> list[tuple[str, str, str]]:
    """Collapse atom lines into one entry per contiguous residue block.

    Args:
        lines : The file's lines.

    Returns:
        Ordered list of (chain, resSeq-plus-icode, resname), one per residue.
    """
    sequence: list[tuple[str, str, str]] = []
    previous: tuple[str, str, str] | None = None
    for line in lines:
        if not _is_atom_record(line):
            continue
        key = (line[_CHAIN_COL], line[_RESSEQ_SLICE], line[17:20].strip())
        if key != previous:
            sequence.append(key)
            previous = key
    return sequence


def blank_chain_column(lines: list[str]) -> tuple[list[str], int]:
    """Replace the chain identifier with a space on every atom line.

    Vacating column 22 leaves Mgltools nothing to absorb into the residue number,
    so four-digit numbering survives intact.

    Args:
        lines : The file's lines.

    Returns:
        Tuple of (rewritten lines, number of lines changed).
    """
    out: list[str] = []
    changed = 0
    for line in lines:
        if _is_atom_record(line) and line[_CHAIN_COL] != " ":
            line = line[:_CHAIN_COL] + " " + line[_CHAIN_COL + 1 :]
            changed += 1
        out.append(line)
    return out, changed


def align_chains(source_lines: list[str], output_lines: list[str]) -> list[str]:
    """Recover the chain of each output residue by position, not by number.

    Args:
        source_lines : Lines of the input receptor, before the chain was removed.
        output_lines : Lines of the PDBQT that MGLTools produced.

    Raises:
        ValueError: If an output residue cannot be found in the remaining input
            sequence, which would mean the order was not preserved.

    Returns:
        Chain identifier for each output residue block, in order.
    """
    source = residue_sequence(source_lines)
    output = residue_sequence(output_lines)

    assigned: list[str] = []
    cursor = 0
    for _, out_key, out_resname in output:
        while cursor < len(source) and source[cursor][1] != out_key:
            cursor += 1
        if cursor >= len(source):
            raise ValueError(
                f"Cannot align output residue {out_key.strip()!r} ({out_resname}) against the input sequence; the order was not preserved."
            )
        if source[cursor][2] != out_resname:
            logger.debug(
                "Residue renamed during preparation: {old} to  {new} at {key}",
                old=source[cursor][2],
                new=out_resname,
                key=out_key.strip(),
            )
        assigned.append(source[cursor][0])
        cursor += 1
    return assigned


def restore_chain_column(
    lines: list[str], assigned: list[str]
) -> tuple[list[str], int]:
    """Write chain identifiers back onto atom lines, block by block.

    Overwrites whatever Mgltools put in the chain column, which is 'U' (the first
    character of its 'UNK' placeholder) when no chain was present in the input.

    Args:
        lines : The output file's lines.
        assigned : Chain identifier per residue block, from :func:`align_chains`.

    Returns:
        Tuple of (rewritten lines, number of atoms restored).
    """
    out: list[str] = []
    block = -1
    previous: tuple[str, str, str] | None = None
    restored = 0

    for line in lines:
        if _is_atom_record(line):
            key = (line[_CHAIN_COL], line[_RESSEQ_SLICE], line[17:20].strip())
            if key != previous:
                block += 1
                previous = key
            if block < len(assigned):
                line = line[:_CHAIN_COL] + assigned[block] + line[_CHAIN_COL + 1 :]
                restored += 1
        out.append(line)

    return out, restored


def _split_chain_and_resseq(field: str) -> tuple[str, str]:
    """Separate a fused chain/residue-number field such as A1086.

    Args:
        field : The single token standing in for chain and residue number.

    Returns:
        Tuple of (chain, resseq). Chain is a space when it cannot be identified.
    """
    resseq_part = field[1:].lstrip("-").rstrip(string.ascii_letters)
    if len(field) >= 2 and field[0].isalpha() and resseq_part.isdigit():
        return field[0], field[1:]
    return " ", field


def _split_resseq_and_icode(resseq: str) -> tuple[str, str]:
    """Peel a trailing insertion code off a residue number.

    Args:
        resseq : Residue number, possibly with a trailing letter.

    Returns:
        Tuple of (resseq, icode). Icode is a space when absent.
    """
    if resseq and resseq[-1].isalpha():
        return resseq[:-1], resseq[-1]
    return resseq, " "


def parse_pqr_atom_line(line: str) -> PqrAtom | None:
    """Parse a PQR atom line whose columns may have run together.

    Works backwards from the five trailing floats (x, y, z, charge, radius),
    then splits the remaining prefix on whitespace. This tolerates the fused
    columns PDB2PQR emits for large negative numbers as well as the extra
    padding added by --whitespace.

    WARNING: PQR only. A PDB line's numeric tail is x, y, z, occupancy,
    B-factor, so this would silently relabel occupancy as charge.

    Args:
        line : A raw ATOM/HETATM line.

    Returns:
        The parsed fields, or None if the line does not look like a PQR atom
        record (caller should then pass the line through untouched).
    """
    floats = _FLOAT_RE.findall(line)
    if len(floats) < 5:
        return None

    x, y, z = float(floats[-5]), float(floats[-4]), float(floats[-3])
    charge, radius = float(floats[-2]), float(floats[-1])

    # Locate where the numeric tail begins so the prefix can be split separately
    tail_parts = [re.escape(f) for f in floats[-5:]]
    tail_re = re.compile(r"[\s\-]*".join(tail_parts).replace(r"\-", "-"))
    tail_match = tail_re.search(line)
    if not tail_match:
        return None

    prefix = line[: tail_match.start()]
    rec_match = _RECORD_RE.match(prefix)
    if not rec_match:
        return None

    remainder = prefix[rec_match.end() :].split()
    if len(remainder) < 3:
        return None

    atom_name, resname = remainder[0], remainder[1]
    if len(remainder) >= 4:
        chain, resseq = remainder[2], remainder[3]
    else:
        chain, resseq = _split_chain_and_resseq(remainder[2])
    resseq, icode = _split_resseq_and_icode(resseq)

    try:
        resseq_int = int(resseq)
    except ValueError:
        logger.debug("Unparseable residue number {resseq!r}", resseq=resseq)
        return None

    return PqrAtom(
        record=rec_match.group(1),
        serial=int(rec_match.group(2)),
        atom_name=atom_name,
        resname=resname,
        chain=chain,
        resseq=resseq_int,
        icode=icode,
        x=x,
        y=y,
        z=z,
        charge=charge,
        radius=radius,
    )


def format_pqr_atom_line(atom: PqrAtom) -> str:
    """Render parsed fields back into strict fixed-width PDB columns.

    Args:
        atom : The parsed atom fields.

    Returns:
        A single line MGLTools' fixed-width parser can read.
    """
    if len(atom.atom_name) >= 4 or atom.atom_name[0].isdigit():
        atom_name_fmt = f"{atom.atom_name:<4s}"
    else:
        atom_name_fmt = f" {atom.atom_name:<3s}"

    numeric_tail = (
        _fmt_float(atom.x, 8, 3, need_space=False)
        + _fmt_float(atom.y, 8, 3)
        + _fmt_float(atom.z, 8, 3)
        + _fmt_float(atom.charge, 8, 4)
        + _fmt_float(atom.radius, 7, 4)
    )
    return (
        f"{atom.record:<6s}"
        f"{atom.serial:>5d} "
        f"{atom_name_fmt}"
        f" "
        f"{atom.resname:>3s}"
        f" "
        f"{atom.chain:1s}"
        f"{atom.resseq:>4d}"
        f"{atom.icode:1s}"
        f"   " + numeric_tail
    )


class MglTools:
    """Class to handle receptor preparation using MGLTools."""

    _SCRIPT = "prepare_receptor4.py"
    """The script name for MGLTools' receptor preparation tool."""

    _ALLOWED_INPUT_EXTENSIONS = MglToolsConfig._ALLOWED_INPUT_EXTENSIONS
    """Allowed input file extensions, kept in sync with the config's schema."""

    _HETATM_CLEANUP_KEYWORDS = ("waters", "nonstdres")
    """Cleanup keywords that imply HETATM records will be dropped."""

    def __init__(
        self, input_file: str, output_basename: str | None, config_obj: MglToolsConfig
    ) -> None:
        """Initialize MglTools with a given configuration object.

        Args:
            input_file : Path to the input receptor file.
            output_basename : Basename used for the generated receptor files. If None, defaults to same as input file without extension.
            config_obj : Configuration object for prepare_receptor4.py.
        """
        self.input_file = os.path.abspath(input_file)
        # check the input file exists
        if not os.path.exists(self.input_file):
            raise FileNotFoundError(f"Input file {input_file} does not exist.")
        extension = os.path.splitext(self.input_file)[1].lower()
        if extension not in self._ALLOWED_INPUT_EXTENSIONS:
            raise ValueError(
                f"Input file {input_file} must be one of {sorted(self._ALLOWED_INPUT_EXTENSIONS)}."
            )
        self.config = config_obj
        if output_basename is None:
            self.output_basename = os.path.splitext(self.input_file)[0]
        else:
            self.output_basename = (
                os.path.splitext(output_basename)[0]
                if output_basename.endswith(".pdbqt")
                else output_basename
            )
        self._write_paths()

    def _write_paths(self) -> None:
        """Write the input and output paths into the configuration and revalidate.

        Warns if the configuration already held a different receptor or outfile.
        """
        input_cfg = self.config.data_dict.get("mgltools", {}).get("input_output", {})

        existing_input = input_cfg.get("receptor")
        if existing_input is not None and existing_input != self.input_file:
            logger.warning(
                "Config 'receptor' is '{existing_input}' but '{input_file}' was provided. Overwriting with provided value.",
                existing_input=existing_input,
                input_file=self.input_file,
            )
        existing_outfile = input_cfg.get("outfile")
        if existing_outfile is not None and existing_outfile != self.pdbqt_file:
            logger.warning(
                "Config 'outfile' is '{existing_outfile}' but '{outfile}' was provided. Overwriting with provided value.",
                existing_outfile=existing_outfile,
                outfile=self.pdbqt_file,
            )

        input_cfg["receptor"] = self.input_file
        input_cfg["outfile"] = self.pdbqt_file

        self.config.validate()

    @property
    def pdbqt_file(self) -> str:
        """Path to the PDBQT that prepare_receptor4.py writes."""
        return f"{self.output_basename}.pdbqt"

    @property
    def _input_extension(self) -> str:
        """Lowercased extension of the input receptor, including the dot."""
        return os.path.splitext(self.input_file)[1].lower()

    @property
    def _cleanup_mode(self) -> str:
        """The configured -U cleanup mode."""
        return (
            self.config.data_dict.get("mgltools", {})
            .get("receptor_perception", {})
            .get("cleanup", "nphs_lps_waters_nonstdres")
        )

    @property
    def _cleanup_drops_hetatm(self) -> bool:
        """Whether the configured cleanup mode removes any HETATM records."""
        return any(kw in self._cleanup_mode for kw in self._HETATM_CLEANUP_KEYWORDS)

    @staticmethod
    def _read_lines(path: str) -> list[str]:
        """Read a file into lines with terminators preserved.

        Args:
            path : File to read.

        Returns:
            The file's lines.
        """
        with open(path) as fh:
            return fh.readlines()

    @staticmethod
    def _write_lines(path: str, lines: list[str]) -> None:
        """Write lines back to a file.

        Args:
            path : File to write.
            lines : Lines to write, terminators included.
        """
        with open(path, "w") as fh:
            fh.writelines(lines)

    @staticmethod
    def _backup_original(file_path: str) -> str:
        """Copy the original receptor aside before it is overwritten.

        Existing backups are never replaced, so repeated passes cannot clobber
        the pristine copy with already-rewritten content.

        Args:
            file_path : Path to the file about to be rewritten.

        Returns:
            Path to the backup file.
        """
        stem, extension = os.path.splitext(file_path)
        backup_path = f"{stem}_org{extension}"
        if not os.path.exists(backup_path):
            shutil.copy2(file_path, backup_path)
            logger.info(
                "Backed up original receptor to {backup_path}",
                backup_path=backup_path,
            )
        return backup_path

    def _delete_chain(self) -> list[str] | None:
        """Remove the chain column when Mgltools PQR parser would corrupt it.

        Returns:
            The input's lines as they were before blanking, needed later to
            restore the chain identifiers, or None if no workaround applied.
        """
        if self._input_extension != ".pqr":
            return None

        lines = self._read_lines(self.input_file)
        if not has_fused_resseq(lines):
            return None

        logger.warning(
            "{path} has four-digit residue numbers adjacent to a chain ID."
            " Mgltools would absorb the chain letter and truncate the residue number, so the chain column is removed and restored afterwards.",
            path=self.input_file,
        )

        original = list(lines)
        self._backup_original(self.input_file)
        removed, changed = blank_chain_column(lines)
        self._write_lines(self.input_file, removed)
        logger.info(
            "Blanked the chain column on {changed} lines of {path}",
            changed=changed,
            path=self.input_file,
        )
        return original

    def _restore_chain(self, source_lines: list[str]) -> None:
        """Restore chain identifiers in the generated PDBQT by residue order.

        Args:
            source_lines : The input's lines from before the chain was removed.
        """
        output_lines = self._read_lines(self.pdbqt_file)
        assigned = align_chains(source_lines, output_lines)
        restored_lines, restored = restore_chain_column(output_lines, assigned)
        self._write_lines(self.pdbqt_file, restored_lines)
        logger.info(
            "Restored the chain identifier on {restored} atoms of {path}",
            restored=restored,
            path=self.pdbqt_file,
        )

    def _check_resseq(self) -> None:
        """Verify the generated PDBQT has usable residue numbering."""
        names_per_key: dict[tuple[str, str], set[str]] = {}
        for line in self._read_lines(self.pdbqt_file):
            if not _is_atom_record(line):
                continue
            resseq = line[22:26].strip()
            if not resseq.lstrip("-").isdigit():
                raise RuntimeError(
                    f"{self.pdbqt_file} has a non-numeric residue number {resseq!r}; MGLTools truncated the resSeq field."
                )
            key = (line[_CHAIN_COL], resseq)
            names_per_key.setdefault(key, set()).add(line[17:20].strip())

        collapsed = {k: v for k, v in names_per_key.items() if len(v) > 1}
        if collapsed:
            sample = list(collapsed.items())[:3]
            raise RuntimeError(
                f"{self.pdbqt_file} has residue numbers shared by different residue names (e.g. {sample}); numbering was corrupted."
            )

    def _strip_hetatm(self, file_path: str, cleanup: str) -> bool:
        """Remove HETATM records without reformatting anything.
        Args:
            file_path : Path to the receptor file to filter.
            cleanup : The cleanup mode passed to prepare_receptor4.py (-U flag).

        Returns:
            True if the file was modified.
        """
        with open(file_path) as fh:
            lines = fh.read().splitlines()

        kept, removed_waters, removed_hetatm = strip_hetatm_lines(lines, cleanup)
        if not (removed_waters or removed_hetatm):
            return False

        logger.warning(
            "Cleanup '{cleanup}' implies HETATM removal, so {removed_waters} water "
            "and {removed_hetatm} other HETATM records are being stripped from "
            "{path} before preparation.",
            cleanup=cleanup,
            removed_waters=removed_waters,
            removed_hetatm=removed_hetatm,
            path=file_path,
        )
        self._backup_original(file_path)
        with open(file_path, "w") as fh:
            fh.write("\n".join(kept) + "\n")
        return True

    def _rewrite_lines(
        self,
        lines: list[str],
        cleanup: str,
        blank_chain: bool = False,
    ) -> tuple[list[str], int, int]:
        """Rewrite atom lines into fixed columns, dropping HETATMs per cleanup.

        Args:
            lines : The file's lines, without terminators.
            cleanup : The cleanup mode passed to prepare_receptor4.py (-U flag).
            blank_chain : Emit a space in the chain column instead of the parsed
                identifier

        Returns:
            Tuple of (rewritten lines, waters removed, other HETATMs removed).
        """
        kept, removed_waters, removed_hetatm = strip_hetatm_lines(lines, cleanup)

        fixed: list[str] = []
        for line in kept:
            if not line.startswith(("ATOM", "HETATM")):
                fixed.append(line)
                continue

            atom = parse_pqr_atom_line(line)
            if atom is None:
                logger.debug(
                    "Could not parse atom line in , passing through: {line!r}",
                    line=line,
                )
                fixed.append(line)
                continue

            if blank_chain:
                atom = replace(atom, chain=" ")

            for label, val in (("x", atom.x), ("y", atom.y), ("z", atom.z)):
                if abs(val) >= 1000.0:
                    logger.warning(
                        "Coordinate {label}={val:.3f} overflows PDB 8.3f column in (serial {serial})",
                        label=label,
                        val=val,
                        serial=atom.serial,
                    )

            fixed.append(format_pqr_atom_line(atom))

        return fixed, removed_waters, removed_hetatm

    def _fix_pqr_spacing(
        self, pqr_path: str, cleanup: str, blank_chain: bool = False
    ) -> str:
        """Rewrite a PQR file with strict PDB column widths.

        PQR files from PDB2PQR can have columns that run together when
        coordinates or charges are large negative numbers (e.g.
        -18.713-100.168). MGLTools uses a strict
        fixed-width PDB parser that fails on both.

        The original file is backed up to *_org.pqr before overwriting.

        Args:
            pqr_path : Path to the PQR file to fix.
            cleanup : The cleanup mode passed to prepare_receptor4.py (-U flag).
                Used to determine which HETATM lines to strip.
            blank_chain : Also vacate the chain column, for the MGLTools
                four-digit residue-number workaround.

        Returns:
            Path to the overwritten fixed PQR file.
        """
        self._backup_original(pqr_path)

        with open(pqr_path) as fh:
            lines = fh.read().splitlines()

        fixed, _, _ = self._rewrite_lines(lines, cleanup, blank_chain=blank_chain)

        with open(pqr_path, "w") as fh:
            fh.write("\n".join(fixed) + "\n")

        logger.info("Fixed PQR spacing written to {path}", path=pqr_path)
        return pqr_path

    def _repair_input(self) -> bool:
        """Attempt a format repair appropriate to the input type.

        Returns:
            True if the input was modified and a retry is worth attempting.
        """
        extension = self._input_extension
        if extension == ".pqr":
            logger.info("Input is a PQR. Fixing column spacing and retrying.")
            self._fix_pqr_spacing(self.input_file, cleanup=self._cleanup_mode)
            return True
        return False

    def run(self) -> str:
        """Run the MGLTools receptor preparation process.

        Returns:
            Path to the generated PDBQT file.
        """
        source_lines = self._delete_chain()
        if self._cleanup_drops_hetatm:
            self._strip_hetatm(self.input_file, cleanup=self._cleanup_mode)
        mgltools_config = self.config.to_cli()
        cmd = ["pythonsh", self._SCRIPT]
        outdir = os.path.dirname(self.output_basename)
        if outdir and not os.path.exists(outdir):
            logger.debug(
                "Output directory {outdir} does not exist. Creating it.", outdir=outdir
            )
            os.makedirs(outdir)
        if mgltools_config:
            cmd.extend(mgltools_config)

        logger.debug("Running MGLTools with command: {cmd}", cmd=" ".join(cmd))
        result = run_mgltools_command(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            logger.warning(
                "{script} failed on first attempt:\n{stderr}",
                script=self._SCRIPT,
                stderr=result.stderr,
            )
            if not self._repair_input():
                logger.error(
                    "{script} stderr:\n{stderr}\nstdout:\n{stdout}",
                    script=self._SCRIPT,
                    stderr=result.stderr,
                    stdout=result.stdout,
                )
                raise RuntimeError(
                    f"{self._SCRIPT} failed with exit code {result.returncode} "
                    f"and no repair is available for a "
                    f"'{self._input_extension}' input: {result.stderr}"
                )

            # a partial output from the failed attempt must not satisfy the
            # existence check below
            if os.path.exists(self.pdbqt_file):
                os.remove(self.pdbqt_file)
            result = run_mgltools_command(cmd, capture_output=True, text=True)
            if result.returncode != 0:
                logger.error(
                    "{script} stderr:\n{stderr}\nstdout:\n{stdout}",
                    script=self._SCRIPT,
                    stderr=result.stderr,
                    stdout=result.stdout,
                )
                raise RuntimeError(
                    f"Even after repairing the input, {self._SCRIPT} failed with "
                    f"exit code {result.returncode}"
                )

        if result.stdout:
            logger.debug(
                "{script} stdout:\n{stdout}",
                script=self._SCRIPT,
                stdout=result.stdout,
            )

        if not os.path.exists(self.pdbqt_file):
            raise RuntimeError(
                f"{self._SCRIPT} completed but output not found at {self.pdbqt_file}"
            )

        if source_lines is not None:
            self._restore_chain(source_lines)
        self._check_resseq()

        logger.info(
            "PDBQT receptor written to {pdbqt_file}", pdbqt_file=self.pdbqt_file
        )
        return self.pdbqt_file
