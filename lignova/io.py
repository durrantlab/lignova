# SPDX-License-Identifier: Apache-2.0
# Copyright 2026 University of Pittsburgh — Of the Commonwealth System of Higher Education
# Source: https://github.com/durrantlab/lignova

r"""Implements methods for reading and writing files."""

import gzip
import os
import subprocess
from contextlib import contextmanager
from tempfile import NamedTemporaryFile

import MDAnalysis as mda
from loguru import logger


def get_file_ext(file_path: str) -> str:
    r"""Get the file extension from a file path.

    Args:
        file_path : Path to file.

    Returns:
        File extension.
    """
    return os.path.splitext(file_path)[-1]


def write_text(
    text: str | mda.Universe,
    write_path: None | str = None,
    file_ext: None | str = None,
) -> str:
    r"""General method to write text files.

    Args:
        text : Text to write to file.
        write_path : Path to write to file. If None, then a NamedTemporaryFile will
            be created instead.
        file_ext : Specify the file extension if write_path is None.

    Returns:
        Path to file that was just written.
    """
    if write_path is None:
        if file_ext is None:
            raise ValueError("file_ext must be specified if write_path is None.")
        with NamedTemporaryFile(
            mode="w+", encoding="utf-8", suffix=file_ext, delete=False
        ) as temp_file:
            write_path = temp_file.name
            logger.info(f"Writing to temporary file: {write_path}")
            if isinstance(text, mda.core.groups.AtomGroup):
                text.write(temp_file.name)
            else:
                temp_file.write(text)
    else:
        logger.info(f"Writing to file: {write_path}")
        with open(write_path, "w", encoding="utf-8") as file:
            file.write(text)

    return write_path


@contextmanager
def decompress(input_file: str):
    """Decompress a gzipped file to a temporary file.

    Yields the path to the decompressed temporary file and cleans it up on
    exit. If the input file is not gzipped, yields the original file path

    Args:
        input_file: Path to any file, possibly gzip-compressed.

    Yields:
        Path to the decompressed file.
    """
    ext = get_file_ext(input_file)
    if ext == ".gz":
        inner_ext = get_file_ext(os.path.splitext(input_file)[0])
    elif ext.endswith("gz"):
        inner_ext = ext[: -len("gz")]
    else:
        yield input_file
        return

    tmp = NamedTemporaryFile(suffix=inner_ext, delete=False)
    try:
        with gzip.open(input_file, "rb") as f_in:
            tmp.write(f_in.read())
        tmp.close()
        logger.info(f"Decompressed {input_file} to temporary {inner_ext} file")
        yield tmp.name
    finally:
        if os.path.exists(tmp.name):
            os.unlink(tmp.name)


def _find_project_root() -> str:
    """Walk up to the directory containing pixi.toml."""
    current = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    for _ in range(10):
        if os.path.exists(os.path.join(current, "pixi.toml")):
            return current
        parent = os.path.dirname(current)
        if parent == current:
            break
        current = parent
    raise RuntimeError(
        "mgltools not found. Run from within the lignova project (pixi.toml must be present)."
    )


def mgltools_env_exists() -> bool:
    """Whether the mgltools pixi environment is installed. Never installs."""
    try:
        root = _find_project_root()
    except RuntimeError:
        return False
    return os.path.isdir(os.path.join(root, ".pixi", "envs", "mgltools"))


def _get_mgltools_prefix(auto_install: bool = False) -> str:
    """Find the isolated mgltools pixi environment, optionally installing it.

    Args:
        auto_install: Install the environment if missing. When False (default),
            a missing environment raises instead of triggering a long install.
    """
    current = _find_project_root()
    mgltools_env = os.path.join(current, ".pixi", "envs", "mgltools")
    if os.path.isdir(mgltools_env):
        return mgltools_env

    if not auto_install:
        raise RuntimeError(
            f"mgltools environment not found at '{mgltools_env}'. Install it with "
            "`pixi install -e mgltools`, or pass auto_install=True."
        )

    logger.info("mgltools environment not found. Installing it.")
    result = subprocess.run(
        ["pixi", "install", "-e", "mgltools"],
        cwd=current,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(f"Failed to install mgltools environment:\n{result.stderr}")
    logger.info("mgltools environment installed successfully.")
    return mgltools_env


def run_mgltools_command(
    cmd: list[str], auto_install: bool = False, **subprocess_kwargs
) -> subprocess.CompletedProcess:
    """Run a command inside the isolated mgltools pixi environment.

    Args:
        cmd: Command and arguments to run
        auto_install: Install the environment if missing. Defaults is False
        **subprocess_kwargs: Additional kwargs passed to subprocess.run


    Returns:
        The CompletedProcess result.
    """
    mgltools_prefix = _get_mgltools_prefix(auto_install=auto_install)
    mgltools_bin = os.path.join(mgltools_prefix, "bin")

    # Resolve command names to full paths within the mgltools env
    resolved_cmd = []
    for arg in cmd:
        candidate = os.path.join(mgltools_bin, arg)
        if os.path.exists(candidate):
            resolved_cmd.append(candidate)
        else:
            resolved_cmd.append(arg)

    # Build a clean env scoped to mgltools only
    env = os.environ.copy()
    env["PATH"] = f"{mgltools_bin}{os.pathsep}{env.get('PATH', '')}"

    for var in ("LD_LIBRARY_PATH", "DYLD_LIBRARY_PATH", "DYLD_FALLBACK_LIBRARY_PATH"):
        env.pop(var, None)

    mgltools_lib = os.path.join(mgltools_prefix, "lib")
    if os.path.isdir(mgltools_lib):
        env["LD_LIBRARY_PATH"] = mgltools_lib
        env["DYLD_FALLBACK_LIBRARY_PATH"] = mgltools_lib

    return subprocess.run(resolved_cmd, env=env, **subprocess_kwargs)
