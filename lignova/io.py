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


def _get_mgltools_prefix() -> str:
    """Find or install the isolated mgltools pixi environment."""
    current = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    for _ in range(10):
        if os.path.exists(os.path.join(current, "pixi.toml")):
            mgltools_env = os.path.join(current, ".pixi", "envs", "mgltools")
            if not os.path.isdir(mgltools_env):
                logger.info("mgltools environment not found. Installing it.")
                result = subprocess.run(
                    ["pixi", "install", "-e", "mgltools"],
                    cwd=current,
                    capture_output=True,
                    text=True,
                )
                if result.returncode != 0:
                    raise RuntimeError(
                        f"Failed to install mgltools environment:\n{result.stderr}"
                    )
                logger.info("mgltools environment installed successfully.")
            return mgltools_env
        parent = os.path.dirname(current)
        if parent == current:
            break
        current = parent
    raise RuntimeError(
        "mgltools not found. Run from within the lignova project (pixi.toml must be present)."
    )


def run_mgltools_command(
    cmd: list[str], **subprocess_kwargs
) -> subprocess.CompletedProcess:
    """Run a command inside the isolated mgltools pixi environment.

    Automatically installs the environment if it doesn't exist.

    Args:
        cmd: Command and arguments to run
        **subprocess_kwargs: Additional kwargs passed to subprocess.run

    Returns:
        The CompletedProcess result.
    """
    mgltools_prefix = _get_mgltools_prefix()
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
