#!/usr/bin/env bash
set -euo pipefail

if ! command -v pixi >/dev/null 2>&1; then
    echo "ERROR: pixi not found in PATH. Install pixi first:" >&2
    echo "       https://pixi.prefix.dev/latest/installation" >&2
    exit 1
fi

if command -v nvidia-smi >/dev/null 2>&1 && nvidia-smi -L >/dev/null 2>&1; then
    echo "GPU detected - installing normally"
    pixi install -e gnina
    pixi run -e gnina bash ./run_scripts/install_gnina.sh
else
    echo "No GPU detected,using CONDA_OVERRIDE_CUDA=12.8 for CPU-only mode"
    CONDA_OVERRIDE_CUDA=12.8 pixi install -e gnina
    CONDA_OVERRIDE_CUDA=12.8 pixi run -e gnina bash ./run_scripts/install_gnina.sh
fi