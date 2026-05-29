#!/usr/bin/env bash
set -euo pipefail


require_tool() {
    if ! command -v "$1" >/dev/null 2>&1; then
        echo "ERROR: required tool '$1' not found in PATH." >&2
        echo "       Activate the gnina pixi environment first:" >&2
        echo "         pixi shell -e gnina   (or run via 'pixi run -e gnina ...')" >&2
        exit 1
    fi
}
for tool in cmake make git gcc g++ nvcc; do
    require_tool "$tool"
done

if [ -z "${CONDA_PREFIX:-}" ]; then
    echo "ERROR: CONDA_PREFIX is not set. Activate the gnina pixi env first." >&2
    exit 1
fi
echo " Using CONDA_PREFIX = ${CONDA_PREFIX}"


export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH:-}"
export LDFLAGS="-L${CONDA_PREFIX}/lib -Wl,-rpath,${CONDA_PREFIX}/lib"
export CPATH="${CONDA_PREFIX}/include:${CPATH:-}"


SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
GNINA_SRC="${GNINA_SRC:-${REPO_ROOT}/build/gnina}"
GNINA_BUILD="${GNINA_SRC}/build"
mkdir -p "$(dirname "$GNINA_SRC")"

if [ ! -d "${GNINA_SRC}/.git" ]; then
    echo " Cloning gnina into ${GNINA_SRC}"
    git clone https://github.com/gnina/gnina.git "${GNINA_SRC}"
else
    echo " Source already present at ${GNINA_SRC} (skipping clone)"
fi


CMAKE_FILE="${GNINA_SRC}/CMakeLists.txt"
PATCH_MARKER="# --- patched by install_gnina.sh ---"

if ! grep -qF "$PATCH_MARKER" "$CMAKE_FILE"; then
    echo " Patching ${CMAKE_FILE}"

    sed -i '/ExternalProject_Add(libmolgrid/,/)/{
      s|CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${EXTERNAL_INSTALL_LOCATION}|CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${EXTERNAL_INSTALL_LOCATION} -DCMAKE_PREFIX_PATH='"${CONDA_PREFIX}"' -DOPENBABEL3_INCLUDE_DIR='"${CONDA_PREFIX}"'/include/openbabel3 -DOPENBABEL3_LIBRARIES='"${CONDA_PREFIX}"'/lib/libopenbabel.so -DCMAKE_CXX_STANDARD=17|
    }' "$CMAKE_FILE"

    sed -i 's|GIT_REPOSITORY https://github.com/gnina/libmolgrid|GIT_REPOSITORY https://github.com/gnina/libmolgrid\n     PATCH_COMMAND sed -i "s/-Werror//g" <SOURCE_DIR>/CMakeLists.txt|' "$CMAKE_FILE"

    echo "$PATCH_MARKER" >> "$CMAKE_FILE"
else
    echo " ${CMAKE_FILE} already patched (skipping)"
fi


mkdir -p "$GNINA_BUILD"
cd "$GNINA_BUILD"

echo " Configuring CMake build in ${GNINA_BUILD}"

CUDA_INC_DIR=""
while IFS= read -r hdr; do
    dir="$(dirname "$hdr")"
    if [ -f "${dir}/cuda_runtime.h" ]; then
        CUDA_INC_DIR="$dir"
        break
    fi
done < <(find "${CONDA_PREFIX}" -name cuda.h 2>/dev/null)

if [ -z "$CUDA_INC_DIR" ]; then
    echo "ERROR: could not locate the CUDA toolkit include dir under ${CONDA_PREFIX}" >&2
    exit 1
fi
CUDA_ROOT="$(dirname "$CUDA_INC_DIR")"
echo " CUDA headers at:  ${CUDA_INC_DIR}"
echo " CUDA toolkit root: ${CUDA_ROOT}"

export CPATH="${CUDA_INC_DIR}:${CPATH}"

export LIBRARY_PATH="${CUDA_ROOT}/lib:${LIBRARY_PATH:-}"
export LD_LIBRARY_PATH="${CUDA_ROOT}/lib:${LD_LIBRARY_PATH}"
cmake .. \
    -DCMAKE_INSTALL_PREFIX="${CONDA_PREFIX}" \
    -DCMAKE_PREFIX_PATH="${CONDA_PREFIX}" \
    -DCUDA_TOOLKIT_ROOT_DIR="${CUDA_ROOT}" \
    -DCUDA_INCLUDE_DIRS="${CUDA_INC_DIR}" \
    -DCUDAToolkit_ROOT="${CUDA_ROOT}" \
    -DCUDAToolkit_INCLUDE_DIR="${CUDA_INC_DIR}" \
    -DCMAKE_CUDA_ARCHITECTURES=all \
    -DCMAKE_CXX_STANDARD=17 \
    -DCMAKE_BUILD_TYPE=Release
NPROC="$(nproc)"
echo " Building with -j${NPROC}"
make -j"$NPROC"

echo " Installing into ${CONDA_PREFIX}"
make install

# Create libcuda.so.1 symlinks so the binary can load on no-GPU hosts.
for stub_dir in "${CONDA_PREFIX}/targets/x86_64-linux/lib/stubs" "${CONDA_PREFIX}/lib/stubs"; do
    if [ -f "${stub_dir}/libcuda.so" ] && [ ! -e "${stub_dir}/libcuda.so.1" ]; then
        echo " Creating libcuda.so.1 symlink in ${stub_dir}"
        ln -s libcuda.so "${stub_dir}/libcuda.so.1"
    fi
done

BIN="${CONDA_PREFIX}/bin/gnina"
if [ ! -x "$BIN" ]; then
    echo "ERROR: gnina binary was not produced at ${BIN}" >&2
    exit 1
fi
echo " Verifying installation"
HELP_OUT="$("$BIN" --help 2>&1 || true)"
VERSION_OUT="$("$BIN" --version 2>&1 || true)"

if echo "$HELP_OUT" | grep -q "^Usage:" || "$BIN" --version >/dev/null 2>&1; then
    echo " Done. Binary at: ${BIN}"
elif echo "$HELP_OUT" | grep -q "libcuda.so.1"; then
    echo " Built successfully, but libcuda.so.1 missing (no-GPU host)."
    echo " Binary at: ${BIN}"
else
    echo "WARNING: gnina --help failed unexpectedly:" >&2
    echo "$HELP_OUT" | head -10 >&2
    exit 1
fi