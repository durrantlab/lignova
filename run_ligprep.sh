#!/bin/bash
#SBATCH --job-name=lig_prep
#SBATCH --partition=htc
#SBATCH --cluster=htc
#SBATCH --account=jdurrant
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=2-12:00:00
#SBATCH --array=111
#SBATCH --output=../logs/lig_prep/%x_%A_%a.out
#SBATCH --error=../logs/lig_prep/%x_%A_%a.err


set -euo pipefail
module purge
mkdir -p ../logs/lig_prep


# Avoid thread oversubscription inside Python/RDKit stacks
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

# paths 
PARQUET="../../lignova_parquets/protein_clustered_data.parquet"
OUTDIR="../prepared_prot"
CONFIG_YAML="${OUTDIR}/gypsum.yaml"
BATCH_SIZE=30
LIG_SOURCE="../../lignova_parquets/new_hdf5_w_smi.parquet"
SMI_NAME="lig.smi"
CHUNK_SIZE=250
PER_RUN_CPUS=8
TASK_ID="${SLURM_ARRAY_TASK_ID}"
START_IDX=$(( TASK_ID * BATCH_SIZE ))
END_IDX=$(( START_IDX + BATCH_SIZE ))

echo "Array task ${TASK_ID}: indices [${START_IDX}, ${END_IDX})"
echo "Job ID: ${SLURM_JOB_ID:-N/A}  Node: $(hostname)  Start: $(date)"
echo "PWD: $(pwd)"

# run data prep for this batch 
srun pixi run -e dev python3 -m run_scripts.data_prep \
  -m ligand \
  -p "${PARQUET}" \
  -o "${OUTDIR}" \
  -l "${LIG_SOURCE}" \
  -osmi "${SMI_NAME}" \
  -ch "${CHUNK_SIZE}" \
  --start-index "${START_IDX}" \
  --end-index "${END_IDX}"

TXT_LIST="${OUTDIR}/smiles_filelist_${START_IDX}_${END_IDX}.txt"
echo "Expecting list: ${TXT_LIST}"
[[ -f "${TXT_LIST}" ]] || { echo "ERROR: missing ${TXT_LIST}" >&2; exit 1; }

mapfile -t SMI_FILES < "${TXT_LIST}"
echo "Got ${#SMI_FILES[@]} chunked .smi files"

if (( ${#SMI_FILES[@]} == 0 )); then
  echo "No .smi files listed; exiting."
  exit 0
fi

# parallel gypsum runs inside the node 
MAX_CONCURRENT=$(( SLURM_CPUS_PER_TASK / PER_RUN_CPUS ))
(( MAX_CONCURRENT < 1 )) && MAX_CONCURRENT=1

echo "Running up to ${MAX_CONCURRENT} gypsum jobs concurrently; ${PER_RUN_CPUS} CPUs each"

run_one() {
  local SMI="$1"

  # Output directory: keep it next to the smi chunks
  local OUT_DIR
  OUT_DIR="$(dirname "$SMI")/gypsum_out_$(basename "${SMI%.smi}")"

  #rename the gypsum output file to include the smi file index
  local base produced target
  base="$(basename "${SMI}" .smi)"
  produced="${OUT_DIR}/gypsum_dl_success.sdf"
  target="${OUT_DIR}/${base}.sdf"
  #check if gypsum output file already exists
  if [[ -f "${target}" ]]; then
    echo "Skipping ${SMI}, target ${target} already exists."
    return 0
  fi


  mkdir -p "$OUT_DIR"
  echo "Processing ${SMI} into ${OUT_DIR} using gypsum"

  if ! GYPSUM_NUM_PROCS="${PER_RUN_CPUS}" \
  pixi run -e dev python3 -m run_scripts.conformer_generation \
    -s "${SMI}" \
    -c "${CONFIG_YAML}" \
    -o "${OUT_DIR}"; then
    echo "Warning: gypsum failed for ${SMI}" >&2
    echo "${SMI}" >> "${OUTDIR}/failed_smi_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}.txt"
    echo "Continuing to next file."
    return 0
    fi


  if [[ -f "${produced}" ]]; then
    mv -f "${produced}" "${target}"
    echo "Produced ${target} for ${SMI}"
  else
    echo "Warning: expected ${produced} not found for ${SMI}" >&2
  fi

}

pids=()
for SMI in "${SMI_FILES[@]}"; do
  run_one "${SMI}" &
  pids+=($!)

  while (( ${#pids[@]} >= MAX_CONCURRENT )); do
    if ! wait -n; then
      echo "Warning: one of the gypsum processes failed." >&2
    fi
    alive=()
    for pid in "${pids[@]}"; do kill -0 "$pid" 2>/dev/null && alive+=("$pid"); done
    pids=("${alive[@]}")
  done
done
wait

echo "Done array task ${TASK_ID} at $(date)"
crc-job-stats || true