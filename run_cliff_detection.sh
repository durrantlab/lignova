#!/bin/bash
#SBATCH --job-name=cliff_2048
#SBATCH -A bio260240
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --time=00:30:00
#SBATCH --output=../logs/cliff_sweep/%x_%A_%a.out
#SBATCH --error=../logs/cliff_sweep/%x_%A_%a.err

set -uo pipefail
module purge
mkdir -p ../logs/cliff_sweep

MODE="${1:-all}"          # all | range | merge

ENRICHED="../lignova_parquets/measurements_merged.parquet"
OUTDIR="../cliff_sweep/fp_2048"

FP_SIZES="2048"
RADIUS=2
FLOOR=0.55
CUTOFF=0.55
MIN_DELTA=2.0
CHUNK=250                 # genes per array task (range mode only)

export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export MKL_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK}"
export NUMEXPR_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

mkdir -p "${OUTDIR}"

echo "Mode:        ${MODE}"
echo "Job ID:      ${SLURM_JOB_ID}"
echo "Array task:  ${SLURM_ARRAY_TASK_ID:-<none>}"
echo "Node:        $(hostname)"
echo "Start time:  $(date)"
echo "Out dir:     ${OUTDIR}"
echo "FP sizes:    ${FP_SIZES}"

if [[ ! -f "${ENRICHED}" ]]; then
    echo "ERROR: enriched parquet not found: ${ENRICHED}" >&2
    exit 1
fi

#Ensure that 'all' and 'merge' modes must be single jobs.
if [[ "${MODE}" != "range" && -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "ERROR: mode '${MODE}' must be a single job, not an array. " \
         "Submit without --array, or use 'range' for the array." >&2
    exit 2
fi

# Args for run_scripts.cliff_detection.py
COMMON=(
--parquet "${ENRICHED}"
--out     "${OUTDIR}"
--fp-sizes ${FP_SIZES}
--radius  "${RADIUS}"
--floor   "${FLOOR}"
--cutoff  "${CUTOFF}"
--min-delta "${MIN_DELTA}"
)

case "${MODE}" in
all)
      ARGS=( "${COMMON[@]}" --select all )
      ;;
range)
      : "${SLURM_ARRAY_TASK_ID:?range mode must be submitted as an array (--array=...)}"
      START=$(( SLURM_ARRAY_TASK_ID * CHUNK + 1 ))
      END=$((   START + CHUNK - 1 ))
      echo "Range task ${SLURM_ARRAY_TASK_ID}: genes ${START}-${END}"
      ARGS=( "${COMMON[@]}" --select range --range "${START}" "${END}" --no-merge )
      ;;
merge)
      ARGS=( --out "${OUTDIR}" --merge-only )
      ;;
  *)
      echo "ERROR: unknown mode '${MODE}' (use all | range | merge)" >&2
      exit 2
      ;;
esac


set +e
pixi run python3 -m run_scripts.cliff_detection "${ARGS[@]}"
rc=$?
set -e

if [[ "${rc}" -eq 0 ]]; then
    echo "OK"
elif [[ "${rc}" -eq 3 && "${MODE}" == "range" ]]; then
    echo "No targets in range ${START}-${END}; nothing to do (clean skip)."
    rc=0
else
    echo "Pipeline failed (exit ${rc})" >&2
fi

echo "Finished at: $(date)"
exit "${rc}"