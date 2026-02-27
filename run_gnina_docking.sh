#!/bin/bash
#SBATCH --job-name=gnina_dock
#SBATCH --partition=htc
#SBATCH --cluster=htc
#SBATCH --account=jdurrant
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=32G
#SBATCH --time=1-22:00:00
#SBATCH --array=51,96,104,163,190,201-300
#SBATCH --output=../logs/dock/%x_%A_%a.out
#SBATCH --error=../logs/dock/%x_%A_%a.err

set -euo pipefail
module purge
mkdir -p ../logs/dock

# Configuration
PREPARED_DIR="../prepared_prot"

# Parallelism settings
PER_RUN_CPUS=12
MAX_CONCURRENT=$(( SLURM_CPUS_PER_TASK / PER_RUN_CPUS ))
(( MAX_CONCURRENT < 1 )) && MAX_CONCURRENT=1

# Avoid thread oversubscription
export OMP_NUM_THREADS="${PER_RUN_CPUS}"
export MKL_NUM_THREADS="${PER_RUN_CPUS}"
export OPENBLAS_NUM_THREADS="${PER_RUN_CPUS}"
export NUMEXPR_NUM_THREADS="${PER_RUN_CPUS}"

# Load GNINA module
module purge
module load gnina

TASK_ID="${SLURM_ARRAY_TASK_ID}"

echo "Job ID: ${SLURM_JOB_ID}"
echo "Array Task ID: ${TASK_ID}"
echo "Node: $(hostname)"
echo "Start time: $(date)"
echo "PWD: $(pwd)"

# Get the Nth directory (0-indexed)
mapfile -t PDB_DIRS < <(find "${PREPARED_DIR}" -mindepth 1 -maxdepth 1 -type d | sort)

if (( TASK_ID >= ${#PDB_DIRS[@]} )); then
    echo "Task ID ${TASK_ID} exceeds number of directories (${#PDB_DIRS[@]}), exiting."
    exit 0
fi

pdb_dir="${PDB_DIRS[${TASK_ID}]}"
pdb_name=$(basename "${pdb_dir}")

echo "Processing directory: ${pdb_dir}"
echo "PDB name: ${pdb_name}"

# Validate required files exist
receptor=$(find "${pdb_dir}" -maxdepth 1 -name "*_protonated.pdbqt" -type f | head -n 1)
if [[ -n "${receptor}" ]]; then
    echo "Receptor (pdbqt): ${receptor}"
else
    receptor=$(find "${pdb_dir}" -maxdepth 1 -name "*_protonated.pqr" -type f | head -n 1)
    if [[ -n "${receptor}" ]]; then
        echo "Receptor (pqr fallback): ${receptor}"
    else
        echo "ERROR: No receptor (*_protonated.pdbqt or *_protonated.pqr) found for ${pdb_name}, exiting."
        exit 1
    fi
fi


box_ligand=$(find "${pdb_dir}" -maxdepth 1 -name "*_ligand.pdb" -type f | head -n 1)
if [[ -z "${box_ligand}" ]]; then
    echo "ERROR: No autobox ligand (*_ligand.pdb) found for ${pdb_name}, exiting."
    exit 1
fi

mapfile -t LIGAND_FILES < <(find "${pdb_dir}" -path "*/gypsum_out_*/*.sdf" -type f)
if (( ${#LIGAND_FILES[@]} == 0 )); then
    echo "ERROR: No ligands (gypsum_out_*/*.sdf) found for ${pdb_name}, exiting."
    exit 1
fi

echo "Receptor: ${receptor}"
echo "Box ligand: ${box_ligand}"
echo "Ligands found: ${#LIGAND_FILES[@]}"
echo "Max concurrent jobs: ${MAX_CONCURRENT}"
echo "CPUs per job: ${PER_RUN_CPUS}"

# Track failed dockings
failed_file="${pdb_dir}/failed_docking_${SLURM_JOB_ID}_${TASK_ID}.txt"

# Function to dock one ligand
dock_one() {
    local ligand="$1"
    local ligand_name
    ligand_name=$(basename "${ligand}" .sdf)
    
    local out_dir
    out_dir=$(dirname "${ligand}")
    local out_file="${out_dir}/${ligand_name}_docked.sdf.gz"
    
    # Skip if already docked
    if [[ -f "${out_file}" ]]; then
        echo "Skipping ${ligand_name}, already docked: ${out_file}"
        return 0
    fi
    
    echo "Docking ${ligand_name}..."
    
    # Create unique config file for this ligand (avoids race conditions)
    local config_file="${out_dir}/${ligand_name}_gnina_config.yaml"
    local cmd_file="${out_dir}/${ligand_name}_cmd.txt"
    
    if ! pixi run -e dev python3 -m run_scripts.dock \
        -r "${receptor}" \
        -l "${ligand}" \
        -b "${box_ligand}" \
        -c "${config_file}" \
        -o "${cmd_file}"; then
        echo "Warning: Failed to generate command for ${ligand_name}" >&2
        echo "${ligand}" >> "${failed_file}"
        return 0
    fi
    
    # run GNINA command
    local gnina_cmd
    gnina_cmd=$(head -n 1 "${cmd_file}")
    
    if [[ -z "${gnina_cmd}" ]]; then
        echo "Warning: Empty command for ${ligand_name}" >&2
        echo "${ligand}" >> "${failed_file}"
        return 0
    fi
    
    echo "Running: ${gnina_cmd}"
    
    if ! eval "${gnina_cmd}"; then
        echo "Warning: GNINA failed for ${ligand_name}" >&2
        echo "${ligand}" >> "${failed_file}"
        return 0
    fi
    
    echo "Completed docking for ${ligand_name}"
}

# Run docking jobs in parallel
pids=()
for ligand in "${LIGAND_FILES[@]}"; do
    dock_one "${ligand}" &
    pids+=($!)
    
    # Wait if we've reached max concurrent jobs
    while (( ${#pids[@]} >= MAX_CONCURRENT )); do
        if ! wait -n; then
            echo "Warning: One of the docking processes failed." >&2
        fi
        # Remove finished PIDs
        alive=()
        for pid in "${pids[@]}"; do
            kill -0 "$pid" 2>/dev/null && alive+=("$pid")
        done
        pids=("${alive[@]}")
    done
done

# Wait for remaining jobs
wait

echo "Completed ${pdb_name} at $(date)"

if [[ -f "${failed_file}" ]]; then
    failed_count=$(wc -l < "${failed_file}")
    echo "Failed dockings: ${failed_count}"
    echo "See: ${failed_file}"
else
    echo "All dockings completed successfully!"
fi

crc-job-stats || true