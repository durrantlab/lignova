#!/bin/bash
#SBATCH --job-name=gnina_dock
#SBATCH --partition=htc
#SBATCH --cluster=htc
#SBATCH --account=jdurrant
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=2-22:00:00
#SBATCH --array=130,141,239,253,266-269,271,275,283,291,298,310,338-343,403-404%20
#SBATCH --output=../logs/dock/%x_%A_%a.out
#SBATCH --error=../logs/dock/%x_%A_%a.err


#to mitigate the 500 job array limit, Added an offset to the task ID 
ARRAY_OFFSET="${ARRAY_OFFSET:-4500}"
TASK_ID=$(( SLURM_ARRAY_TASK_ID + ARRAY_OFFSET ))

set -euo pipefail
module purge
mkdir -p ../logs/dock

# Configuration
PREPARED_DIR="../prepared_prot"

# Affinity filtering
REQUIRE_AFFINITY="${REQUIRE_AFFINITY:-1}"
AFFINITY_PATH="../../lignova_parquets/data_w_rotatable_bonds.parquet"
CLUSTERED_PATH="../../lignova_parquets/protein_clustered_data.parquet"

# Parallelism settings
PER_RUN_CPUS=8
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


#Rename log files to include the offset-adjusted task ID
ORIG_OUT="../logs/dock/${SLURM_JOB_NAME}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out"
ORIG_ERR="../logs/dock/${SLURM_JOB_NAME}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err"
NEW_OUT="../logs/dock/${SLURM_JOB_NAME}_${SLURM_ARRAY_JOB_ID}_${TASK_ID}.out"
NEW_ERR="../logs/dock/${SLURM_JOB_NAME}_${SLURM_ARRAY_JOB_ID}_${TASK_ID}.err"


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

ORIGINAL_LIGAND_COUNT=${#LIGAND_FILES[@]}
echo "Receptor: ${receptor}"
echo "Box ligand: ${box_ligand}"
echo "Ligands found: ${ORIGINAL_LIGAND_COUNT}"
echo "Max concurrent jobs: ${MAX_CONCURRENT}"
echo "CPUs per job: ${PER_RUN_CPUS}"

# Track failed dockings
failed_file="${pdb_dir}/failed_docking_${SLURM_JOB_ID}_${TASK_ID}.txt"

# filter SDFs to only compounds with usable affinity
if [[ "${REQUIRE_AFFINITY}" == "1" ]]; then
    if [[ -z "${AFFINITY_PATH:-}" || -z "${CLUSTERED_PATH:-}" ]]; then
        echo "ERROR: REQUIRE_AFFINITY=1 but AFFINITY_PATH/CLUSTERED_PATH not set"
        exit 1
    fi
    echo "Filtering SDFs by affinity availability."
    if ! pixi run -e dev python3 -m run_scripts.affinity_filter \
        -d "${pdb_dir}" \
        -a "${AFFINITY_PATH}" \
        -c "${CLUSTERED_PATH}"; then
        echo "ERROR: Filter step failed — aborting (REQUIRE_AFFINITY=1)" >&2
        exit 1
    fi
    NEW_LIGANDS=()
    for lig in "${LIGAND_FILES[@]}"; do
        filtered="${lig%.sdf}_filtered.sdf"
        if [[ -f "${filtered}" ]]; then
            NEW_LIGANDS+=("${filtered}")
        else
            echo "Skipping ${lig}: no filtered SDF (either excluded or filter skipped)"
        fi
    done
    LIGAND_FILES=("${NEW_LIGANDS[@]}")
    if (( ${#LIGAND_FILES[@]} == 0 )); then
        echo "No ligands remain after affinity filtering.Skipping."
        exit 0
    fi
    echo "After filtering: ${#LIGAND_FILES[@]} SDFs remain (from ${ORIGINAL_LIGAND_COUNT} originally)"
fi

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
    
    # Create unique config file for this ligand
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
        alive=()
        for pid in "${pids[@]}"; do
            kill -0 "$pid" 2>/dev/null && alive+=("$pid")
        done
        pids=("${alive[@]}")
    done
done

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

# Rename logs to include offset-adjusted task ID
if [[ "${ARRAY_OFFSET}" -ne 0 ]]; then
    mv "${ORIG_OUT}" "${NEW_OUT}" 2>/dev/null || true
    mv "${ORIG_ERR}" "${NEW_ERR}" 2>/dev/null || true
fi
