#!/usr/bin/env bash
#SBATCH --account=m3239
#SBATCH --constraint=cpu
#SBATCH --qos=regular
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=04:00:00
#SBATCH --job-name=laser-bp5
#SBATCH --output=laser-bp5-%j.out
#SBATCH --error=laser-bp5-%j.err

# Convert each HDF5 acquisition in a separate Python process so its memory is
# returned to the operating system before the next acquisition is loaded.

set -uo pipefail

SUBMIT_DIR=${SLURM_SUBMIT_DIR:-$PWD}
LASER_DATA_ROOT=${LASER_DATA_ROOT:-/global/cfs/cdirs/m3239/2026_FES-AmSC/data/laser}
LASER_BP_OUT_DIR=${LASER_BP_OUT_DIR:-$LASER_DATA_ROOT/bp_output}
LASER_CONDA_ENV=${LASER_CONDA_ENV:-IFE_AmSC}

# sbatch executes a copy of this file from /var/spool/slurmd, so BASH_SOURCE
# cannot be used to locate repository files. Resolve them from the directory
# where sbatch was invoked, while supporting both repository layouts in use.
if [[ -n "${LASER_SHIM_DIR:-}" ]]; then
    : # Respect an explicit override.
elif [[ -f "$SUBMIT_DIR/Laser/shim/laserWrite_multi.py" ]]; then
    LASER_SHIM_DIR="$SUBMIT_DIR/Laser/shim"
elif [[ -f "$SUBMIT_DIR/LASER/laserWrite_multi.py" ]]; then
    LASER_SHIM_DIR="$SUBMIT_DIR/LASER"
elif [[ -f "$SUBMIT_DIR/laserWrite_multi.py" ]]; then
    LASER_SHIM_DIR="$SUBMIT_DIR"
else
    LASER_SHIM_DIR="$SUBMIT_DIR/Laser/shim"
fi

if [[ ! -d "$LASER_DATA_ROOT" ]]; then
    echo "ERROR: LASER_DATA_ROOT does not exist: $LASER_DATA_ROOT" >&2
    exit 1
fi
if [[ ! -f "$LASER_SHIM_DIR/laserWrite_multi.py" ]]; then
    echo "ERROR: laserWrite_multi.py not found under LASER_SHIM_DIR: $LASER_SHIM_DIR" >&2
    echo "Set LASER_SHIM_DIR explicitly with sbatch --export if the shim is elsewhere." >&2
    exit 1
fi

module load python
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$LASER_CONDA_ENV"

mkdir -p "$LASER_BP_OUT_DIR"

mapfile -d '' h5_files < <(
    find "$LASER_DATA_ROOT" \
        -type f \
        -name '*.h5' \
        ! -path "$LASER_BP_OUT_DIR/*" \
        ! -path '*/documentation/*' \
        -print0
)

if (( ${#h5_files[@]} == 0 )); then
    echo "ERROR: No HDF5 inputs found under $LASER_DATA_ROOT" >&2
    exit 1
fi

echo "LASER source root : $LASER_DATA_ROOT"
echo "BP5 output        : $LASER_BP_OUT_DIR"
echo "Conda environment : $LASER_CONDA_ENV"
echo "Discovered        : ${#h5_files[@]} HDF5 file(s)"

converted_or_skipped=0
failed=0

for h5_path in "${h5_files[@]}"; do
    echo
    echo "Processing: $h5_path"
    if python "$LASER_SHIM_DIR/laserWrite_multi.py" \
        "$h5_path" \
        --out-dir "$LASER_BP_OUT_DIR" \
        --skip-existing \
        --stop-on-error; then
        ((converted_or_skipped += 1))
    else
        ((failed += 1))
        echo "ERROR: Conversion failed: $h5_path" >&2
    fi
done

echo
echo "Processed successfully or already present: $converted_or_skipped"
echo "Failed: $failed"

if (( failed > 0 )); then
    exit 1
fi
