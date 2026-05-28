#!/usr/bin/env bash

# -----------------------------------------------------------------------------
# RHINO Campaign Archive Creation
#
# Creates HPC Campaign archives from RHINO ADIOS/openPMD outputs.
#
# Workflow:
#   1. Discover RHINO .bp5 outputs
#   2. Partition datasets into campaign archives
#   3. Create campaign archives
#   4. Add datasets with RHINO run identifiers
#
# Archives are split to avoid oversized campaign archives and to support
# scalable indexing/querying workflows.
#
# Archive naming:
#   rhino1.aca
#   rhino2.aca
#   ...
#
# Dataset naming:
#   Derived from RHINO run identifiers embedded in filenames.
#
# Requires:
#   - hpc_campaign
#   - RHINO ADIOS/openPMD outputs
#
# RHINO / IFE_AmSC
# -----------------------------------------------------------------------------

set -euo pipefail

# -----------------------------------------------------------------------------
# Load workflow configuration
# -----------------------------------------------------------------------------

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$SCRIPT_DIR/config_campaign.sh"

# -----------------------------------------------------------------------------
# Move to RHINO data root
# -----------------------------------------------------------------------------

cd "$RHINO_DATA_ROOT"

echo "--------------------------------------------------"
echo "RHINO Campaign Archive Creation"
echo "--------------------------------------------------"
echo "Data root          : $RHINO_DATA_ROOT"
echo "Campaign store     : $CAMPAIGN_STORE"
echo "Archive prefix     : $ARCHIVE_PREFIX"
echo "Datasets/archive   : $DATASETS_PER_ARCHIVE"
echo "--------------------------------------------------"

# -----------------------------------------------------------------------------
# Discover input datasets
# -----------------------------------------------------------------------------

shopt -s nullglob

files=()

for d in "${INPUT_DIRS[@]}"; do
    echo "Scanning directory: $d"

    for f in "$d"/*.bp5; do
        files+=("$f")
    done
done

shopt -u nullglob

IFS=$'\n' files=($(sort <<<"${files[*]}"))
unset IFS

# -----------------------------------------------------------------------------
# Validate dataset discovery
# -----------------------------------------------------------------------------

if (( ${#files[@]} == 0 )); then
    echo "ERROR: No .bp5 files discovered."
    exit 1
fi

echo "Discovered ${#files[@]} dataset(s)."

estimated_archives=$(( (${#files[@]} + DATASETS_PER_ARCHIVE - 1) / DATASETS_PER_ARCHIVE ))

echo "Estimated archive count : $estimated_archives"

# -----------------------------------------------------------------------------
# Create campaign archives
# -----------------------------------------------------------------------------

archive_idx=1
count=0

for f in "${files[@]}"; do

    name=$(basename "$f" .bp5)

    # Extract RHINO run identifier from filename
    run_id=$(echo "$name" | grep -oP '\d{2}-\d{2}-\d{2}.*' || true)

    # Fallback if regex extraction fails
    if [[ -z "$run_id" ]]; then
        run_id="$name"
    fi

    archive="${ARCHIVE_PREFIX}${archive_idx}.aca"

    echo
    echo "Adding dataset:"
    echo "  File    : $f"
    echo "  Run ID  : $run_id"
    echo "  Archive : $archive"

    # Create new archive when count == 0
    if (( count == 0 )); then

        hpc_campaign manager \
            --campaign_store "$CAMPAIGN_STORE" \
            "$archive" \
            --truncate \
            data "$f" \
            --name "$run_id"

    else

        hpc_campaign manager \
            "${CAMPAIGN_NAMESPACE}/$archive" \
            data "$f" \
            --name "$run_id"

    fi

    ((count++))

    # Start next archive group
    if (( count == DATASETS_PER_ARCHIVE )); then
        ((archive_idx++))
        count=0
    fi

done

# -----------------------------------------------------------------------------
# Summary
# -----------------------------------------------------------------------------

echo
echo "--------------------------------------------------"
echo "Campaign archive creation complete."
echo "Created archive groups with up to"
echo "$DATASETS_PER_ARCHIVE datasets each."
echo "--------------------------------------------------"