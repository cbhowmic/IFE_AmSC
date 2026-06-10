#!/usr/bin/env bash

# -----------------------------------------------------------------------------
# RHINO Campaign Index Creation
#
# Creates a campaign index from RHINO campaign archives.
#
# The campaign index is stored as an .acx file and can be inspected using
# hpc_campaign index commands or queried directly with sqlite3.
#
# RHINO / IFE_AmSC
# -----------------------------------------------------------------------------

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$SCRIPT_DIR/config_campaign.sh"

cd "$CAMPAIGN_STORE"

echo "--------------------------------------------------"
echo "RHINO Campaign Index Creation"
echo "--------------------------------------------------"
echo "Campaign store : $CAMPAIGN_STORE"
echo "Archive prefix : $ARCHIVE_PREFIX"
echo "Index          : $CAMPAIGN_INDEX"
echo "--------------------------------------------------"

shopt -s nullglob
# archive_files=("${CAMPAIGN_NAMESPACE}/${ARCHIVE_PREFIX}"*.aca)
archive_files=("${ARCHIVE_PREFIX}"*.aca)
shopt -u nullglob

if (( ${#archive_files[@]} == 0 )); then
    echo "ERROR: No archives found to index."
    echo "Expected archives like: ${CAMPAIGN_NAMESPACE}/${ARCHIVE_PREFIX}*.aca"
    exit 1
fi

echo "Found ${#archive_files[@]} archive(s):"
for archive in "${archive_files[@]}"; do
    echo "  $archive"
done

echo
echo "Creating index: $CAMPAIGN_INDEX"

rm -f "$CAMPAIGN_INDEX"

hpc_campaign index "$CAMPAIGN_INDEX" add "${archive_files[@]}"

echo
echo "Inspecting index:"
hpc_campaign index "$CAMPAIGN_INDEX" ls

echo
echo "--------------------------------------------------"
echo "Campaign index creation complete."
echo "Index: $CAMPAIGN_INDEX"
echo "--------------------------------------------------"