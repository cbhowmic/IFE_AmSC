# -----------------------------------------------------------------------------
# RHINO Campaign Workflow Configuration
#
# Central configuration file for RHINO campaign archive and index creation.
#
# This file defines:
#   - RHINO data locations
#   - campaign storage locations
#   - archive naming conventions
#   - archive sizing parameters
#   - input directories containing ADIOS/openPMD outputs
#
# These parameters are separated from the workflow logic so that:
#   - campaign scripts remain reusable
#   - dataset locations can change without modifying workflow code
#   - multiple RHINO campaign configurations can be supported easily
#
# Used by:
#   create_campaign_archives.sh
#   create_campaign_index.sh
#   validate_campaign.sh
#
# RHINO / IFE_AmSC
# -----------------------------------------------------------------------------

RHINO_DATA_ROOT="/global/cfs/cdirs/m3239/2026_FES-AmSC/data/rhino"

CAMPAIGN_STORE="/global/homes/b/bhowmic/campaign-store/IFE"     #needs to be updated as per user's campaign store location

CAMPAIGN_NAMESPACE="IFE"

ARCHIVE_PREFIX="rhino"

DATASETS_PER_ARCHIVE=80

CAMPAIGN_INDEX="IFE/rhino.acx"

INPUT_DIRS=(
  "surrogate_bp_output/2026-04-29"
  "surrogate_bp_output/2026-04-30"
  "surrogate_bp_output/2026-05-01"
)