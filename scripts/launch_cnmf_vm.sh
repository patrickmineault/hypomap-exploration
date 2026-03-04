#!/bin/bash
# Launch a GCP Spot VM to run cNMF on hypothalamus scRNA-seq.
#
# Run upload_cnmf_data.sh first to upload the expression parquet to GCS.
#
# Prerequisites (one-time):
#   gcloud auth login
#   gcloud config set project <PROJECT_ID>
#   gcloud services enable compute.googleapis.com
#
# Usage:
#   bash scripts/upload_cnmf_data.sh   # upload data first
#   bash scripts/launch_cnmf_vm.sh     # full run
#   bash scripts/launch_cnmf_vm.sh --trial-run  # quick test (~5 min)
#
# Estimated cost: ~$1.50-3 on Spot (n2-highmem-32, 2-4 hrs)

set -euo pipefail

# --- Parse flags ---
TRIAL_RUN="false"
for arg in "$@"; do
    case "$arg" in
        --trial-run) TRIAL_RUN="true" ;;
    esac
done

# --- Configuration ---
PROJECT=$(gcloud config get-value project 2>/dev/null)
ZONE="us-central1-a"
INSTANCE_NAME="cnmf-runner"
MACHINE_TYPE="n2-highmem-32"  # 32 vCPU, 256 GB RAM
BOOT_DISK_SIZE="100GB"
GCS_BUCKET="gs://neuroai-abc"

if [ "$TRIAL_RUN" = "true" ]; then
    MACHINE_TYPE="e2-highmem-4"  # 4 vCPU, 32 GB — plenty for 500 cells
    INSTANCE_NAME="cnmf-trial"
fi

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
STARTUP_SCRIPT="${SCRIPT_DIR}/cnmf_startup.sh"

echo "=== cNMF VM Launcher ==="
echo "Project: $PROJECT"
echo "Zone: $ZONE"
echo "Machine: $MACHINE_TYPE"
[ "$TRIAL_RUN" = "true" ] && echo "Mode: TRIAL RUN"
echo ""

# --- 1. Verify data is in GCS ---
GCS_PARQUET="${GCS_BUCKET}/cnmf_input/scrna_expression_log2cpm.parquet"
if ! gsutil -q stat "$GCS_PARQUET" 2>/dev/null; then
    echo "ERROR: Expression parquet not found in GCS."
    echo "Run first: bash scripts/upload_cnmf_data.sh"
    exit 1
fi
echo "Expression parquet found in GCS."

# --- 2. Create Spot VM ---
echo ""
echo "Creating Spot VM: $INSTANCE_NAME..."
gcloud compute instances create "$INSTANCE_NAME" \
    --zone="$ZONE" \
    --machine-type="$MACHINE_TYPE" \
    --provisioning-model=SPOT \
    --instance-termination-action=DELETE \
    --boot-disk-size="$BOOT_DISK_SIZE" \
    --image-family=debian-12 \
    --image-project=debian-cloud \
    --scopes=storage-full,compute-rw \
    --metadata=TRIAL_RUN="$TRIAL_RUN" \
    --metadata-from-file=startup-script="$STARTUP_SCRIPT"

echo ""
echo "=== VM created successfully ==="
echo ""
echo "Monitor commands:"
echo "  # View startup log (serial console):"
echo "  gcloud compute instances get-serial-port-output $INSTANCE_NAME --zone=$ZONE"
echo ""
echo "  # SSH into VM and tail log:"
echo "  gcloud compute ssh $INSTANCE_NAME --zone=$ZONE -- tail -f /var/log/cnmf.log"
echo ""
echo "  # Check if VM is still running:"
echo "  gcloud compute instances describe $INSTANCE_NAME --zone=$ZONE --format='value(status)'"
echo ""
echo "  # Download results when done:"
echo "  gsutil -m cp -r ${GCS_BUCKET}/cnmf_output/ data/processed/mouse_abc/cnmf/"
echo ""
echo "  # Delete VM manually (if needed):"
echo "  gcloud compute instances delete $INSTANCE_NAME --zone=$ZONE --quiet"
