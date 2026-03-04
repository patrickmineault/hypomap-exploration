#!/bin/bash
# Launch a GCP Spot VM to run cNMF on hypothalamus scRNA-seq.
#
# The VM downloads expression data directly from Allen ABC Atlas,
# runs cNMF, and uploads results to GCS.
#
# Prerequisites (one-time):
#   gcloud auth login
#   gcloud config set project <PROJECT_ID>
#   gcloud services enable compute.googleapis.com
#
# Usage:
#   bash scripts/launch_cnmf_vm.sh              # full run
#   bash scripts/launch_cnmf_vm.sh --trial-run  # quick test (~10 min)
#
# Estimated cost: ~$1.50-3 on Spot (n2-highmem-32, 2-4 hrs)

set -euo pipefail

# --- Parse flags ---
TRIAL_RUN="false"
ZONE="us-central1-a"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --trial-run) TRIAL_RUN="true"; shift ;;
        --zone) ZONE="$2"; shift 2 ;;
        *) echo "Unknown flag: $1"; exit 1 ;;
    esac
done

# --- Configuration ---
PROJECT=$(gcloud config get-value project 2>/dev/null)
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

# --- 1. Create GCS bucket if needed ---
if ! gsutil ls "$GCS_BUCKET" &>/dev/null; then
    echo "Creating GCS bucket ${GCS_BUCKET}..."
    gsutil mb -l us-central1 "$GCS_BUCKET"
fi

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
