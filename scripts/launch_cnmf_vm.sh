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
#   bash scripts/launch_cnmf_vm.sh                  # full run (all cells)
#   bash scripts/launch_cnmf_vm.sh --trial-run      # quick test (~10 min)
#   bash scripts/launch_cnmf_vm.sh --neurons-only   # neurons only (~96k cells)
#
# Estimated cost: ~$1.50-3 on Spot (n2-highmem-32, 2-4 hrs)

set -euo pipefail

# --- Parse flags ---
TRIAL_RUN="false"
NEURONS_ONLY="false"
ZONE="us-central1-a"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --trial-run) TRIAL_RUN="true"; shift ;;
        --neurons-only) NEURONS_ONLY="true"; shift ;;
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
    MACHINE_TYPE="e2-highmem-8"  # 8 vCPU, 64 GB — ABC metadata loading needs headroom
    INSTANCE_NAME="cnmf-trial"
fi
if [ "$NEURONS_ONLY" = "true" ]; then
    INSTANCE_NAME="cnmf-neurons"
fi

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
STARTUP_SCRIPT="${SCRIPT_DIR}/cnmf_startup.sh"

echo "=== cNMF VM Launcher ==="
echo "Project: $PROJECT"
echo "Zone: $ZONE"
echo "Machine: $MACHINE_TYPE"
[ "$TRIAL_RUN" = "true" ] && echo "Mode: TRIAL RUN"
[ "$NEURONS_ONLY" = "true" ] && echo "Mode: NEURONS ONLY"
echo ""

# --- 1. Create GCS bucket if needed ---
if ! gcloud storage ls "$GCS_BUCKET" &>/dev/null; then
    echo "Creating GCS bucket ${GCS_BUCKET}..."
    gcloud storage buckets create "$GCS_BUCKET" --location=us-central1
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
    --metadata=TRIAL_RUN="$TRIAL_RUN",NEURONS_ONLY="$NEURONS_ONLY" \
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
echo "  gcloud storage cp -r ${GCS_BUCKET}/cnmf_output/ data/processed/mouse_abc/cnmf/"
echo ""
echo "  # Delete VM manually (if needed):"
echo "  gcloud compute instances delete $INSTANCE_NAME --zone=$ZONE --quiet"
