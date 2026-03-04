#!/bin/bash
# Upload expression h5ad to GCS for cNMF processing.
#
# Run this BEFORE launch_cnmf_vm.sh. It uploads the CPM expression h5ad
# that the VM will download and feed to cNMF.
#
# Prerequisites:
#   gcloud auth login
#   uv run snakemake extract_scrna_expression --cores 1  (generate the h5ad first)
#
# Usage:
#   bash scripts/upload_cnmf_data.sh

set -euo pipefail

GCS_BUCKET="gs://neuroai-abc"

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

# --- 1. Create GCS bucket if needed ---
if ! gsutil ls "$GCS_BUCKET" &>/dev/null; then
    echo "Creating GCS bucket ${GCS_BUCKET}..."
    gsutil mb -l us-central1 "$GCS_BUCKET"
fi

# --- 2. Upload expression h5ad ---
H5AD_FILE="${PROJECT_ROOT}/data/processed/mouse_abc/scrna_expression_cpm.h5ad"
GCS_H5AD="${GCS_BUCKET}/cnmf_input/scrna_expression_cpm.h5ad"

if [ ! -f "$H5AD_FILE" ]; then
    echo "ERROR: $H5AD_FILE not found."
    echo "Generate it first: uv run snakemake extract_scrna_expression --cores 1"
    exit 1
fi

if gsutil -q stat "$GCS_H5AD" 2>/dev/null; then
    LOCAL_SIZE=$(stat -f%z "$H5AD_FILE" 2>/dev/null || stat -c%s "$H5AD_FILE")
    REMOTE_SIZE=$(gsutil stat "$GCS_H5AD" 2>/dev/null | grep "Content-Length" | awk '{print $2}')
    if [ "$LOCAL_SIZE" = "$REMOTE_SIZE" ]; then
        echo "Expression h5ad already in GCS (same size), skipping upload."
        exit 0
    else
        echo "Expression h5ad in GCS has different size, re-uploading..."
    fi
fi

echo "Uploading expression h5ad to GCS..."
echo "  Source: $H5AD_FILE ($(du -h "$H5AD_FILE" | cut -f1))"
gsutil -q cp "$H5AD_FILE" "$GCS_H5AD"
echo "  Uploaded to $GCS_H5AD"
