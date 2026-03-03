#!/bin/bash
# Upload expression data to GCS for cNMF processing.
#
# Run this BEFORE launch_cnmf_vm.sh. It uploads the expression parquet
# that the VM will download and process.
#
# Prerequisites:
#   gcloud auth login
#   uv run snakemake extract_scrna_expression --cores 1  (generate the parquet first)
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

# --- 2. Upload expression data ---
PARQUET_FILE="${PROJECT_ROOT}/data/processed/mouse_abc/scrna_expression_log2cpm.parquet"
GCS_PARQUET="${GCS_BUCKET}/cnmf_input/scrna_expression_log2cpm.parquet"

if [ ! -f "$PARQUET_FILE" ]; then
    echo "ERROR: $PARQUET_FILE not found."
    echo "Generate it first: uv run snakemake extract_scrna_expression --cores 1"
    exit 1
fi

if gsutil -q stat "$GCS_PARQUET" 2>/dev/null; then
    LOCAL_SIZE=$(stat -f%z "$PARQUET_FILE" 2>/dev/null || stat -c%s "$PARQUET_FILE")
    REMOTE_SIZE=$(gsutil stat "$GCS_PARQUET" 2>/dev/null | grep "Content-Length" | awk '{print $2}')
    if [ "$LOCAL_SIZE" = "$REMOTE_SIZE" ]; then
        echo "Expression parquet already in GCS (same size), skipping upload."
        exit 0
    else
        echo "Expression parquet in GCS has different size, re-uploading..."
    fi
fi

echo "Uploading expression parquet to GCS..."
echo "  Source: $PARQUET_FILE ($(du -h "$PARQUET_FILE" | cut -f1))"
gsutil -q cp "$PARQUET_FILE" "$GCS_PARQUET"
echo "  Uploaded to $GCS_PARQUET"
