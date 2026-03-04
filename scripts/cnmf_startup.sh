#!/bin/bash
# GCP VM startup script for running cNMF on hypothalamus scRNA-seq.
#
# This script runs automatically when the VM boots. It:
# 1. Installs system deps + uv
# 2. Clones the repo from GitHub and downloads data from GCS
# 3. Runs the full cNMF pipeline (prepare → factorize → combine → consensus)
# 4. Uploads results to GCS
# 5. Self-deletes the VM
#
# All output is logged to /var/log/cnmf.log (also visible via serial console).

set -uo pipefail
export HOME="${HOME:-/root}"

GCS_BUCKET="gs://neuroai-abc"
LOG=/var/log/cnmf.log
exec > >(tee -a "$LOG") 2>&1

# Upload logs to GCS on any error or exit
upload_logs() {
    echo ""
    echo "=== Uploading logs (exit code: $?) ==="
    echo "Time: $(date -u)"
    gsutil -q cp /var/log/cnmf*.log "${GCS_BUCKET}/cnmf_output/logs/" 2>/dev/null || true
}
trap upload_logs EXIT
set -e

echo "=== cNMF startup script ==="
echo "Time: $(date -u)"
echo "Instance: $(hostname)"

# --- Configuration ---
REPO_URL="https://github.com/patrickmineault/hypomap-exploration.git"
GIT_REF="main"  # branch/tag/commit to checkout
WORK_DIR="/opt/hypomap"
TRIAL_RUN=$(curl -sf -H "Metadata-Flavor: Google" \
    http://metadata.google.internal/computeMetadata/v1/instance/attributes/TRIAL_RUN 2>/dev/null || echo "false")
N_WORKERS=28  # n2-highmem-32 has 32 vCPUs; leave 4 for OS
K_VALUES="50"
N_ITER=100
SELECTED_K=50

if [ "$TRIAL_RUN" = "true" ]; then
    echo "*** TRIAL RUN MODE ***"
    N_WORKERS=1
fi

# --- 1. Install system dependencies ---
echo ""
echo "=== Installing system deps ==="
apt-get update -qq
apt-get install -y -qq python3-dev build-essential git curl

# Install uv
if ! command -v uv &>/dev/null; then
    echo "Installing uv..."
    curl -LsSf https://astral.sh/uv/install.sh | sh
    export PATH="$HOME/.local/bin:$PATH"
fi
echo "uv version: $(uv --version)"

# --- 2. Clone repo from GitHub ---
echo ""
echo "=== Cloning repo ==="
git clone --depth 1 --branch "$GIT_REF" "$REPO_URL" "$WORK_DIR"
cd "$WORK_DIR"
echo "  Checked out $(git rev-parse --short HEAD)"

# --- 3. Install Python deps ---
echo ""
echo "=== Installing Python environment ==="
cd "$WORK_DIR"
uv sync --no-dev 2>&1 | tail -5

# --- 4. Extract expression data from Allen ABC Atlas ---
echo ""
echo "=== Extracting scRNA-seq expression (downloads from Allen ABC Atlas) ==="
echo "Time: $(date -u)"
uv run python -m hypomap.preprocessing.extract_scrna_expression
echo "Time: $(date -u)"

# --- 5. Run cNMF pipeline ---
if [ "$TRIAL_RUN" = "true" ]; then
    echo ""
    echo "=== cNMF trial run (500 cells, 500 genes, K=10, 5 iter) ==="
    echo "Time: $(date -u)"
    uv run python scripts/run_cnmf.py --trial-run
else
    # 5a. Prepare
    echo ""
    echo "=== cNMF prepare ==="
    echo "Time: $(date -u)"
    uv run python scripts/run_cnmf.py \
        --step prepare \
        --k-values $K_VALUES \
        --n-iter $N_ITER

    # 5b. Factorize (parallel workers)
    echo ""
    echo "=== cNMF factorize ($N_WORKERS workers) ==="
    echo "Time: $(date -u)"

    pids=()
    for i in $(seq 0 $((N_WORKERS - 1))); do
        uv run python scripts/run_cnmf.py \
            --step factorize \
            --worker-index "$i" \
            --total-workers "$N_WORKERS" \
            > /var/log/cnmf_worker_${i}.log 2>&1 &
        pids+=($!)
    done

    echo "  Launched ${#pids[@]} workers: ${pids[*]}"

    # Wait for all workers, track failures
    failed=0
    for pid in "${pids[@]}"; do
        if ! wait "$pid"; then
            echo "  Worker PID $pid failed!"
            ((failed++))
        fi
    done

    if [ "$failed" -gt 0 ]; then
        echo "ERROR: $failed workers failed. Check /var/log/cnmf_worker_*.log"
        # Upload logs even on failure
        gsutil -q cp /var/log/cnmf*.log "${GCS_BUCKET}/cnmf_output/logs/"
        exit 1
    fi

    echo "  All $N_WORKERS workers completed successfully."

    # 5c. Combine
    echo ""
    echo "=== cNMF combine ==="
    echo "Time: $(date -u)"
    uv run python scripts/run_cnmf.py --step combine

    # 5d. Consensus
    echo ""
    echo "=== cNMF consensus (K=$SELECTED_K) ==="
    echo "Time: $(date -u)"
    uv run python scripts/run_cnmf.py \
        --step consensus \
        --selected-k $SELECTED_K
fi

# --- 6. Upload results to GCS ---
echo ""
echo "=== Uploading results ==="
echo "Time: $(date -u)"
gsutil -q -m cp -r "data/processed/mouse_abc/cnmf/hypo_cnmf/" "${GCS_BUCKET}/cnmf_output/" || true
echo "  Uploaded to ${GCS_BUCKET}/cnmf_output/"

echo ""
echo "=== cNMF pipeline complete ==="
echo "Time: $(date -u)"

# --- 7. Self-delete the VM (skip for trial runs so you can SSH in) ---
if [ "$TRIAL_RUN" = "true" ]; then
    echo "Trial run complete. VM left running for inspection."
    echo "Delete manually: gcloud compute instances delete $(hostname) --zone=\$(curl -sf -H 'Metadata-Flavor: Google' http://metadata.google.internal/computeMetadata/v1/instance/zone | awk -F/ '{print \$NF}') --quiet"
else
    echo "Self-deleting VM in 60 seconds..."
    sleep 60
    ZONE=$(curl -s -H "Metadata-Flavor: Google" http://metadata.google.internal/computeMetadata/v1/instance/zone | awk -F/ '{print $NF}')
    INSTANCE=$(hostname)
    gcloud compute instances delete "$INSTANCE" --zone="$ZONE" --quiet
fi
