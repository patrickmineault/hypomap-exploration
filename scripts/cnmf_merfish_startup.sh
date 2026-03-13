#!/bin/bash
# GCP VM startup script for running cNMF on MERFISH expression within a bounding box.
#
# This script runs automatically when the VM boots. It:
# 1. Installs system deps + uv
# 2. Clones the repo and downloads MERFISH data from Allen ABC Atlas
# 3. Extracts expression for cells within the bounding box
# 4. Runs the full cNMF pipeline (prepare → factorize → combine → consensus)
# 5. Uploads results to GCS
# 6. Self-deletes the VM
#
# All output is logged to /var/log/cnmf.log (also visible via serial console).

set -uo pipefail
export HOME="${HOME:-/root}"

GCS_BUCKET="gs://neuroai-abc"
LOG=/var/log/cnmf.log

exec >> "$LOG" 2>&1

upload_logs() {
    echo ""
    echo "=== Uploading logs (exit code: $?) ==="
    echo "Time: $(date -u)"
    gcloud storage cp /var/log/cnmf*.log "${GCS_BUCKET}/cnmf_output/logs/" 2>/dev/null || true
}
trap upload_logs EXIT
set -e

heartbeat() {
    echo "[heartbeat] $1 — $(date -u)"
    gcloud storage cp "$LOG" "${GCS_BUCKET}/cnmf_output/logs/cnmf_merfish.log" 2>/dev/null || true
}

echo "=== cNMF MERFISH startup script ==="
echo "Time: $(date -u)"
echo "Instance: $(hostname)"

# --- Configuration from metadata ---
REPO_URL="https://github.com/patrickmineault/hypomap-exploration.git"
GIT_REF="main"
WORK_DIR="/opt/hypomap"

TRIAL_RUN=$(curl -sf -H "Metadata-Flavor: Google" \
    http://metadata.google.internal/computeMetadata/v1/instance/attributes/TRIAL_RUN 2>/dev/null || echo "false")
NEURONS_ONLY=$(curl -sf -H "Metadata-Flavor: Google" \
    http://metadata.google.internal/computeMetadata/v1/instance/attributes/NEURONS_ONLY 2>/dev/null || echo "false")
BBOX=$(curl -sf -H "Metadata-Flavor: Google" \
    http://metadata.google.internal/computeMetadata/v1/instance/attributes/BBOX 2>/dev/null || echo "-0.1 0.8 1.2 2.05 -2.5 -1.55")
CNMF_RUN_NAME=$(curl -sf -H "Metadata-Flavor: Google" \
    http://metadata.google.internal/computeMetadata/v1/instance/attributes/CNMF_RUN_NAME 2>/dev/null || echo "nmf_arh_me_vmh")

N_WORKERS=28  # n2-highmem-32 has 32 vCPUs; leave 4 for OS
K_VALUES="10 20 30 40 50"
N_ITER=100

if [ "$TRIAL_RUN" = "true" ]; then
    echo "*** TRIAL RUN MODE ***"
    N_WORKERS=1
fi

if [ "$NEURONS_ONLY" = "true" ]; then
    echo "*** NEURONS ONLY MODE ***"
fi

echo "BBOX: $BBOX"
echo "Run name: $CNMF_RUN_NAME"

# --- 1. Install system dependencies ---
echo ""
echo "=== Installing system deps ==="
apt-get update -qq
apt-get install -y -qq python3-dev build-essential git curl
heartbeat "System deps installed"

if ! command -v uv &>/dev/null; then
    echo "Installing uv..."
    curl -LsSf https://astral.sh/uv/install.sh | sh
    export PATH="$HOME/.local/bin:$PATH"
fi
echo "uv version: $(uv --version)"

# --- 2. Clone repo ---
echo ""
echo "=== Cloning repo ==="
git clone --depth 1 --branch "$GIT_REF" "$REPO_URL" "$WORK_DIR"
cd "$WORK_DIR"
echo "  Checked out $(git rev-parse --short HEAD)"
heartbeat "Repo cloned"

# --- 3. Install Python deps ---
echo ""
echo "=== Installing Python environment ==="
uv sync --no-dev 2>&1 | tail -5
heartbeat "Python deps installed"

# --- 4. Generate cells_with_coords.parquet (extract metadata + downsample) ---
echo ""
echo "=== Extracting ABC MERFISH cell metadata ==="
echo "Time: $(date -u)"
uv run python -m hypomap.datasets.mouse_abc \
    --divisions HY TH STR PAL P MY MB \
    --output-dir data/processed/mouse_abc_extended
heartbeat "ABC metadata extracted"

echo ""
echo "=== Downsampling to cells_with_coords ==="
echo "Time: $(date -u)"
uv run python -m hypomap.preprocessing.downsample --dataset mouse_abc_extended
heartbeat "Coordinate data generated"

# --- 5. Download MERFISH expression h5ad ---
# The h5ad is downloaded as a side effect of step 4 (abc_atlas_cache).
# Verify it exists; if not, download explicitly.
echo ""
echo "=== Ensuring MERFISH expression h5ad is cached ==="
echo "Time: $(date -u)"
uv run python -c "
from abc_atlas_access.abc_atlas_cache.abc_project_cache import AbcProjectCache
c = AbcProjectCache.from_cache_dir('data/raw/abc_atlas_cache')
c.load_latest_manifest()
p = c.get_file_path(directory='MERFISH-C57BL6J-638850', file_name='C57BL6J-638850/log2')
print(f'MERFISH h5ad: {p}')
"
heartbeat "MERFISH h5ad ready"

# --- 6. Extract bounding box expression ---
echo ""
echo "=== Extracting MERFISH expression for bounding box ==="
echo "Time: $(date -u)"
EXTRACT_FLAGS="--bbox $BBOX --output-name ${CNMF_RUN_NAME}_expression_cpm.h5ad"
if [ "$TRIAL_RUN" = "true" ]; then
    EXTRACT_FLAGS="$EXTRACT_FLAGS --max-cells 500"
fi
if [ "$NEURONS_ONLY" = "true" ]; then
    EXTRACT_FLAGS="$EXTRACT_FLAGS --neurons-only"
fi
uv run python -m hypomap.preprocessing.extract_merfish_expression $EXTRACT_FLAGS
echo "Time: $(date -u)"
heartbeat "MERFISH expression extraction complete"

INPUT_H5AD="data/processed/mouse_abc/${CNMF_RUN_NAME}_expression_cpm.h5ad"

# --- 7. Run cNMF pipeline ---
if [ "$TRIAL_RUN" = "true" ]; then
    echo ""
    echo "=== cNMF trial run ==="
    echo "Time: $(date -u)"
    uv run python scripts/run_cnmf.py --trial-run --run-name "$CNMF_RUN_NAME" --input-h5ad "$INPUT_H5AD"
    heartbeat "cNMF trial run complete"
else
    echo ""
    echo "=== cNMF prepare ==="
    echo "Time: $(date -u)"
    uv run python scripts/run_cnmf.py \
        --step prepare \
        --k-values $K_VALUES \
        --n-iter $N_ITER \
        --run-name "$CNMF_RUN_NAME" \
        --input-h5ad "$INPUT_H5AD"
    heartbeat "cNMF prepare complete"

    echo ""
    echo "=== cNMF factorize ($N_WORKERS workers) ==="
    echo "Time: $(date -u)"

    pids=()
    for i in $(seq 0 $((N_WORKERS - 1))); do
        uv run python scripts/run_cnmf.py \
            --step factorize \
            --worker-index "$i" \
            --total-workers "$N_WORKERS" \
            --run-name "$CNMF_RUN_NAME" \
            > /var/log/cnmf_worker_${i}.log 2>&1 &
        pids+=($!)
    done

    echo "  Launched ${#pids[@]} workers: ${pids[*]}"

    failed=0
    for pid in "${pids[@]}"; do
        if ! wait "$pid"; then
            echo "  Worker PID $pid failed!"
            ((failed++))
        fi
    done

    if [ "$failed" -gt 0 ]; then
        echo "ERROR: $failed workers failed. Check /var/log/cnmf_worker_*.log"
        gcloud storage cp /var/log/cnmf*.log "${GCS_BUCKET}/cnmf_output/logs/"
        exit 1
    fi

    echo "  All $N_WORKERS workers completed successfully."
    heartbeat "cNMF factorize complete"

    echo ""
    echo "=== cNMF combine ==="
    echo "Time: $(date -u)"
    uv run python scripts/run_cnmf.py --step combine --run-name "$CNMF_RUN_NAME"
    heartbeat "cNMF combine complete"

    echo ""
    echo "=== cNMF consensus (all K values) ==="
    echo "Time: $(date -u)"
    for K in $K_VALUES; do
        echo "  Running consensus for K=$K..."
        uv run python scripts/run_cnmf.py \
            --step consensus \
            --selected-k "$K" \
            --run-name "$CNMF_RUN_NAME"
        heartbeat "cNMF consensus K=$K complete"
    done
fi

# --- 8. Upload results to GCS ---
echo ""
echo "=== Uploading results ==="
echo "Time: $(date -u)"
gcloud storage cp -r "data/processed/mouse_abc/cnmf/${CNMF_RUN_NAME}/" "${GCS_BUCKET}/cnmf_output/" || true
gcloud storage cp "$INPUT_H5AD" "${GCS_BUCKET}/cnmf_output/${CNMF_RUN_NAME}/" || true
echo "  Uploaded to ${GCS_BUCKET}/cnmf_output/${CNMF_RUN_NAME}/"

echo ""
echo "=== cNMF MERFISH pipeline complete ==="
echo "Time: $(date -u)"

# --- 9. Self-delete the VM (skip for trial runs) ---
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
