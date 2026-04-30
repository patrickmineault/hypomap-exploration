#!/usr/bin/env bash
# download_lewislab.sh
# ====================
# Downloads the 10 mouse ligand-receptor pair databases from the Lewis Lab
# GitHub repository (LewisLabUCSD/Ligand-Receptor-Pairs).
#
# These are the raw inputs to build_decomposition.py (PART A).
#
# Usage (from project root):
#   bash scripts/download_lewislab_np.sh
#
# Requires: curl (or wget)

set -euo pipefail

DEST="data/raw/lewislab"
REPO="https://raw.githubusercontent.com/LewisLabUCSD/Ligand-Receptor-Pairs/master/Mouse"

FILES=(
    "Mouse-2023-Zhao-LR-pairs.tsv"
    "Mouse-2022-Dimitrov-LR-pairs.csv"
    "Mouse-2020-Jin-LR-pairs.csv"
    "Mouse-2020-Shao-LR-pairs.txt"
    "Mouse-2020-Baccin-LR-pairs.xlsx"
    "Mouse-2020-Cain-LR-pairs.xlsx"
    "Mouse-2019-Sheikh-LR-pairs.xlsx"
    "Mouse-2018-Skelly-LR-pairs.xlsx"
    "Mouse-2016-Yuzwa-LR-pairs.xlsx"
    "Mouse-2016-Ding-LR-pairs.xlsx"
)

mkdir -p "$DEST"

for f in "${FILES[@]}"; do
    if [ -f "$DEST/$f" ]; then
        echo "  SKIP  $f (already exists)"
    else
        echo "  FETCH $f"
        curl -sL "$REPO/$f" -o "$DEST/$f"
    fi
done

echo ""
echo "Downloaded ${#FILES[@]} files to $DEST/"
ls -lh "$DEST/"
