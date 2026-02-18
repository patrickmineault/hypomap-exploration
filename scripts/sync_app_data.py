"""Copy exactly the data files referenced by app/app.py from data/ to app/data/.

Usage:
    uv run python scripts/sync_app_data.py
"""

import shutil
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SRC = REPO_ROOT / "data"
DST = REPO_ROOT / "app" / "data"

# Shared files
SHARED_FILES = [
    "generated/mouse_common/np_map.csv",
    "generated/mouse_common/np_system_blacklist.csv",
    "generated/mouse_common/region_descriptions.csv",
    "generated/mouse_common/hormone_map.csv",
    "raw/mouse_abc/abc_cluster_annotations.csv/cluster_annotation-Table 1.csv",
    "raw/abc_atlas_cache/metadata/Allen-CCF-2020/20230630/views/parcellation_to_parcellation_term_membership_color.csv",
    "raw/abc_atlas_cache/metadata/Allen-CCF-2020/20230630/views/parcellation_to_parcellation_term_membership_acronym.csv",
]

# Per-dataset files
DATASETS = ["mouse_abc", "mouse_abc_extended"]
PER_DATASET_FILES = [
    "processed/{dataset}/cells_with_coords.parquet",
    "processed/{dataset}/coronal_atlas_regions.json",
    "processed/{dataset}/cluster_ligand_receptor_profile.parquet",
    "processed/{dataset}/cluster_np_expression.parquet",
]


def main():
    all_files = list(SHARED_FILES)
    for dataset in DATASETS:
        for template in PER_DATASET_FILES:
            all_files.append(template.format(dataset=dataset))

    copied = 0
    total_bytes = 0

    for rel_path in all_files:
        src = SRC / rel_path
        dst = DST / rel_path

        if not src.exists():
            print(f"  MISSING  {rel_path}")
            continue

        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
        size = src.stat().st_size
        total_bytes += size
        copied += 1
        print(f"  COPIED   {rel_path}  ({size / 1024 / 1024:.1f} MB)")

    print(f"\nDone: {copied}/{len(all_files)} files, {total_bytes / 1024 / 1024:.1f} MB total")


if __name__ == "__main__":
    main()
