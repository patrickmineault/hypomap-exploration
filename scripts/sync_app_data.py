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
    "processed/mouse_common/np_map.csv",
    "generated/mouse_common/np_system_blacklist.csv",
    "generated/mouse_common/region_descriptions.csv",
    "generated/mouse_common/hormone_map.csv",
    "raw/mouse_abc/abc_cluster_annotations.csv/cluster_annotation-Table 1.csv",
    "raw/abc_atlas_cache/metadata/Allen-CCF-2020/20230630/views/parcellation_to_parcellation_term_membership_color.csv",
    "raw/abc_atlas_cache/metadata/Allen-CCF-2020/20230630/views/parcellation_to_parcellation_term_membership_acronym.csv",
]

# Per-dataset files
DATASETS = ["mouse_abc", "mouse_abc_extended", "mouse_abc_whole_brain", "mouse_langlieb"]
PER_DATASET_FILES = [
    "processed/{dataset}/cells_with_coords.parquet",
    "processed/{dataset}/coronal_atlas_regions.msgpack",
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

    print(f"\nCopied: {copied}/{len(all_files)} files, {total_bytes / 1024 / 1024:.1f} MB total")

    # Delete files in DST that are not in the manifest
    expected = {Path(rel) for rel in all_files}
    removed = 0
    for dst_file in DST.rglob("*"):
        if not dst_file.is_file():
            continue
        rel = dst_file.relative_to(DST)
        if rel not in expected:
            dst_file.unlink()
            removed += 1
            print(f"  DELETED  {rel}")

    # Remove empty directories left behind
    for dst_dir in sorted(DST.rglob("*"), reverse=True):
        if dst_dir.is_dir() and not any(dst_dir.iterdir()):
            dst_dir.rmdir()

    print(f"Deleted: {removed} stale files")

    # Check total app/ directory size
    app_dir = REPO_ROOT / "app"
    app_bytes = sum(f.stat().st_size for f in app_dir.rglob("*") if f.is_file())
    app_mb = app_bytes / 1024 / 1024
    print(f"\nTotal app/ size: {app_mb:.1f} MB")
    if app_mb > 200:
        print(f"\033[91mERROR: app/ is {app_mb:.0f} MB — too large for Plotly Cloud (200 MB limit)\033[0m")
    elif app_mb > 150:
        print(f"\033[93mWARNING: app/ is {app_mb:.0f} MB — approaching Plotly Cloud 200 MB limit\033[0m")


if __name__ == "__main__":
    main()
