"""Extract and process Langlieb et al. 2023 Slide-seq data.

Reads per-puck metadata and RCTD mapping matrices from Mapping_Matrices.tar.gz,
assigns dominant cell type per bead, converts CCF voxel coords to RAS mm, and
outputs a single cells_with_coords.parquet.
"""

import gzip
import io
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import mmread

DATA_DIR = Path(__file__).parent.parent.parent / "data"
RAW_DIR = DATA_DIR / "raw" / "mouse_langlieb" / "cellspatial_data"
TAR_PATH = RAW_DIR / "Mapping_Matrices.tar.gz"
OUTPUT_DIR = DATA_DIR / "processed" / "mouse_langlieb"
OUTPUT_PATH = OUTPUT_DIR / "cells_with_coords.parquet"

# CCF voxel size in mm (25 µm resolution)
VOXEL_MM = 0.025

# Neurotransmitter prefixes in Langlieb cell type names
NT_PREFIXES = ("Ex_", "Inh_", "Chol_", "Dop_", "Nor_", "Ser_")


def parse_neurotransmitter(name: str) -> str:
    """Extract neurotransmitter class from a Langlieb cell type name."""
    for prefix in NT_PREFIXES:
        if name.startswith(prefix):
            return prefix.rstrip("_")
    return "Non-neuronal"


def load_puck_from_tar(tar: tarfile.TarFile, puck_id: str) -> pd.DataFrame | None:
    """Load one puck's metadata + dominant cell type from an open tar.

    Returns a DataFrame with columns:
        cell_id, CCF_X, CCF_Y, CCF_Z, TopStruct, DeepCCF, CCF_acronym,
        IsOutsideCCF, puck_id, cluster, neurotransmitter
    or None if required files are missing.
    """
    prefix = f"Mapping_Matrices/{puck_id}"
    meta_name = f"{prefix}.mapping.metadata.tsv.gz"
    mtx_name = f"{prefix}.mapping.mtx.gz"
    ct_name = f"{prefix}.mapping.MappedCellTypes.txt.gz"

    # Check all three files exist
    try:
        meta_info = tar.getmember(meta_name)
        mtx_info = tar.getmember(mtx_name)
        ct_info = tar.getmember(ct_name)
    except KeyError:
        return None

    # Load metadata
    with tar.extractfile(meta_info) as f:
        raw = gzip.decompress(f.read())
    meta = pd.read_csv(
        io.BytesIO(raw),
        sep="\t",
        index_col=0,
        dtype={"CCF_X": np.int32, "CCF_Y": np.int32, "CCF_Z": np.int32},
    )

    # Load cell type labels
    with tar.extractfile(ct_info) as f:
        raw = gzip.decompress(f.read())
    ct_lines = raw.decode("utf-8").strip().split("\n")
    ct_keys = []
    ct_names = []
    for line in ct_lines:
        key, name = line.split("=", 1)
        ct_keys.append(key)
        ct_names.append(name)
    ct_names_arr = np.array(ct_names)

    # Load sparse RCTD weight matrix (beads x cell types)
    with tar.extractfile(mtx_info) as f:
        raw = gzip.decompress(f.read())
    mat = mmread(io.BytesIO(raw)).tocsr()  # (n_beads, n_cell_types)

    # Assign dominant cell type per bead (argmax of RCTD weights)
    # Beads with all-zero weights get index 0 — we'll filter these via IsOutsideCCF
    n_beads = mat.shape[0]
    dominant_idx = np.zeros(n_beads, dtype=np.int32)
    for i in range(n_beads):
        row = mat.getrow(i)
        if row.nnz > 0:
            dominant_idx[i] = row.indices[row.data.argmax()]

    # Align: mtx rows correspond 1:1 to metadata rows (same order)
    if len(meta) != n_beads:
        print(
            f"  Warning: {puck_id} metadata rows ({len(meta)}) != "
            f"mtx rows ({n_beads}), skipping"
        )
        return None

    meta = meta.copy()
    meta["cluster"] = ct_names_arr[dominant_idx]
    meta["neurotransmitter"] = meta["cluster"].map(parse_neurotransmitter)
    meta["puck_id"] = puck_id
    meta["cell_id"] = meta.index.astype(str)

    keep_cols = [
        "cell_id", "CCF_X", "CCF_Y", "CCF_Z",
        "TopStruct", "DeepCCF", "CCF_acronym",
        "IsOutsideCCF", "puck_id", "cluster", "neurotransmitter",
    ]
    return meta[keep_cols].reset_index(drop=True)


def extract_langlieb() -> pd.DataFrame:
    """Extract all pucks and produce cells_with_coords.parquet."""
    if not TAR_PATH.exists():
        raise FileNotFoundError(f"Mapping_Matrices.tar.gz not found at {TAR_PATH}")

    print(f"Opening {TAR_PATH} ...")
    puck_dfs = []

    with tarfile.open(TAR_PATH, "r:gz") as tar:
        # Discover puck numbers from tar members
        puck_ids = sorted(
            {
                m.name.split("/")[1].split(".")[0]
                for m in tar.getmembers()
                if m.name.startswith("Mapping_Matrices/Puck_Num_")
            }
        )
        print(f"Found {len(puck_ids)} pucks")

        for i, puck_id in enumerate(puck_ids):
            df = load_puck_from_tar(tar, puck_id)
            if df is not None:
                puck_dfs.append(df)
                if (i + 1) % 20 == 0:
                    print(f"  Processed {i + 1}/{len(puck_ids)} pucks ...")

    print(f"Concatenating {len(puck_dfs)} pucks ...")
    cells = pd.concat(puck_dfs, ignore_index=True)
    print(f"Total beads: {len(cells):,}")

    # Filter out beads outside CCF
    n_before = len(cells)
    cells = cells[cells["IsOutsideCCF"] != True].copy()  # noqa: E712
    # Handle string "TRUE"/"FALSE" values
    cells = cells[cells["IsOutsideCCF"].astype(str).str.upper() != "TRUE"].copy()
    print(f"After IsOutsideCCF filter: {len(cells):,} (removed {n_before - len(cells):,})")

    # Filter beads with empty TopStruct
    cells = cells[cells["TopStruct"].notna() & (cells["TopStruct"] != "")].copy()
    print(f"After TopStruct filter: {len(cells):,}")

    # Convert CCF voxel indices to mm then to RAS convention
    # Langlieb CCF coords are 25 µm voxel indices into the Allen CCF volume
    # Allen CCF axes: X = anterior-posterior, Y = dorsal-ventral, Z = left-right
    # RAS convention (matching mouse_abc.py):
    #   x_ras = Z_ccf * 0.025 - 5.7   (mediolateral, centered at midline)
    #   z_ras = Y_ccf * 0.025 - 5.4   (dorsoventral)
    #   y_ras = -(X_ccf * 0.025 - 6.78) - 1.77  (anteroposterior)
    cells["x"] = (cells["CCF_Z"] * VOXEL_MM - 5.7).astype(np.float32)
    cells["y"] = (cells["CCF_Y"] * VOXEL_MM - 5.4).astype(np.float32)

    # Each Slide-seq puck is a physical coronal section, but it's slightly
    # tilted relative to the CCF coronal plane, so individual bead CCF_X
    # values span ~0.3-0.5mm within a single puck.  Using per-bead z would
    # split each puck into many thin slices (streak artifacts).  Instead,
    # assign every bead in a puck a single AP (z) value derived from the
    # bead closest to the origin (x=0, y=0), matching the ABC pipeline.
    _z_map = {}
    for puck_id, grp in cells.groupby("puck_id"):
        dist = np.sqrt(grp["x"].values ** 2 + grp["y"].values ** 2)
        origin_ccf_x = grp["CCF_X"].values[dist.argmin()]
        _z_map[puck_id] = round(-(origin_ccf_x * VOXEL_MM - 6.78) - 1.77, 2)
    cells["z"] = cells["puck_id"].map(_z_map).astype(np.float32)

    # Rename TopStruct to region for consistency with load_cell_data()
    cells = cells.rename(columns={"TopStruct": "region"})

    # Select final columns
    out_cols = [
        "cell_id", "x", "y", "z",
        "cluster", "neurotransmitter", "region",
        "DeepCCF", "CCF_acronym", "puck_id",
    ]
    cells = cells[out_cols]

    # Save
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    cells.to_parquet(OUTPUT_PATH, index=False)
    print(f"\nSaved {len(cells):,} beads to {OUTPUT_PATH}")

    # Summary
    print(f"\n=== Langlieb Dataset Summary ===")
    print(f"Total beads: {len(cells):,}")
    print(f"Pucks: {cells['puck_id'].nunique()}")
    print(f"Slices (unique z): {cells['z'].nunique()}")
    print(f"Cell types (cluster): {cells['cluster'].nunique()}")
    print(f"Neurotransmitter classes: {cells['neurotransmitter'].nunique()}")
    print(f"\nRegions:")
    print(cells["region"].value_counts().to_string())
    print(f"\nx range: [{cells['x'].min():.2f}, {cells['x'].max():.2f}] mm")
    print(f"y range: [{cells['y'].min():.2f}, {cells['y'].max():.2f}] mm")
    print(f"z range: [{cells['z'].min():.2f}, {cells['z'].max():.2f}] mm")

    return cells


if __name__ == "__main__":
    extract_langlieb()
