"""Extract MERFISH expression for a spatial bounding box from Allen ABC Atlas.

Uses the pre-processed cells_with_coords.parquet for spatial filtering, then
reads expression values from the cached MERFISH h5ad (550 genes, log2(CPM+1)).
Outputs CPM-valued h5ad suitable for cNMF.

Output:
    data/processed/mouse_abc/merfish_expression_cpm.h5ad   (cells x 550 genes, CPM)

Usage:
    python -m hypomap.preprocessing.extract_merfish_expression \\
        --bbox -0.1 0.8 1.2 2.05 -2.5 -1.55
    python -m hypomap.preprocessing.extract_merfish_expression \\
        --bbox -0.1 0.8 1.2 2.05 -2.5 -1.55 --neurons-only --max-cells 500
"""

import argparse
from pathlib import Path

import anndata
import numpy as np
import pandas as pd

DATA_DIR = Path(__file__).parent.parent.parent / "data"
COORDS_PATH = DATA_DIR / "processed" / "mouse_abc_extended" / "cells_with_coords.parquet"
MERFISH_H5AD = (
    DATA_DIR / "raw" / "abc_atlas_cache" / "expression_matrices"
    / "MERFISH-C57BL6J-638850" / "20230830" / "C57BL6J-638850-log2.h5ad"
)
OUTPUT_DIR = DATA_DIR / "processed" / "mouse_abc"

CHUNK_SIZE = 50_000


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--bbox",
        type=float,
        nargs=6,
        required=True,
        metavar=("X_MIN", "X_MAX", "Y_MIN", "Y_MAX", "Z_MIN", "Z_MAX"),
        help="Bounding box in RAS coordinates: x_min x_max y_min y_max z_min z_max",
    )
    parser.add_argument(
        "--neurons-only",
        action="store_true",
        help="Filter to neuronal cells only (is_neuron == True)",
    )
    parser.add_argument(
        "--max-cells",
        type=int,
        default=None,
        help="Limit to N cells (for trial runs)",
    )
    parser.add_argument(
        "--output-name",
        type=str,
        default="merfish_expression_cpm.h5ad",
        help="Output filename (default: merfish_expression_cpm.h5ad)",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    x_min, x_max, y_min, y_max, z_min, z_max = args.bbox

    print("=== Extracting MERFISH Expression for Bounding Box ===\n")
    print(f"  Bbox: x=[{x_min}, {x_max}], y=[{y_min}, {y_max}], z=[{z_min}, {z_max}]")
    if args.neurons_only:
        print("  Mode: neurons only")
    if args.max_cells:
        print(f"  Max cells: {args.max_cells}")

    # 1. Load coordinate data and filter to bounding box
    print(f"\nLoading coordinates from {COORDS_PATH}...")
    if not COORDS_PATH.exists():
        raise FileNotFoundError(
            f"Coordinate file not found: {COORDS_PATH}\n"
            "Run the Snakemake pipeline first to generate cells_with_coords.parquet"
        )
    cells = pd.read_parquet(COORDS_PATH)
    cells = cells[~cells["region"].str.contains("unassigned", case=False)]
    print(f"  Total cells (excluding unassigned): {len(cells):,}")

    mask = (
        (cells["x"] >= x_min) & (cells["x"] <= x_max)
        & (cells["y"] >= y_min) & (cells["y"] <= y_max)
        & (cells["z"] >= z_min) & (cells["z"] <= z_max)
    )
    cells = cells[mask]
    print(f"  Cells in bounding box: {len(cells):,}")

    if args.neurons_only:
        n_before = len(cells)
        cells = cells[cells["is_neuron"] == True]  # noqa: E712
        print(f"  Neurons-only filter: {n_before} → {len(cells)} cells")

    if len(cells) == 0:
        raise ValueError("No cells found in bounding box")

    # 2. Subsample if requested
    if args.max_cells is not None and len(cells) > args.max_cells:
        rng = np.random.default_rng(42)
        cells = cells.iloc[rng.choice(len(cells), size=args.max_cells, replace=False)]
        print(f"  Subsampled to {len(cells)} cells (--max-cells {args.max_cells})")

    keep_ids = set(cells["cell_id"])
    print(f"\n  Regions: {sorted(cells['region'].unique())}")
    print(f"  Classes: {cells['class'].nunique()} unique")

    # 3. Open MERFISH expression h5ad
    print(f"\nOpening MERFISH h5ad: {MERFISH_H5AD}")
    if not MERFISH_H5AD.exists():
        raise FileNotFoundError(
            f"MERFISH h5ad not found: {MERFISH_H5AD}\n"
            "Download it via: uv run python -c \"from abc_atlas_access.abc_atlas_cache.abc_project_cache import AbcProjectCache; "
            "c = AbcProjectCache.from_cache_dir('data/raw/abc_atlas_cache'); c.load_latest_manifest(); "
            "c.get_file_path(directory='MERFISH-C57BL6J-638850', file_name='C57BL6J-638850/log2')\""
        )
    adata = anndata.read_h5ad(MERFISH_H5AD, backed="r")
    print(f"  Shape: {adata.shape}")

    # 4. Find matching cell indices
    h5ad_labels = pd.Series(adata.obs.index)
    cell_mask = h5ad_labels.isin(keep_ids)
    cell_indices = np.where(cell_mask.values)[0]
    print(f"  Matched {len(cell_indices):,} / {len(keep_ids):,} cells in h5ad")

    if len(cell_indices) == 0:
        raise ValueError("No bounding box cells found in MERFISH h5ad")

    # 5. Get gene symbols
    if "gene_symbol" in adata.var.columns:
        gene_symbols = adata.var["gene_symbol"].values
    else:
        gene_symbols = adata.var.index.values
    print(f"  Genes: {len(gene_symbols)}")

    # 6. Read expression in chunks
    print(f"\nReading expression in chunks of {CHUNK_SIZE}...")
    expr_chunks = []
    cell_ids = []

    for start in range(0, len(cell_indices), CHUNK_SIZE):
        end = min(start + CHUNK_SIZE, len(cell_indices))
        chunk_idx = cell_indices[start:end]

        chunk = adata.X[chunk_idx, :]
        if hasattr(chunk, "toarray"):
            chunk = chunk.toarray()
        expr_chunks.append(chunk)
        cell_ids.extend(adata.obs.index[chunk_idx].tolist())
        print(f"  Processed {end}/{len(cell_indices)} cells...")

    adata.file.close()

    # 7. Combine and convert to CPM
    print("Combining chunks...")
    expr_matrix = np.vstack(expr_chunks)
    del expr_chunks

    # Deduplicate gene symbols (keep highest mean expression)
    n_unique = len(set(gene_symbols))
    if n_unique < len(gene_symbols):
        print(f"  Deduplicating {len(gene_symbols) - n_unique} duplicate gene symbols...")
        expr_df = pd.DataFrame(expr_matrix, index=cell_ids, columns=gene_symbols)
        mean_expr = expr_df.mean(axis=0)
        keep_mask = []
        seen = {}
        for i, col in enumerate(expr_df.columns):
            if col not in seen:
                seen[col] = i
                keep_mask.append(True)
            else:
                if mean_expr.iloc[i] > mean_expr.iloc[seen[col]]:
                    keep_mask[seen[col]] = False
                    seen[col] = i
                    keep_mask.append(True)
                else:
                    keep_mask.append(False)
        expr_df = expr_df.iloc[:, keep_mask]
        gene_symbols = expr_df.columns.values
        expr_matrix = expr_df.values
        del expr_df

    print(f"  Expression shape: {expr_matrix.shape}")
    print(f"  Value range (log2): [{expr_matrix.min():.2f}, {expr_matrix.max():.2f}]")

    # Convert log2(CPM+1) → CPM
    print("Converting log2(CPM+1) → CPM (float32)...")
    cpm_matrix = (np.power(2, expr_matrix.astype(np.float32)) - 1)
    del expr_matrix

    # 8. Save as h5ad
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    adata_out = anndata.AnnData(
        X=cpm_matrix,
        obs=pd.DataFrame(index=cell_ids),
        var=pd.DataFrame(index=gene_symbols),
    )
    adata_out.obs.index.name = "cell_label"
    adata_out.var.index.name = "gene_symbol"

    out_path = OUTPUT_DIR / args.output_name
    print(f"Saving h5ad to {out_path}...")
    adata_out.write_h5ad(out_path, compression="gzip")
    print(f"  File size: {out_path.stat().st_size / 1e6:.1f} MB")

    # 9. Save cell metadata sidecar
    meta_out = cells[cells["cell_id"].isin(cell_ids)].copy()
    meta_path = out_path.with_suffix("").with_name(out_path.stem + "_metadata.parquet")
    meta_out.to_parquet(meta_path, index=False)
    print(f"  Metadata: {len(meta_out)} cells → {meta_path}")

    print(f"\n=== Done ===")
    print(f"  {len(cell_ids)} cells x {len(gene_symbols)} genes")
    print(f"  Output: {out_path}")


if __name__ == "__main__":
    main()
