"""Download 10Xv3 scRNA-seq expression for hypothalamus cells from Allen ABC Atlas.

Downloads the pre-partitioned hypothalamus expression matrix (WMB-10Xv3-HY/log2)
which contains log2(CPM+1) values for all QC-passing cells. Also saves a metadata
sidecar with cell annotations (cluster, class, subclass, supertype).

Output:
    data/processed/mouse_abc/scrna_expression_cpm.h5ad     (~100-300k cells x ~30k genes, CPM values)
    data/processed/mouse_abc/scrna_cell_metadata.parquet    (cell annotations)

Usage:
    python -m hypomap.preprocessing.extract_scrna_expression
"""

from pathlib import Path

import anndata
import numpy as np
import pandas as pd
from abc_atlas_access.abc_atlas_cache.abc_project_cache import AbcProjectCache

DATA_DIR = Path(__file__).parent.parent.parent / "data"
CACHE_DIR = DATA_DIR / "raw" / "abc_atlas_cache"
OUTPUT_DIR = DATA_DIR / "processed" / "mouse_abc"

# ABC Atlas directory/file names for 10Xv3 scRNA-seq
CELL_METADATA_DIR = "WMB-10X"
CELL_METADATA_FILE = "cell_metadata"
GENE_METADATA_DIR = "WMB-10X"
GENE_METADATA_FILE = "gene"
EXPRESSION_DIR = "WMB-10Xv3"
EXPRESSION_FILE = "WMB-10Xv3-HY/log2"

# Feature matrix label for hypothalamus partition
HY_LABEL = "WMB-10Xv3-HY"

CHUNK_SIZE = 50_000

# Metadata columns to keep in the sidecar
METADATA_COLS = [
    "cell_label",
    "class",
    "subclass",
    "supertype",
    "cluster",
    "feature_matrix_label",
    "donor_label",
    "sex",
    "dataset_label",
    "region_of_interest_acronym",
]


def main():
    print("=== Extracting 10Xv3 scRNA-seq Expression for Hypothalamus ===\n")

    # 1. Initialize cache
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    print(f"Initializing ABC cache at: {CACHE_DIR}")
    cache = AbcProjectCache.from_cache_dir(CACHE_DIR)
    cache.load_latest_manifest()

    # 2. Load cell metadata
    print(f"Loading cell metadata from {CELL_METADATA_DIR}/{CELL_METADATA_FILE}...")
    cell_meta = cache.get_metadata_dataframe(
        directory=CELL_METADATA_DIR,
        file_name=CELL_METADATA_FILE,
        dtype={"cell_label": str},
    ).set_index("cell_label")
    print(f"  Total 10X cells: {len(cell_meta)}")

    # 3. Filter to hypothalamus partition
    hy_mask = cell_meta["feature_matrix_label"] == HY_LABEL
    hy_cells = cell_meta[hy_mask]
    print(f"  Hypothalamus cells ({HY_LABEL}): {len(hy_cells)}")

    if len(hy_cells) == 0:
        raise ValueError(
            f"No cells found with feature_matrix_label == '{HY_LABEL}'. "
            f"Available labels: {cell_meta['feature_matrix_label'].unique().tolist()[:10]}"
        )

    # 4. Load gene metadata
    print(f"Loading gene metadata from {GENE_METADATA_DIR}/{GENE_METADATA_FILE}...")
    gene_meta = cache.get_metadata_dataframe(
        directory=GENE_METADATA_DIR,
        file_name=GENE_METADATA_FILE,
    )
    print(f"  Total genes: {len(gene_meta)}")

    # 5. Download/open expression h5ad (backed for memory efficiency)
    print(f"Getting expression file: {EXPRESSION_DIR}/{EXPRESSION_FILE}")
    print("  This may take a while if downloading...")
    expr_path = cache.get_file_path(
        directory=EXPRESSION_DIR,
        file_name=EXPRESSION_FILE,
    )
    print(f"  Loading h5ad (backed) from {expr_path}")
    adata = anndata.read_h5ad(expr_path, backed="r")
    print(f"  Expression shape: {adata.shape}")

    # 6. Build gene symbol mapping
    # The h5ad var index has gene identifiers (Ensembl IDs); map to symbols.
    # Some Ensembl IDs share a gene symbol (pseudogenes, readthrough transcripts),
    # producing ~40 duplicate column names. We'll deduplicate after reading.
    if "gene_symbol" in adata.var.columns:
        gene_symbols = adata.var["gene_symbol"].values
    else:
        gene_symbols = adata.var.index.values

    n_unique = len(set(gene_symbols))
    n_total = len(gene_symbols)
    print(f"  Gene columns: {n_total} total, {n_unique} unique symbols")
    if n_total != n_unique:
        print(f"  ({n_total - n_unique} duplicate symbols will be deduplicated by max mean expression)")

    # 7. Identify cell indices in the h5ad that match our hypothalamus cells
    hy_cell_labels = set(hy_cells.index)
    h5ad_cell_labels = pd.Series(adata.obs.index)
    cell_mask = h5ad_cell_labels.isin(hy_cell_labels)
    cell_indices = np.where(cell_mask.values)[0]
    print(f"  Cells matched in h5ad: {len(cell_indices)} / {len(hy_cell_labels)}")

    if len(cell_indices) == 0:
        raise ValueError(
            "No hypothalamus cells found in the expression h5ad. "
            "Check that cell_label indexing is consistent."
        )

    # 8. Read expression in chunks
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

        n_done = end
        print(f"  Processed {n_done}/{len(cell_indices)} cells...")

    adata.file.close()

    # 9. Combine into DataFrame
    print("Combining chunks...")
    expr_matrix = np.vstack(expr_chunks)
    del expr_chunks  # free memory

    expr_df = pd.DataFrame(
        expr_matrix,
        index=cell_ids,
        columns=gene_symbols,
    )
    expr_df.index.name = "cell_label"
    print(f"  Expression DataFrame shape: {expr_df.shape}")

    # Deduplicate gene symbols: keep the column with highest mean expression
    if expr_df.columns.duplicated().any():
        n_dupes = expr_df.columns.duplicated().sum()
        print(f"  Deduplicating {n_dupes} duplicate gene symbols (keeping highest mean expression)...")
        mean_expr = expr_df.mean(axis=0)
        keep_mask = []
        seen = {}
        for i, col in enumerate(expr_df.columns):
            if col not in seen:
                seen[col] = i
                keep_mask.append(True)
            else:
                # Compare means: keep the better one
                if mean_expr.iloc[i] > mean_expr.iloc[seen[col]]:
                    keep_mask[seen[col]] = False
                    seen[col] = i
                    keep_mask.append(True)
                else:
                    keep_mask.append(False)
        expr_df = expr_df.iloc[:, keep_mask]
        print(f"  After dedup: {expr_df.shape[1]} unique genes")

    print(f"  Value range: [{expr_df.values.min():.2f}, {expr_df.values.max():.2f}]")
    print(f"  Memory usage: {expr_df.memory_usage(deep=True).sum() / 1e9:.1f} GB")

    # 10. Save expression as h5ad with CPM values
    # The source data is log2(CPM+1); reverse to CPM for downstream tools (e.g. cNMF)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    print("\nConverting log2(CPM+1) → CPM (float32)...")
    cpm_matrix = (np.power(2, expr_df.values) - 1).astype(np.float32)

    adata_out = anndata.AnnData(
        X=cpm_matrix,
        obs=pd.DataFrame(index=expr_df.index),
        var=pd.DataFrame(index=expr_df.columns),
    )
    adata_out.obs.index.name = "cell_label"
    adata_out.var.index.name = "gene_symbol"

    h5ad_path = OUTPUT_DIR / "scrna_expression_cpm.h5ad"
    print(f"Saving h5ad to {h5ad_path}...")
    adata_out.write_h5ad(h5ad_path, compression="gzip")
    print(f"  File size: {h5ad_path.stat().st_size / 1e9:.2f} GB")
    del cpm_matrix, adata_out

    # 11. Save cell metadata sidecar
    available_cols = [c for c in METADATA_COLS if c in hy_cells.columns or c == "cell_label"]
    meta_out = hy_cells.reset_index()
    meta_out = meta_out[[c for c in available_cols if c in meta_out.columns]]

    # Only keep cells that were actually in the expression data
    meta_out = meta_out[meta_out["cell_label"].isin(cell_ids)]

    meta_path = OUTPUT_DIR / "scrna_cell_metadata.parquet"
    print(f"Saving cell metadata to {meta_path}...")
    meta_out.to_parquet(meta_path, index=False)
    print(f"  {len(meta_out)} cells, {len(meta_out.columns)} columns")

    # Summary
    print(f"\n=== Done ===")
    print(f"  Expression: {expr_df.shape[0]} cells x {expr_df.shape[1]} genes")
    print(f"  Output: {h5ad_path}")
    print(f"  Metadata: {meta_path}")


if __name__ == "__main__":
    main()
