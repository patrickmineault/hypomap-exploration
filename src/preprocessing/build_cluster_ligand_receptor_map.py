"""Map neuropeptide and hormone ligand/receptor expression to ABC cell clusters.

Fetches expression data from the Allen Brain Cell Atlas for neuropeptide
and hormone ligands/receptors, then aggregates by cluster to create expression profiles.

Usage:
    # Using measured MERFISH panel (subset of genes):
    python -m src.preprocessing.build_cluster_ligand_receptor_map

    # Using imputed data (all genes, ~50GB download):
    python -m src.preprocessing.build_cluster_ligand_receptor_map --use-imputed
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

# Paths
DATA_DIR = Path(__file__).parent.parent.parent / "data"
CACHE_DIR = DATA_DIR / "raw" / "abc_atlas_cache"
OUTPUT_DIR = DATA_DIR / "processed" / "mouse_abc"
NP_MAP_PATH = DATA_DIR / "generated" / "mouse_common" / "np_map.csv"
HORMONE_MAP_PATH = DATA_DIR / "generated" / "mouse_common" / "hormone_map.csv"

# Imputed data paths (external drive for large files)
IMPUTED_CACHE_DIR = Path("/Volumes/ExtDrive/mouse_abc_data")
IMPUTED_DIRECTORY = "MERFISH-C57BL6J-638850-imputed"
IMPUTED_EXPR_FILE = "C57BL6J-638850-imputed/log2"

# Expression threshold for "expressing" (log2 scale)
# Similar to the the Allen Brain Atlas MERFISH documentation
EXPRESSION_THRESHOLD = 3.0


def get_ligand_receptor_genes() -> tuple[list[str], set[str], set[str]]:
    """Load ligand/receptor genes from np_map.csv and hormone_map.csv.

    Handles semicolon-separated gene names (e.g., "Calcrl;Ramp2" for heterodimer
    receptors) by splitting them into individual genes for expression lookup.

    Returns:
        Tuple of (all_genes list, ligand_genes set, receptor_genes set)
    """
    ligand_genes = set()
    receptor_genes = set()

    # Load neuropeptide genes
    np_map = pd.read_csv(NP_MAP_PATH)
    for lg in np_map["Ligand_Gene"].dropna().unique():
        for g in lg.split(";"):
            ligand_genes.add(g.strip())
    for rg in np_map["Receptor_Gene"].dropna().unique():
        for g in rg.split(";"):
            receptor_genes.add(g.strip())

    print(f"  Loaded {len(np_map)} ligand-receptor pairs from np_map.csv")
    for lclass in np_map["Ligand_Class"].unique():
        n = len(np_map[np_map["Ligand_Class"] == lclass])
        print(f"    {lclass}: {n} interactions")

    # Load hormone genes
    if HORMONE_MAP_PATH.exists():
        hormone_map = pd.read_csv(HORMONE_MAP_PATH)
        for lg in hormone_map["Ligand_Gene"].dropna().unique():
            for g in lg.split(";"):
                ligand_genes.add(g.strip())
        for rg in hormone_map["Receptor_Gene"].dropna().unique():
            for g in rg.split(";"):
                receptor_genes.add(g.strip())

        print(f"  Loaded {len(hormone_map)} ligand-receptor pairs from hormone_map.csv")
        if "hormone_class" in hormone_map.columns:
            for hclass in hormone_map["hormone_class"].dropna().unique():
                n = len(hormone_map[hormone_map["hormone_class"] == hclass])
                print(f"    {hclass}: {n} interactions")
    else:
        print(f"  Warning: hormone_map.csv not found at {HORMONE_MAP_PATH}")

    all_genes = sorted(ligand_genes | receptor_genes)

    return all_genes, ligand_genes, receptor_genes


def fetch_expression_data(
    cells: pd.DataFrame, genes: list[str], use_imputed: bool = False
) -> tuple[pd.DataFrame, list[str], list[str]]:
    """Fetch expression data from ABC Atlas for specified genes.

    Args:
        cells: Cell metadata DataFrame with cell_id as index
        genes: List of gene symbols to fetch
        use_imputed: If True, use imputed dataset (all genes) from external drive

    Returns:
        Tuple of (expression DataFrame, genes_found, genes_missing)
    """
    import anndata
    from abc_atlas_access.abc_atlas_cache.abc_project_cache import AbcProjectCache

    if use_imputed:
        # Use imputed dataset from external drive
        cache_dir = IMPUTED_CACHE_DIR
        directory = IMPUTED_DIRECTORY
        expr_file = IMPUTED_EXPR_FILE
        print(f"Using IMPUTED dataset (all genes)")
        print(f"Cache directory: {cache_dir}")
    else:
        # Use measured MERFISH panel
        cache_dir = CACHE_DIR
        directory = "MERFISH-C57BL6J-638850"
        expr_file = "C57BL6J-638850/log2"
        print(f"Using MEASURED MERFISH panel")

    print(f"Initializing ABC Atlas cache at {cache_dir}")
    cache = AbcProjectCache.from_cache_dir(cache_dir)
    cache.load_latest_manifest()

    # Get path to expression file (will download if not present)
    print(f"Getting expression file: {directory}/{expr_file}")
    print("This may take a while if downloading (~50GB for imputed data)...")
    expr_path = cache.get_file_path(directory=directory, file_name=expr_file)

    print(f"Loading expression data from {expr_path}")
    adata = anndata.read_h5ad(expr_path, backed="r")

    # Check which genes are available
    measured_genes = set(adata.var["gene_symbol"].tolist())
    genes_found = [g for g in genes if g in measured_genes]
    genes_missing = [g for g in genes if g not in measured_genes]

    print(f"Genes requested: {len(genes)}")
    print(f"Genes available: {len(genes_found)}")
    if genes_missing:
        print(f"Genes missing: {len(genes_missing)}")
        print(
            f"  Missing: {genes_missing[:10]}{'...' if len(genes_missing) > 10 else ''}"
        )

    if not genes_found:
        raise ValueError("No requested genes found in expression data!")

    # Create gene mask
    gene_mask = adata.var["gene_symbol"].isin(genes_found)
    gene_indices = np.where(gene_mask)[0]

    # Filter to hypothalamus cells
    hy_cell_ids = set(cells.index)
    cell_mask = pd.Series(adata.obs.index).isin(hy_cell_ids)
    cell_indices = np.where(cell_mask)[0]

    print(
        f"Extracting expression for {len(genes_found)} genes across {len(cell_indices)} cells..."
    )

    # Extract expression data in chunks to avoid memory issues
    chunk_size = 50000
    expr_data = []
    cell_ids = []

    for start in range(0, len(cell_indices), chunk_size):
        end = min(start + chunk_size, len(cell_indices))
        chunk_idx = cell_indices[start:end]

        # Get chunk of expression data
        chunk = adata.X[chunk_idx, :][:, gene_indices]
        # Handle both sparse and dense arrays
        if hasattr(chunk, "toarray"):
            chunk = chunk.toarray()
        expr_data.append(chunk)
        cell_ids.extend(adata.obs.index[chunk_idx].tolist())

        if (start // chunk_size + 1) % 3 == 0:
            print(f"  Processed {end}/{len(cell_indices)} cells...")

    adata.file.close()

    # Combine into DataFrame
    expr_matrix = np.vstack(expr_data)
    gene_names = adata.var.loc[gene_mask, "gene_symbol"].tolist()
    expr_df = pd.DataFrame(expr_matrix, index=cell_ids, columns=gene_names)

    return expr_df, genes_found, genes_missing


def aggregate_by_cluster(
    cells: pd.DataFrame,
    expr_df: pd.DataFrame,
    genes: list[str],
    ligand_genes: set[str],
    receptor_genes: set[str],
) -> pd.DataFrame:
    """Aggregate expression statistics by cluster.

    Args:
        cells: Cell metadata with 'cluster' column
        expr_df: Expression DataFrame (cells x genes)
        genes: List of gene names
        ligand_genes: Set of ligand gene names
        receptor_genes: Set of receptor gene names

    Returns:
        DataFrame with cluster-level statistics
    """
    print("Aggregating expression by cluster...")

    # Join cluster info with expression
    cluster_col = cells[["cluster"]].copy()
    cells_with_expr = cluster_col.join(expr_df, how="inner")

    # Compute stats per cluster
    cluster_stats = []
    clusters = cells_with_expr["cluster"].unique()

    for i, cluster in enumerate(clusters):
        if (i + 1) % 500 == 0:
            print(f"  Processed {i+1}/{len(clusters)} clusters...")

        grp = cells_with_expr[cells_with_expr["cluster"] == cluster]
        n_cells = len(grp)

        for gene in genes:
            if gene not in grp.columns:
                continue

            vals = grp[gene].dropna()
            if len(vals) == 0:
                continue

            cluster_stats.append(
                {
                    "cluster": cluster,
                    "gene": gene,
                    "mean_expr": vals.mean(),
                    "median_expr": vals.median(),
                    "max_expr": vals.max(),
                    "pct_expressing": (vals > EXPRESSION_THRESHOLD).mean() * 100,
                    "n_cells": n_cells,
                    "n_expressing": (vals > EXPRESSION_THRESHOLD).sum(),
                    "is_ligand": gene in ligand_genes,
                    "is_receptor": gene in receptor_genes,
                }
            )

    return pd.DataFrame(cluster_stats)


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Map neuropeptide ligand/receptor expression to ABC cell clusters."
    )
    parser.add_argument(
        "--use-imputed",
        action="store_true",
        help="Use imputed expression data (all genes, ~50GB download to external drive)",
    )
    parser.add_argument(
        "--metadata-dir",
        type=Path,
        default=None,
        help="Directory containing cell_metadata.parquet (default: data/processed/mouse_abc)",
    )
    return parser.parse_args()


def main():
    """Main function to build cluster ligand-receptor map."""

    args = parse_args()

    print("=== Building Cluster Ligand-Receptor Expression Map ===\n")

    if args.use_imputed:
        print("MODE: Using IMPUTED expression data (all genes)")
        print(f"Data will be downloaded to: {IMPUTED_CACHE_DIR}\n")
    else:
        print("MODE: Using MEASURED MERFISH panel (subset of genes)")
        print("Use --use-imputed for all genes\n")

    # 1. Load genes
    print("Loading neuropeptide and hormone ligand/receptor genes...")
    genes, ligand_genes, receptor_genes = get_ligand_receptor_genes()
    print(f"  Ligands: {len(ligand_genes)}")
    print(f"  Receptors: {len(receptor_genes)}")
    print(f"  Total unique genes: {len(genes)}")

    # 2. Load cell metadata
    metadata_dir = args.metadata_dir if args.metadata_dir else OUTPUT_DIR
    print(f"\nLoading cell metadata from {metadata_dir}...")
    cells = pd.read_parquet(metadata_dir / "cell_metadata.parquet")
    cells = cells.set_index("cell_id") if "cell_id" in cells.columns else cells
    print(f"  Cells: {len(cells)}")
    print(f"  Clusters: {cells['cluster'].nunique()}")

    # 3. Fetch expression data
    print("\nFetching expression data from ABC Atlas...")
    expr_df, genes_found, genes_missing = fetch_expression_data(
        cells, genes, use_imputed=args.use_imputed
    )

    # 4. Save per-cell expression
    expr_path = metadata_dir / "neuropeptide_expression.parquet"
    expr_df.to_parquet(expr_path)
    print(f"\nSaved per-cell expression to {expr_path}")
    print(f"  Shape: {expr_df.shape}")

    # Save list of missing genes for reference
    if genes_missing:
        missing_path = metadata_dir / "ligand_receptor_genes_missing.txt"
        with open(missing_path, "w") as f:
            f.write(
                "# Ligand/receptor genes not in MERFISH panel (would need imputed data)\n"
            )
            for g in sorted(genes_missing):
                f.write(f"{g}\n")
        print(f"  Missing genes saved to {missing_path}")

    # 5. Aggregate by cluster (only for genes we have)
    cluster_df = aggregate_by_cluster(
        cells, expr_df, genes_found, ligand_genes, receptor_genes
    )

    # 6. Save cluster profiles
    profile_path = metadata_dir / "cluster_ligand_receptor_profile.parquet"
    cluster_df.to_parquet(profile_path)
    print(f"\nSaved cluster profiles to {profile_path}")
    print(f"  Shape: {cluster_df.shape}")

    # 7. Summary stats
    print("\n=== Summary ===")
    expressing_clusters = cluster_df[cluster_df["pct_expressing"] > 10]
    print(f"Cluster-gene pairs with >10% expressing: {len(expressing_clusters)}")

    # Top ligand-expressing clusters
    top_ligands = (
        cluster_df[cluster_df["is_ligand"]]
        .groupby("gene")
        .apply(lambda x: x.nlargest(3, "pct_expressing")[["cluster", "pct_expressing"]])
        .reset_index(level=0)
    )

    print("\nTop clusters per ligand (sample):")
    for gene in list(ligand_genes)[:5]:
        gene_data = top_ligands[top_ligands["gene"] == gene]
        if len(gene_data) > 0:
            top = gene_data.iloc[0]
            print(
                f"  {gene}: {top['cluster'][:50]}... ({top['pct_expressing']:.1f}% expressing)"
            )


if __name__ == "__main__":
    main()
