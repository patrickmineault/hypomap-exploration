"""Build cluster ligand/receptor expression profiles for Langlieb et al. data.

Uses snRNA-seq cluster summary files (avg expression + nonzero counts) instead
of per-cell expression. This avoids needing the full counts h5ad.

Input files (from Google Drive singlenuclei_data folder):
  - Single_Nuc_Cluster_Avg_Expression.csv.gz  (clusters × genes, mean expression)
  - Single_Nuc_Cluster_NonZero_Counts.csv.gz   (clusters × genes, nonzero cell count)
  - CellType_Metadata.tsv                       (cluster metadata with num_cells_postQC)

Usage:
    python -m hypomap.preprocessing.build_cluster_ligand_receptor_map_langlieb
"""

import argparse
import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

DATA_DIR = Path(__file__).parent.parent.parent / "data"
SNRNASEQ_DIR = DATA_DIR / "raw" / "mouse_langlieb" / "singlenuclei_data"
OUTPUT_DIR = DATA_DIR / "processed" / "mouse_langlieb"
NP_MAP_PATH = DATA_DIR / "generated" / "mouse_common" / "np_map.csv"
HORMONE_MAP_PATH = DATA_DIR / "generated" / "mouse_common" / "hormone_map.csv"
ABC_PROFILE_PATH = DATA_DIR / "processed" / "mouse_abc" / "cluster_ligand_receptor_profile.parquet"


def strip_key(name: str) -> str:
    """Strip numeric key prefix from Langlieb cluster names.

    E.g. '300-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-0-1-0-0=Ex_Rorb_Foxp2' -> 'Ex_Rorb_Foxp2'
    """
    return re.sub(r".*=", "", name)


def strip_gene_key(name: str) -> str:
    """Strip Ensembl ID suffix from gene column names.

    E.g. 'Xkr4=ENSMUSG00000051951' -> 'Xkr4'
    """
    return name.split("=")[0] if "=" in name else name


def get_ligand_receptor_genes() -> tuple[list[str], set[str], set[str]]:
    """Load ligand/receptor genes from np_map.csv and hormone_map.csv."""
    ligand_genes = set()
    receptor_genes = set()

    np_map = pd.read_csv(NP_MAP_PATH)
    for lg in np_map["Ligand_Gene"].dropna().unique():
        for g in lg.split(";"):
            ligand_genes.add(g.strip())
    for rg in np_map["Receptor_Gene"].dropna().unique():
        for g in rg.split(";"):
            receptor_genes.add(g.strip())

    if HORMONE_MAP_PATH.exists():
        hormone_map = pd.read_csv(HORMONE_MAP_PATH)
        for lg in hormone_map["Ligand_Gene"].dropna().unique():
            for g in lg.split(";"):
                ligand_genes.add(g.strip())
        for rg in hormone_map["Receptor_Gene"].dropna().unique():
            for g in rg.split(";"):
                receptor_genes.add(g.strip())

    all_genes = sorted(ligand_genes | receptor_genes)
    print(f"  Ligands: {len(ligand_genes)}, Receptors: {len(receptor_genes)}, Total: {len(all_genes)}")
    return all_genes, ligand_genes, receptor_genes


def compute_gene_scale_factors(profile_df: pd.DataFrame, abc_profile_path: Path) -> pd.DataFrame:
    """Normalize Langlieb mean_expr gene-by-gene using ABC as reference.

    For each gene present in both datasets, compute scale_factor = lang_95 / abc_95
    (95th percentile ratio). Divide Langlieb mean_expr by this factor so that
    expression scales match ABC MERFISH. Genes missing from ABC use the median
    scale factor as fallback.
    """
    if not abc_profile_path.exists():
        warnings.warn(f"ABC profile not found at {abc_profile_path}, skipping normalization")
        return profile_df

    abc_profile = pd.read_parquet(abc_profile_path)
    print(f"\nNormalizing to ABC scale using {abc_profile_path}...")

    common_genes = sorted(
        set(profile_df["gene"].unique()) & set(abc_profile["gene"].unique())
    )
    print(f"  Common genes: {len(common_genes)}")

    gene_scales = {}
    for gene in common_genes:
        lang_vals = profile_df.loc[profile_df["gene"] == gene, "mean_expr"].values.astype(np.float64)
        lang_vals = lang_vals[np.isfinite(lang_vals)]

        abc_vals = abc_profile.loc[abc_profile["gene"] == gene, "mean_expr"].values.astype(np.float64)
        abc_vals = abc_vals[np.isfinite(abc_vals)]

        if len(lang_vals) < 10 or len(abc_vals) < 10:
            continue

        lang_95 = np.percentile(lang_vals, 95)
        abc_95 = np.percentile(abc_vals, 95)

        if lang_95 < 0.01 or abc_95 < 0.01:
            continue

        gene_scales[gene] = lang_95 / abc_95

    if not gene_scales:
        warnings.warn("No valid gene scale factors computed, skipping normalization")
        return profile_df

    median_scale = np.median(list(gene_scales.values()))
    print(f"  Computed scale factors for {len(gene_scales)} genes")
    print(f"  Median scale factor: {median_scale:.3f}")

    scale_series = profile_df["gene"].map(gene_scales).fillna(median_scale)
    profile_df = profile_df.copy()
    profile_df["mean_expr"] = profile_df["mean_expr"] / scale_series

    return profile_df


def main():
    print("=== Building Langlieb Cluster Ligand-Receptor Expression Map ===\n")

    # 1. Load gene lists
    print("Loading ligand/receptor gene lists...")
    all_genes, ligand_genes, receptor_genes = get_ligand_receptor_genes()

    # 2. Load cluster metadata for num_cells_postQC
    meta_path = SNRNASEQ_DIR / "CellType_Metadata.tsv"
    print(f"Loading cluster metadata from {meta_path}...")
    meta = pd.read_csv(meta_path, sep="\t")
    # Build mapping: cluster_name -> num_cells_postQC
    meta_clean = meta[meta["Annotation"] != "REMOVE"].copy()
    num_cells = dict(zip(meta_clean["Annotation"], meta_clean["num_cells_postQC"]))
    print(f"  {len(num_cells)} clusters with cell counts")

    # 3. Load average expression (clusters × genes)
    avg_path = SNRNASEQ_DIR / "Single_Nuc_Cluster_Avg_Expression.csv.gz"
    print(f"Loading average expression from {avg_path}...")
    avg_df = pd.read_csv(avg_path, index_col=0)
    # Strip keys from index and columns
    avg_df.index = [strip_key(idx) for idx in avg_df.index]
    avg_df.columns = [strip_gene_key(col) for col in avg_df.columns]
    print(f"  Shape: {avg_df.shape}")

    # 4. Load nonzero counts (clusters × genes)
    nz_path = SNRNASEQ_DIR / "Single_Nuc_Cluster_NonZero_Counts.csv.gz"
    print(f"Loading nonzero counts from {nz_path}...")
    nz_df = pd.read_csv(nz_path, index_col=0)
    nz_df.index = [strip_key(idx) for idx in nz_df.index]
    nz_df.columns = [strip_gene_key(col) for col in nz_df.columns]
    print(f"  Shape: {nz_df.shape}")

    # 5. Find which target genes are available
    available_genes = set(avg_df.columns)
    genes_found = [g for g in all_genes if g in available_genes]
    genes_missing = [g for g in all_genes if g not in available_genes]
    print(f"\nGenes requested: {len(all_genes)}, found: {len(genes_found)}, missing: {len(genes_missing)}")
    if genes_missing:
        print(f"  Missing (first 10): {genes_missing[:10]}")

    # 6. Build cluster profiles
    print("\nBuilding cluster-level profiles...")
    rows = []
    clusters = avg_df.index.tolist()

    for cluster in clusters:
        n_cells = num_cells.get(cluster, 0)
        if n_cells == 0:
            continue

        for gene in genes_found:
            mean_expr = float(avg_df.loc[cluster, gene])
            nonzero = float(nz_df.loc[cluster, gene]) if gene in nz_df.columns else 0
            pct_expressing = nonzero / n_cells * 100

            rows.append({
                "cluster": cluster,
                "gene": gene,
                "mean_expr": mean_expr,
                "pct_expressing": pct_expressing,
                "n_cells": int(n_cells),
                "n_expressing": int(nonzero),
                "is_ligand": gene in ligand_genes,
                "is_receptor": gene in receptor_genes,
            })

    profile_df = pd.DataFrame(rows)
    print(f"  Generated {len(profile_df)} cluster-gene pairs")

    # 7. Normalize mean_expr to ABC scale (gene-by-gene)
    parser = argparse.ArgumentParser()
    parser.add_argument("--abc-profile", type=Path, default=ABC_PROFILE_PATH,
                        help="Path to ABC cluster_ligand_receptor_profile.parquet for normalization")
    args = parser.parse_args()

    profile_df = compute_gene_scale_factors(profile_df, args.abc_profile)

    # 8. Save
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    output_path = OUTPUT_DIR / "cluster_ligand_receptor_profile.parquet"
    profile_df.to_parquet(output_path, index=False)
    print(f"\nSaved to {output_path}")

    # Summary
    expressing = profile_df[profile_df["pct_expressing"] > 10]
    print(f"Cluster-gene pairs with >10% expressing: {len(expressing)}")


if __name__ == "__main__":
    main()
