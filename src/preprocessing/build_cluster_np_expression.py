"""Precompute cluster-level neuropeptide system expression for fast app lookups.

For each (cluster, NP system) pair, precomputes:
- max_ligand_expr: Maximum expression across all ligand genes
- max_receptor_expr: Maximum expression across receptor complexes (with AND logic for heterodimers)

This reduces per-request lookup from O(clusters × genes) to O(clusters).

Usage:
    python -m src.preprocessing.build_cluster_np_expression
    python -m src.preprocessing.build_cluster_np_expression --metadata-dir data/processed/mouse_abc_subcortical
"""

import argparse
from pathlib import Path
import pandas as pd
import numpy as np

DATA_DIR = Path(__file__).parent.parent.parent / "data"
NP_MAP_PATH = DATA_DIR / "generated" / "mouse_common" / "np_map.csv"
DEFAULT_METADATA_DIR = DATA_DIR / "processed" / "mouse_abc"


def load_np_systems():
    """Load NP system definitions with ligand/receptor genes."""
    df = pd.read_csv(NP_MAP_PATH)

    systems = {}
    for system in df['System'].unique():
        system_df = df[df['System'] == system]

        # Ligands: simple set of genes
        ligand_genes = set()
        for lg in system_df['Ligand_Gene'].dropna().unique():
            for g in lg.split(';'):
                ligand_genes.add(g.strip())

        # Receptors: list of tuples for AND logic (heterodimer complexes)
        receptor_complexes = []
        for rg in system_df['Receptor_Gene'].dropna().unique():
            genes = tuple(g.strip() for g in rg.split(';'))
            if genes not in receptor_complexes:
                receptor_complexes.append(genes)

        systems[system] = {
            'ligands': list(ligand_genes),
            'receptors': receptor_complexes,
        }

    return systems


def load_cluster_expression(profile_path):
    """Load cluster-level gene expression data."""
    df = pd.read_parquet(profile_path)

    # Pivot to cluster -> gene -> mean_expr
    cluster_expr = {}
    for _, row in df.iterrows():
        cluster = row['cluster']
        gene = row['gene']
        if cluster not in cluster_expr:
            cluster_expr[cluster] = {}
        cluster_expr[cluster][gene] = float(row['mean_expr'])

    return cluster_expr


def compute_cluster_system_expression(clusters, systems, cluster_expr):
    """Compute max ligand/receptor expression for each cluster-system pair."""
    rows = []

    for cluster in clusters:
        expr = cluster_expr.get(cluster, {})

        for system_name, system_info in systems.items():
            ligand_genes = system_info['ligands']
            receptor_complexes = system_info['receptors']

            # Max ligand expression (any ligand gene above 0)
            max_ligand = 0.0
            for gene in ligand_genes:
                if gene in expr:
                    max_ligand = max(max_ligand, expr[gene])

            # Max receptor expression with AND logic for heterodimers
            # Each receptor_complex is a tuple of genes that must ALL be expressed
            max_receptor = 0.0
            for receptor_complex in receptor_complexes:
                complex_exprs = []
                all_present = True
                for gene in receptor_complex:
                    if gene in expr and expr[gene] > 0:
                        complex_exprs.append(expr[gene])
                    else:
                        all_present = False
                        break
                # If all genes in complex are expressed, use minimum as complex expression
                if all_present and complex_exprs:
                    max_receptor = max(max_receptor, min(complex_exprs))

            rows.append({
                'cluster': cluster,
                'system': system_name,
                'max_ligand_expr': max_ligand,
                'max_receptor_expr': max_receptor,
            })

    return pd.DataFrame(rows)


def main(metadata_dir=None):
    if metadata_dir is None:
        metadata_dir = DEFAULT_METADATA_DIR

    profile_path = metadata_dir / "cluster_ligand_receptor_profile.parquet"
    output_path = metadata_dir / "cluster_np_expression.parquet"

    print("=== Building Cluster NP Expression Lookup ===\n")

    # Load data
    print("Loading NP systems...")
    systems = load_np_systems()
    print(f"  {len(systems)} NP systems")

    print(f"Loading cluster expression from {profile_path}...")
    cluster_expr = load_cluster_expression(profile_path)
    clusters = sorted(cluster_expr.keys())
    print(f"  {len(clusters)} clusters")

    # Compute lookup table
    print("\nComputing cluster-system expression matrix...")
    df = compute_cluster_system_expression(clusters, systems, cluster_expr)
    print(f"  Generated {len(df)} rows ({len(clusters)} clusters × {len(systems)} systems)")

    # Summary stats
    expressing_ligand = df[df['max_ligand_expr'] > 0]
    expressing_receptor = df[df['max_receptor_expr'] > 0]
    print(f"  Cluster-system pairs with ligand expression: {len(expressing_ligand)}")
    print(f"  Cluster-system pairs with receptor expression: {len(expressing_receptor)}")

    # Save
    df.to_parquet(output_path, index=False)
    print(f"\nSaved to {output_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Precompute cluster-level NP system expression for fast app lookups."
    )
    parser.add_argument(
        "--metadata-dir",
        type=Path,
        default=None,
        help="Directory with cluster_ligand_receptor_profile.parquet (default: data/processed/mouse_abc)",
    )
    args = parser.parse_args()
    main(metadata_dir=args.metadata_dir)
