# Snakefile for HypoMap preprocessing pipeline
#
# Usage:
#   snakemake -n              # Dry run
#   snakemake --cores 1       # Run full pipeline
#   snakemake <target>        # Run specific target
#   snakemake --dag | dot -Tpng > dag.png  # Visualize DAG

# =============================================================================
# Default target: process all datasets
# =============================================================================

rule all:
    input:
        "data/processed/mouse_abc/cells_with_coords.parquet",
        "data/processed/mouse_common/gene_descriptions.csv",
        "data/processed/mouse_abc/cluster_ligand_receptor_profile.parquet",
        "data/processed/mouse_abc/cluster_np_expression.parquet",
        "data/processed/mouse_abc/coronal_atlas_regions.json",
        "data/processed/mouse_common/hypothalamus_connectivity.csv"

# =============================================================================
# Mouse ABC (Allen Brain Cell Census) Pipeline
# =============================================================================

rule extract_mouse_abc_metadata:
    """Extract cell metadata from ABC MERFISH data, filtered to hypothalamus."""
    output:
        metadata="data/processed/mouse_abc/cell_metadata.parquet",
        genes="data/processed/mouse_abc/genes.parquet"
    shell:
        "python -m src.datasets.mouse_abc"

rule downsample_mouse_abc:
    """Downsample mouse ABC cells (already have CCF coordinates)."""
    input:
        "data/processed/mouse_abc/cell_metadata.parquet"
    output:
        "data/processed/mouse_abc/cells_with_coords.parquet"
    shell:
        "python -m src.preprocessing.downsample --dataset mouse_abc"

rule build_cluster_ligand_receptor_map:
    """Map neuropeptide ligand/receptor expression to ABC clusters."""
    input:
        metadata="data/processed/mouse_abc/cell_metadata.parquet",
        np_map="data/generated/mouse_common/np_map.csv"
    output:
        expression="data/processed/mouse_abc/neuropeptide_expression.parquet",
        profiles="data/processed/mouse_abc/cluster_ligand_receptor_profile.parquet"
    shell:
        "python -m src.preprocessing.build_cluster_ligand_receptor_map --use-imputed"

rule build_lateralized_regions:
    """Precompute lateralized region boundaries for coronal atlas app."""
    input:
        "data/processed/mouse_abc/cells_with_coords.parquet"
    output:
        "data/processed/mouse_abc/coronal_atlas_regions.json"
    shell:
        "python -m src.preprocessing.build_lateralized_regions"

rule build_cluster_np_expression:
    """Precompute cluster-system expression lookup for fast NP mode."""
    input:
        profiles="data/processed/mouse_abc/cluster_ligand_receptor_profile.parquet",
        np_map="data/generated/mouse_common/np_map.csv"
    output:
        "data/processed/mouse_abc/cluster_np_expression.parquet"
    shell:
        "python -m src.preprocessing.build_cluster_np_expression"

rule extract_hypothalamus_connectivity:
    """Download hypothalamus connectivity from Allen Brain Connectivity Atlas."""
    output:
        connectivity="data/processed/mouse_common/hypothalamus_connectivity.csv",
        summary="data/processed/mouse_common/hypothalamus_connectivity_summary.csv",
        structures="data/processed/mouse_common/hypothalamus_structures.csv"
    shell:
        "python -m src.preprocessing.extract_hypothalamus_connectivity"

# =============================================================================
# Convenience targets
# =============================================================================

rule mouse_abc:
    """Process mouse ABC dataset only."""
    input:
        "data/processed/mouse_abc/cells_with_coords.parquet"

rule gene_descriptions:
    """Build gene descriptions from cluster names."""
    input:
        "data/processed/mouse_common/gene_descriptions.csv"

rule ligand_receptor_map:
    """Build cluster ligand-receptor expression profiles."""
    input:
        "data/processed/mouse_abc/cluster_ligand_receptor_profile.parquet"
