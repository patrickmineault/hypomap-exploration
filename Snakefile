# Snakefile for HypoMap preprocessing pipeline
#
# Usage:
#   snakemake -n              # Dry run
#   snakemake --cores 1       # Run full pipeline
#   snakemake <target>        # Run specific target
#   snakemake --dag | dot -Tpng > dag.png  # Visualize DAG

# =============================================================================
# Configuration
# =============================================================================

MOUSE_H5AD = "data/raw/mouse_hypomap/d3be7423-d664-4913-89a9-a506cae4c28f.h5ad"
HUMAN_H5AD = "data/raw/human_hypomap/480e89e7-84ad-4fa8-adc3-f7c562a77a78.h5ad"
VISIUM_TAR = "data/raw/human_hypomap/GSE278848_tissue_positions_lists.tar.gz"
VISIUM_H5_DIR = "data/raw/human_hypomap/visium_h5"

# =============================================================================
# Default target: process all datasets
# =============================================================================

rule all:
    input:
        "data/processed/mouse_hypomap/cells_with_coords.parquet",
        "data/processed/human_hypomap/cells_with_coords.parquet",
        "data/processed/mouse_abc/cells_with_coords.parquet",
        "data/processed/mouse_common/gene_descriptions.csv",
        "data/processed/mouse_abc/cluster_ligand_receptor_profile.parquet"

# =============================================================================
# Mouse HypoMap Pipeline
# =============================================================================

rule extract_mouse_hypomap_metadata:
    """Extract cell metadata and gene list from mouse HypoMap h5ad."""
    input:
        MOUSE_H5AD
    output:
        metadata="data/processed/mouse_hypomap/cell_metadata.parquet",
        genes="data/processed/mouse_hypomap/genes.parquet"
    shell:
        "python -m src.datasets.mouse_hypomap"

rule downsample_mouse_hypomap:
    """Downsample mouse cells and assign 3D coordinates (Allen CCF)."""
    input:
        "data/processed/mouse_hypomap/cell_metadata.parquet"
    output:
        "data/processed/mouse_hypomap/cells_with_coords.parquet"
    shell:
        "python -m src.preprocessing.downsample --dataset mouse_hypomap"

# =============================================================================
# Human HypoMap Pipeline
# =============================================================================

rule extract_human_hypomap_metadata:
    """Extract cell metadata and gene list from human hypothalamus h5ad."""
    input:
        HUMAN_H5AD
    output:
        metadata="data/processed/human_hypomap/cell_metadata.parquet",
        genes="data/processed/human_hypomap/genes.parquet"
    shell:
        "python -m src.datasets.human_hypomap"

rule extract_human_visium:
    """Extract Visium spatial data and create 3D coordinates."""
    input:
        tar=VISIUM_TAR,
        h5_dir=VISIUM_H5_DIR
    output:
        "data/processed/human_hypomap/cells_with_coords.parquet"
    shell:
        "python -m src.preprocessing.extract_human_visium"

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

rule download_neuronchat_db:
    """Download NeuronChat mouse interaction database."""
    output:
        "data/raw/mouse_common/interactionDB_mouse.txt"
    shell:
        """
        mkdir -p data/raw/mouse_common
        curl -L -o {output} https://raw.githubusercontent.com/Wei-BioMath/NeuronChatAnalysis2022/main/NeuronChatDB_table/interactionDB_mouse.txt
        """

rule build_ligand_receptor_list:
    """Process NeuronChat database and add manual neuropeptide entries."""
    input:
        "data/raw/mouse_common/interactionDB_mouse.txt"
    output:
        "data/processed/mouse_common/ligand_receptor_mouse.csv"
    shell:
        "python -m src.preprocessing.build_ligand_receptor_list"

rule build_cluster_ligand_receptor_map:
    """Map neuropeptide ligand/receptor expression to ABC clusters."""
    input:
        metadata="data/processed/mouse_abc/cell_metadata.parquet",
        lr_genes="data/processed/mouse_common/ligand_receptor_mouse.csv"
    output:
        expression="data/processed/mouse_abc/neuropeptide_expression.parquet",
        profiles="data/processed/mouse_abc/cluster_ligand_receptor_profile.parquet"
    shell:
        "python -m src.preprocessing.build_cluster_ligand_receptor_map --use-imputed"

# =============================================================================
# Optional: Connectivity data
# =============================================================================

rule extract_connectivity:
    """Download hypothalamus connectivity from Allen Brain Atlas API."""
    output:
        connectivity="data/raw/mouse_hypomap/hypothalamus_connectivity.csv",
        structures="data/raw/mouse_hypomap/hypothalamus_structures.csv"
    shell:
        "python -m src.preprocessing.extract_hypothalamus_connectivity"

# =============================================================================
# Common/Cross-dataset Processing
# =============================================================================

rule build_gene_descriptions:
    """Extract gene names from cluster hierarchies and fetch UniProt descriptions."""
    input:
        hypomap="data/processed/mouse_hypomap/cells_with_coords.parquet",
        abc="data/processed/mouse_abc/cell_metadata.parquet"
    output:
        "data/processed/mouse_common/gene_descriptions.csv"
    shell:
        "python -m src.preprocessing.build_gene_descriptions"

# =============================================================================
# Convenience targets
# =============================================================================

rule mouse_hypomap:
    """Process mouse HypoMap dataset only."""
    input:
        "data/processed/mouse_hypomap/cells_with_coords.parquet"

rule human_hypomap:
    """Process human HypoMap dataset only."""
    input:
        "data/processed/human_hypomap/cells_with_coords.parquet"

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

rule clean:
    """Remove all processed data (keeps raw data)."""
    shell:
        "rm -rf data/processed/mouse_hypomap/*.parquet data/processed/human_hypomap/*.parquet data/processed/mouse_abc/*.parquet"
