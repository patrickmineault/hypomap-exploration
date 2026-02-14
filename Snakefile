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
        "data/processed/mouse_abc/cluster_ligand_receptor_profile.parquet",
        "data/processed/mouse_abc/cluster_np_expression.parquet",
        "data/processed/mouse_abc/coronal_atlas_regions.json",
        "data/processed/mouse_common/hypothalamus_connectivity.csv",
        "data/processed/mouse_abc_subcortical/cells_with_coords.parquet",
        "data/processed/mouse_abc_subcortical/cluster_ligand_receptor_profile.parquet",
        "data/processed/mouse_abc_subcortical/cluster_np_expression.parquet",
        "data/processed/mouse_abc_subcortical/coronal_atlas_regions.json",

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
    """Map neuropeptide and hormone ligand/receptor expression to ABC clusters."""
    input:
        metadata="data/processed/mouse_abc/cell_metadata.parquet",
        np_map="data/generated/mouse_common/np_map.csv",
        hormone_map="data/generated/mouse_common/hormone_map.csv"
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
# Mouse ABC Subcortical (HY + TH + STR + PAL) Pipeline
# =============================================================================

rule extract_mouse_abc_subcortical_metadata:
    """Extract cell metadata from ABC MERFISH data for subcortical divisions."""
    output:
        metadata="data/processed/mouse_abc_subcortical/cell_metadata.parquet",
        genes="data/processed/mouse_abc_subcortical/genes.parquet"
    shell:
        "python -m src.datasets.mouse_abc --divisions HY TH STR PAL --output-dir data/processed/mouse_abc_subcortical"

rule downsample_mouse_abc_subcortical:
    """Downsample mouse ABC subcortical cells."""
    input:
        "data/processed/mouse_abc_subcortical/cell_metadata.parquet"
    output:
        "data/processed/mouse_abc_subcortical/cells_with_coords.parquet"
    shell:
        "python -m src.preprocessing.downsample --dataset mouse_abc_subcortical"

rule build_lateralized_regions_subcortical:
    """Precompute lateralized region boundaries for subcortical coronal atlas."""
    input:
        "data/processed/mouse_abc_subcortical/cells_with_coords.parquet"
    output:
        "data/processed/mouse_abc_subcortical/coronal_atlas_regions.json"
    shell:
        "python -m src.preprocessing.build_lateralized_regions --input data/processed/mouse_abc_subcortical/cells_with_coords.parquet --output data/processed/mouse_abc_subcortical/coronal_atlas_regions.json"

rule build_cluster_ligand_receptor_map_subcortical:
    """Map neuropeptide and hormone ligand/receptor expression for subcortical clusters."""
    input:
        metadata="data/processed/mouse_abc_subcortical/cell_metadata.parquet",
        np_map="data/generated/mouse_common/np_map.csv",
        hormone_map="data/generated/mouse_common/hormone_map.csv"
    output:
        expression="data/processed/mouse_abc_subcortical/neuropeptide_expression.parquet",
        profiles="data/processed/mouse_abc_subcortical/cluster_ligand_receptor_profile.parquet"
    shell:
        "python -m src.preprocessing.build_cluster_ligand_receptor_map --use-imputed --metadata-dir data/processed/mouse_abc_subcortical"

rule build_cluster_np_expression_subcortical:
    """Precompute cluster-system expression lookup for subcortical NP mode."""
    input:
        profiles="data/processed/mouse_abc_subcortical/cluster_ligand_receptor_profile.parquet",
        np_map="data/generated/mouse_common/np_map.csv"
    output:
        "data/processed/mouse_abc_subcortical/cluster_np_expression.parquet"
    shell:
        "python -m src.preprocessing.build_cluster_np_expression --metadata-dir data/processed/mouse_abc_subcortical"

# =============================================================================
# Dash app
# =============================================================================

rule sync_app_data:
    """Copy data files needed by the Dash app into app/data/."""
    input:
        cells="data/processed/mouse_abc/cells_with_coords.parquet",
        regions="data/processed/mouse_abc/coronal_atlas_regions.json",
        lr_profile="data/processed/mouse_abc/cluster_ligand_receptor_profile.parquet",
        np_expr="data/processed/mouse_abc/cluster_np_expression.parquet",
        sub_cells="data/processed/mouse_abc_subcortical/cells_with_coords.parquet",
        sub_regions="data/processed/mouse_abc_subcortical/coronal_atlas_regions.json",
        sub_lr_profile="data/processed/mouse_abc_subcortical/cluster_ligand_receptor_profile.parquet",
        sub_np_expr="data/processed/mouse_abc_subcortical/cluster_np_expression.parquet",
        np_map="data/generated/mouse_common/np_map.csv",
        np_blacklist="data/generated/mouse_common/np_system_blacklist.csv",
        hormone_map="data/generated/mouse_common/hormone_map.csv",
        annotations="data/raw/mouse_abc/abc_cluster_annotations.csv/cluster_annotation-Table 1.csv",
        region_desc="data/generated/mouse_common/region_descriptions.csv",
    output:
        cells="app/data/processed/mouse_abc/cells_with_coords.parquet",
        regions="app/data/processed/mouse_abc/coronal_atlas_regions.json",
        lr_profile="app/data/processed/mouse_abc/cluster_ligand_receptor_profile.parquet",
        np_expr="app/data/processed/mouse_abc/cluster_np_expression.parquet",
        sub_cells="app/data/processed/mouse_abc_subcortical/cells_with_coords.parquet",
        sub_regions="app/data/processed/mouse_abc_subcortical/coronal_atlas_regions.json",
        sub_lr_profile="app/data/processed/mouse_abc_subcortical/cluster_ligand_receptor_profile.parquet",
        sub_np_expr="app/data/processed/mouse_abc_subcortical/cluster_np_expression.parquet",
        np_map="app/data/generated/mouse_common/np_map.csv",
        np_blacklist="app/data/generated/mouse_common/np_system_blacklist.csv",
        hormone_map="app/data/generated/mouse_common/hormone_map.csv",
        annotations="app/data/raw/mouse_abc/abc_cluster_annotations.csv/cluster_annotation-Table 1.csv",
        region_desc="app/data/generated/mouse_common/region_descriptions.csv",
    run:
        import shutil
        from pathlib import Path
        for src, dst in zip(input, output):
            Path(dst).parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(src, dst)

rule app:
    """Sync data and start the Dash app."""
    input:
        rules.sync_app_data.output
    shell:
        "python -m app.app"

# =============================================================================
# Convenience targets
# =============================================================================

rule mouse_abc:
    """Process mouse ABC dataset only."""
    input:
        "data/processed/mouse_abc/cells_with_coords.parquet"

rule ligand_receptor_map:
    """Build cluster ligand-receptor expression profiles."""
    input:
        "data/processed/mouse_abc/cluster_ligand_receptor_profile.parquet"
