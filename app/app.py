"""Main Dash application for HypoMap Coronal Atlas Viewer."""

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

import dash_bootstrap_components as dbc
from dash import Dash

from app.callbacks import register_callbacks
from app.layouts import create_layout

# Feature flags
ENABLE_QUANTILE_TOGGLE = True  # Expression metric toggle (log2CPM vs quantile)
ENABLE_REGION_HIGHLIGHT = True  # Region click-to-select highlight
ENABLE_HORMONE_MODE = True  # Hormone receptor visualization mode

# Data paths (use app/data symlink for Plotly Cloud compatibility)
DATA_DIR = Path(__file__).parent / "data"
CLUSTER_ANNOTATIONS_PATH = (
    DATA_DIR
    / "raw"
    / "mouse_abc"
    / "abc_cluster_annotations.csv"
    / "cluster_annotation-Table 1.csv"
)
NP_MAP_PATH = DATA_DIR / "generated" / "mouse_common" / "np_map.csv"
NP_BLACKLIST_PATH = DATA_DIR / "generated" / "mouse_common" / "np_system_blacklist.csv"
REGION_DESCRIPTIONS_PATH = (
    DATA_DIR / "generated" / "mouse_common" / "region_descriptions.csv"
)
HORMONE_MAP_PATH = DATA_DIR / "generated" / "mouse_common" / "hormone_map.csv"

# Per-dataset paths
DATASET_PATHS = {
    "mouse_abc": {
        "cells": DATA_DIR / "processed" / "mouse_abc" / "cells_with_coords.parquet",
        "regions": DATA_DIR / "processed" / "mouse_abc" / "coronal_atlas_regions.json",
        "lr_profile": DATA_DIR
        / "processed"
        / "mouse_abc"
        / "cluster_ligand_receptor_profile.parquet",
        "np_expr": DATA_DIR
        / "processed"
        / "mouse_abc"
        / "cluster_np_expression.parquet",
    },
    "mouse_abc_extended": {
        "cells": DATA_DIR
        / "processed"
        / "mouse_abc_extended"
        / "cells_with_coords.parquet",
        "regions": DATA_DIR
        / "processed"
        / "mouse_abc_extended"
        / "coronal_atlas_regions.json",
        "lr_profile": DATA_DIR
        / "processed"
        / "mouse_abc_extended"
        / "cluster_ligand_receptor_profile.parquet",
        "np_expr": DATA_DIR
        / "processed"
        / "mouse_abc_extended"
        / "cluster_np_expression.parquet",
    },
}


def load_ccf_region_colors():
    """Load official Allen CCF colors for brain regions.

    Merges the parcellation color and acronym tables from the ABC atlas cache
    to produce a mapping from region acronym to hex color.
    """
    cache_base = (
        DATA_DIR
        / "raw"
        / "abc_atlas_cache"
        / "metadata"
        / "Allen-CCF-2020"
        / "20230630"
        / "views"
    )
    color_path = cache_base / "parcellation_to_parcellation_term_membership_color.csv"
    acronym_path = (
        cache_base / "parcellation_to_parcellation_term_membership_acronym.csv"
    )

    if not color_path.exists() or not acronym_path.exists():
        print(
            f"Warning: CCF color tables not found at {cache_base}, using empty color map"
        )
        return {}

    colors = pd.read_csv(color_path)[["parcellation_index", "structure_color"]]
    acronyms = pd.read_csv(acronym_path)[["parcellation_index", "structure"]]
    merged = pd.merge(colors, acronyms, on="parcellation_index")
    mapping = (
        merged.drop_duplicates(subset="structure")
        .set_index("structure")["structure_color"]
        .to_dict()
    )
    # Remove the 'unassigned' entry if present
    mapping.pop("unassigned", None)
    print(f"Loaded {len(mapping)} CCF region colors")
    return mapping


# Load CCF colors once at module level — used for all datasets
CCF_REGION_COLORS = load_ccf_region_colors()

# Map dataset name to region colors (same CCF colors for all datasets)
DATASET_REGION_COLORS = {
    "mouse_abc": CCF_REGION_COLORS,
    "mouse_abc_extended": CCF_REGION_COLORS,
}


def load_cell_data(dataset_name="mouse_abc"):
    """Load cell data with coordinates and precomputed region assignments."""
    paths = DATASET_PATHS[dataset_name]
    cells_path = paths["cells"]
    regions_path = paths["regions"]

    if not cells_path.exists():
        raise FileNotFoundError(f"Cell data not found at {cells_path}")
    if not regions_path.exists():
        raise FileNotFoundError(
            f"Coronal regions not found at {regions_path}\n"
            "Run: uv run python -m hypomap.preprocessing.build_lateralized_regions"
        )

    df = pd.read_parquet(cells_path)
    print(f"[{dataset_name}] Loaded {len(df):,} cells with coordinates")

    # Filter out *-unassigned regions (e.g. HY-unassigned, TH-unassigned)
    unassigned_mask = df["region"].str.endswith("-unassigned")
    df = df[~unassigned_mask].copy()
    print(f"[{dataset_name}] After removing *-unassigned: {len(df):,} cells")

    # Load precomputed region data
    with open(regions_path, "r") as f:
        region_data = json.load(f)

    slices = region_data["slices"]
    print(f"[{dataset_name}] Found {len(slices)} discrete Z slices")

    # Add z_slice and region_display from precomputed data
    df["z_slice"] = df["z"].round(1)
    df["region_display"] = (
        df["cell_id"].map(region_data["cell_regions"]).fillna(df["region"])
    )

    # Convert boundaries/centroids keys back to float
    boundaries = {float(k): v for k, v in region_data["boundaries"].items()}
    centroids = {
        float(k): {r: tuple(c) for r, c in regions.items()}
        for k, regions in region_data["centroids"].items()
    }

    return df, slices, boundaries, centroids


def load_cluster_nt_mapping():
    """Load cluster to neurotransmitter type mapping."""
    if not CLUSTER_ANNOTATIONS_PATH.exists():
        print(f"Warning: Cluster annotations not found at {CLUSTER_ANNOTATIONS_PATH}")
        return {}

    df = pd.read_csv(CLUSTER_ANNOTATIONS_PATH)

    # Create mapping from cluster_id_label to nt_type_label
    # The cluster_id_label format matches what's in our cell data
    nt_mapping = {}
    for _, row in df.iterrows():
        cluster_label = row.get("cluster_id_label", "")
        nt_type = row.get("nt_type_label", "NA")
        if cluster_label and pd.notna(nt_type):
            nt_mapping[cluster_label] = nt_type

    print(f"Loaded NT mapping for {len(nt_mapping)} clusters")
    nt_types = sorted(set(nt_mapping.values()))
    print(f"NT types: {nt_types}")

    return nt_mapping, nt_types


def load_np_systems():
    """Load neuropeptide system definitions."""
    if not NP_MAP_PATH.exists():
        print(f"Warning: NP map not found at {NP_MAP_PATH}")
        return {}, [], {}

    df = pd.read_csv(NP_MAP_PATH)

    # Load blacklist of systems to exclude (both ligand AND receptor dead)
    blacklisted = set()
    if NP_BLACKLIST_PATH.exists():
        blacklist_df = pd.read_csv(NP_BLACKLIST_PATH)
        blacklisted = set(blacklist_df["system"].tolist())
        print(f"Loaded blacklist: {len(blacklisted)} systems excluded")

    # Group by system to get ligand and receptor genes
    # Receptor genes can be semicolon-separated for heterodimers (AND logic)
    # e.g., "Calcrl;Ramp2" means both must be expressed
    systems = {}
    for system in df["System"].unique():
        # Skip blacklisted systems
        if system in blacklisted:
            continue

        system_df = df[df["System"] == system]

        # Ligands: simple set of individual genes
        ligand_genes = set()
        for lg in system_df["Ligand_Gene"].dropna().unique():
            for g in lg.split(";"):
                ligand_genes.add(g.strip())

        # Receptors: list of tuples for AND logic
        # Each tuple contains genes that must ALL be expressed
        receptor_complexes = []
        for rg in system_df["Receptor_Gene"].dropna().unique():
            genes = tuple(g.strip() for g in rg.split(";"))
            if genes not in receptor_complexes:
                receptor_complexes.append(genes)

        systems[system] = {
            "ligands": ligand_genes,
            "receptors": receptor_complexes,  # List of tuples for AND logic
        }

    system_names = sorted(systems.keys())
    print(f"Loaded {len(systems)} neuropeptide systems (after blacklist)")

    # Build gene-to-system info mapping for tooltips
    gene_info = {}
    for system in df["System"].unique():
        # Skip blacklisted systems
        if system in blacklisted:
            continue

        system_df = df[df["System"] == system]
        ligand_genes = list(system_df["Ligand_Gene"].dropna().unique())
        receptor_genes = list(system_df["Receptor_Gene"].dropna().unique())
        functional_role = (
            system_df["Functional_Role"].dropna().iloc[0]
            if len(system_df["Functional_Role"].dropna()) > 0
            else ""
        )

        # Map each ligand gene to system info
        for gene in ligand_genes:
            if gene not in gene_info:
                gene_info[gene] = []
            gene_info[gene].append(
                {
                    "system": system,
                    "role": "ligand",
                    "ligand_genes": ligand_genes,
                    "receptor_genes": receptor_genes,
                    "functional_role": functional_role,
                }
            )

        # Map each receptor gene to system info
        for gene in receptor_genes:
            if gene not in gene_info:
                gene_info[gene] = []
            gene_info[gene].append(
                {
                    "system": system,
                    "role": "receptor",
                    "ligand_genes": ligand_genes,
                    "receptor_genes": receptor_genes,
                    "functional_role": functional_role,
                }
            )

    print(f"Built gene info for {len(gene_info)} genes")

    return systems, system_names, gene_info


def load_hormone_systems():
    """Load hormone system definitions (receptor genes only)."""
    if not HORMONE_MAP_PATH.exists():
        print(f"Warning: Hormone map not found at {HORMONE_MAP_PATH}")
        return {}, []

    df = pd.read_csv(HORMONE_MAP_PATH)

    # Group by System to get unique receptor genes per hormone
    systems = {}
    for system in df["System"].unique():
        system_df = df[df["System"] == system]
        receptor_genes = set(system_df["Receptor_Gene"].dropna().unique())
        functional_role = (
            system_df["Functional_Role"].dropna().iloc[0]
            if len(system_df["Functional_Role"].dropna()) > 0
            else ""
        )
        systems[system] = {
            "receptors": receptor_genes,
            "functional_role": functional_role,
        }

    hormone_names = sorted(systems.keys())
    print(f"Loaded {len(systems)} hormone systems")

    return systems, hormone_names


def load_cluster_expression(dataset_name="mouse_abc"):
    """Load cluster-level ligand/receptor expression profiles."""
    lr_profile_path = DATASET_PATHS[dataset_name]["lr_profile"]
    if not lr_profile_path.exists():
        print(f"Warning: Cluster expression not found at {lr_profile_path}")
        return {}

    df = pd.read_parquet(lr_profile_path)

    # Create nested dict: cluster -> gene -> expression info
    cluster_expr = {}
    for _, row in df.iterrows():
        cluster = row["cluster"]
        gene = row["gene"]
        if cluster not in cluster_expr:
            cluster_expr[cluster] = {}
        cluster_expr[cluster][gene] = {
            "mean_expr": float(row["mean_expr"]),
            "pct_expressing": float(row["pct_expressing"]),
            "is_ligand": bool(row["is_ligand"]),
            "is_receptor": bool(row["is_receptor"]),
        }

    print(f"Loaded expression profiles for {len(cluster_expr)} clusters")
    return cluster_expr


def compute_gene_quantiles(cluster_expression):
    """Compute quantile ranks for each gene's expression across all clusters.

    Returns dict: gene -> cluster -> quantile (0-100)
    """
    from scipy import stats

    # First, collect all expression values per gene
    gene_values = {}  # gene -> list of (cluster, expr)
    for cluster, genes in cluster_expression.items():
        for gene, info in genes.items():
            if gene not in gene_values:
                gene_values[gene] = []
            gene_values[gene].append((cluster, info["mean_expr"]))

    # Compute quantile rank for each cluster's expression
    gene_quantiles = {}
    for gene, values in gene_values.items():
        if len(values) < 2:
            continue

        clusters = [v[0] for v in values]
        exprs = np.array([v[1] for v in values])

        # Compute percentile ranks (0-100)
        ranks = stats.rankdata(exprs, method="average")
        percentiles = (ranks - 1) / (len(ranks) - 1) * 100

        gene_quantiles[gene] = {
            cluster: float(pct) for cluster, pct in zip(clusters, percentiles)
        }

    print(f"Computed quantiles for {len(gene_quantiles)} genes")
    return gene_quantiles


def load_region_descriptions():
    """Load region descriptions."""
    if not REGION_DESCRIPTIONS_PATH.exists():
        return {}

    df = pd.read_csv(REGION_DESCRIPTIONS_PATH)
    region_info = {}
    for _, row in df.iterrows():
        region_info[row["acronym"]] = {
            "full_name": row.get("full_name", ""),
            "division": row.get("division", ""),
            "description": row.get("description", ""),
        }
    return region_info


def load_cluster_np_expression(dataset_name="mouse_abc"):
    """Load precomputed cluster-system expression lookup for fast NP mode.

    Returns a dict: system -> cluster -> (max_ligand_expr, max_receptor_expr)
    """
    np_expr_path = DATASET_PATHS[dataset_name]["np_expr"]
    if not np_expr_path.exists():
        print(f"Warning: Cluster NP expression not found at {np_expr_path}")
        print("Run: uv run python -m hypomap.preprocessing.build_cluster_np_expression")
        return {}

    df = pd.read_parquet(np_expr_path)

    # Build lookup: system -> cluster -> (ligand, receptor)
    lookup = {}
    for _, row in df.iterrows():
        system = row["system"]
        cluster = row["cluster"]
        if system not in lookup:
            lookup[system] = {}
        lookup[system][cluster] = (
            float(row["max_ligand_expr"]),
            float(row["max_receptor_expr"]),
        )

    print(f"Loaded NP expression lookup for {len(lookup)} systems")
    return lookup


def load_dataset_bundle(dataset_name):
    """Load all data for a single dataset variant."""
    cells_df, slices, region_boundaries, region_centroids = load_cell_data(dataset_name)
    cluster_expression = load_cluster_expression(dataset_name)
    cluster_np_expression = load_cluster_np_expression(dataset_name)
    gene_quantiles = compute_gene_quantiles(cluster_expression)
    region_colors = DATASET_REGION_COLORS.get(dataset_name, CCF_REGION_COLORS)

    # Auto-assign colors for any regions not in the color map
    import plotly.express as px

    extended_palette = (
        px.colors.qualitative.Dark24
        + px.colors.qualitative.Light24
        + px.colors.qualitative.Alphabet
    )
    all_regions = set(cells_df["region"].unique())
    missing_regions = all_regions - set(region_colors.keys())
    if missing_regions:
        existing_count = len(region_colors)
        for i, region in enumerate(sorted(missing_regions)):
            region_colors[region] = extended_palette[
                (existing_count + i) % len(extended_palette)
            ]
        print(
            f"[{dataset_name}] Auto-assigned colors for {len(missing_regions)} new regions"
        )

    return {
        "cells_df": cells_df,
        "slices": slices,
        "region_boundaries": region_boundaries,
        "region_centroids": region_centroids,
        "cluster_expression": cluster_expression,
        "cluster_np_expression": cluster_np_expression,
        "gene_quantiles": gene_quantiles,
        "region_colors": region_colors,
    }


def create_app():
    """Create and configure the Dash application."""
    print("\n" + "=" * 50)
    print("Loading HypoMap Coronal Atlas data...")
    print("=" * 50 + "\n")

    # Load shared data (independent of dataset)
    nt_mapping, nt_types = load_cluster_nt_mapping()
    np_systems, np_system_names, gene_info = load_np_systems()
    hormone_systems, hormone_names = (
        load_hormone_systems() if ENABLE_HORMONE_MODE else ({}, [])
    )
    region_descriptions = load_region_descriptions()

    # Load per-dataset data bundles
    datasets = {}
    for ds_name in DATASET_PATHS:
        try:
            datasets[ds_name] = load_dataset_bundle(ds_name)
            print(
                f"  Loaded {ds_name}: {len(datasets[ds_name]['cells_df']):,} cells, {len(datasets[ds_name]['slices'])} slices"
            )
        except FileNotFoundError as e:
            print(f"  Skipping {ds_name}: {e}")

    if not datasets:
        raise FileNotFoundError("No datasets found. Run preprocessing first.")

    # Default dataset is extended if available
    default_dataset = (
        "mouse_abc_extended"
        if "mouse_abc_extended" in datasets
        else list(datasets.keys())[0]
    )

    # Collect all data into app_data dict
    app_data = {
        "datasets": datasets,
        "default_dataset": default_dataset,
        # Shared data
        "nt_mapping": nt_mapping,
        "nt_types": nt_types,
        "np_systems": np_systems,
        "np_system_names": np_system_names,
        "hormone_systems": hormone_systems,
        "hormone_names": hormone_names,
        "gene_info": gene_info,
        "region_descriptions": region_descriptions,
    }

    # Get hierarchy levels
    cell_type_levels = ["class", "subclass", "supertype", "cluster"]

    # Get region list from default dataset (sorted alphabetically)
    default_data = datasets[default_dataset]
    region_list = sorted(default_data["cells_df"]["region"].unique())

    # Dataset names for UI toggle
    dataset_names = list(datasets.keys())

    print(f"\nData summary:")
    for ds_name, ds_data in datasets.items():
        print(
            f"  [{ds_name}] Cells: {len(ds_data['cells_df']):,}, Slices: {len(ds_data['slices'])}"
        )
    print(f"  NT types: {len(nt_types)}")
    print(f"  NP systems: {len(np_system_names)}")
    print(f"  Hormone systems: {len(hormone_names)}")
    print(f"  Regions (default): {len(region_list)}")

    # Create Dash app
    assets_path = Path(__file__).parent / "assets"
    print(f"  Assets path: {assets_path} (exists: {assets_path.exists()})")
    app = Dash(
        __name__,
        assets_folder=str(assets_path),
        external_stylesheets=[dbc.themes.BOOTSTRAP],
        title="HypoMap Coronal Atlas",
        suppress_callback_exceptions=True,
    )

    # Set layout
    default_n_cells = len(default_data["cells_df"])
    default_subsample = max(5, min(30, int(60000 / default_n_cells * 100)))
    app.layout = create_layout(
        cell_type_levels=cell_type_levels,
        nt_types=nt_types,
        np_system_names=np_system_names,
        hormone_names=hormone_names if ENABLE_HORMONE_MODE else None,
        region_list=region_list if ENABLE_REGION_HIGHLIGHT else None,
        region_descriptions=region_descriptions if ENABLE_REGION_HIGHLIGHT else None,
        enable_region_highlight=ENABLE_REGION_HIGHLIGHT,
        enable_quantile_toggle=ENABLE_QUANTILE_TOGGLE,
        dataset_names=dataset_names,
        default_dataset=default_dataset,
        default_slices=default_data["slices"],
        default_subsample=default_subsample,
    )

    # Register callbacks
    register_callbacks(
        app,
        app_data,
        enable_region_highlight=ENABLE_REGION_HIGHLIGHT,
        enable_quantile_toggle=ENABLE_QUANTILE_TOGGLE,
    )

    return app


# Create app instance at module level (used by both __main__ and WSGI imports)
app = create_app()


def main():
    """Entry point for the fireitup console script."""
    print("\n" + "=" * 50)
    print("HypoMap Coronal Atlas Viewer")
    print("=" * 50)
    print("\nStarting server at http://localhost:8050")
    print("Press Ctrl+C to stop\n")

    app.run(debug=True, port=8050, use_reloader=False)


if __name__ == "__main__":
    main()
