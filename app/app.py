"""Main Dash application for HypoMap Coronal Atlas Viewer."""

import sys
from pathlib import Path
import pandas as pd
import numpy as np
import json

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from dash import Dash
import dash_bootstrap_components as dbc

from app.layouts import create_layout
from app.callbacks import register_callbacks

# Data paths
DATA_DIR = Path(__file__).parent.parent / "data"
CELLS_PATH = DATA_DIR / "processed" / "mouse_abc" / "cells_with_coords.parquet"
CORONAL_REGIONS_PATH = DATA_DIR / "processed" / "mouse_abc" / "coronal_atlas_regions.json"
CLUSTER_ANNOTATIONS_PATH = DATA_DIR / "raw" / "mouse_abc" / "abc_cluster_annotations.csv" / "cluster_annotation-Table 1.csv"
NP_MAP_PATH = DATA_DIR / "generated" / "mouse_common" / "np_map.csv"
NP_BLACKLIST_PATH = DATA_DIR / "generated" / "mouse_common" / "np_system_blacklist.csv"
CLUSTER_LR_PROFILE_PATH = DATA_DIR / "processed" / "mouse_abc" / "cluster_ligand_receptor_profile.parquet"
CLUSTER_NP_EXPR_PATH = DATA_DIR / "processed" / "mouse_abc" / "cluster_np_expression.parquet"
REGION_DESCRIPTIONS_PATH = DATA_DIR / "raw" / "mouse_common.csv"

# Region colors for hypothalamic regions
REGION_COLORS = {
    'ARH': '#DC143C',  # Crimson
    'VMH': '#0066CC',  # Strong blue
    'PVH': '#228B22',  # Forest green
    'DMH': '#8B008B',  # Dark magenta
    'LHA': '#FF6600',  # Vivid orange
    'SCH': '#FFD700',  # Gold
    'MPO': '#8B4513',  # Saddle brown
    'LPO': '#FF1493',  # Deep pink
    'TU': '#696969',   # Dim gray
    'ZI': '#008B8B',   # Dark cyan
    'AHN': '#FF4500',  # Orange red
    'PH': '#4169E1',   # Royal blue
    'SO': '#C71585',   # Medium violet red
    'PVZ': '#32CD32',  # Lime green
    'PVpo': '#00CC66', # Spring green
    'PVi': '#66CDAA',  # Medium aquamarine
    'MBO': '#DAA520',  # Goldenrod
    'AVPV': '#9932CC', # Dark orchid
    'MM': '#B8860B',   # Dark goldenrod
    'MPN': '#DB7093',  # Pale violet red
    'SFO': '#20B2AA',  # Light sea green
    'MEPO': '#778899', # Light slate gray
    'STN': '#6B8E23',  # Olive drab
    'PS': '#BC8F8F',   # Rosy brown
    'PMd': '#CD5C5C',  # Indian red
    'PMv': '#F08080',  # Light coral
    'LM': '#DEB887',   # Burlywood
    'SUMl': '#D2691E', # Chocolate
    'SUMm': '#A0522D', # Sienna
    'TMd': '#8FBC8F',  # Dark sea green
    'TMv': '#3CB371',  # Medium sea green
    'VLPO': '#BA55D3', # Medium orchid
}


def load_cell_data():
    """Load cell data with coordinates and precomputed region assignments."""
    if not CELLS_PATH.exists():
        raise FileNotFoundError(f"Cell data not found at {CELLS_PATH}")
    if not CORONAL_REGIONS_PATH.exists():
        raise FileNotFoundError(
            f"Coronal regions not found at {CORONAL_REGIONS_PATH}\n"
            "Run: uv run python -m src.preprocessing.build_lateralized_regions"
        )

    df = pd.read_parquet(CELLS_PATH)
    print(f"Loaded {len(df):,} cells with coordinates")

    # Filter out HY-unassigned region
    df = df[df['region'] != 'HY-unassigned'].copy()
    print(f"After removing HY-unassigned: {len(df):,} cells")

    # Load precomputed region data
    with open(CORONAL_REGIONS_PATH, 'r') as f:
        region_data = json.load(f)

    slices = region_data['slices']
    print(f"Found {len(slices)} discrete Z slices")

    # Add z_slice and region_display from precomputed data
    df['z_slice'] = df['z'].round(1)
    df['region_display'] = df['cell_id'].map(region_data['cell_regions']).fillna(df['region'])

    # Convert boundaries/centroids keys back to float
    boundaries = {float(k): v for k, v in region_data['boundaries'].items()}
    centroids = {float(k): {r: tuple(c) for r, c in regions.items()}
                 for k, regions in region_data['centroids'].items()}

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
        cluster_label = row.get('cluster_id_label', '')
        nt_type = row.get('nt_type_label', 'NA')
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
        blacklisted = set(blacklist_df['system'].tolist())
        print(f"Loaded blacklist: {len(blacklisted)} systems excluded")

    # Group by system to get ligand and receptor genes
    # Receptor genes can be semicolon-separated for heterodimers (AND logic)
    # e.g., "Calcrl;Ramp2" means both must be expressed
    systems = {}
    for system in df['System'].unique():
        # Skip blacklisted systems
        if system in blacklisted:
            continue

        system_df = df[df['System'] == system]

        # Ligands: simple set of individual genes
        ligand_genes = set()
        for lg in system_df['Ligand_Gene'].dropna().unique():
            for g in lg.split(';'):
                ligand_genes.add(g.strip())

        # Receptors: list of tuples for AND logic
        # Each tuple contains genes that must ALL be expressed
        receptor_complexes = []
        for rg in system_df['Receptor_Gene'].dropna().unique():
            genes = tuple(g.strip() for g in rg.split(';'))
            if genes not in receptor_complexes:
                receptor_complexes.append(genes)

        systems[system] = {
            'ligands': ligand_genes,
            'receptors': receptor_complexes,  # List of tuples for AND logic
        }

    system_names = sorted(systems.keys())
    print(f"Loaded {len(systems)} neuropeptide systems (after blacklist)")

    # Build gene-to-system info mapping for tooltips
    gene_info = {}
    for system in df['System'].unique():
        # Skip blacklisted systems
        if system in blacklisted:
            continue

        system_df = df[df['System'] == system]
        ligand_genes = list(system_df['Ligand_Gene'].dropna().unique())
        receptor_genes = list(system_df['Receptor_Gene'].dropna().unique())
        functional_role = system_df['Functional_Role'].dropna().iloc[0] if len(system_df['Functional_Role'].dropna()) > 0 else ''

        # Map each ligand gene to system info
        for gene in ligand_genes:
            if gene not in gene_info:
                gene_info[gene] = []
            gene_info[gene].append({
                'system': system,
                'role': 'ligand',
                'ligand_genes': ligand_genes,
                'receptor_genes': receptor_genes,
                'functional_role': functional_role,
            })

        # Map each receptor gene to system info
        for gene in receptor_genes:
            if gene not in gene_info:
                gene_info[gene] = []
            gene_info[gene].append({
                'system': system,
                'role': 'receptor',
                'ligand_genes': ligand_genes,
                'receptor_genes': receptor_genes,
                'functional_role': functional_role,
            })

    print(f"Built gene info for {len(gene_info)} genes")

    return systems, system_names, gene_info


def load_cluster_expression():
    """Load cluster-level ligand/receptor expression profiles."""
    if not CLUSTER_LR_PROFILE_PATH.exists():
        print(f"Warning: Cluster expression not found at {CLUSTER_LR_PROFILE_PATH}")
        return {}

    df = pd.read_parquet(CLUSTER_LR_PROFILE_PATH)

    # Create nested dict: cluster -> gene -> expression info
    cluster_expr = {}
    for _, row in df.iterrows():
        cluster = row['cluster']
        gene = row['gene']
        if cluster not in cluster_expr:
            cluster_expr[cluster] = {}
        cluster_expr[cluster][gene] = {
            'mean_expr': float(row['mean_expr']),
            'pct_expressing': float(row['pct_expressing']),
            'is_ligand': bool(row['is_ligand']),
            'is_receptor': bool(row['is_receptor']),
        }

    print(f"Loaded expression profiles for {len(cluster_expr)} clusters")
    return cluster_expr


def load_region_descriptions():
    """Load region descriptions."""
    if not REGION_DESCRIPTIONS_PATH.exists():
        return {}

    df = pd.read_csv(REGION_DESCRIPTIONS_PATH)
    region_info = {}
    for _, row in df.iterrows():
        region_info[row['acronym']] = {
            'full_name': row.get('full_name', ''),
            'description': row.get('description', ''),
        }
    return region_info


def load_cluster_np_expression():
    """Load precomputed cluster-system expression lookup for fast NP mode.

    Returns a dict: system -> cluster -> (max_ligand_expr, max_receptor_expr)
    """
    if not CLUSTER_NP_EXPR_PATH.exists():
        print(f"Warning: Cluster NP expression not found at {CLUSTER_NP_EXPR_PATH}")
        print("Run: uv run python -m src.preprocessing.build_cluster_np_expression")
        return {}

    df = pd.read_parquet(CLUSTER_NP_EXPR_PATH)

    # Build lookup: system -> cluster -> (ligand, receptor)
    lookup = {}
    for _, row in df.iterrows():
        system = row['system']
        cluster = row['cluster']
        if system not in lookup:
            lookup[system] = {}
        lookup[system][cluster] = (
            float(row['max_ligand_expr']),
            float(row['max_receptor_expr']),
        )

    print(f"Loaded NP expression lookup for {len(lookup)} systems")
    return lookup


def create_app():
    """Create and configure the Dash application."""
    print("\n" + "=" * 50)
    print("Loading HypoMap Coronal Atlas data...")
    print("=" * 50 + "\n")

    # Load all data (includes precomputed region boundaries/centroids)
    cells_df, slices, region_boundaries, region_centroids = load_cell_data()
    nt_mapping, nt_types = load_cluster_nt_mapping()
    np_systems, np_system_names, gene_info = load_np_systems()
    cluster_expression = load_cluster_expression()
    cluster_np_expression = load_cluster_np_expression()
    region_descriptions = load_region_descriptions()

    # Collect all data into app_data dict
    app_data = {
        'cells_df': cells_df,
        'slices': slices,
        'nt_mapping': nt_mapping,
        'nt_types': nt_types,
        'np_systems': np_systems,
        'np_system_names': np_system_names,
        'gene_info': gene_info,
        'cluster_expression': cluster_expression,
        'cluster_np_expression': cluster_np_expression,
        'region_boundaries': region_boundaries,
        'region_centroids': region_centroids,
        'region_colors': REGION_COLORS,
        'region_descriptions': region_descriptions,
    }

    # Get hierarchy levels
    cell_type_levels = ['class', 'subclass', 'supertype', 'cluster']

    print(f"\nData summary:")
    print(f"  Cells: {len(cells_df):,}")
    print(f"  Slices: {len(slices)}")
    print(f"  NT types: {len(nt_types)}")
    print(f"  NP systems: {len(np_system_names)}")
    print(f"  Regions: {len(cells_df['region'].unique())}")

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
    app.layout = create_layout(
        cell_type_levels=cell_type_levels,
        nt_types=nt_types,
        np_system_names=np_system_names,
    )

    # Register callbacks
    register_callbacks(app, app_data)

    return app


# Create app instance
app = create_app()

if __name__ == "__main__":
    print("\n" + "=" * 50)
    print("HypoMap Coronal Atlas Viewer")
    print("=" * 50)
    print("\nStarting server at http://localhost:8050")
    print("Press Ctrl+C to stop\n")

    app.run(debug=True, port=8050)
