"""Main Dash application for HypoMap 3D Viewer."""

import sys
from pathlib import Path
import pandas as pd

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from dash import Dash
import dash_bootstrap_components as dbc

from app.layouts import create_layout
from app.callbacks import register_callbacks

# Import datasets module
from src.datasets import (
    get_processed_datasets,
    load_dataset,
    get_marker_genes,
    get_key_receptors,
    get_region_colors,
    detect_cell_type_levels,
    detect_region_column,
)

# Path to gene descriptions
GENE_DESCRIPTIONS_PATH = Path(__file__).parent.parent / "data" / "processed" / "mouse_common" / "gene_descriptions.csv"

# Path to region descriptions
REGION_DESCRIPTIONS_PATH = Path(__file__).parent.parent / "data" / "raw" / "mouse_common.csv"

# Path to ligand-receptor data
CLUSTER_LR_PROFILE_PATH = Path(__file__).parent.parent / "data" / "processed" / "mouse_abc" / "cluster_ligand_receptor_profile.parquet"
LR_PAIRS_PATH = Path(__file__).parent.parent / "data" / "processed" / "mouse_common" / "ligand_receptor_mouse.csv"


def load_gene_descriptions():
    """Load gene descriptions lookup table."""
    if GENE_DESCRIPTIONS_PATH.exists():
        df = pd.read_csv(GENE_DESCRIPTIONS_PATH)
        # Create lookup dict: gene_symbol -> {protein_name, description, ...}
        gene_info = {}
        for _, row in df.iterrows():
            gene_info[row['gene_symbol']] = {
                'protein_name': row.get('protein_name', ''),
                'description': row.get('description', ''),
                'go_biological_process': row.get('go_biological_process', ''),
            }
        print(f"Loaded {len(gene_info)} gene descriptions")
        return gene_info
    else:
        print(f"Warning: Gene descriptions not found at {GENE_DESCRIPTIONS_PATH}")
        return {}


def load_region_descriptions():
    """Load region descriptions lookup table."""
    if REGION_DESCRIPTIONS_PATH.exists():
        df = pd.read_csv(REGION_DESCRIPTIONS_PATH)
        # Create lookup dict: acronym -> {full_name, description}
        region_info = {}
        for _, row in df.iterrows():
            region_info[row['acronym']] = {
                'full_name': row.get('full_name', ''),
                'description': row.get('description', ''),
            }
        print(f"Loaded {len(region_info)} region descriptions")
        return region_info
    else:
        print(f"Warning: Region descriptions not found at {REGION_DESCRIPTIONS_PATH}")
        return {}


def load_signaling_data():
    """Load ligand-receptor signaling data for cluster-level expression.

    Returns:
        Tuple of (ligand_to_receptors, cluster_expression, ligand_genes, receptor_genes)
    """
    if not CLUSTER_LR_PROFILE_PATH.exists() or not LR_PAIRS_PATH.exists():
        print("Warning: Signaling data not found")
        return {}, {}, set(), set()

    # Load ligand-receptor pairs
    lr_pairs = pd.read_csv(LR_PAIRS_PATH)

    # Filter to neuropeptides/hormones/amines (exclude GABA/Glu etc)
    lr_pairs = lr_pairs[lr_pairs['ligand_type'].isin(['Neuropeptide', 'Hormone', 'Amine'])]

    # Create lookup: ligand -> [receptors]
    ligand_to_receptors = lr_pairs.groupby('gene_ligand')['gene_receptor'].apply(list).to_dict()

    # Get sets of ligand and receptor genes
    ligand_genes = set(lr_pairs['gene_ligand'].unique())
    receptor_genes = set(lr_pairs['gene_receptor'].unique())

    # Load cluster expression profiles
    cluster_lr = pd.read_parquet(CLUSTER_LR_PROFILE_PATH)

    # Create lookup: cluster -> {gene: pct_expressing}
    cluster_expression = cluster_lr.pivot_table(
        index='cluster', columns='gene', values='pct_expressing', aggfunc='first'
    ).to_dict('index')

    print(f"Loaded signaling data: {len(ligand_genes)} ligands, {len(receptor_genes)} receptors, {len(cluster_expression)} clusters")

    return ligand_to_receptors, cluster_expression, ligand_genes, receptor_genes


def load_all_datasets():
    """Load all available processed datasets.

    Returns:
        Dictionary mapping dataset name to (cells_df, config) tuple
    """
    # Temporarily disabled datasets
    DISABLED_DATASETS = {'human_hypomap'}

    datasets = {}
    available = [d for d in get_processed_datasets() if d not in DISABLED_DATASETS]

    if not available:
        raise FileNotFoundError(
            "No processed datasets found.\n"
            "Run preprocessing first:\n"
            "  python -m src.datasets.mouse_hypomap\n"
            "  python -m src.preprocessing.downsample --dataset mouse_hypomap\n"
            "  python -m src.datasets.human_hypomap\n"
            "  python -m src.preprocessing.downsample --dataset human_hypomap\n"
            "  python -m src.datasets.mouse_abc"
        )

    for name in available:
        try:
            cells_df, config = load_dataset(name)
            datasets[name] = {
                'cells_df': cells_df,
                'config': config,
                'region_col': detect_region_column(cells_df),
                'cell_type_levels': detect_cell_type_levels(cells_df),
                'regions': cells_df[detect_region_column(cells_df)].unique().tolist()
                    if detect_region_column(cells_df) else [],
                'marker_genes': get_marker_genes(name),
                'key_receptors': get_key_receptors(name),
                'region_colors': get_region_colors(name),
            }
            print(f"Loaded {name}: {len(cells_df)} cells")
        except Exception as e:
            print(f"Warning: Could not load {name}: {e}")

    return datasets


def create_app():
    """Create and configure the Dash application."""
    # Load all datasets
    datasets = load_all_datasets()
    available_datasets = list(datasets.keys())

    if not available_datasets:
        raise ValueError("No datasets could be loaded")

    # Load gene descriptions
    gene_descriptions = load_gene_descriptions()

    # Load region descriptions
    region_descriptions = load_region_descriptions()

    # Load signaling data
    ligand_to_receptors, cluster_expression, ligand_genes, receptor_genes = load_signaling_data()
    signaling_data = {
        'ligand_to_receptors': ligand_to_receptors,
        'cluster_expression': cluster_expression,
        'ligand_genes': ligand_genes,
        'receptor_genes': receptor_genes,
    }

    # Use mouse_abc as default for now
    default_dataset = 'mouse_abc' if 'mouse_abc' in available_datasets else available_datasets[0]
    default_data = datasets[default_dataset]

    print(f"\nAvailable datasets: {available_datasets}")
    print(f"Default: {default_dataset}")
    print(f"  Region column: {default_data['region_col']}")
    print(f"  Cell type levels: {default_data['cell_type_levels']}")
    print(f"  Regions: {len(default_data['regions'])}")

    # Create Dash app
    app = Dash(
        __name__,
        external_stylesheets=[dbc.themes.BOOTSTRAP],
        title="HypoMap 3D Viewer",
        suppress_callback_exceptions=True,
    )

    # Set layout with dataset options
    app.layout = create_layout(
        datasets=datasets,
        available_datasets=available_datasets,
        default_dataset=default_dataset,
        ligand_options=sorted(ligand_to_receptors.keys()),
    )

    # Register callbacks
    register_callbacks(app, datasets, gene_descriptions, region_descriptions, signaling_data)

    return app


# Create app instance
app = create_app()

if __name__ == "__main__":
    print("\n" + "=" * 50)
    print("HypoMap 3D Atlas Viewer")
    print("=" * 50)
    print("\nStarting server at http://localhost:8050")
    print("Press Ctrl+C to stop\n")

    app.run(debug=True, port=8050)
