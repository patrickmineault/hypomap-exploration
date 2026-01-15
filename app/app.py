"""Main Dash application for HypoMap 3D Viewer."""

import sys
from pathlib import Path

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


def load_all_datasets():
    """Load all available processed datasets.

    Returns:
        Dictionary mapping dataset name to (cells_df, config) tuple
    """
    datasets = {}
    available = get_processed_datasets()

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

    # Use human as default (more spatially precise), fall back to first available
    default_dataset = 'human' if 'human' in available_datasets else available_datasets[0]
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
    )

    # Register callbacks
    register_callbacks(app, datasets)

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
