"""Unified data loading interface for HypoMap datasets."""

from typing import Dict, List, Optional, Tuple

import pandas as pd

from .base import DatasetConfig
from .mouse_abc import get_mouse_abc_config, get_mouse_abc_extended_config

# Registry of available datasets
DATASETS = {
    "mouse_abc": get_mouse_abc_config,
    "mouse_abc_extended": get_mouse_abc_extended_config,
}

# Optional datasets (modules may not exist)
try:
    from .mouse_hypomap import get_mouse_hypomap_config

    DATASETS["mouse_hypomap"] = get_mouse_hypomap_config
except ImportError:
    pass

try:
    from .human_hypomap import get_human_hypomap_config

    DATASETS["human_hypomap"] = get_human_hypomap_config
except ImportError:
    pass


def get_available_datasets() -> List[str]:
    """Return list of all registered dataset names."""
    return list(DATASETS.keys())


def get_processed_datasets() -> List[str]:
    """Return list of datasets that have been processed."""
    processed = []
    for name in DATASETS:
        config = DATASETS[name]()
        if config.is_processed():
            processed.append(name)
    return processed


def get_config(name: str) -> DatasetConfig:
    """Get configuration for a dataset.

    Args:
        name: Dataset name ("mouse_hypomap", "human_hypomap", or "mouse_abc")

    Returns:
        DatasetConfig instance
    """
    if name not in DATASETS:
        raise ValueError(f"Unknown dataset: {name}. Available: {list(DATASETS.keys())}")
    return DATASETS[name]()


def load_dataset(name: str) -> Tuple[pd.DataFrame, DatasetConfig]:
    """Load a processed dataset.

    Args:
        name: Dataset name ("mouse_hypomap", "human_hypomap", or "mouse_abc")

    Returns:
        Tuple of (cells_df, config)
    """
    config = get_config(name)

    if not config.cells_with_coords_path.exists():
        raise FileNotFoundError(
            f"Processed data not found for {name} at {config.cells_with_coords_path}\n"
            f"Run preprocessing first: python -m hypomap.datasets.{name}"
        )

    cells_df = pd.read_parquet(config.cells_with_coords_path)
    return cells_df, config


def load_cell_metadata(name: str) -> Tuple[pd.DataFrame, DatasetConfig]:
    """Load cell metadata (before downsampling/coordinate assignment).

    Args:
        name: Dataset name

    Returns:
        Tuple of (metadata_df, config)
    """
    config = get_config(name)

    if not config.cell_metadata_path.exists():
        raise FileNotFoundError(
            f"Cell metadata not found for {name}. Run extraction first."
        )

    metadata_df = pd.read_parquet(config.cell_metadata_path)
    return metadata_df, config


def get_gene_list(name: str) -> List[str]:
    """Get list of genes for a dataset.

    Args:
        name: Dataset name

    Returns:
        List of gene symbols
    """
    config = get_config(name)

    if not config.genes_path.exists():
        raise FileNotFoundError(f"Gene list not found for {name}")

    genes_df = pd.read_parquet(config.genes_path)
    return genes_df["gene"].tolist()


def get_marker_genes(name: str) -> List[str]:
    """Get marker genes for a dataset.

    Args:
        name: Dataset name

    Returns:
        List of marker gene symbols
    """
    config = get_config(name)
    return config.marker_genes


def get_key_receptors(name: str) -> List[str]:
    """Get key receptors for a dataset.

    Args:
        name: Dataset name

    Returns:
        List of receptor gene symbols
    """
    config = get_config(name)
    return config.key_receptors


def get_region_colors(name: str) -> Dict[str, str]:
    """Get region color mapping for a dataset.

    Args:
        name: Dataset name

    Returns:
        Dictionary mapping region names to hex colors
    """
    config = get_config(name)
    return config.region_colors


def detect_cell_type_levels(cells_df: pd.DataFrame) -> List[str]:
    """Detect cell type hierarchy columns in a dataframe.

    Args:
        cells_df: DataFrame with cell data

    Returns:
        List of cell type column names, sorted by hierarchy level
    """
    # Check for ABC-style columns first (class, subclass, supertype, cluster)
    abc_cols = ["class", "subclass", "supertype", "cluster"]
    abc_found = [col for col in abc_cols if col in cells_df.columns]
    if abc_found:
        return abc_found

    # Fall back to HypoMap-style columns (C7_named, C66_named, etc.)
    cell_type_cols = [
        col for col in cells_df.columns if col.startswith("C") and "_named" in col
    ]
    cell_type_cols = sorted(cell_type_cols, key=lambda x: int(x.split("_")[0][1:]))
    return cell_type_cols


def detect_region_column(cells_df: pd.DataFrame) -> Optional[str]:
    """Detect the region column in a dataframe.

    Args:
        cells_df: DataFrame with cell data

    Returns:
        Name of region column, or None if not found
    """
    # Check for exact matches first
    if "region" in cells_df.columns:
        return "region"
    if "Region" in cells_df.columns:
        return "Region"

    # Fall back to pattern matching
    region_cols = [col for col in cells_df.columns if "region" in col.lower()]
    return region_cols[0] if region_cols else None
    return region_cols[0] if region_cols else None
