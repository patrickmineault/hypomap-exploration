"""Dataset handling for HypoMap viewer.

This module provides a unified interface for loading and processing
single-cell hypothalamus atlases from different species.
"""

from .base import DatasetConfig, STANDARD_COLUMNS
from .loader import (
    get_available_datasets,
    get_processed_datasets,
    get_config,
    load_dataset,
    load_cell_metadata,
    get_gene_list,
    get_marker_genes,
    get_key_receptors,
    get_region_colors,
    detect_cell_type_levels,
    detect_region_column,
)

__all__ = [
    "DatasetConfig",
    "STANDARD_COLUMNS",
    "get_available_datasets",
    "get_processed_datasets",
    "get_config",
    "load_dataset",
    "load_cell_metadata",
    "get_gene_list",
    "get_marker_genes",
    "get_key_receptors",
    "get_region_colors",
    "detect_cell_type_levels",
    "detect_region_column",
]
