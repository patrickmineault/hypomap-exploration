"""Base classes and types for dataset handling."""

from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional, Dict


@dataclass
class DatasetConfig:
    """Configuration for a single-cell dataset.

    Attributes:
        name: Short identifier (e.g., "mouse_hypomap", "human_hypomap", "mouse_abc")
        species: Species name (e.g., "Mus musculus", "Homo sapiens")
        h5ad_path: Path to the source h5ad file
        processed_dir: Directory for processed outputs
        cell_type_columns: List of cell type hierarchy columns (e.g., ["C7_named", "C66_named"])
        region_column: Name of the region column in obs
        gene_column: var column containing gene symbols (None if var_names are already symbols)
        spatial_data_path: Optional path to spatial transcriptomics data
    """
    name: str
    species: str
    h5ad_path: Path
    processed_dir: Path
    cell_type_columns: List[str]
    region_column: str
    gene_column: Optional[str] = None
    spatial_data_path: Optional[Path] = None

    # Additional metadata
    marker_genes: List[str] = field(default_factory=list)
    key_receptors: List[str] = field(default_factory=list)
    region_colors: Dict[str, str] = field(default_factory=dict)

    @property
    def cell_metadata_path(self) -> Path:
        """Path to cell metadata parquet file."""
        return self.processed_dir / "cell_metadata.parquet"

    @property
    def cells_with_coords_path(self) -> Path:
        """Path to downsampled cells with coordinates."""
        return self.processed_dir / "cells_with_coords.parquet"

    @property
    def genes_path(self) -> Path:
        """Path to gene list parquet file."""
        return self.processed_dir / "genes.parquet"

    def is_processed(self) -> bool:
        """Check if dataset has been processed."""
        return (
            self.cell_metadata_path.exists() and
            self.cells_with_coords_path.exists() and
            self.genes_path.exists()
        )


# Standard output column names for the common intermediate format
STANDARD_COLUMNS = {
    "cell_id": "Unique cell identifier",
    "region": "Anatomical region (standardized name)",
    "x": "X coordinate (species-specific units)",
    "y": "Y coordinate (species-specific units)",
    "z": "Z coordinate (species-specific units)",
    "umap_1": "UMAP dimension 1",
    "umap_2": "UMAP dimension 2",
    "sample_id": "Sample/section identifier",
    "donor_id": "Donor/subject identifier",
}
