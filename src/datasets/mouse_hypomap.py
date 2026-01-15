"""Mouse HypoMap dataset adapter.

Handles extraction and processing of the mouse hypothalamus single-cell atlas.
"""

import scanpy as sc
import pandas as pd
from pathlib import Path
from typing import List, Optional

from .base import DatasetConfig

# Paths
DATA_DIR = Path(__file__).parent.parent.parent / "data"
RAW_DIR = DATA_DIR / "raw" / "mouse_hypomap"
PROCESSED_DIR = DATA_DIR / "processed" / "mouse_hypomap"

# The mouse h5ad file
H5AD_FILE = RAW_DIR / "d3be7423-d664-4913-89a9-a506cae4c28f.h5ad"

# Fallback to old location for backward compatibility
if not H5AD_FILE.exists():
    _old_path = DATA_DIR / "d3be7423-d664-4913-89a9-a506cae4c28f.h5ad"
    if _old_path.exists():
        H5AD_FILE = _old_path

# Mouse-specific marker genes and receptors
MOUSE_MARKER_GENES = [
    'Agrp', 'Pomc', 'Npy', 'Sf1', 'Oxt', 'Avp', 'Crh', 'Trh', 'Ghrh', 'Sst',
]

MOUSE_KEY_RECEPTORS = [
    'Lepr', 'Mc4r', 'Glp1r', 'Ghsr', 'Insr', 'Npy1r', 'Npy2r', 'Crhr1', 'Oxtr', 'Avpr1a',
]

# Mouse region colors (from existing region_mapping.py)
MOUSE_REGION_COLORS = {
    'Arcuate hypothalamic nucleus': '#DC143C',
    'Ventromedial hypothalamic nucleus': '#0066CC',
    'Paraventricular hypothalamic nucleus': '#228B22',
    'Dorsomedial nucleus of the hypothalamus': '#8B008B',
    'Lateral hypothalamic area': '#FF6600',
    'Suprachiasmatic nucleus': '#FFD700',
    'Medial preoptic area': '#8B4513',
    'Lateral preoptic area': '#FF1493',
    'Tuberal nucleus': '#696969',
    'Zona incerta': '#008B8B',
    'Anterior hypothalamic nucleus': '#FF4500',
    'Posterior hypothalamic nucleus': '#4169E1',
    'Supraoptic nucleus': '#C71585',
    'Periventricular hypothalamic nucleus': '#32CD32',
    'Periventricular hypothalamic nucleus, posterior part': '#00CC66',
    'Periventricular hypothalamic nucleus, intermediate part': '#66CDAA',
    '(Pre)Mammillary region': '#DAA520',
    '(Anterior/Preoptic)Periventricular region': '#9932CC',
    'Mammillary body': '#DAA520',
    'Preoptic area': '#CD853F',
    'NA': '#A0A0A0',
}


def get_mouse_hypomap_config() -> DatasetConfig:
    """Get configuration for the mouse HypoMap dataset."""
    return DatasetConfig(
        name="mouse_hypomap",
        species="Mus musculus",
        h5ad_path=H5AD_FILE,
        processed_dir=PROCESSED_DIR,
        cell_type_columns=[],  # Will be detected dynamically
        region_column="Region_summarized",  # Detected from data
        gene_column="feature_name",  # Gene symbols are in var['feature_name']
        marker_genes=MOUSE_MARKER_GENES,
        key_receptors=MOUSE_KEY_RECEPTORS,
        region_colors=MOUSE_REGION_COLORS,
    )


def detect_cell_type_columns(obs_columns: List[str]) -> List[str]:
    """Detect cell type hierarchy columns (C7_named, C66_named, etc.)."""
    cell_type_cols = [col for col in obs_columns if col.startswith('C') and '_named' in col]
    # Sort by the number in the column name
    cell_type_cols = sorted(cell_type_cols, key=lambda x: int(x.split('_')[0][1:]))
    return cell_type_cols


def detect_region_column(obs_columns: List[str]) -> Optional[str]:
    """Detect the region column."""
    region_cols = [col for col in obs_columns if 'region' in col.lower() or 'Region' in col]
    return region_cols[0] if region_cols else None


def extract_mouse_metadata(h5ad_path: Optional[Path] = None) -> pd.DataFrame:
    """Extract cell metadata from mouse h5ad file.

    Args:
        h5ad_path: Path to h5ad file. If None, uses default.

    Returns:
        DataFrame with cell metadata
    """
    if h5ad_path is None:
        h5ad_path = H5AD_FILE

    print(f"Loading mouse h5ad file: {h5ad_path}")
    adata = sc.read_h5ad(h5ad_path, backed='r')

    print(f"Loaded {adata.n_obs} cells, {adata.n_vars} genes")
    print(f"Available obs columns: {list(adata.obs.columns)}")

    # Extract metadata
    metadata = adata.obs.copy()

    # Detect columns
    cell_type_cols = detect_cell_type_columns(metadata.columns)
    region_col = detect_region_column(metadata.columns)

    print(f"Cell type columns found: {cell_type_cols}")
    print(f"Region column found: {region_col}")

    # Add UMAP coordinates if available
    if 'X_umap' in adata.obsm:
        umap_coords = adata.obsm['X_umap']
        metadata['umap_1'] = umap_coords[:, 0]
        metadata['umap_2'] = umap_coords[:, 1]
        print("Added UMAP coordinates")

    # Ensure output directory exists
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

    # Save gene names - use feature_name if available (gene symbols)
    if 'feature_name' in adata.var.columns:
        gene_symbols = adata.var['feature_name'].tolist()
        ensembl_ids = adata.var_names.tolist()
        genes_df = pd.DataFrame({
            'gene': gene_symbols,
            'ensembl_id': ensembl_ids,
        })
        # Remove duplicates, keeping first
        genes_df = genes_df.drop_duplicates(subset=['gene'], keep='first')
        print("Using feature_name column for gene symbols")
    else:
        genes_df = pd.DataFrame({'gene': adata.var_names.tolist()})

    genes_path = PROCESSED_DIR / "genes.parquet"
    genes_df.to_parquet(genes_path)
    print(f"Saved {len(genes_df)} gene names to {genes_path}")

    # Save full metadata
    metadata_path = PROCESSED_DIR / "cell_metadata.parquet"
    metadata.to_parquet(metadata_path)
    print(f"Saved metadata for {len(metadata)} cells to {metadata_path}")

    return metadata


def print_mouse_summary(metadata: pd.DataFrame):
    """Print summary statistics of the mouse dataset."""
    print("\n=== Mouse Dataset Summary ===")
    print(f"Total cells: {len(metadata)}")

    # Region distribution
    region_col = detect_region_column(metadata.columns)
    if region_col:
        print(f"\nRegion distribution ({region_col}):")
        print(metadata[region_col].value_counts())

    # Cell type levels
    cell_type_cols = detect_cell_type_columns(metadata.columns)
    if cell_type_cols:
        print(f"\nCell type hierarchy levels: {cell_type_cols}")
        for col in cell_type_cols[:3]:
            n_types = metadata[col].nunique()
            print(f"  {col}: {n_types} unique types")


if __name__ == "__main__":
    metadata = extract_mouse_metadata()
    print_mouse_summary(metadata)
