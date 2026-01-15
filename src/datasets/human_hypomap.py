"""Human HypoMap dataset adapter.

Handles extraction and processing of the human hypothalamus single-cell atlas.
"""

import tarfile
import scanpy as sc
import pandas as pd
from pathlib import Path
from typing import Tuple, List, Optional, Dict
import re

from .base import DatasetConfig

# Paths
DATA_DIR = Path(__file__).parent.parent.parent / "data"
RAW_DIR = DATA_DIR / "raw" / "human_hypomap"
PROCESSED_DIR = DATA_DIR / "processed" / "human_hypomap"

# Human h5ad file
H5AD_FILE = RAW_DIR / "480e89e7-84ad-4fa8-adc3-f7c562a77a78.h5ad"
SPATIAL_POSITIONS_TAR = RAW_DIR / "GSE278848_tissue_positions_lists.tar.gz"

# Human marker genes (uppercase)
HUMAN_MARKER_GENES = [
    'AGRP', 'POMC', 'NPY', 'SF1', 'OXT', 'AVP', 'CRH', 'TRH', 'GHRH', 'SST',
]

HUMAN_KEY_RECEPTORS = [
    'LEPR', 'MC4R', 'GLP1R', 'GHSR', 'INSR', 'NPY1R', 'NPY2R', 'CRHR1', 'OXTR', 'AVPR1A',
]

# Human region mapping (abbreviation -> full name) and colors
HUMAN_REGIONS = {
    'ARC': ('Arcuate nucleus', '#DC143C'),
    'VMH': ('Ventromedial hypothalamus', '#0066CC'),
    'PVN': ('Paraventricular nucleus', '#228B22'),
    'DMH': ('Dorsomedial hypothalamus', '#8B008B'),
    'LH': ('Lateral hypothalamus', '#FF6600'),
    'SCN': ('Suprachiasmatic nucleus', '#FFD700'),
    'MPOA': ('Medial preoptic area', '#8B4513'),
    'POA': ('Preoptic area', '#CD853F'),
    'LPOA': ('Lateral preoptic area', '#FF1493'),
    'TMN': ('Tuberomammillary nucleus', '#FF4500'),
    'SON': ('Supraoptic nucleus', '#C71585'),
    'ME': ('Median eminence', '#4169E1'),
    'Fx/OT/ac': ('Fornix/Optic tract/Anterior commissure', '#696969'),
    'MAM': ('Mammillary region', '#DAA520'),
    'Perivent': ('Periventricular region', '#32CD32'),
    'LTN': ('Lateral tuberal nucleus', '#008B8B'),
    'Vent': ('Ventricle', '#A9A9A9'),
    'Vascular': ('Vascular', '#B22222'),
    'Thalamus': ('Thalamus', '#9932CC'),
    'Thalamaus': ('Thalamus', '#9932CC'),  # Handle typo in dataset
    'NA': ('Unknown', '#808080'),
}

HUMAN_REGION_COLORS = {abbr: color for abbr, (name, color) in HUMAN_REGIONS.items()}


def get_human_hypomap_config() -> DatasetConfig:
    """Get configuration for the human HypoMap dataset."""
    return DatasetConfig(
        name="human_hypomap",
        species="Homo sapiens",
        h5ad_path=H5AD_FILE,
        processed_dir=PROCESSED_DIR,
        cell_type_columns=[],  # Will be detected: C0_named through C4_named
        region_column="region",
        gene_column="feature_name",  # Gene symbols are in var['feature_name']
        spatial_data_path=SPATIAL_POSITIONS_TAR,
        marker_genes=HUMAN_MARKER_GENES,
        key_receptors=HUMAN_KEY_RECEPTORS,
        region_colors=HUMAN_REGION_COLORS,
    )


def detect_cell_type_columns(obs_columns: List[str]) -> List[str]:
    """Detect cell type hierarchy columns (C0_named, C1_named, etc.)."""
    cell_type_cols = [col for col in obs_columns if col.startswith('C') and '_named' in col]
    # Sort by the number in the column name
    cell_type_cols = sorted(cell_type_cols, key=lambda x: int(x.split('_')[0][1:]))
    return cell_type_cols


def load_spatial_positions() -> Dict[str, pd.DataFrame]:
    """Load spatial positions from 10x Visium tissue positions files.

    Returns:
        Dictionary mapping sample ID to DataFrame with columns:
        barcode, in_tissue, array_row, array_col, pixel_x, pixel_y
    """
    if not SPATIAL_POSITIONS_TAR.exists():
        print(f"Warning: Spatial positions file not found: {SPATIAL_POSITIONS_TAR}")
        return {}

    spatial_data = {}

    with tarfile.open(SPATIAL_POSITIONS_TAR, 'r:gz') as tar:
        for member in tar.getmembers():
            if member.name.endswith('_tissue_positions_list.csv'):
                # Extract sample ID from filename (e.g., "bT5r4_b8_tissue_positions_list.csv")
                sample_id = member.name.replace('_tissue_positions_list.csv', '')

                f = tar.extractfile(member)
                if f is not None:
                    df = pd.read_csv(f, header=None, names=[
                        'barcode', 'in_tissue', 'array_row', 'array_col', 'pixel_x', 'pixel_y'
                    ])
                    spatial_data[sample_id] = df

    print(f"Loaded spatial positions for {len(spatial_data)} samples")
    return spatial_data


def parse_sample_id(sample_id: str) -> Tuple[str, Optional[str], Optional[int]]:
    """Parse sample ID to extract donor, section type, and slice number.

    Args:
        sample_id: Sample ID like "3u5kk_A.01" or "buFpQ_A3"

    Returns:
        Tuple of (donor_id, section_type, slice_num)
    """
    # Try pattern like "3u5kk_A.01" or "3u5kk_R.11"
    match = re.match(r'^([^_]+)_([AR])\.(\d+)$', sample_id)
    if match:
        return match.group(1), match.group(2), int(match.group(3))

    # Try pattern like "buFpQ_A3" or "1C41i_R2"
    match = re.match(r'^([^_]+)_([AR])(\d+)$', sample_id)
    if match:
        return match.group(1), match.group(2), int(match.group(3))

    # Try pattern like "TweV8_N1" or "znZv1_S5"
    match = re.match(r'^([^_]+)_([A-Z]+)(\d*)$', sample_id)
    if match:
        num = int(match.group(3)) if match.group(3) else None
        return match.group(1), match.group(2), num

    # Fallback: just return the whole thing as donor
    return sample_id, None, None


def extract_human_metadata(h5ad_path: Optional[Path] = None) -> pd.DataFrame:
    """Extract cell metadata from human h5ad file.

    Args:
        h5ad_path: Path to h5ad file. If None, uses default.

    Returns:
        DataFrame with cell metadata
    """
    if h5ad_path is None:
        h5ad_path = H5AD_FILE

    print(f"Loading human h5ad file: {h5ad_path}")
    adata = sc.read_h5ad(h5ad_path, backed='r')

    print(f"Loaded {adata.n_obs} cells, {adata.n_vars} genes")
    print(f"Available obs columns: {list(adata.obs.columns)}")

    # Extract metadata
    metadata = adata.obs.copy()

    # Detect cell type columns
    cell_type_cols = detect_cell_type_columns(metadata.columns)
    print(f"Cell type columns found: {cell_type_cols}")

    # Add UMAP coordinates if available
    if 'X_umap' in adata.obsm:
        umap_coords = adata.obsm['X_umap']
        metadata['umap_1'] = umap_coords[:, 0]
        metadata['umap_2'] = umap_coords[:, 1]
        print("Added UMAP coordinates")

    # Ensure output directory exists
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

    # Save gene names - map Ensembl IDs to symbols
    print("Creating gene name mapping (Ensembl -> symbols)...")
    gene_symbols = adata.var['feature_name'].tolist()
    ensembl_ids = adata.var_names.tolist()

    genes_df = pd.DataFrame({
        'gene': gene_symbols,
        'ensembl_id': ensembl_ids,
    })
    # Remove duplicates, keeping first occurrence
    genes_df = genes_df.drop_duplicates(subset=['gene'], keep='first')

    genes_path = PROCESSED_DIR / "genes.parquet"
    genes_df.to_parquet(genes_path)
    print(f"Saved {len(genes_df)} unique gene symbols to {genes_path}")

    # Save full metadata
    metadata_path = PROCESSED_DIR / "cell_metadata.parquet"
    metadata.to_parquet(metadata_path)
    print(f"Saved metadata for {len(metadata)} cells to {metadata_path}")

    return metadata


def print_human_summary(metadata: pd.DataFrame):
    """Print summary statistics of the human dataset."""
    print("\n=== Human Dataset Summary ===")
    print(f"Total cells: {len(metadata)}")

    # Region distribution
    if 'region' in metadata.columns:
        print("\nRegion distribution:")
        print(metadata['region'].value_counts())

    # Cell type levels
    cell_type_cols = detect_cell_type_columns(metadata.columns)
    if cell_type_cols:
        print(f"\nCell type hierarchy levels: {cell_type_cols}")
        for col in cell_type_cols:
            n_types = metadata[col].nunique()
            print(f"  {col}: {n_types} unique types")

    # Sample/donor info
    if 'Sample_ID' in metadata.columns:
        print(f"\nSamples: {metadata['Sample_ID'].nunique()}")
    if 'donor_id' in metadata.columns:
        print(f"Donors: {metadata['donor_id'].nunique()}")
        print(metadata['donor_id'].value_counts())


if __name__ == "__main__":
    metadata = extract_human_metadata()
    print_human_summary(metadata)

    # Also test spatial loading
    print("\n=== Testing Spatial Data Loading ===")
    spatial = load_spatial_positions()
    for sample_id, df in list(spatial.items())[:3]:
        print(f"\n{sample_id}: {len(df)} spots")
        print(df.head(3))
