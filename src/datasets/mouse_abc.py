"""Allen Brain Cell Census (ABC) Hypothalamus dataset adapter.

Handles extraction and processing of the ABC MERFISH hypothalamus data.
Uses AbcProjectCache to load data from the Allen Institute cloud storage.
"""

import pandas as pd
from pathlib import Path
from typing import Optional

from .base import DatasetConfig

# Paths
DATA_DIR = Path(__file__).parent.parent.parent / "data"
DOWNLOAD_BASE = DATA_DIR / "raw" / "abc_atlas_cache"
PROCESSED_DIR = DATA_DIR / "processed" / "mouse_abc"

# ABC hypothalamus marker genes (mouse gene symbols)
ABC_MARKER_GENES = [
    'Agrp', 'Pomc', 'Npy', 'Sf1', 'Oxt', 'Avp', 'Crh', 'Trh', 'Ghrh', 'Sst',
]

ABC_KEY_RECEPTORS = [
    'Lepr', 'Mc4r', 'Glp1r', 'Ghsr', 'Insr', 'Npy1r', 'Npy2r', 'Crhr1', 'Oxtr', 'Avpr1a',
]

# Region colors for hypothalamic subregions (matching mouse hypomap where possible)
ABC_REGION_COLORS = {
    'AHN': '#FF4500',     # Anterior hypothalamic nucleus
    'ARH': '#DC143C',     # Arcuate nucleus
    'AVPV': '#9932CC',    # Anteroventral periventricular nucleus
    'DMH': '#8B008B',     # Dorsomedial hypothalamus
    'LHA': '#FF6600',     # Lateral hypothalamic area
    'LPO': '#FF1493',     # Lateral preoptic area
    'ME': '#DAA520',      # Median eminence
    'MEA': '#CD853F',     # Medial amygdalar nucleus
    'MPN': '#8B4513',     # Medial preoptic nucleus
    'MPO': '#8B4513',     # Medial preoptic area
    'PH': '#4169E1',      # Posterior hypothalamus
    'PMd': '#DAA520',     # Premammillary nucleus, dorsal
    'PMv': '#DAA520',     # Premammillary nucleus, ventral
    'PVH': '#228B22',     # Paraventricular hypothalamic nucleus
    'PVHd': '#228B22',
    'PVi': '#32CD32',     # Periventricular hypothalamic nucleus, intermediate
    'PVp': '#00CC66',     # Periventricular hypothalamic nucleus, posterior
    'PVpo': '#66CDAA',    # Periventricular hypothalamic nucleus, preoptic
    'SCH': '#FFD700',     # Suprachiasmatic nucleus
    'SFO': '#008B8B',     # Subfornical organ
    'SO': '#C71585',      # Supraoptic nucleus
    'STN': '#696969',     # Subthalamic nucleus
    'TU': '#696969',      # Tuberal nucleus
    'VMH': '#0066CC',     # Ventromedial hypothalamus
    'ZI': '#008B8B',      # Zona incerta
    'HY-unassigned': '#666666',
    'NA': '#A0A0A0',
}


def get_mouse_abc_config() -> DatasetConfig:
    """Get configuration for the mouse ABC hypothalamus dataset."""
    return DatasetConfig(
        name="mouse_abc",
        species="Mus musculus",
        h5ad_path=DOWNLOAD_BASE,  # Not an h5ad, but points to cache dir
        processed_dir=PROCESSED_DIR,
        cell_type_columns=['class', 'subclass', 'supertype', 'cluster'],
        region_column="region",
        gene_column=None,  # Genes loaded separately
        marker_genes=ABC_MARKER_GENES,
        key_receptors=ABC_KEY_RECEPTORS,
        region_colors=ABC_REGION_COLORS,
    )


def extract_mouse_abc_metadata(cache_dir: Optional[Path] = None) -> pd.DataFrame:
    """Extract cell metadata from ABC MERFISH data, filtered to hypothalamus.

    Args:
        cache_dir: Path to ABC cache directory. If None, uses default.

    Returns:
        DataFrame with cell metadata for hypothalamic cells
    """
    from abc_atlas_access.abc_atlas_cache.abc_project_cache import AbcProjectCache

    if cache_dir is None:
        cache_dir = DOWNLOAD_BASE

    # Ensure cache directory exists
    cache_dir.mkdir(parents=True, exist_ok=True)

    print(f"Initializing ABC cache at: {cache_dir}")
    cache = AbcProjectCache.from_cache_dir(cache_dir)
    cache.load_latest_manifest()

    # 1. Load MERFISH cell metadata with cluster annotations
    print("Loading cell metadata with cluster annotations...")
    cell = cache.get_metadata_dataframe(
        directory="MERFISH-C57BL6J-638850",
        file_name="cell_metadata_with_cluster_annotation",
        dtype={"cell_label": str}
    ).set_index("cell_label", drop=False)
    print(f"Loaded {len(cell)} total cells")

    # 2. Load reconstructed coordinates (CCF registration)
    print("Loading reconstructed coordinates...")
    recon = cache.get_metadata_dataframe(
        directory="MERFISH-C57BL6J-638850-CCF",
        file_name="reconstructed_coordinates",
        dtype={"cell_label": str}
    ).set_index("cell_label", drop=False)

    # Rename coordinate columns to avoid confusion
    recon = recon.rename(columns={
        'x': 'x_reconstructed',
        'y': 'y_reconstructed',
        'z': 'z_reconstructed'
    })

    # 3. Load parcellation terms (anatomy labels)
    print("Loading parcellation terms...")
    parc_terms = cache.get_metadata_dataframe(
        directory="Allen-CCF-2020",
        file_name="parcellation_to_parcellation_term_membership_acronym"
    ).set_index("parcellation_index")
    parc_terms.columns = [f"parcellation_{c}" for c in parc_terms.columns]

    # 4. Join tables
    print("Joining tables...")
    # Join parcellation_index from reconstructed coordinates
    cell = cell.join(
        recon.set_index("cell_label")[["parcellation_index", "x_reconstructed", "y_reconstructed", "z_reconstructed"]],
        on="cell_label",
        how="inner"
    )

    # Join parcellation terms
    cell = cell.join(parc_terms, on="parcellation_index", how="left")

    print(f"After joining: {len(cell)} cells")

    # 5. Filter to hypothalamus (HY)
    print("Filtering to hypothalamus (HY)...")
    hy_mask = cell['parcellation_division'] == 'HY'
    cell = cell[hy_mask].copy()
    print(f"HY cells: {len(cell)}")

    # 6. Rename/standardize columns for hypomap compatibility
    cell = cell.rename(columns={
        'cell_label': 'cell_id',
        'x_reconstructed': 'x',
        'y_reconstructed': 'y',
        'z_reconstructed': 'z',
        'parcellation_structure': 'region',
        'donor_label': 'donor_id',
        'brain_section_label': 'sample_id',
    })

    # Ensure output directory exists
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

    # 7. Save cell metadata
    metadata_path = PROCESSED_DIR / "cell_metadata.parquet"
    cell.to_parquet(metadata_path)
    print(f"Saved metadata for {len(cell)} cells to {metadata_path}")

    # 8. Save gene list (from ABC gene table)
    print("Loading gene table...")
    gene_table = cache.get_metadata_dataframe(
        directory="WMB-10X",
        file_name="gene"
    )
    genes_df = pd.DataFrame({
        'gene': gene_table['gene_symbol'].tolist(),
        'gene_identifier': gene_table['gene_identifier'].tolist() if 'gene_identifier' in gene_table.columns else None,
    })
    # Remove duplicates
    genes_df = genes_df.drop_duplicates(subset=['gene'], keep='first')

    genes_path = PROCESSED_DIR / "genes.parquet"
    genes_df.to_parquet(genes_path)
    print(f"Saved {len(genes_df)} genes to {genes_path}")

    return cell


def print_mouse_abc_summary(metadata: pd.DataFrame):
    """Print summary statistics of the mouse ABC hypothalamus dataset."""
    print("\n=== Mouse ABC Hypothalamus Dataset Summary ===")
    print(f"Total cells: {len(metadata)}")

    # Region distribution
    if 'region' in metadata.columns:
        print("\nRegion distribution:")
        print(metadata['region'].value_counts().head(20))

    # Cell type levels
    cell_type_cols = ['class', 'subclass', 'supertype', 'cluster']
    available_cols = [c for c in cell_type_cols if c in metadata.columns]
    if available_cols:
        print(f"\nCell type hierarchy levels: {available_cols}")
        for col in available_cols:
            n_types = metadata[col].nunique()
            print(f"  {col}: {n_types} unique types")

    # Neuron vs non-neuron
    if 'class' in metadata.columns:
        is_neuron = (
            metadata["class"].str.contains("Glut", case=False, na=False) |
            metadata["class"].str.contains("GABA", case=False, na=False) |
            metadata["class"].str.contains("Dopa", case=False, na=False) |
            metadata["class"].str.contains("Sero", case=False, na=False)
        )
        print(f"\nNeuronal cells: {is_neuron.sum()} ({100*is_neuron.mean():.1f}%)")
        print(f"Non-neuronal cells: {(~is_neuron).sum()} ({100*(~is_neuron).mean():.1f}%)")


if __name__ == "__main__":
    metadata = extract_mouse_abc_metadata()
    print_mouse_abc_summary(metadata)
