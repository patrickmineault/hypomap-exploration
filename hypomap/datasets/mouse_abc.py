"""Allen Brain Cell Census (ABC) dataset adapter.

Handles extraction and processing of the ABC MERFISH data.
Uses AbcProjectCache to load data from the Allen Institute cloud storage.
Supports filtering to hypothalamus-only or expanded subcortical divisions.
"""

import argparse
import pandas as pd
from pathlib import Path
from typing import Optional

from .base import DatasetConfig

# Paths
DATA_DIR = Path(__file__).parent.parent.parent / "data"
DOWNLOAD_BASE = DATA_DIR / "raw" / "abc_atlas_cache"
PROCESSED_DIR = DATA_DIR / "processed" / "mouse_abc"
EXTENDED_PROCESSED_DIR = DATA_DIR / "processed" / "mouse_abc_extended"

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


# Additional region colors for extended structures outside hypothalamus
ABC_EXTENDED_REGION_COLORS = {
    **ABC_REGION_COLORS,
    # Thalamus (TH) regions
    'AD': '#1E90FF',      # Anterodorsal nucleus
    'AM': '#4682B4',      # Anteromedial nucleus
    'AV': '#5F9EA0',      # Anteroventral nucleus
    'CM': '#6495ED',      # Central medial nucleus
    'IAD': '#7B68EE',     # Interanterodorsal nucleus
    'IMD': '#6A5ACD',     # Intermediodorsal nucleus
    'LD': '#4169E1',      # Lateral dorsal nucleus
    'LGd': '#0000CD',     # Lateral geniculate, dorsal
    'LGv': '#0000FF',     # Lateral geniculate, ventral
    'LP': '#1E90FF',      # Lateral posterior nucleus
    'MD': '#00BFFF',      # Mediodorsal nucleus
    'MG': '#87CEEB',      # Medial geniculate
    'PCN': '#4682B4',     # Paracentral nucleus
    'PF': '#B0C4DE',      # Parafascicular nucleus
    'PO': '#5B9BD5',      # Posterior thalamic nuclear group
    'PT': '#2F4F4F',      # Parataenial nucleus
    'PVT': '#008080',     # Paraventricular thalamic nucleus
    'RE': '#20B2AA',      # Nucleus reuniens
    'RH': '#48D1CC',      # Rhomboid nucleus
    'RT': '#40E0D0',      # Reticular nucleus
    'SGN': '#00CED1',     # Suprageniculate nucleus
    'SMT': '#5F9EA0',     # Submedial thalamic nucleus
    'SPFm': '#7FFFD4',    # Subparafascicular nucleus, magnocellular
    'SPFp': '#66CDAA',    # Subparafascicular nucleus, parvicellular
    'VAL': '#3CB371',     # Ventral anterior-lateral
    'VM': '#2E8B57',      # Ventral medial nucleus
    'VPL': '#006400',     # Ventral posterolateral
    'VPLpc': '#228B22',   # VPL, parvicellular
    'VPM': '#32CD32',     # Ventral posteromedial
    'VPMpc': '#00FF00',   # VPM, parvicellular
    'Xi': '#98FB98',      # Xiphoid thalamic nucleus
    # Striatum (STR) regions
    'ACB': '#FF7F50',     # Nucleus accumbens
    'CP': '#FF6347',      # Caudoputamen
    'FS': '#FF4500',      # Fundus of striatum
    'OT': '#FFD700',      # Olfactory tubercle
    'LSr': '#FFA07A',     # Lateral septal nucleus, rostral
    'LSv': '#FF8C69',     # Lateral septal nucleus, ventral
    'LSc': '#E9967A',     # Lateral septal nucleus, caudal
    'CEA': '#CD5C5C',     # Central amygdalar nucleus
    'AAA': '#BC8F8F',     # Anterior amygdalar area
    'IA': '#F08080',      # Intercalated amygdalar nucleus
    'SF': '#DEB887',      # Septofimbrial nucleus
    'SH': '#D2B48C',      # Septohippocampal nucleus
    'BA': '#C4A882',      # Bed nucleus of anterior commissure
    # Pallidum (PAL) regions
    'BST': '#DDA0DD',     # Bed nuclei of the stria terminalis
    'GPe': '#DA70D6',     # Globus pallidus, external
    'GPi': '#BA55D3',     # Globus pallidus, internal
    'MS': '#9370DB',      # Medial septal nucleus
    'NDB': '#8A2BE2',     # Diagonal band nucleus
    'SI': '#9400D3',      # Substantia innominata
    'MA': '#800080',      # Magnocellular nucleus
    'TRS': '#C71585',     # Triangular nucleus of septum
    # Pons (P) and Medulla (MY) - generic unassigned
    'P-unassigned': '#666666',
    'MY-unassigned': '#666666',
    # Midbrain (MB) regions
    'MB-unassigned': '#666666',
    'SCm': '#2E8B57',     # Superior colliculus, motor
    'SCs': '#3CB371',     # Superior colliculus, sensory
    'IC': '#008B8B',      # Inferior colliculus
    'PAG': '#9932CC',     # Periaqueductal gray
    'MRN': '#708090',     # Midbrain reticular nucleus
    'APN': '#4682B4',     # Anterior pretectal nucleus
    'SNr': '#8B0000',     # Substantia nigra, reticular part
    'SNc': '#B22222',     # Substantia nigra, compact part
    'VTA': '#DC143C',     # Ventral tegmental area
    'PPN': '#D2691E',     # Pedunculopontine nucleus
    'CUN': '#A0522D',     # Cuneiform nucleus
    'RN': '#CD5C5C',      # Red nucleus
    'IPN': '#9370DB',     # Interpeduncular nucleus
    'NPC': '#6A5ACD',     # Nucleus of posterior commissure
    'NOT': '#483D8B',     # Nucleus of the optic tract
    'DR': '#FF8C00',      # Dorsal raphe nucleus
    'PPT': '#5F9EA0',     # Posterior pretectal nucleus
    'RR': '#BC8F8F',      # Retrorubral area
    'NB': '#4169E1',      # Nucleus of brachium of IC
    'SAG': '#7B68EE',     # Nucleus sagulum
    'IF': '#FFA500',      # Interfascicular nucleus raphe
    'CLI': '#FFD700',     # Central linear nucleus raphe
    'PBG': '#32CD32',     # Parabigeminal nucleus
    'AT': '#66CDAA',      # Anterior tegmental nucleus
    'VTN': '#20B2AA',     # Ventral tegmental nucleus
    'RL': '#DAA520',      # Rostral linear nucleus raphe
    'MPT': '#6495ED',     # Medial pretectal area
    'MT': '#00CED1',      # Medial terminal nucleus of accessory optic tract
    'OP': '#1E90FF',      # Olivary pretectal nucleus
    'LT': '#48D1CC',      # Lateral terminal nucleus of accessory optic tract
    'III': '#FF6347',     # Oculomotor nucleus
    'RPF': '#8FBC8F',     # Retroparafascicular nucleus
    'SCO': '#778899',     # Subcommissural organ
    'PN': '#DB7093',      # Paranigral nucleus
    'Pa4': '#F4A460',     # Paratrochlear nucleus
    'EW': '#FA8072',      # Edinger-Westphal nucleus
    'MA3': '#E9967A',     # Medial accessory oculomotor nucleus
    'DT': '#87CEEB',      # Dorsal terminal nucleus of accessory optic tract
    'IV': '#FF7F50',      # Trochlear nucleus
    # Generic division-unassigned colors
    'TH-unassigned': '#666666',
    'STR-unassigned': '#666666',
    'PAL-unassigned': '#666666',
}


def get_mouse_abc_extended_config() -> DatasetConfig:
    """Get configuration for the mouse ABC extended dataset (HY + TH + STR + PAL + P + MY + MB)."""
    return DatasetConfig(
        name="mouse_abc_extended",
        species="Mus musculus",
        h5ad_path=DOWNLOAD_BASE,
        processed_dir=EXTENDED_PROCESSED_DIR,
        cell_type_columns=['class', 'subclass', 'supertype', 'cluster'],
        region_column="region",
        gene_column=None,
        marker_genes=ABC_MARKER_GENES,
        key_receptors=ABC_KEY_RECEPTORS,
        region_colors=ABC_EXTENDED_REGION_COLORS,
    )


def classify_neurons(metadata: pd.DataFrame) -> pd.Series:
    """Classify cells as neuronal or non-neuronal.

    For ABC Atlas, neurons are identified by 'class' column containing
    neurotransmitter-related terms (Glutamatergic, GABAergic, Dopaminergic, Serotonergic).

    Args:
        metadata: Cell metadata DataFrame

    Returns:
        Boolean Series where True indicates neuronal cells
    """
    if 'class' in metadata.columns:
        is_neuron = (
            metadata['class'].str.contains('Glut', case=False, na=False) |
            metadata['class'].str.contains('GABA', case=False, na=False) |
            metadata['class'].str.contains('Dopa', case=False, na=False) |
            metadata['class'].str.contains('Sero', case=False, na=False) |
            metadata['class'].str.contains('Chol', case=False, na=False)  # Cholinergic
        )
        return is_neuron

    # Fallback: check neurotransmitter column
    if 'neurotransmitter' in metadata.columns:
        is_neuron = metadata['neurotransmitter'].notna() & (metadata['neurotransmitter'] != '')
        return is_neuron

    # Default: unknown (all False)
    return pd.Series(False, index=metadata.index)


def extract_mouse_abc_metadata(
    cache_dir: Optional[Path] = None,
    divisions: Optional[list[str]] = None,
    output_dir: Optional[Path] = None,
) -> pd.DataFrame:
    """Extract cell metadata from ABC MERFISH data, filtered to specified divisions.

    Args:
        cache_dir: Path to ABC cache directory. If None, uses default.
        divisions: List of parcellation_division values to include (default: ['HY']).
        output_dir: Output directory for processed files. If None, uses default.

    Returns:
        DataFrame with cell metadata for filtered cells
    """
    from abc_atlas_access.abc_atlas_cache.abc_project_cache import AbcProjectCache

    if cache_dir is None:
        cache_dir = DOWNLOAD_BASE
    if divisions is None:
        divisions = ['HY']
    if output_dir is None:
        output_dir = PROCESSED_DIR

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
    # Both dataframes have cell_label as index, so join on index directly
    cell = cell.join(
        recon[["parcellation_index", "x_reconstructed", "y_reconstructed", "z_reconstructed"]],
        how="inner"
    )

    # Join parcellation terms
    cell = cell.join(parc_terms, on="parcellation_index", how="left")

    print(f"After joining: {len(cell)} cells")

    # 5. Filter to requested divisions
    print(f"Filtering to divisions: {divisions}...")
    div_mask = cell['parcellation_division'].isin(divisions)
    cell = cell[div_mask].copy()
    print(f"Filtered cells: {len(cell)}")

    # 6. Rename/standardize columns for hypomap compatibility
    # Drop original section-local coordinates (we use CCF-registered reconstructed coords)
    cols_to_drop = [c for c in ['x', 'y', 'z'] if c in cell.columns]
    if cols_to_drop:
        cell = cell.drop(columns=cols_to_drop)

    cell = cell.rename(columns={
        'cell_label': 'cell_id',
        'x_reconstructed': 'x',
        'y_reconstructed': 'y',
        'z_reconstructed': 'z',
        'parcellation_structure': 'region',
        'donor_label': 'donor_id',
        'brain_section_label': 'sample_id',
    })

    # 7. Classify neurons vs non-neurons
    cell['is_neuron'] = classify_neurons(cell)
    n_neurons = cell['is_neuron'].sum()
    print(f"Classified {n_neurons} neurons ({100*n_neurons/len(cell):.1f}%)")

    # Ensure output directory exists
    output_dir.mkdir(parents=True, exist_ok=True)

    # 8. Save cell metadata
    metadata_path = output_dir / "cell_metadata.parquet"
    cell.to_parquet(metadata_path)
    print(f"Saved metadata for {len(cell)} cells to {metadata_path}")

    # 9. Save gene list (from ABC gene table)
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

    genes_path = output_dir / "genes.parquet"
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
    if 'is_neuron' in metadata.columns:
        is_neuron = metadata['is_neuron']
        print(f"\nNeuronal cells: {is_neuron.sum()} ({100*is_neuron.mean():.1f}%)")
        print(f"Non-neuronal cells: {(~is_neuron).sum()} ({100*(~is_neuron).mean():.1f}%)")
    elif 'class' in metadata.columns:
        is_neuron = classify_neurons(metadata)
        print(f"\nNeuronal cells: {is_neuron.sum()} ({100*is_neuron.mean():.1f}%)")
        print(f"Non-neuronal cells: {(~is_neuron).sum()} ({100*(~is_neuron).mean():.1f}%)")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract cell metadata from ABC MERFISH data."
    )
    parser.add_argument(
        "--divisions",
        nargs="+",
        default=["HY"],
        help="Parcellation divisions to include (default: HY)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Output directory for processed files (default: data/processed/mouse_abc)",
    )
    args = parser.parse_args()

    metadata = extract_mouse_abc_metadata(
        divisions=args.divisions,
        output_dir=args.output_dir,
    )
    print_mouse_abc_summary(metadata)
