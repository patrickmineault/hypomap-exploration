"""Mapping between HypoMap regions and Allen CCF structures.

This module provides region mappings, colors, and gene lists for both
mouse and human hypothalamus datasets.
"""

from typing import Dict, List, Optional

# =============================================================================
# Mouse HypoMap Mappings
# =============================================================================

# Mapping from HypoMap region names to Allen CCF structure IDs
# Format: 'HypoMap region name': (CCF acronym, CCF structure ID)
HYPOMAP_TO_CCF = {
    'Arcuate hypothalamic nucleus': ('ARH', 223),
    'Ventromedial hypothalamic nucleus': ('VMH', 693),
    'Paraventricular hypothalamic nucleus': ('PVH', 38),
    'Dorsomedial nucleus of the hypothalamus': ('DMH', 830),
    'Lateral hypothalamic area': ('LHA', 194),
    'Suprachiasmatic nucleus': ('SCH', 286),
    'Medial preoptic area': ('MPO', 523),
    'Lateral preoptic area': ('LPO', 226),
    'Tuberal nucleus': ('TU', 557),
    'Zona incerta': ('ZI', 797),
    'Anterior hypothalamic nucleus': ('AHN', 88),
    'Posterior hypothalamic nucleus': ('PH', 946),
    'Supraoptic nucleus': ('SO', 390),
    'Periventricular hypothalamic nucleus': ('PVZ', 141),
    'Periventricular hypothalamic nucleus, posterior part': ('PVpo', 149),
    'Periventricular hypothalamic nucleus, intermediate part': ('PVi', 156),
    '(Pre)Mammillary region': ('MBO', 331),
    '(Anterior/Preoptic)Periventricular region': ('AVPV', 272),
    'Mammillary body': ('MBO', 331),
    'Preoptic area': ('POA', 803),
}

# Consistent color scheme for mouse regions (saturated colors for visibility)
MOUSE_REGION_COLORS = {
    'Arcuate hypothalamic nucleus': '#DC143C',  # Crimson
    'Ventromedial hypothalamic nucleus': '#0066CC',  # Strong blue
    'Paraventricular hypothalamic nucleus': '#228B22',  # Forest green
    'Dorsomedial nucleus of the hypothalamus': '#8B008B',  # Dark magenta
    'Lateral hypothalamic area': '#FF6600',  # Vivid orange
    'Suprachiasmatic nucleus': '#FFD700',  # Gold
    'Medial preoptic area': '#8B4513',  # Saddle brown
    'Lateral preoptic area': '#FF1493',  # Deep pink
    'Tuberal nucleus': '#696969',  # Dim gray
    'Zona incerta': '#008B8B',  # Dark cyan
    'Anterior hypothalamic nucleus': '#FF4500',  # Orange red
    'Posterior hypothalamic nucleus': '#4169E1',  # Royal blue
    'Supraoptic nucleus': '#C71585',  # Medium violet red
    'Periventricular hypothalamic nucleus': '#32CD32',  # Lime green
    'Periventricular hypothalamic nucleus, posterior part': '#00CC66',  # Spring green
    'Periventricular hypothalamic nucleus, intermediate part': '#66CDAA',  # Medium aquamarine
    '(Pre)Mammillary region': '#DAA520',  # Goldenrod
    '(Anterior/Preoptic)Periventricular region': '#9932CC',  # Dark orchid
    'Mammillary body': '#DAA520',  # Goldenrod
    'Preoptic area': '#CD853F',  # Peru
    'NA': '#A0A0A0',  # Gray for unknown
}

# Backward compatibility alias
REGION_COLORS = MOUSE_REGION_COLORS

# =============================================================================
# Human HypoMap Mappings
# =============================================================================

# Human region abbreviations to full names
HUMAN_REGION_NAMES = {
    'ARC': 'Arcuate nucleus',
    'VMH': 'Ventromedial hypothalamus',
    'PVN': 'Paraventricular nucleus',
    'DMH': 'Dorsomedial hypothalamus',
    'LH': 'Lateral hypothalamus',
    'SCN': 'Suprachiasmatic nucleus',
    'MPOA': 'Medial preoptic area',
    'POA': 'Preoptic area',
    'LPOA': 'Lateral preoptic area',
    'TMN': 'Tuberomammillary nucleus',
    'SON': 'Supraoptic nucleus',
    'ME': 'Median eminence',
    'Fx/OT/ac': 'Fornix/Optic tract/Anterior commissure',
    'MAM': 'Mammillary region',
    'Perivent': 'Periventricular region',
    'LTN': 'Lateral tuberal nucleus',
    'Vent': 'Ventricle',
    'Vascular': 'Vascular',
    'Thalamus': 'Thalamus',
    'Thalamaus': 'Thalamus',  # Handle typo in dataset
    'NA': 'Unknown',
}

# Color scheme for human regions
HUMAN_REGION_COLORS = {
    'ARC': '#DC143C',           # Crimson
    'VMH': '#0066CC',           # Strong blue
    'PVN': '#228B22',           # Forest green
    'DMH': '#8B008B',           # Dark magenta
    'LH': '#FF6600',            # Vivid orange
    'SCN': '#FFD700',           # Gold
    'MPOA': '#8B4513',          # Saddle brown
    'POA': '#CD853F',           # Peru
    'LPOA': '#FF1493',          # Deep pink
    'TMN': '#FF4500',           # Orange red
    'SON': '#C71585',           # Medium violet red
    'ME': '#4169E1',            # Royal blue
    'Fx/OT/ac': '#696969',      # Dim gray
    'MAM': '#DAA520',           # Goldenrod
    'Perivent': '#32CD32',      # Lime green
    'LTN': '#008B8B',           # Dark cyan
    'Vent': '#A9A9A9',          # Dark gray
    'Vascular': '#B22222',      # Firebrick
    'Thalamus': '#9932CC',      # Dark orchid
    'Thalamaus': '#9932CC',     # Dark orchid (typo variant)
    'NA': '#808080',            # Gray for unknown
}

# =============================================================================
# Mouse Gene Lists
# =============================================================================

# Key receptors for hypothalamic function (mouse gene symbols)
MOUSE_KEY_RECEPTORS = [
    'Lepr',   # Leptin receptor
    'Mc4r',   # Melanocortin 4 receptor
    'Glp1r',  # GLP-1 receptor
    'Ghsr',   # Ghrelin receptor
    'Insr',   # Insulin receptor
    'Npy1r',  # NPY receptor Y1
    'Npy2r',  # NPY receptor Y2
    'Crhr1',  # CRH receptor 1
    'Oxtr',   # Oxytocin receptor
    'Avpr1a', # Vasopressin receptor 1A
]

# Marker genes commonly used in hypothalamus studies (mouse)
MOUSE_MARKER_GENES = [
    'Agrp',   # AgRP neurons
    'Pomc',   # POMC neurons
    'Npy',    # NPY neurons
    'Sf1',    # VMH marker (Nr5a1)
    'Oxt',    # Oxytocin
    'Avp',    # Vasopressin
    'Crh',    # CRH
    'Trh',    # TRH
    'Ghrh',   # GHRH
    'Sst',    # Somatostatin
]

# Backward compatibility aliases
KEY_RECEPTORS = MOUSE_KEY_RECEPTORS
MARKER_GENES = MOUSE_MARKER_GENES

# =============================================================================
# Human Gene Lists (uppercase symbols)
# =============================================================================

HUMAN_KEY_RECEPTORS = [
    'LEPR',   # Leptin receptor
    'MC4R',   # Melanocortin 4 receptor
    'GLP1R',  # GLP-1 receptor
    'GHSR',   # Ghrelin receptor
    'INSR',   # Insulin receptor
    'NPY1R',  # NPY receptor Y1
    'NPY2R',  # NPY receptor Y2
    'CRHR1',  # CRH receptor 1
    'OXTR',   # Oxytocin receptor
    'AVPR1A', # Vasopressin receptor 1A
]

HUMAN_MARKER_GENES = [
    'AGRP',   # AgRP neurons
    'POMC',   # POMC neurons
    'NPY',    # NPY neurons
    'NR5A1',  # VMH marker (SF1)
    'OXT',    # Oxytocin
    'AVP',    # Vasopressin
    'CRH',    # CRH
    'TRH',    # TRH
    'GHRH',   # GHRH
    'SST',    # Somatostatin
]


# =============================================================================
# Utility Functions
# =============================================================================

def get_ccf_id(hypomap_region: str) -> Optional[int]:
    """Get Allen CCF structure ID for a mouse HypoMap region."""
    if hypomap_region in HYPOMAP_TO_CCF:
        return HYPOMAP_TO_CCF[hypomap_region][1]
    return None


def get_ccf_acronym(hypomap_region: str) -> Optional[str]:
    """Get Allen CCF acronym for a mouse HypoMap region."""
    if hypomap_region in HYPOMAP_TO_CCF:
        return HYPOMAP_TO_CCF[hypomap_region][0]
    return None


def get_region_color(region_name: str, dataset: str = "mouse_hypomap") -> str:
    """Get color for a region, with fallback.

    Args:
        region_name: Name of the region
        dataset: Dataset name ("mouse_hypomap", "human_hypomap", "mouse_abc", etc.)

    Returns:
        Hex color string
    """
    if dataset.startswith("human"):
        return HUMAN_REGION_COLORS.get(region_name, '#CCCCCC')
    return MOUSE_REGION_COLORS.get(region_name, '#CCCCCC')


def get_region_colors_for_dataset(dataset: str) -> Dict[str, str]:
    """Get all region colors for a dataset.

    Args:
        dataset: Dataset name ("mouse_hypomap", "human_hypomap", "mouse_abc", etc.)

    Returns:
        Dictionary mapping region names to colors
    """
    if dataset.startswith("human"):
        return HUMAN_REGION_COLORS.copy()
    return MOUSE_REGION_COLORS.copy()


def get_marker_genes_for_dataset(dataset: str) -> List[str]:
    """Get marker genes for a dataset.

    Args:
        dataset: Dataset name ("mouse_hypomap", "human_hypomap", "mouse_abc", etc.)

    Returns:
        List of marker gene symbols
    """
    if dataset.startswith("human"):
        return HUMAN_MARKER_GENES.copy()
    return MOUSE_MARKER_GENES.copy()


def get_key_receptors_for_dataset(dataset: str) -> List[str]:
    """Get key receptors for a dataset.

    Args:
        dataset: Dataset name ("mouse_hypomap", "human_hypomap", "mouse_abc", etc.)

    Returns:
        List of receptor gene symbols
    """
    if dataset.startswith("human"):
        return HUMAN_KEY_RECEPTORS.copy()
    return MOUSE_KEY_RECEPTORS.copy()


def normalize_region_name(name: str, dataset: str = "mouse_hypomap") -> str:
    """Normalize region name for matching.

    Args:
        name: Region name or abbreviation
        dataset: Dataset name ("mouse_hypomap", "human_hypomap", "mouse_abc", etc.)

    Returns:
        Normalized region name
    """
    name = name.strip()

    if dataset.startswith("human"):
        # Human regions are already abbreviated, return full name if available
        return HUMAN_REGION_NAMES.get(name, name)

    # Mouse: map abbreviations to full names
    mouse_variations = {
        'ARH': 'Arcuate hypothalamic nucleus',
        'VMH': 'Ventromedial hypothalamic nucleus',
        'PVH': 'Paraventricular hypothalamic nucleus',
        'PVN': 'Paraventricular hypothalamic nucleus',
        'DMH': 'Dorsomedial nucleus of the hypothalamus',
        'LHA': 'Lateral hypothalamic area',
        'LH': 'Lateral hypothalamic area',
        'SCN': 'Suprachiasmatic nucleus',
        'SCH': 'Suprachiasmatic nucleus',
        'MPOA': 'Medial preoptic area',
        'MPO': 'Medial preoptic area',
        'ARC': 'Arcuate hypothalamic nucleus',
    }

    return mouse_variations.get(name.upper(), name)


def get_human_region_full_name(abbrev: str) -> str:
    """Get full name for a human region abbreviation.

    Args:
        abbrev: Region abbreviation (e.g., "ARC")

    Returns:
        Full region name
    """
    return HUMAN_REGION_NAMES.get(abbrev, abbrev)
