"""Assign 3D coordinates to cells based on dataset type.

For mouse: Uses Allen CCF (Common Coordinate Framework) or synthetic centroids
For human: Uses spatial transcriptomics positions + slice ordering for z-axis
"""

import numpy as np
import pandas as pd
from pathlib import Path
from typing import Optional, Tuple
import re

import sys
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

try:
    from src.atlas.allen_api import get_atlas, ALLEN_SDK_AVAILABLE
    from src.atlas.region_mapping import HYPOMAP_TO_CCF
except ImportError:
    ALLEN_SDK_AVAILABLE = False
    HYPOMAP_TO_CCF = {}



def assign_coordinates(
    cells_df: pd.DataFrame,
    dataset: str,
    region_col: Optional[str] = None,
) -> pd.DataFrame:
    """Assign 3D coordinates based on dataset type.

    Args:
        cells_df: DataFrame with cell metadata
        dataset: Dataset name ("mouse_hypomap", "human_hypomap", or "mouse_abc")
        region_col: Name of region column

    Returns:
        DataFrame with added x, y, z columns
    """
    if dataset == "mouse_hypomap":
        return assign_mouse_coordinates(cells_df, region_col)
    elif dataset == "human_hypomap":
        return assign_human_coordinates(cells_df, region_col)
    elif dataset in ("mouse_abc", "mouse_abc_subcortical"):
        # mouse_abc datasets already have x, y, z from reconstructed_coordinates
        if 'x' in cells_df.columns and 'y' in cells_df.columns and 'z' in cells_df.columns:
            print(f"  {dataset}: Using existing coordinates from ABC reconstructed_coordinates")
            return cells_df
        else:
            raise ValueError(f"{dataset} dataset should already have x, y, z coordinates")
    else:
        raise ValueError(f"Unknown dataset: {dataset}")


# =============================================================================
# Mouse coordinate assignment (Allen CCF)
# =============================================================================

# Approximate centroids in CCF coordinates (in um) for fallback
MOUSE_REGION_CENTROIDS = {
    'Arcuate hypothalamic nucleus': (6500, 5000, 5700),
    'Ventromedial hypothalamic nucleus': (6200, 4800, 5700),
    'Paraventricular hypothalamic nucleus': (5500, 4200, 5700),
    'Dorsomedial nucleus of the hypothalamus': (6000, 4500, 5700),
    'Lateral hypothalamic area': (6000, 4500, 5300),
    'Suprachiasmatic nucleus': (5700, 5500, 5700),
    'Medial preoptic area': (5200, 4800, 5700),
    'Lateral preoptic area': (5200, 4800, 5300),
    'Tuberal nucleus': (6300, 5200, 5700),
    'Zona incerta': (6000, 4000, 5200),
    'Anterior hypothalamic nucleus': (5600, 4500, 5700),
    'Posterior hypothalamic nucleus': (6800, 4500, 5700),
    'Supraoptic nucleus': (5600, 5200, 5400),
    'Periventricular hypothalamic nucleus': (5800, 4500, 5700),
    'Periventricular hypothalamic nucleus, posterior part': (6200, 4500, 5700),
    'Periventricular hypothalamic nucleus, intermediate part': (6000, 4500, 5700),
    '(Pre)Mammillary region': (7000, 5000, 5700),
    '(Anterior/Preoptic)Periventricular region': (5300, 4500, 5700),
    'Mammillary body': (7000, 5000, 5700),
    'Preoptic area': (5100, 4800, 5700),
}

MOUSE_REGION_SIZES = {
    'Arcuate hypothalamic nucleus': 200,
    'Ventromedial hypothalamic nucleus': 300,
    'Paraventricular hypothalamic nucleus': 250,
    'Dorsomedial nucleus of the hypothalamus': 250,
    'Lateral hypothalamic area': 500,
    'Suprachiasmatic nucleus': 150,
    'Medial preoptic area': 300,
    'Lateral preoptic area': 250,
    'Tuberal nucleus': 200,
    'Zona incerta': 350,
    'Anterior hypothalamic nucleus': 200,
    'Posterior hypothalamic nucleus': 200,
    'Supraoptic nucleus': 150,
    'Periventricular hypothalamic nucleus': 200,
    'Periventricular hypothalamic nucleus, posterior part': 200,
    'Periventricular hypothalamic nucleus, intermediate part': 150,
    '(Pre)Mammillary region': 300,
    '(Anterior/Preoptic)Periventricular region': 200,
    'Mammillary body': 250,
    'Preoptic area': 350,
}


def assign_mouse_coordinates(
    cells_df: pd.DataFrame,
    region_col: Optional[str] = 'Region',
) -> pd.DataFrame:
    """Assign 3D CCF coordinates to mouse cells."""
    if not ALLEN_SDK_AVAILABLE or not HYPOMAP_TO_CCF:
        print("Allen SDK not available, using synthetic coordinates")
        return assign_mouse_synthetic_coordinates(cells_df, region_col)

    atlas = get_atlas(resolution=25)

    # Pre-compute region coordinates
    region_coords = {}
    for region_name, (acronym, ccf_id) in HYPOMAP_TO_CCF.items():
        coords = atlas.get_structure_coords(ccf_id)
        if coords is not None and len(coords) > 0:
            region_coords[region_name] = coords
            print(f"  {region_name}: {len(coords)} voxels")
        else:
            print(f"  Warning: No coords for {region_name} (ID: {ccf_id})")

    # Assign coordinates
    x_coords, y_coords, z_coords = [], [], []

    for idx, row in cells_df.iterrows():
        region = row[region_col] if region_col else None

        if region in region_coords:
            coords = region_coords[region]
            voxel_idx = np.random.randint(len(coords))
            voxel = coords[voxel_idx]
            jitter = np.random.uniform(-0.5, 0.5, 3)
            coord = (voxel + jitter) * atlas.resolution

            z_coords.append(coord[0])  # AP
            y_coords.append(coord[1])  # DV
            x_coords.append(coord[2])  # ML
        else:
            x_coords.append(np.nan)
            y_coords.append(np.nan)
            z_coords.append(np.nan)

    cells_df = cells_df.copy()
    cells_df['x'] = x_coords
    cells_df['y'] = y_coords
    cells_df['z'] = z_coords

    return cells_df


def assign_mouse_synthetic_coordinates(
    cells_df: pd.DataFrame,
    region_col: Optional[str] = 'Region',
) -> pd.DataFrame:
    """Assign synthetic coordinates when Allen SDK is not available."""
    x_coords, y_coords, z_coords = [], [], []

    for idx, row in cells_df.iterrows():
        region = row[region_col] if region_col else None

        if region in MOUSE_REGION_CENTROIDS:
            centroid = MOUSE_REGION_CENTROIDS[region]
            size = MOUSE_REGION_SIZES.get(region, 200)
            jitter = np.random.normal(0, size / 2, 3)
            coord = np.array(centroid) + jitter

            z_coords.append(coord[0])  # AP
            y_coords.append(coord[1])  # DV
            x_coords.append(coord[2])  # ML
        else:
            x_coords.append(np.nan)
            y_coords.append(np.nan)
            z_coords.append(np.nan)

    cells_df = cells_df.copy()
    cells_df['x'] = x_coords
    cells_df['y'] = y_coords
    cells_df['z'] = z_coords

    return cells_df


# =============================================================================
# Human coordinate assignment (UMAP scaled to micrometers + section-based z)
# =============================================================================

# Section thickness in micrometers (FFPE sections are 5 µm per Tadross et al. 2025 Methods)
SECTION_THICKNESS_UM = 5.0

# Inter-section spacing in micrometers
# NOTE: The actual spacing between captured sections is not specified in the paper.
# Using 200 µm as a reasonable estimate for visualization purposes.
# The paper captured 9 sections spanning anterior-posterior hypothalamus (~10-15mm),
# so actual spacing could be ~1-2mm between sections.
INTER_SECTION_SPACING_UM = 200.0

# Target tissue dimensions for UMAP scaling (human hypothalamus is ~5mm in xy)
TISSUE_SIZE_UM = 5000.0  # 5mm = 5000 micrometers


def extract_section_number(sample_id: str) -> Tuple[str, Optional[int]]:
    """Extract donor ID and section number from sample ID.

    Args:
        sample_id: Sample ID like "3u5kk_A.01", "znZv1_S10", "siletti_10X389_2"

    Returns:
        Tuple of (donor_id, section_number or None)
    """
    # Pattern: "donor_type.number" (e.g., "3u5kk_A.01")
    match = re.match(r'^([^_]+)_[AR]\.(\d+)$', sample_id)
    if match:
        return match.group(1), int(match.group(2))

    # Pattern: "donor_typeNumber" (e.g., "znZv1_S10", "1C41i_A4")
    match = re.match(r'^([^_]+)_[A-Z]+(\d+)$', sample_id)
    if match:
        return match.group(1), int(match.group(2))

    # Pattern: "donor_code_number" (e.g., "siletti_10X389_2")
    match = re.match(r'^(.+)_(\d+)$', sample_id)
    if match:
        return match.group(1), int(match.group(2))

    # Fallback: no section number detected
    return sample_id, None


def assign_human_coordinates(
    cells_df: pd.DataFrame,
    region_col: Optional[str] = 'region',
) -> pd.DataFrame:
    """Assign 3D coordinates to human cells in micrometers.

    Strategy:
    - x, y: UMAP coordinates scaled to anatomically realistic tissue dimensions (~5mm)
    - z: Based on section number extracted from Sample_ID (with inter-section spacing)

    All coordinates are in micrometers (µm).
    """
    cells_df = cells_df.copy()

    # Check required columns
    if 'umap_1' not in cells_df.columns or 'umap_2' not in cells_df.columns:
        raise ValueError("UMAP coordinates (umap_1, umap_2) required for human dataset")

    # Calculate UMAP scaling to convert to micrometers
    # Scale UMAP range to target tissue size
    umap_range_x = cells_df['umap_1'].max() - cells_df['umap_1'].min()
    umap_range_y = cells_df['umap_2'].max() - cells_df['umap_2'].min()
    umap_range = max(umap_range_x, umap_range_y)

    # Scale factor: UMAP units -> micrometers
    scale_factor = TISSUE_SIZE_UM / umap_range if umap_range > 0 else 1.0

    print(f"  Scaling UMAP coordinates to micrometers (scale: {scale_factor:.2f} µm/unit)")
    print(f"  Target tissue size: {TISSUE_SIZE_UM:.0f} µm ({TISSUE_SIZE_UM/1000:.1f} mm)")

    # Assign x, y from scaled UMAP (centered at 0)
    umap_center_x = (cells_df['umap_1'].max() + cells_df['umap_1'].min()) / 2
    umap_center_y = (cells_df['umap_2'].max() + cells_df['umap_2'].min()) / 2

    cells_df['x'] = (cells_df['umap_1'] - umap_center_x) * scale_factor
    cells_df['y'] = (cells_df['umap_2'] - umap_center_y) * scale_factor

    # Assign z based on section numbers from Sample_ID
    cells_df['z'] = 0.0

    if 'Sample_ID' in cells_df.columns:
        # Extract donor and section info
        donors = []
        sections = []
        for sample_id in cells_df['Sample_ID']:
            donor, section = extract_section_number(sample_id)
            donors.append(donor)
            sections.append(section)
        cells_df['_donor'] = donors
        cells_df['_section'] = sections

        # For each donor, assign z based on section number
        # Sections are spaced by INTER_SECTION_SPACING_UM
        for donor in cells_df['_donor'].unique():
            donor_mask = cells_df['_donor'] == donor
            donor_sections = cells_df.loc[donor_mask, '_section']

            # Get section numbers (use 0 for cells without section info)
            valid_sections = donor_sections.dropna()

            if len(valid_sections) > 0:
                # Get unique sections for this donor and map to z positions
                unique_sections = sorted(valid_sections.unique())
                section_to_z = {s: i * INTER_SECTION_SPACING_UM for i, s in enumerate(unique_sections)}

                # Assign z coordinates
                for idx in cells_df[donor_mask].index:
                    section = cells_df.at[idx, '_section']
                    if section is not None and section in section_to_z:
                        z_base = section_to_z[section]
                    else:
                        z_base = 0.0

                    # Add jitter within section thickness
                    cells_df.at[idx, 'z'] = z_base + np.random.uniform(0, SECTION_THICKNESS_UM)

        # Report section statistics
        n_with_section = cells_df['_section'].notna().sum()
        n_sections = cells_df['_section'].nunique()
        n_donors = cells_df['_donor'].nunique()
        print(f"  Found {n_sections} sections across {n_donors} donors")
        print(f"  Cells with section info: {n_with_section}/{len(cells_df)} ({100*n_with_section/len(cells_df):.1f}%)")

        # Clean up temporary columns
        cells_df = cells_df.drop(columns=['_donor', '_section'])
    else:
        print("  Warning: No Sample_ID column, using flat z=0")

    # Add small jitter to x, y to avoid overlapping points
    jitter_um = 5.0  # 5 micrometer jitter
    cells_df['x'] += np.random.uniform(-jitter_um, jitter_um, len(cells_df))
    cells_df['y'] += np.random.uniform(-jitter_um, jitter_um, len(cells_df))

    # Report final coordinate ranges
    print("  Final coordinate ranges (µm):")
    print(f"    x: {cells_df['x'].min():.0f} to {cells_df['x'].max():.0f}")
    print(f"    y: {cells_df['y'].min():.0f} to {cells_df['y'].max():.0f}")
    print(f"    z: {cells_df['z'].min():.0f} to {cells_df['z'].max():.0f}")

    return cells_df


if __name__ == "__main__":
    # Test coordinate assignment
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", "-d", default="mouse_hypomap", choices=["mouse_hypomap", "human_hypomap", "mouse_abc"])
    args = parser.parse_args()

    # Create test dataframe
    if args.dataset == "mouse_hypomap":
        test_df = pd.DataFrame({
            'Region': ['Arcuate hypothalamic nucleus', 'Ventromedial hypothalamic nucleus'] * 5
        })
    else:
        test_df = pd.DataFrame({
            'region': ['ARC', 'VMH', 'NA'] * 5,
            'umap_1': np.random.randn(15),
            'umap_2': np.random.randn(15),
            'Sample_ID': ['test_sample'] * 15,
        })

    result = assign_coordinates(test_df, args.dataset)
    print("\nTest result:")
    print(result[['x', 'y', 'z']].describe())
