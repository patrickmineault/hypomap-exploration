"""Extract human Visium spatial transcriptomics data.

This creates the spatial dataset for human hypothalamus visualization using
true spatial coordinates from 10x Visium data (not UMAP from snRNA-seq).
"""

import scanpy as sc
import pandas as pd
from pathlib import Path

DATA_DIR = Path(__file__).parent.parent.parent / "data"
RAW_DIR = DATA_DIR / "raw" / "human_hypomap"
PROCESSED_DIR = DATA_DIR / "processed" / "human_hypomap"

VISIUM_H5_DIR = RAW_DIR / "visium_h5"
SPATIAL_POSITIONS_FILE = PROCESSED_DIR / "visium_spots.parquet"


def extract_visium_data():
    """Extract and combine all Visium sections with spatial coordinates."""

    # Load spatial positions
    spots_df = pd.read_parquet(SPATIAL_POSITIONS_FILE)
    print(f"Loaded spatial positions: {len(spots_df)} spots from {spots_df['sample_id'].nunique()} samples")

    all_spots = []
    all_genes = None

    for h5_file in sorted(VISIUM_H5_DIR.glob("*.h5")):
        # Extract sample ID from filename
        # e.g., GSM8555882_OgTiF_b6_filtered_feature_bc_matrix.h5 -> OgTiF_b6
        sample_id = h5_file.stem.split('_', 1)[1].replace('_filtered_feature_bc_matrix', '')

        print(f"\nProcessing {sample_id}...")

        # Load expression data
        adata = sc.read_10x_h5(h5_file)
        adata.var_names_make_unique()

        # Get sample spatial coords
        sample_spots = spots_df[spots_df['sample_id'] == sample_id].copy()
        if len(sample_spots) == 0:
            print(f"  WARNING: No spatial coords for {sample_id}, skipping")
            continue

        # Create spot dataframe with coords
        sample_spots = sample_spots.set_index('barcode')

        # Match barcodes
        common_barcodes = sample_spots.index.intersection(adata.obs_names)
        print(f"  {len(common_barcodes)} spots with both expression and coords")

        # Filter to common barcodes
        adata = adata[common_barcodes].copy()
        sample_spots = sample_spots.loc[common_barcodes]

        # Build spot metadata
        spot_data = pd.DataFrame({
            'spot_id': [f"{sample_id}_{bc}" for bc in common_barcodes],
            'sample_id': sample_id,
            'barcode': common_barcodes.tolist(),
            'x': sample_spots['x'].values,
            'y': sample_spots['y'].values,
            'z': sample_spots['z'].values,
        })

        all_spots.append(spot_data)

        # Track genes (use first sample as reference)
        if all_genes is None:
            all_genes = adata.var_names.tolist()

    # Combine all spots
    combined_spots = pd.concat(all_spots, ignore_index=True)
    combined_spots = combined_spots.set_index('spot_id')

    print("\n=== Combined Data ===")
    print(f"Total spots: {len(combined_spots)}")
    print(f"Samples: {combined_spots['sample_id'].nunique()}")
    print("Coordinate ranges:")
    print(f"  x: {combined_spots['x'].min():.1f} to {combined_spots['x'].max():.1f} µm")
    print(f"  y: {combined_spots['y'].min():.1f} to {combined_spots['y'].max():.1f} µm")
    print(f"  z: {combined_spots['z'].min():.1f} to {combined_spots['z'].max():.1f} µm")

    # Save
    output_path = PROCESSED_DIR / "visium_spatial.parquet"
    combined_spots.to_parquet(output_path)
    print(f"\nSaved to {output_path}")

    # Save gene list
    genes_df = pd.DataFrame({'gene': all_genes})
    genes_path = PROCESSED_DIR / "visium_genes.parquet"
    genes_df.to_parquet(genes_path)
    print(f"Saved {len(all_genes)} genes to {genes_path}")

    return combined_spots


def create_human_spatial_cells():
    """Create cells_with_coords.parquet using Visium spatial data.

    This replaces the UMAP-based coordinates with true spatial coordinates.
    """

    # Load Visium spatial data
    visium_path = PROCESSED_DIR / "visium_spatial.parquet"
    if not visium_path.exists():
        print("Extracting Visium data first...")
        extract_visium_data()

    spots_df = pd.read_parquet(visium_path)
    spots_df = spots_df.reset_index()

    # Rename columns to match expected format
    cells_df = spots_df.rename(columns={
        'spot_id': 'cell_id',
    })

    # Add placeholder columns for compatibility
    cells_df['region'] = 'Visium'  # Spots don't have pre-assigned regions
    cells_df['cell_type'] = 'Spatial spot'

    # Save as the main human spatial data
    output_path = PROCESSED_DIR / "cells_with_coords.parquet"
    cells_df.to_parquet(output_path)
    print(f"Saved {len(cells_df)} spatial spots to {output_path}")

    return cells_df


if __name__ == "__main__":
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

    print("=== Extracting Human Visium Spatial Data ===\n")
    extract_visium_data()

    print("\n=== Creating Spatial Cells File ===\n")
    create_human_spatial_cells()
