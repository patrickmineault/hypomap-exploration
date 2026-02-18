"""Precompute lateralized regions and boundaries for the coronal atlas app."""

import argparse
import pandas as pd
import numpy as np
import alphashape
from shapely.geometry import Polygon, MultiPolygon
import json
from pathlib import Path

# Paths
DATA_DIR = Path(__file__).parent.parent.parent / "data"
DEFAULT_INPUT = DATA_DIR / "processed" / "mouse_abc" / "cells_with_coords.parquet"
DEFAULT_OUTPUT = DATA_DIR / "processed" / "mouse_abc" / "coronal_atlas_regions.json"

# Midline X coordinate
MIDLINE_X = 0.0

# Alphashape parameter (higher = tighter fit, lower = more convex)
ALPHA = 6.0

# Laterality threshold (fraction of cells needed on each side to split L/R)
LATERALITY_THRESHOLD = 0.01


def main(input_path=None, output_path=None):
    if input_path is None:
        input_path = DEFAULT_INPUT
    if output_path is None:
        output_path = DEFAULT_OUTPUT

    print("Loading cell data...")
    df = pd.read_parquet(input_path)
    print(f"Loaded {len(df):,} cells")

    # Filter out *-unassigned regions (e.g. HY-unassigned, TH-unassigned)
    unassigned_mask = df['region'].str.endswith('-unassigned')
    df = df[~unassigned_mask].copy()
    print(f"After removing *-unassigned: {len(df):,} cells")

    # Add z_slice column
    df['z_slice'] = df['z'].round(2)
    slices = sorted(df['z_slice'].unique())
    print(f"Found {len(slices)} slices: {slices}")

    # Compute lateralized region assignments
    print("\nComputing lateralized regions...")
    df['region_display'] = df['region']

    for z_slice in slices:
        slice_mask = df['z_slice'] == z_slice

        for region in df.loc[slice_mask, 'region'].unique():
            region_mask = slice_mask & (df['region'] == region)
            region_cells = df.loc[region_mask]

            if len(region_cells) < 5:
                continue

            # Check if region spans both sides of midline
            left_count = (region_cells['x'] < MIDLINE_X).sum()
            right_count = (region_cells['x'] >= MIDLINE_X).sum()

            # Split if both sides have significant presence
            total = len(region_cells)
            if left_count > LATERALITY_THRESHOLD * total and right_count > LATERALITY_THRESHOLD * total:
                left_mask = region_mask & (df['x'] < MIDLINE_X)
                right_mask = region_mask & (df['x'] >= MIDLINE_X)
                df.loc[left_mask, 'region_display'] = f"{region}-L"
                df.loc[right_mask, 'region_display'] = f"{region}-R"

    n_lateralized = (df['region_display'] != df['region']).sum()
    print(f"Lateralized {n_lateralized:,} cells into L/R regions")

    # Compute boundaries using alphashape
    print(f"\nComputing region boundaries (alphashape, alpha={ALPHA})...")
    boundaries = {}
    for z_slice in slices:
        slice_df = df[df['z_slice'] == z_slice]
        boundaries[str(z_slice)] = {}

        for region in slice_df['region_display'].unique():
            region_df = slice_df[slice_df['region_display'] == region]
            if len(region_df) >= 3:
                points = region_df[['x', 'y']].values
                try:
                    shape = alphashape.alphashape(points, ALPHA)

                    # Extract boundary coordinates from the shape
                    if shape.is_empty:
                        continue

                    # Handle both Polygon and MultiPolygon cases
                    if isinstance(shape, Polygon):
                        coords = list(shape.exterior.coords)
                        boundaries[str(z_slice)][region] = [list(c) for c in coords]
                    elif isinstance(shape, MultiPolygon):
                        # Use the largest polygon
                        largest = max(shape.geoms, key=lambda p: p.area)
                        coords = list(largest.exterior.coords)
                        boundaries[str(z_slice)][region] = [list(c) for c in coords]
                except Exception:
                    pass

    # Compute centroids
    print("Computing region centroids...")
    centroids = {}
    for z_slice in slices:
        slice_df = df[df['z_slice'] == z_slice]
        centroids[str(z_slice)] = {}

        for region in slice_df['region_display'].unique():
            region_df = slice_df[slice_df['region_display'] == region]
            if len(region_df) >= 1:
                centroids[str(z_slice)][region] = [
                    float(region_df['x'].mean()),
                    float(region_df['y'].mean())
                ]

    # Create cell_id -> region_display mapping
    print("Creating cell region mapping...")
    cell_regions = df.set_index('cell_id')['region_display'].to_dict()

    # Save everything
    output = {
        'slices': [float(s) for s in slices],
        'boundaries': boundaries,
        'centroids': centroids,
        'cell_regions': cell_regions,
    }

    print(f"\nSaving to {output_path}...")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w') as f:
        json.dump(output, f)

    print(f"Done! Output size: {output_path.stat().st_size / 1024 / 1024:.1f} MB")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Precompute lateralized regions and boundaries for the coronal atlas app."
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=None,
        help="Input cells_with_coords.parquet path",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output coronal_atlas_regions.json path",
    )
    args = parser.parse_args()
    main(input_path=args.input, output_path=args.output)
