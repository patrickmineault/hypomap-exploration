"""Precompute lateralized regions and boundaries for the coronal atlas app."""

import pandas as pd
import numpy as np
from scipy.spatial import ConvexHull
import json
from pathlib import Path

# Paths
DATA_DIR = Path(__file__).parent.parent.parent / "data"
CELLS_PATH = DATA_DIR / "processed" / "mouse_abc" / "cells_with_coords.parquet"
OUTPUT_PATH = DATA_DIR / "processed" / "mouse_abc" / "coronal_atlas_regions.json"

# Midline X coordinate
MIDLINE_X = 5.5


def main():
    print("Loading cell data...")
    df = pd.read_parquet(CELLS_PATH)
    print(f"Loaded {len(df):,} cells")

    # Filter out HY-unassigned
    df = df[df['region'] != 'HY-unassigned'].copy()
    print(f"After removing HY-unassigned: {len(df):,} cells")

    # Add z_slice column
    df['z_slice'] = df['z'].round(1)
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

            # Split if both sides have significant presence (>10% each)
            total = len(region_cells)
            if left_count > 0.1 * total and right_count > 0.1 * total:
                left_mask = region_mask & (df['x'] < MIDLINE_X)
                right_mask = region_mask & (df['x'] >= MIDLINE_X)
                df.loc[left_mask, 'region_display'] = f"{region}-L"
                df.loc[right_mask, 'region_display'] = f"{region}-R"

    n_lateralized = (df['region_display'] != df['region']).sum()
    print(f"Lateralized {n_lateralized:,} cells into L/R regions")

    # Compute boundaries (convex hulls)
    print("\nComputing region boundaries...")
    boundaries = {}
    for z_slice in slices:
        slice_df = df[df['z_slice'] == z_slice]
        boundaries[str(z_slice)] = {}

        for region in slice_df['region_display'].unique():
            region_df = slice_df[slice_df['region_display'] == region]
            if len(region_df) >= 3:
                points = region_df[['x', 'y']].values
                try:
                    hull = ConvexHull(points)
                    hull_points = points[hull.vertices].tolist()
                    hull_points.append(hull_points[0])  # Close polygon
                    boundaries[str(z_slice)][region] = hull_points
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

    print(f"\nSaving to {OUTPUT_PATH}...")
    with open(OUTPUT_PATH, 'w') as f:
        json.dump(output, f)

    print(f"Done! Output size: {OUTPUT_PATH.stat().st_size / 1024 / 1024:.1f} MB")


if __name__ == "__main__":
    main()
