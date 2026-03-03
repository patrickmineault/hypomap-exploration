"""Precompute lateralized regions and boundaries for the coronal atlas app."""

import argparse
from pathlib import Path

import alphashape
import msgpack
import numpy as np
import pandas as pd
from shapely.geometry import MultiPolygon, Polygon

# Paths
DATA_DIR = Path(__file__).parent.parent.parent / "data"
DEFAULT_INPUT = DATA_DIR / "processed" / "mouse_abc" / "cells_with_coords.parquet"
DEFAULT_OUTPUT = DATA_DIR / "processed" / "mouse_abc" / "coronal_atlas_regions.msgpack"

# Midline X coordinate
MIDLINE_X = 0.0

# Alphashape parameter (higher = tighter fit, lower = more convex)
ALPHA = 6.0

# Laterality threshold (fraction of cells needed on each side to split L/R)
LATERALITY_THRESHOLD = 0.01


def main(input_path=None, output_path=None, region_col='region'):
    if input_path is None:
        input_path = DEFAULT_INPUT
    if output_path is None:
        output_path = DEFAULT_OUTPUT

    print("Loading cell data...")
    df = pd.read_parquet(input_path)
    print(f"Loaded {len(df):,} cells")

    if region_col != 'region':
        df['region'] = df[region_col]

    # Filter out *-unassigned regions (e.g. HY-unassigned, TH-unassigned)
    unassigned_mask = df['region'].str.endswith('-unassigned')
    df = df[~unassigned_mask].copy()
    print(f"After removing *-unassigned: {len(df):,} cells")

    # Add z_slice column (promote to float64 before rounding to avoid float32 artifacts)
    df['z_slice'] = df['z'].astype(np.float64).round(2)
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
                        coords = np.array(shape.exterior.coords, dtype=np.float32)
                        boundaries[str(z_slice)][region] = coords.tobytes()
                    elif isinstance(shape, MultiPolygon):
                        # Use the largest polygon
                        largest = max(shape.geoms, key=lambda p: p.area)
                        coords = np.array(largest.exterior.coords, dtype=np.float32)
                        boundaries[str(z_slice)][region] = coords.tobytes()
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
                centroids[str(z_slice)][region] = np.array(
                    [region_df['x'].mean(), region_df['y'].mean()], dtype=np.float32
                ).tobytes()

    # Save region_display + z_slice back into cells parquet (avoids recomputing at load time)
    print(f"\nSaving region_display + z_slice into {input_path}...")
    full_df = pd.read_parquet(input_path)
    # Map region_display from the filtered df back to the full df
    full_df["z_slice"] = full_df["z"].round(2)
    if "region_display" in full_df.columns:
        full_df = full_df.drop(columns=["region_display"])
    display_map = df.set_index(df.index)["region_display"]
    full_df["region_display"] = display_map.reindex(full_df.index)
    # Unassigned rows keep their base region
    full_df["region_display"] = full_df["region_display"].fillna(full_df[region_col])
    full_df.to_parquet(input_path, index=False)
    print(f"Updated {input_path}")

    # Save everything
    output = {
        'slices': [float(s) for s in slices],
        'boundaries': boundaries,
        'centroids': centroids,
    }

    print(f"\nSaving to {output_path}...")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'wb') as f:
        f.write(msgpack.packb(output, use_bin_type=True))

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
        help="Output coronal_atlas_regions.msgpack path",
    )
    parser.add_argument(
        "--region-col",
        default="region",
        help="Column name to use for region labels (default: region)",
    )
    args = parser.parse_args()
    main(input_path=args.input, output_path=args.output, region_col=args.region_col)
