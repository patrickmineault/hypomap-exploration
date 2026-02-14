"""Stratified downsampling of HypoMap cell data."""

import argparse
import sys
from pathlib import Path
from typing import Optional

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent.parent))

from hypomap.datasets import detect_cell_type_levels, detect_region_column, get_config

DATA_DIR = Path(__file__).parent.parent.parent / "data"
PROCESSED_DIR = DATA_DIR / "processed"


def stratified_downsample(df, stratify_cols, target_n=50000, min_per_group=10):
    """
    Perform stratified downsampling to preserve cell type representation.

    Args:
        df: DataFrame with cell metadata
        stratify_cols: List of columns to stratify by (e.g., ['Region', 'C66'])
        target_n: Target total number of cells
        min_per_group: Minimum cells to keep per stratum

    Returns:
        Downsampled DataFrame
    """
    if len(df) <= target_n:
        print(f"Dataset already has {len(df)} cells, no downsampling needed")
        return df

    # Create combined stratification key
    df = df.copy()

    # Handle missing values in stratification columns
    for col in stratify_cols:
        if col not in df.columns:
            print(f"Warning: Column {col} not found, skipping")
            stratify_cols = [c for c in stratify_cols if c != col]

    if not stratify_cols:
        # No stratification columns, random sample
        return df.sample(n=target_n, random_state=42)

    df['_strat_key'] = df[stratify_cols].astype(str).agg('_'.join, axis=1)

    # Count cells per stratum
    strat_counts = df['_strat_key'].value_counts()
    n_strata = len(strat_counts)

    print(f"Total cells: {len(df)}")
    print(f"Target cells: {target_n}")
    print(f"Number of strata: {n_strata}")

    # Calculate sampling fractions
    # Base fraction to get approximately target_n cells
    base_frac = target_n / len(df)

    sampled_indices = []

    for stratum, count in strat_counts.items():
        stratum_df = df[df['_strat_key'] == stratum]

        # Calculate number to sample from this stratum
        n_sample = max(min_per_group, int(count * base_frac))
        n_sample = min(n_sample, count)  # Can't sample more than available

        # Random sample from stratum
        sampled = stratum_df.sample(n=n_sample, random_state=42)
        sampled_indices.extend(sampled.index.tolist())

    result = df.loc[sampled_indices].drop(columns=['_strat_key'])

    print(f"Sampled {len(result)} cells from {n_strata} strata")

    return result


def run_downsampling(
    dataset: str = "mouse_hypomap",
    cell_type_level: Optional[str] = None,
    target_n: int = 50000,
    assign_coords: bool = True,
):
    """
    Run the full downsampling pipeline for a dataset.

    Args:
        dataset: Dataset name ("mouse_hypomap", "human_hypomap", or "mouse_abc")
        cell_type_level: Which cell type hierarchy level to stratify by.
                        If None, auto-selects based on dataset.
        target_n: Target number of cells after downsampling
        assign_coords: Whether to assign 3D coordinates
    """
    # Get dataset configuration
    config = get_config(dataset)

    # Load metadata
    metadata_path = config.cell_metadata_path
    if not metadata_path.exists():
        raise FileNotFoundError(
            f"Run extraction first for {dataset}:\n"
            f"  python -m hypomap.datasets.{dataset}"
        )

    print(f"Loading {dataset} metadata from {metadata_path}")
    df = pd.read_parquet(metadata_path)
    print(f"Loaded {len(df)} cells")

    # Find region column
    region_col = detect_region_column(df)
    print(f"Using region column: {region_col}")

    # Find cell type columns and select level
    cell_type_cols = detect_cell_type_levels(df)
    print(f"Available cell type levels: {cell_type_cols}")

    if cell_type_level is None:
        # Auto-select: use second level if available (moderate granularity)
        if len(cell_type_cols) >= 2:
            cell_type_col = cell_type_cols[1]
        elif cell_type_cols:
            cell_type_col = cell_type_cols[0]
        else:
            cell_type_col = None
    else:
        cell_type_col = cell_type_level if cell_type_level in df.columns else None

    print(f"Using cell type column: {cell_type_col}")

    # Build stratification columns
    stratify_cols = []
    if region_col:
        stratify_cols.append(region_col)
    if cell_type_col:
        stratify_cols.append(cell_type_col)

    print(f"Stratifying by: {stratify_cols}")

    # Skip downsampling for mouse_abc datasets — display-time subsampling handles performance
    if dataset.startswith("mouse_abc"):
        print(f"Skipping downsampling for {dataset} (display-time subsampling is sufficient)")
        downsampled = df
    else:
        # Perform downsampling
        downsampled = stratified_downsample(df, stratify_cols, target_n=target_n)

    # Add 3D coordinates
    if assign_coords:
        print("\nAssigning 3D coordinates...")
        from hypomap.preprocessing.assign_coordinates import assign_coordinates
        downsampled = assign_coordinates(downsampled, dataset=dataset, region_col=region_col)
    else:
        print("Skipping coordinate assignment")

    # Save result
    output_path = config.cells_with_coords_path
    output_path.parent.mkdir(parents=True, exist_ok=True)
    downsampled.to_parquet(output_path)
    print(f"\nSaved {len(downsampled)} cells to {output_path}")

    # Print summary
    print(f"\n=== {dataset.title()} Downsampling Summary ===")
    if region_col and region_col in downsampled.columns:
        print("\nCells per region:")
        print(downsampled[region_col].value_counts())

    if cell_type_col and cell_type_col in downsampled.columns:
        print(f"\nUnique {cell_type_col} types: {downsampled[cell_type_col].nunique()}")

    # Check coordinate coverage
    if 'x' in downsampled.columns:
        n_with_coords = downsampled['x'].notna().sum()
        print(f"\nCells with 3D coordinates: {n_with_coords} ({100*n_with_coords/len(downsampled):.1f}%)")

    return downsampled


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Downsample HypoMap cell data")
    parser.add_argument(
        "--dataset", "-d",
        default="mouse_hypomap",
        choices=["mouse_hypomap", "human_hypomap", "mouse_abc", "mouse_abc_subcortical"],
        help="Dataset to process (default: mouse_hypomap)"
    )
    parser.add_argument(
        "--cell-type-level", "-c",
        default=None,
        help="Cell type column to stratify by (auto-detected if not specified)"
    )
    parser.add_argument(
        "--target", "-n",
        type=int,
        default=50000,
        help="Target number of cells (default: 50000)"
    )
    parser.add_argument(
        "--no-coords",
        action="store_true",
        help="Skip coordinate assignment"
    )

    args = parser.parse_args()
    run_downsampling(
        dataset=args.dataset,
        cell_type_level=args.cell_type_level,
        target_n=args.target,
        assign_coords=not args.no_coords,
    )
    )
