# HypoMap 3D Atlas Viewer

Interactive 3D visualization of the human and mouse hypothalamus single-cell atlases.

## Overview

This application visualizes spatial transcriptomics data from:
- **Mouse ABC**: Allen Brain Cell Census MERFISH hypothalamus data (~133k cells)

## Installation

```bash
# Clone the repository
git clone <repository-url>
cd hypomap

# Install dependencies with uv
uv sync
```
### Preprocessing

Use Snakemake to run the preprocessing pipeline:

```bash
# Process all datasets
snakemake --cores 1

# Or process individual datasets
snakemake mouse_abc --cores 1       # Mouse ABC Atlas (MERFISH)

# Dry run to see what will be executed
snakemake -n
```

## Running the App

```bash
uv run python -m app.app
```

The app will start at **http://localhost:8050**

## Features

- **3D Scatter Plot**: Interactive visualization of cells/spots in 3D space
- **Dataset Selection**: Switch between human and mouse datasets
- **Color By**: Region, cell type, sample/section
- **Gene Expression**: Query expression of marker genes and receptors
- **Region Filtering**: Focus on specific anatomical regions

## Data Sources

### Human HypoMap (Tadross et al. 2025)
- 433,369 nuclei from snRNA-seq (11 donors)
- 87,353 Visium spots from spatial transcriptomics (9 sections, 7 donors)
- True spatial coordinates from tissue sections

### Mouse HypoMap (Steuernagel et al. 2022)
- Unified mouse hypothalamic single-cell atlas
- Coordinates mapped to Allen Common Coordinate Framework (CCF)

### Mouse ABC (Allen Brain Cell Census)
- ~133,000 MERFISH cells from hypothalamus (parcellation_division = HY)
- 3D reconstructed coordinates registered to CCF
- Hierarchical cell type taxonomy (class → subclass → supertype → cluster)

## Project Structure

```
hypomap/
├── app/
│   ├── app.py          # Main Dash application
│   ├── callbacks.py    # Interactive callbacks
│   └── layouts.py      # UI layout components
├── src/
│   ├── datasets/       # Dataset loaders
│   ├── preprocessing/  # Data extraction & processing
│   └── atlas/          # Allen Brain Atlas integration
├── data/
│   ├── raw/            # Original data files
│   └── processed/      # Processed parquet files
├── references/         # Reference papers
└── requirements.txt
```

## References

1. Tadross JA, Steuernagel L, et al. (2025). A comprehensive spatio-cellular map of the human hypothalamus. *Nature* 639, 708-716.

2. Steuernagel L, et al. (2022). HypoMap—a unified single-cell gene expression atlas of the murine hypothalamus. *Nature Metabolism* 4, 1402-1419.

3. Yao Z, et al. (2023). A high-resolution transcriptomic and spatial atlas of cell types in the whole mouse brain. *Nature* 624, 317-332.

4. Oh SW, et al. (2014). A mesoscale connectome of the mouse brain. *Nature* 508, 207-214.

## License

See LICENSE file for details.
