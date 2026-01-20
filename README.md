# HypoMap 3D Atlas Viewer

Interactive 3D visualization of the human and mouse hypothalamus single-cell atlases.

## Overview

This application visualizes spatial transcriptomics data from:
- **Human HypoMap**: 87,353 Visium spots from 9 tissue sections (Tadross et al. 2025, Nature)
- **Mouse HypoMap**: Cells mapped to Allen CCF coordinates (Steuernagel et al. 2022)
- **Mouse ABC**: Allen Brain Cell Census MERFISH hypothalamus data (~133k cells)

## Installation

```bash
# Clone the repository
git clone <repository-url>
cd hypomap

# Install dependencies with uv
uv sync
```

## Data Setup

### Required Raw Data

Download the following files to `data/raw/`:

**Human Hypothalamus** (`data/raw/human_hypomap/`):
- `480e89e7-84ad-4fa8-adc3-f7c562a77a78.h5ad` - snRNA-seq data from [CZ CELLxGENE](https://cellxgene.cziscience.com/collections/d0941303-7ce3-4422-9249-cf31eb98c480)
- `GSE278848_tissue_positions_lists.tar.gz` - Visium spatial positions from [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE278848)
- `visium_h5/*.h5` - Visium gene expression matrices from GEO

**Mouse Hypothalamus** (`data/raw/mouse_hypomap/`):
- `d3be7423-d664-4913-89a9-a506cae4c28f.h5ad` - Mouse HypoMap from CZ CELLxGENE
- `hypothalamus_connectivity.csv` - Intra-hypothalamic connectivity from Allen Mouse Brain Connectivity Atlas ([Oh et al. 2014](https://www.nature.com/articles/nature13186))
- `hypothalamus_structures.csv` - Allen CCF structure IDs for hypothalamic regions

**Curated Data** (`data/generated/mouse_common/`):
- `np_map.csv` - Curated neuropeptide ligand-receptor pairs with hypothalamic nuclei annotations

### Preprocessing

Use Snakemake to run the preprocessing pipeline:

```bash
# Process all datasets
snakemake --cores 1

# Or process individual datasets
snakemake mouse_hypomap --cores 1   # Mouse HypoMap (CCF coordinates)
snakemake human_hypomap --cores 1   # Human HypoMap (Visium spatial)
snakemake mouse_abc --cores 1       # Mouse ABC Atlas (MERFISH)

# Dry run to see what will be executed
snakemake -n
```

This creates processed files in `data/processed/`:
- `mouse_hypomap/cells_with_coords.parquet` - Mouse cells with CCF coordinates
- `human_hypomap/cells_with_coords.parquet` - Human Visium spots with spatial coordinates
- `mouse_abc/cell_metadata.parquet` - Mouse ABC hypothalamus cells (MERFISH)

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
