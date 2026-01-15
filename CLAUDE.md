# CLAUDE.md

## Project Overview

HypoMap is a 3D atlas viewer for mouse and human hypothalamus single-cell data. It visualizes spatial transcriptomics from:
- **Mouse**: ~385k cells from Steuernagel et al. 2022, mapped to Allen CCF coordinates
- **Human**: ~87k Visium spots from Tadross et al. 2025, with true tissue spatial coordinates

The app is built with Dash/Plotly and allows interactive exploration of cell types, gene expression, and anatomical regions.

## Key Directories

- `data/raw/` - Original h5ad files and Visium data (not in git)
- `data/processed/` - Parquet files ready for visualization
- `src/datasets/` - Dataset loaders (mouse.py, human.py)
- `src/preprocessing/` - Data processing scripts
- `app/` - Dash application

## Reproducibility

**IMPORTANT**: All data processing must go through the Snakemake pipeline.

When modifying any preprocessing code:
1. Update the relevant rule in `Snakefile`
2. Run `snakemake -n` to verify the DAG
3. Run `snakemake --cores 1` to regenerate outputs

```bash
# Check what needs to run
snakemake -n

# Run full pipeline
snakemake --cores 1

# Run specific target
snakemake data/processed/mouse/cells_with_coords.parquet
```

If you add a new preprocessing step, add a corresponding rule to the Snakefile with proper input/output declarations.

## Running the App

```bash
python -m app.app
# Opens at http://localhost:8050
```

## Common Tasks

- **Add new marker genes**: Edit `src/datasets/mouse.py` or `human.py`
- **Change downsampling**: Edit `src/preprocessing/downsample.py`, then `snakemake`
- **Update coordinates**: Edit `src/preprocessing/assign_coordinates.py`, then `snakemake`
