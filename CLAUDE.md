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

## Environment

**IMPORTANT**: All commands must be run within the uv environment using `uv run`.

```bash
# Example: run any python command
uv run python -m app.app
```

## Reproducibility

**IMPORTANT**: All data processing must go through the Snakemake pipeline.

When modifying any preprocessing code:
1. Update the relevant rule in `Snakefile`
2. Run `uv run snakemake -n` to verify the DAG
3. Run `uv run snakemake --cores 1` to regenerate outputs

```bash
# Check what needs to run
uv run snakemake -n

# Run full pipeline
uv run snakemake --cores 1

# Run specific target
uv run snakemake data/processed/mouse/cells_with_coords.parquet
```

If you add a new preprocessing step, add a corresponding rule to the Snakefile with proper input/output declarations.

## Running the App

```bash
uv run python -m app.app
# Opens at http://localhost:8050
```

## Plotly Performance with Subplots

**CRITICAL**: Plotly's `fig.add_trace()` is extremely slow with subplots (~200ms per trace). For the 18-subplot coronal atlas, naive implementation took 8-10 seconds per update.

**Solution** (see `app/callbacks.py:create_slice_figure`):
1. Build traces as plain Python dicts, not `go.Scatter()` objects
2. Use `fig.add_traces(list_of_dicts)` to batch-add all traces at once
3. Use `scattergl` (WebGL) instead of `scatter` for large point sets
4. Combine multiple shapes into single traces (e.g., all region boundaries per slice)
5. Add annotations via `layout.annotations` list, not `fig.add_annotation()`

This reduced render time from 8-10 seconds to <1 second.

# Marimo

The notebooks in `notebooks/` use Marimo.  Marimo notebooks are reactive computational notebooks written in standard Python with annotations. Unlike traditional notebooks, cells form a directed acyclic graph (DAG) and automatically re-execute when their dependencies change.

## Critical Rules

1. **No variable redeclaration** — Each variable can only be defined in one cell
2. **No circular dependencies** — The dependency graph must be acyclic
3. **UI values require separate cells** — Access `.value` in a different cell than where the UI element is defined
4. **Underscore prefix = cell-local** — Variables like `_temp` won't be visible to other cells

## Cell Structure

Only edit code inside `@app.cell`. Marimo handles function parameters and returns:

```python
@app.cell
def _():
    # your code here
    return
```

## Quick Reference

- **Display**: Last expression auto-displays (like Jupyter)
- **Markdown**: `mo.md("# Title")`
- **Layout**: `mo.hstack([a, b])`, `mo.vstack([a, b])`, `mo.tabs({"Tab1": content})`
- **SQL**: `df = mo.sql(f"""SELECT * FROM table""")`
- **Plots**: Return figure directly; for matplotlib use `plt.gca()` not `plt.show()`
- **Data**: Prefer polars over pandas

## Example: Reactive UI

```python
@app.cell
def _():
    import marimo as mo
    import altair as alt
    import polars as pl
    return

@app.cell
def _():
    n_points = mo.ui.slider(10, 100, value=50, label="Number of points")
    n_points  # display the slider
    return

@app.cell
def _():
    # This cell re-runs automatically when slider changes
    df = pl.DataFrame({
        "x": np.random.rand(n_points.value),
        "y": np.random.rand(n_points.value)
    })
    alt.Chart(df).mark_circle().encode(x="x", y="y")
    return
```

## Common Mistakes

| Problem | Solution |
|---------|----------|
| "Variable already defined" | Move definition to a single cell, reference elsewhere |
| "Cycle detected" | Reorganize so cell A doesn't depend on B while B depends on A |
| UI value is stale/None | Access `.value` in a downstream cell, not where UI is created |

## After Editing

Run `marimo check --fix` to catch formatting issues and common pitfalls.