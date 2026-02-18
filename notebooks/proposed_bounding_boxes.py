# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "auto-mix-prep",
#     "marimo>=0.19.11",
#     "matplotlib",
#     "numpy",
#     "pandas",
#     "pyarrow",
#     "scipy",
# ]
# ///

import marimo

__generated_with = "0.19.11"
app = marimo.App(width="medium")


@app.cell
def _():
    import json
    import math
    from pathlib import Path

    import marimo as mo
    import matplotlib
    import matplotlib.patches as mpatches
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd

    matplotlib.rcParams["figure.dpi"] = 150

    return Path, json, math, mo, mpatches, np, pd, plt


@app.cell
def _(Path, json, np, pd):
    # Load cell data
    _data_dir = Path("../data/processed/mouse_abc_subcortical")
    cells_df = pd.read_parquet(_data_dir / "cells_with_coords.parquet")

    # Filter unassigned regions
    cells_df = cells_df[
        ~cells_df["region"].str.contains("unassigned", case=False)
    ].copy()
    cells_df["z_slice"] = cells_df["z"].round(2)

    # Compute dominant neuropeptide per cell
    _np_expr = pd.read_parquet(_data_dir / "neuropeptide_expression.parquet")
    _np_vals = _np_expr.loc[cells_df["cell_id"]].values.astype(np.float32)
    _np_genes = np.array(_np_expr.columns.tolist())

    # Filter: keep genes with SD above median (removes low-variance genes)
    _sds = np.std(_np_vals, axis=0)
    _sd_thresh = np.median(_sds)
    _keep = _sds > _sd_thresh
    _np_vals = _np_vals[:, _keep]
    _kept_genes = _np_genes[_keep]
    _kept_sds = _sds[_keep]
    _kept_means = np.mean(_np_vals, axis=0)

    # Z-score and assign dominant neuropeptide (argmax z-score per cell)
    _zscored = (_np_vals - _kept_means) / _kept_sds
    cells_df["dominant_neuropeptide"] = _kept_genes[np.argmax(_zscored, axis=1)]

    del _np_expr, _np_vals, _zscored  # free memory

    # Load region boundaries
    with open(_data_dir / "coronal_atlas_regions.json", "r") as _f:
        _region_data = json.load(_f)

    all_slices = sorted(cells_df["z_slice"].unique())
    boundaries = {float(k): v for k, v in _region_data["boundaries"].items()}
    centroids = {
        float(k): {r: tuple(c) for r, c in regions.items()}
        for k, regions in _region_data["centroids"].items()
    }

    print(f"Loaded {len(cells_df):,} cells across {len(all_slices)} slices")
    print(
        f"Granularity levels: class ({cells_df['class'].nunique()}), "
        f"subclass ({cells_df['subclass'].nunique()}), "
        f"supertype ({cells_df['supertype'].nunique()}), "
        f"cluster ({cells_df['cluster'].nunique()}), "
        f"dominant_neuropeptide ({cells_df['dominant_neuropeptide'].nunique()})"
    )
    print(
        f"Neuropeptide SD filter: kept {len(_kept_genes)}/{len(_np_genes)} genes "
        f"(SD > {_sd_thresh:.3f})"
    )
    return boundaries, cells_df, centroids


@app.cell(hide_code=True)
def _(math, mpatches, plt):
    def make_bbox_figure(z_slices, x_range, y_range, boundaries, centroids, title):
        """Create a multi-slice subplot figure with a bounding box rectangle."""
        x0, x1 = x_range
        y0, y1 = y_range
        pad = 0.3

        # Find best tiling: minimize empty cells, prefer fewer rows
        n = len(z_slices)
        ncols = min(n, 4)
        nrows = math.ceil(n / ncols)
        best_empty = ncols * nrows - n
        for try_ncols in range(3, 6):
            if try_ncols > n:
                break
            try_nrows = math.ceil(n / try_ncols)
            empty = try_ncols * try_nrows - n
            if empty < best_empty or (empty == best_empty and try_nrows < nrows):
                best_empty = empty
                ncols = try_ncols
                nrows = try_nrows

        x_extent = (x1 - x0) + 2 * pad
        y_extent = (y1 - y0) + 2 * pad
        data_aspect = y_extent / x_extent
        subplot_w = 4
        subplot_h = subplot_w * data_aspect

        fig, axes = plt.subplots(
            nrows,
            ncols,
            figsize=(ncols * subplot_w, nrows * subplot_h + 0.8),
            squeeze=False,
            sharex=True,
            sharey=True,
        )

        for i, z in enumerate(z_slices):
            ax = axes[i // ncols][i % ncols]

            # Bounding box rectangle
            rect = mpatches.Rectangle(
                (x0, y0),
                x1 - x0,
                y1 - y0,
                linewidth=2,
                edgecolor="red",
                facecolor="red",
                alpha=0.08,
            )
            ax.add_patch(rect)
            # Red border on top (full opacity)
            rect_border = mpatches.Rectangle(
                (x0, y0),
                x1 - x0,
                y1 - y0,
                linewidth=2,
                edgecolor="red",
                facecolor="none",
            )
            ax.add_patch(rect_border)

            # Region boundaries
            zr = round(z, 2)
            if zr in boundaries:
                for rgn, coords in boundaries[zr].items():
                    xs = [c[0] for c in coords] + [coords[0][0]]
                    ys = [c[1] for c in coords] + [coords[0][1]]
                    ax.plot(xs, ys, color="gray", linewidth=0.5, alpha=0.6, clip_on=True)

            # Region labels (only if centroid is inside the visible area)
            if zr in centroids:
                for rgn, (cx, cy) in centroids[zr].items():
                    if not (x0 - pad <= cx <= x1 + pad and y0 - pad <= cy <= y1 + pad):
                        continue
                    ax.text(
                        cx,
                        cy,
                        rgn,
                        fontsize=5,
                        color="white",
                        ha="center",
                        va="center",
                        clip_on=True,
                        bbox=dict(
                            boxstyle="round,pad=0.15", fc="black", alpha=0.5, lw=0
                        ),
                    )

            ax.set_title(f"y={z:.2f}", fontsize=9)
            ax.tick_params(labelsize=6)
            ax.set_facecolor("#f8f8f8")
            for spine in ax.spines.values():
                spine.set_visible(False)

        # Set shared limits (y inverted)
        axes[0][0].set_xlim(x0 - pad, x1 + pad)
        axes[0][0].set_ylim(y1 + pad, y0 - pad)
        axes[0][0].set_aspect("equal")

        # Hide unused axes
        for i in range(n, nrows * ncols):
            axes[i // ncols][i % ncols].set_visible(False)

        fig.suptitle(title, fontsize=12)
        fig.tight_layout(rect=[0, 0, 1.0, 0.96])
        return plt.gca()

    return (make_bbox_figure,)


@app.cell
def _(mo):
    mo.md(r"""
    # Proposed bounding boxes
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## HY microcircuit
    """)
    return


@app.cell
def _(boundaries, cells_df, centroids, make_bbox_figure, mo):
    def plot_regions(_selected):
        if _selected and len(_selected) > 0:
            _sel_df = cells_df[cells_df["region"].isin(_selected)]
            if len(_sel_df) > 0:
                _x0 = float(_sel_df["x"].min())
                _x1 = float(_sel_df["x"].max())
                _y0 = float(_sel_df["y"].min())
                _y1 = float(_sel_df["y"].max())
                _sel_slices = sorted(_sel_df["z_slice"].unique())
                _md = mo.md(
                    f"""
        ### {len(_selected)} region(s) selected: {', '.join(_selected)}
        | Property | Value |
        |----------|-------|
        | **X range** | {_x0:.2f} – {_x1:.2f} mm |
        | **Y range** | {_y0:.2f} – {_y1:.2f} mm |
        | **Z slices** | {len(_sel_slices)} ({_sel_slices[0]:.2f} – {_sel_slices[-1]:.2f}) |
        | **Box size** | {_x1 - _x0:.2f} x {_y1 - _y0:.2f} x {(_sel_slices[-1] - _sel_slices[0] + 0.2):.2f} mm |
        | **Volume** | {(_x1 - _x0) * (_y1 - _y0) * (_sel_slices[-1] - _sel_slices[0] + 0.2):.3f} mm³ |
        """
                )
                _region_fig = make_bbox_figure(
                    _sel_slices,
                    (_x0, _x1),
                    (_y0, _y1),
                    boundaries,
                    centroids,
                    f"Bounding box for: {', '.join(_selected)}",
                )
                return mo.vstack([_md, _region_fig])

    plot_regions(["ARH", "ME", "RCH"])
    return (plot_regions,)


@app.cell
def _():
    ## Medium volume
    return


@app.cell
def _(plot_regions):
    plot_regions(["ARH", "ME", "PVH", "DMH", "VMH", "AHN", "RCH"])
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## Whole volume
    """)
    return


@app.cell
def _(plot_regions):
    plot_regions(["ARH", "ME", "PVH", "DMH", "VMH", "AHN", "RCH", "ZI"])
    return


@app.cell
def _():
    ## Whole volume +
    return


@app.cell
def _(plot_regions):
    plot_regions(["ARH", "ME", "PVH", "DMH", "VMH", "AHN", "RCH", "ZI", "BST", "VP"])
    return


if __name__ == "__main__":
    app.run()
