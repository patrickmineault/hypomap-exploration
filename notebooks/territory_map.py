# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "marimo>=0.19.11",
#     "matplotlib",
#     "numpy",
#     "pandas",
#     "pyarrow",
#     "scipy",
# ]
# ///

import marimo

__generated_with = "0.19.4"
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
    from scipy.ndimage import gaussian_filter

    matplotlib.rcParams["figure.dpi"] = 150
    return Path, gaussian_filter, json, math, mo, mpatches, np, pd, plt


@app.cell
def _(Path, json, np, pd):
    _data_dir = Path("../data/processed/mouse_abc_extended")

    cells_df = pd.read_parquet(_data_dir / "cells_with_coords.parquet")
    cells_df = cells_df[
        ~cells_df["region"].str.contains("unassigned", case=False)
    ].copy()
    cells_df["z_slice"] = cells_df["z"].round(2)

    np_expr = pd.read_parquet(_data_dir / "neuropeptide_expression.parquet")
    # Replace inf with NaN (float16 overflow)
    np_expr = np_expr.replace([np.inf, -np.inf], np.nan)

    with open(_data_dir / "coronal_atlas_regions.json", "r") as _f:
        _region_data = json.load(_f)

    boundaries = {float(k): v for k, v in _region_data["boundaries"].items()}
    centroids = {
        float(k): {r: tuple(c) for r, c in regions.items()}
        for k, regions in _region_data["centroids"].items()
    }

    print(f"Loaded {len(cells_df):,} cells, {np_expr.shape[1]} neuropeptide genes")
    return boundaries, cells_df, centroids, np_expr


@app.cell
def _(mo):
    mo.md(r"""
    # ARH/ME Territory Map
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## Bounding box controls
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    x_range_slider = mo.ui.range_slider(
        start=-2.0, stop=2.0, step=0.05, value=[-0.1, 0.8],
        label="X range (mm)",
    )
    y_range_slider = mo.ui.range_slider(
        start=0.0, stop=4.0, step=0.05, value=[1.2, 2.05],
        label="Y range (mm)",
    )
    z_range_slider = mo.ui.range_slider(
        start=-8.0, stop=0.0, step=0.05, value=[-2.5, -1.55],
        label="Z range (mm)",
    )
    bin_size_slider = mo.ui.slider(
        start=0.01, stop=0.1, step=0.005, value=0.03,
        label="Bin size (mm)",
    )
    smooth_sigma_slider = mo.ui.slider(
        start=0.0, stop=3.0, step=0.1, value=1.0,
        label="Smooth sigma (bins)",
    )
    mo.vstack([
        x_range_slider,
        y_range_slider,
        z_range_slider,
        mo.hstack([bin_size_slider, smooth_sigma_slider]),
    ])
    return (
        bin_size_slider,
        smooth_sigma_slider,
        x_range_slider,
        y_range_slider,
        z_range_slider,
    )


@app.cell(hide_code=True)
def _(
    bin_size_slider,
    cells_df,
    np,
    np_expr,
    x_range_slider,
    y_range_slider,
    z_range_slider,
):
    _x0, _x1 = x_range_slider.value
    _y0, _y1 = y_range_slider.value
    _z0, _z1 = z_range_slider.value

    _mask = (
        (cells_df["x"] >= _x0) & (cells_df["x"] <= _x1) &
        (cells_df["y"] >= _y0) & (cells_df["y"] <= _y1) &
        (cells_df["z"] >= _z0) & (cells_df["z"] <= _z1)
    )
    bbox_df = cells_df[_mask].copy()
    bbox_np = np_expr.loc[bbox_df["cell_id"]].copy()

    bbox_slices = sorted(bbox_df["z_slice"].unique())

    # Shared grid edges
    bin_size = bin_size_slider.value
    x_edges = np.arange(_x0, _x1 + bin_size, bin_size)
    y_edges = np.arange(_y0, _y1 + bin_size, bin_size)
    x_centers = 0.5 * (x_edges[:-1] + x_edges[1:])
    y_centers = 0.5 * (y_edges[:-1] + y_edges[1:])

    print(
        f"Bounding box: {len(bbox_df):,} cells, {len(bbox_slices)} slices, "
        f"grid {len(x_centers)}×{len(y_centers)}"
    )
    return (
        bbox_df,
        bbox_np,
        bbox_slices,
        x_centers,
        x_edges,
        y_centers,
        y_edges,
    )


@app.cell(hide_code=True)
def _(gaussian_filter, math, mpatches, np, plt):
    def make_territory_axes(n_slices, x_range, y_range):
        """Create multi-slice subplot figure with smart tiling."""
        x0, x1 = x_range
        y0, y1 = y_range
        n = n_slices
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

        x_extent = x1 - x0
        y_extent = y1 - y0
        data_aspect = y_extent / x_extent
        subplot_w = 2.5
        subplot_h = subplot_w * data_aspect

        fig, axes = plt.subplots(
            nrows, ncols,
            figsize=(ncols * subplot_w, nrows * subplot_h + 0.8),
            squeeze=False, sharex=True, sharey=True,
        )
        for row in axes:
            for ax in row:
                ax.set_xlim(x0, x1)
                ax.set_ylim(y1, y0)  # inverted y
                ax.set_aspect("equal")
        for i in range(n, nrows * ncols):
            axes[i // ncols][i % ncols].set_visible(False)
        return fig, axes, nrows, ncols

    def draw_region_overlays(ax, z_val, boundaries, centroids, x_range, y_range):
        """Draw gray boundary lines and white-on-black region labels."""
        zr = round(z_val, 2)
        x0, x1 = x_range
        y0, y1 = y_range
        if zr in boundaries:
            for _rgn, coords in boundaries[zr].items():
                xs = [c[0] for c in coords] + [coords[0][0]]
                ys = [c[1] for c in coords] + [coords[0][1]]
                ax.plot(xs, ys, color="gray", linewidth=0.5, alpha=0.6, clip_on=True)
        if zr in centroids:
            for rgn, (cx, cy) in centroids[zr].items():
                if not (x0 <= cx <= x1 and y0 <= cy <= y1):
                    continue
                ax.text(
                    cx, cy, rgn, fontsize=5, color="white", ha="center", va="center",
                    clip_on=True,
                    bbox=dict(boxstyle="round,pad=0.15", fc="black", alpha=0.5, lw=0),
                )

    def compute_class_density(bbox_df, label_col, color_col, x_edges, y_edges, slices, sigma, min_pct=2.0):
        """Compute smoothed density per class label per slice.

        Returns dict[z -> list of (label, color, density_2d)].
        """
        counts = bbox_df[label_col].value_counts()
        total = len(bbox_df)
        keep_labels = counts[counts / total * 100 >= min_pct].index.tolist()

        # Get color mapping
        color_map = bbox_df.drop_duplicates(label_col).set_index(label_col)[color_col].to_dict()

        result = {}
        for z in slices:
            slice_df = bbox_df[bbox_df["z_slice"] == z]
            entries = []
            for label in keep_labels:
                ldf = slice_df[slice_df[label_col] == label]
                if len(ldf) < 3:
                    continue
                h, _, _ = np.histogram2d(
                    ldf["x"].values, ldf["y"].values,
                    bins=[x_edges, y_edges],
                )
                if sigma > 0:
                    h = gaussian_filter(h, sigma=sigma)
                entries.append((label, color_map.get(label, "#888888"), h))
            result[z] = entries
        return result

    def compute_np_density(bbox_df, bbox_np, genes, x_edges, y_edges, slices, sigma):
        """Compute weighted density per neuropeptide gene per slice.

        Returns dict[gene -> dict[z -> density_2d]].
        """
        result = {}
        for gene in genes:
            gene_vals = bbox_np[gene].values.astype(np.float32)
            nan_mask = ~np.isnan(gene_vals)
            gene_result = {}
            for z in slices:
                slice_mask = (bbox_df["z_slice"].values == z) & nan_mask
                if slice_mask.sum() < 3:
                    gene_result[z] = np.zeros((len(x_edges) - 1, len(y_edges) - 1))
                    continue
                h, _, _ = np.histogram2d(
                    bbox_df["x"].values[slice_mask],
                    bbox_df["y"].values[slice_mask],
                    bins=[x_edges, y_edges],
                    weights=gene_vals[slice_mask],
                )
                if sigma > 0:
                    h = gaussian_filter(h, sigma=sigma)
                gene_result[z] = h
            result[gene] = gene_result
        return result

    def compute_rgb_image(np_grids, gene_r, gene_g, gene_b, slices, gamma=0.7):
        """Build RGB image from three gene density grids.

        Returns dict[z -> rgb_array (ny, nx, 3)].
        """
        result = {}
        # Collect all values per channel for global normalization
        all_r = np.concatenate([np_grids[gene_r][z].ravel() for z in slices])
        all_g = np.concatenate([np_grids[gene_g][z].ravel() for z in slices])
        all_b = np.concatenate([np_grids[gene_b][z].ravel() for z in slices])
        r_max = np.percentile(all_r[all_r > 0], 98) if (all_r > 0).any() else 1.0
        g_max = np.percentile(all_g[all_g > 0], 98) if (all_g > 0).any() else 1.0
        b_max = np.percentile(all_b[all_b > 0], 98) if (all_b > 0).any() else 1.0

        for z in slices:
            r = np.clip(np_grids[gene_r][z] / r_max, 0, 1)
            g = np.clip(np_grids[gene_g][z] / g_max, 0, 1)
            b = np.clip(np_grids[gene_b][z] / b_max, 0, 1)
            # Apply gamma correction
            r = np.power(r, gamma)
            g = np.power(g, gamma)
            b = np.power(b, gamma)
            # Stack as (ny, nx, 3) — transpose from histogram (nx, ny)
            rgb = np.stack([r.T, g.T, b.T], axis=-1)
            result[z] = rgb
        return result

    def add_legend_patches(fig, labels_colors, ncol=4, fontsize=6):
        """Add a legend with colored patches to the figure."""
        patches = [mpatches.Patch(color=c, label=l) for l, c in labels_colors]
        fig.legend(
            handles=patches, loc="lower center", ncol=ncol, fontsize=fontsize,
            frameon=False, bbox_to_anchor=(0.5, -0.02),
        )
    return (
        add_legend_patches,
        compute_class_density,
        compute_np_density,
        compute_rgb_image,
        draw_region_overlays,
        make_territory_axes,
    )


@app.cell
def _(mo):
    mo.md(r"""
    ## Cell class territory map
    """)
    return


@app.cell
def _(mo):
    granularity_dropdown = mo.ui.dropdown(
        options=["class", "subclass", "supertype"],
        value="subclass",
        label="Granularity",
    )
    min_pct_slider = mo.ui.slider(
        start=0.5, stop=10.0, step=0.5, value=2.0,
        label="Min % of cells",
    )
    mo.hstack([granularity_dropdown, min_pct_slider])
    return granularity_dropdown, min_pct_slider


@app.cell
def _(
    bbox_df,
    bbox_slices,
    compute_class_density,
    granularity_dropdown,
    min_pct_slider,
    smooth_sigma_slider,
    x_edges,
    y_edges,
):
    _label_col = granularity_dropdown.value
    _color_col = _label_col + "_color"
    class_densities = compute_class_density(
        bbox_df, _label_col, _color_col,
        x_edges, y_edges, bbox_slices,
        sigma=smooth_sigma_slider.value,
        min_pct=min_pct_slider.value,
    )
    return (class_densities,)


@app.cell(hide_code=True)
def _(
    add_legend_patches,
    bbox_slices,
    boundaries,
    centroids,
    class_densities,
    draw_region_overlays,
    granularity_dropdown,
    make_territory_axes,
    np,
    plt,
    x_centers,
    x_range_slider,
    y_centers,
    y_range_slider,
):
    _x_range = tuple(x_range_slider.value)
    _y_range = tuple(y_range_slider.value)
    _fig, _axes, _nrows, _ncols = make_territory_axes(
        len(bbox_slices), _x_range, _y_range,
    )

    _all_labels_colors = {}
    for _i, _z in enumerate(bbox_slices):
        _ax = _axes[_i // _ncols][_i % _ncols]
        _ax.set_facecolor("#f0f0f0")
        _ax.set_title(f"z={_z:.2f}", fontsize=8)
        _ax.tick_params(labelsize=5)

        for _label, _color, _density in class_densities.get(_z, []):
            _nonzero = _density[_density > 0]
            if len(_nonzero) == 0:
                continue
            _threshold = np.percentile(_nonzero, 30)
            _ax.contourf(
                x_centers, y_centers, _density.T,
                levels=[_threshold, _density.max() * 10],
                colors=[_color], alpha=0.25,
            )
            _ax.contour(
                x_centers, y_centers, _density.T,
                levels=[_threshold], colors=[_color], linewidths=0.8,
            )
            _all_labels_colors[_label] = _color

        draw_region_overlays(_ax, _z, boundaries, centroids, _x_range, _y_range)

    _fig.suptitle(f"Cell territory map — {granularity_dropdown.value}", fontsize=11)
    _fig.tight_layout(rect=[0, 0.04, 1.0, 0.96])
    _sorted_lc = sorted(_all_labels_colors.items())
    if _sorted_lc:
        add_legend_patches(_fig, _sorted_lc, ncol=min(5, len(_sorted_lc)), fontsize=5)
    plt.gca()
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## Neuropeptide territory map
    """)
    return


@app.cell
def _(mo):
    min_expr_slider = mo.ui.slider(
        start=0.0, stop=2.0, step=0.1, value=0.5,
        label="Min mean expression",
    )
    min_expr_slider
    return (min_expr_slider,)


@app.cell
def _(
    bbox_df,
    bbox_np,
    bbox_slices,
    compute_np_density,
    min_expr_slider,
    np,
    smooth_sigma_slider,
    x_edges,
    y_edges,
):
    # Filter genes by mean expression in bounding box
    _means = bbox_np.mean().astype(np.float32)
    _means = _means.fillna(0)
    np_genes_filtered = _means[_means >= min_expr_slider.value].index.tolist()

    np_densities = compute_np_density(
        bbox_df, bbox_np, np_genes_filtered,
        x_edges, y_edges, bbox_slices,
        sigma=smooth_sigma_slider.value,
    )
    print(f"Computing territories for {len(np_genes_filtered)} neuropeptide genes")
    return np_densities, np_genes_filtered


@app.cell
def _(
    add_legend_patches,
    bbox_slices,
    boundaries,
    centroids,
    draw_region_overlays,
    make_territory_axes,
    np,
    np_densities,
    np_genes_filtered,
    plt,
    x_centers,
    x_range_slider,
    y_centers,
    y_range_slider,
):
    _x_range = tuple(x_range_slider.value)
    _y_range = tuple(y_range_slider.value)
    _fig, _axes, _nrows, _ncols = make_territory_axes(
        len(bbox_slices), _x_range, _y_range,
    )

    # Assign colors from tab20
    _cmap = plt.cm.tab20
    _gene_colors = {
        g: _cmap(i % 20) for i, g in enumerate(np_genes_filtered)
    }

    _shown_genes = set()
    for _i, _z in enumerate(bbox_slices):
        _ax = _axes[_i // _ncols][_i % _ncols]
        _ax.set_facecolor("#f0f0f0")
        _ax.set_title(f"z={_z:.2f}", fontsize=8)
        _ax.tick_params(labelsize=5)

        for _gene in np_genes_filtered:
            _density = np_densities[_gene].get(_z)
            if _density is None:
                continue
            _nonzero = _density[_density > 0]
            if len(_nonzero) == 0:
                continue
            _threshold = np.percentile(_nonzero, 50)
            _ax.contour(
                x_centers, y_centers, _density.T,
                levels=[_threshold], colors=[_gene_colors[_gene]],
                linewidths=0.6, alpha=0.8,
            )
            _shown_genes.add(_gene)

        draw_region_overlays(_ax, _z, boundaries, centroids, _x_range, _y_range)

    _fig.suptitle("Neuropeptide territory map", fontsize=11)
    _fig.tight_layout(rect=[0, 0.04, 1.0, 0.96])
    _sorted_lc = sorted(
        [(g, _gene_colors[g]) for g in _shown_genes],
    )
    if _sorted_lc:
        add_legend_patches(_fig, _sorted_lc, ncol=min(8, len(_sorted_lc)), fontsize=4)
    plt.gca()
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## RGB neuropeptide overlay
    """)
    return


@app.cell
def _(mo, np_expr):
    _all_genes = sorted(np_expr.columns.tolist())
    _tanycyte_supertypes = [
        "Tanycyte NN_1", "Tanycyte NN_2", "Tanycyte NN_3",
    ]
    _all_options = _tanycyte_supertypes + _all_genes
    gene_r_dropdown = mo.ui.dropdown(
        options=_all_options, value="Agrp", label="Red channel",
    )
    gene_g_dropdown = mo.ui.dropdown(
        options=_all_options, value="Pomc", label="Green channel",
    )
    gene_b_dropdown = mo.ui.dropdown(
        options=_all_options, value="Npy", label="Blue channel",
    )
    gamma_slider = mo.ui.slider(
        start=0.2, stop=2.0, step=0.1, value=0.7,
        label="Gamma",
    )
    r_enabled = mo.ui.checkbox(value=True, label="R")
    g_enabled = mo.ui.checkbox(value=True, label="G")
    b_enabled = mo.ui.checkbox(value=True, label="B")
    mo.hstack([
        gene_r_dropdown, r_enabled,
        gene_g_dropdown, g_enabled,
        gene_b_dropdown, b_enabled,
        gamma_slider,
    ])
    return (
        b_enabled,
        g_enabled,
        gamma_slider,
        gene_b_dropdown,
        gene_g_dropdown,
        gene_r_dropdown,
        r_enabled,
    )


@app.cell(hide_code=True)
def _(
    b_enabled,
    bbox_df,
    bbox_np,
    bbox_slices,
    boundaries,
    centroids,
    compute_np_density,
    compute_rgb_image,
    draw_region_overlays,
    g_enabled,
    gamma_slider,
    gaussian_filter,
    gene_b_dropdown,
    gene_g_dropdown,
    gene_r_dropdown,
    make_territory_axes,
    np,
    plt,
    r_enabled,
    smooth_sigma_slider,
    x_edges,
    x_range_slider,
    y_edges,
    y_range_slider,
):
    _x_range = tuple(x_range_slider.value)
    _y_range = tuple(y_range_slider.value)
    _rgb_channels = [gene_r_dropdown.value, gene_g_dropdown.value, gene_b_dropdown.value]
    _enabled = [r_enabled.value, g_enabled.value, b_enabled.value]
    _sigma = smooth_sigma_slider.value

    _tanycyte_prefix = "Tanycyte NN_"

    # Separate gene channels from supertype channels
    _gene_channels = [ch for ch in _rgb_channels if not ch.startswith(_tanycyte_prefix)]
    _gene_densities = {}
    if _gene_channels:
        _gene_densities = compute_np_density(
            bbox_df, bbox_np, _gene_channels,
            x_edges, y_edges, bbox_slices,
            sigma=_sigma,
        )

    # Compute cell-count density for supertype channels
    _st_densities = {}
    for _ch in _rgb_channels:
        if not _ch.startswith(_tanycyte_prefix):
            continue
        # Match supertype ending with this suffix
        _st_match = bbox_df["supertype"].str.endswith(_ch)
        _st_result = {}
        for _z in bbox_slices:
            _slice_mask = (bbox_df["z_slice"].values == _z) & _st_match.values
            if _slice_mask.sum() < 1:
                _st_result[_z] = np.zeros((len(x_edges) - 1, len(y_edges) - 1))
                continue
            _h, _, _ = np.histogram2d(
                bbox_df["x"].values[_slice_mask],
                bbox_df["y"].values[_slice_mask],
                bins=[x_edges, y_edges],
            )
            if _sigma > 0:
                _h = gaussian_filter(_h, sigma=_sigma)
            _st_result[_z] = _h
        _st_densities[_ch] = _st_result

    # Merge into single dict for compute_rgb_image
    _all_densities = {**_gene_densities, **_st_densities}
    _rgb_images = compute_rgb_image(
        _all_densities,
        _rgb_channels[0], _rgb_channels[1], _rgb_channels[2],
        bbox_slices, gamma=gamma_slider.value,
    )

    _fig, _axes, _nrows, _ncols = make_territory_axes(
        len(bbox_slices), _x_range, _y_range,
    )

    for _i, _z in enumerate(bbox_slices):
        _ax = _axes[_i // _ncols][_i % _ncols]
        _ax.set_facecolor("black")
        _ax.set_title(f"z={_z:.2f}", fontsize=8, color="white")
        _ax.tick_params(labelsize=5, colors="white")

        _rgb = _rgb_images[_z].copy()
        for _ch in range(3):
            if not _enabled[_ch]:
                _rgb[:, :, _ch] = 0
        _ax.imshow(
            _rgb,
            extent=[_x_range[0], _x_range[1], _y_range[1], _y_range[0]],
            aspect="auto",
            interpolation="bilinear",
            origin="upper",
        )
        draw_region_overlays(_ax, _z, boundaries, centroids, _x_range, _y_range)

    _labels = [g for g, e in zip(_rgb_channels, _enabled) if e]
    _fig.patch.set_facecolor("black")
    _fig.suptitle(
        f"RGB: {' / '.join(_labels) if _labels else '(none)'}",
        fontsize=11, color="white",
    )
    _fig.tight_layout(rect=[0, 0, 1.0, 0.96])
    plt.gca()
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## Tanycyte territory map
    """)
    return


@app.cell
def _(
    bbox_df,
    bbox_slices,
    boundaries,
    centroids,
    draw_region_overlays,
    gaussian_filter,
    make_territory_axes,
    mpatches,
    np,
    plt,
    smooth_sigma_slider,
    x_centers,
    x_edges,
    x_range_slider,
    y_centers,
    y_edges,
    y_range_slider,
):
    _x_range = tuple(x_range_slider.value)
    _y_range = tuple(y_range_slider.value)
    _sigma = smooth_sigma_slider.value

    # Filter to tanycytes
    _tany = bbox_df[bbox_df["subclass"].str.contains("Tanycyte", case=False)]
    _supertypes = sorted(_tany["supertype"].unique())
    _color_map = _tany.drop_duplicates("supertype").set_index("supertype")["supertype_color"].to_dict()

    _fig, _axes, _nrows, _ncols = make_territory_axes(
        len(bbox_slices), _x_range, _y_range,
    )

    _shown = {}
    for _i, _z in enumerate(bbox_slices):
        _ax = _axes[_i // _ncols][_i % _ncols]
        _ax.set_facecolor("#f0f0f0")
        _ax.set_title(f"z={_z:.2f}", fontsize=8)
        _ax.tick_params(labelsize=5)

        _slice_tany = _tany[_tany["z_slice"] == _z]

        for _st in _supertypes:
            _st_df = _slice_tany[_slice_tany["supertype"] == _st]
            if len(_st_df) < 3:
                continue
            _h, _, _ = np.histogram2d(
                _st_df["x"].values, _st_df["y"].values,
                bins=[x_edges, y_edges],
            )
            if _sigma > 0:
                _h = gaussian_filter(_h, sigma=_sigma)
            _nonzero = _h[_h > 0]
            if len(_nonzero) == 0:
                continue
            _thresh = np.percentile(_nonzero, 20)
            _color = _color_map.get(_st, "#888888")
            _ax.contourf(
                x_centers, y_centers, _h.T,
                levels=[_thresh, _h.max() * 10],
                colors=[_color], alpha=0.3,
            )
            _ax.contour(
                x_centers, y_centers, _h.T,
                levels=[_thresh], colors=[_color], linewidths=1.0,
            )
            _shown[_st] = _color

        # Scatter individual tanycytes on top
        if len(_slice_tany) > 0:
            _colors = [_color_map.get(st, "#888888") for st in _slice_tany["supertype"]]
            _ax.scatter(
                _slice_tany["x"].values, _slice_tany["y"].values,
                c=_colors, s=1, alpha=0.5, linewidths=0,
            )

        draw_region_overlays(_ax, _z, boundaries, centroids, _x_range, _y_range)

    _fig.suptitle(
        f"Tanycyte territories ({len(_tany):,} cells in bbox)", fontsize=11,
    )
    _fig.tight_layout(rect=[0, 0.04, 1.0, 0.96])
    _patches = [mpatches.Patch(color=c, label=l.split(" ", 1)[-1]) for l, c in sorted(_shown.items())]
    if _patches:
        _fig.legend(
            handles=_patches, loc="lower center", ncol=len(_patches),
            fontsize=6, frameon=False, bbox_to_anchor=(0.5, -0.02),
        )
    plt.gca()
    return


@app.cell
def _(mo):
    mo.md(r"""
    ## cNMF program territory map
    """)
    return


@app.cell
def _(mo):
    cnmf_decomposition = mo.ui.dropdown(
        options={
            "All cells (162k)": "hypo_cnmf",
            "Neurons only (96k)": "hypo_cnmf_neurons",
            "ARH/ME/VMH MERFISH (7k)": "nmf_arh_me_vmh",
        },
        value="ARH/ME/VMH MERFISH (7k)",
        label="Decomposition",
    )
    cnmf_decomposition
    return (cnmf_decomposition,)


@app.cell
def _(Path, cnmf_decomposition, mo):
    import re as _re
    _run_name = cnmf_decomposition.value
    _cnmf_dir = Path(mo.notebook_dir()) / ".." / "data" / "processed" / "mouse_abc" / "cnmf" / _run_name
    _available_ks = sorted({
        int(m.group(1))
        for f in _cnmf_dir.glob(f"{_run_name}.usages.k_*.dt_0_1.consensus.txt")
        if (m := _re.search(r"\.k_(\d+)\.", f.name))
    })
    cnmf_k_selector = mo.ui.dropdown(
        options={str(k): k for k in _available_ks},
        value=str(_available_ks[-1]) if _available_ks else None,
        label="K",
    )
    cnmf_k_selector
    return (cnmf_k_selector,)


@app.cell
def _(Path, cells_df, cnmf_decomposition, cnmf_k_selector, mo, pd):
    _root = Path(mo.notebook_dir()) / ".."
    _data = _root / "data" / "processed" / "mouse_abc"
    _run_name = cnmf_decomposition.value
    _k = cnmf_k_selector.value
    _cnmf_dir = _data / "cnmf" / _run_name

    # Load usages
    _usages = pd.read_csv(
        _cnmf_dir / f"{_run_name}.usages.k_{_k}.dt_0_1.consensus.txt",
        sep="\t", index_col=0,
    )
    _usages.columns = [f"P{c}" for c in _usages.columns]
    _usages_norm = _usages.div(_usages.sum(axis=1), axis=0)
    _programs = [c for c in _usages_norm.columns if c.startswith("P")]

    # Detect MERFISH vs 10X
    _sample_ids = set(_usages_norm.index[:100])
    _cell_ids = set(cells_df["cell_id"])
    _is_merfish = len(_sample_ids & _cell_ids) > len(_sample_ids) * 0.5

    if _is_merfish:
        _coords = cells_df[["cell_id", "cluster"]].set_index("cell_id")
        _merged = _usages_norm[_programs].join(_coords[["cluster"]], how="inner")
        _cluster_usage = _merged.groupby("cluster")[_programs].mean()
    else:
        from abc_atlas_access.abc_atlas_cache.abc_project_cache import AbcProjectCache
        _cache = AbcProjectCache.from_cache_dir(_root / "data" / "raw" / "abc_atlas_cache")
        _cache.load_latest_manifest()
        _cell_meta = _cache.get_metadata_dataframe(
            directory="WMB-10X", file_name="cell_metadata", dtype={"cell_label": str},
        ).set_index("cell_label")
        _membership = _cache.get_metadata_dataframe(
            directory="WMB-taxonomy", file_name="cluster_to_cluster_annotation_membership",
        )
        _taxonomy = _membership.pivot_table(
            index="cluster_alias", columns="cluster_annotation_term_set_name",
            values="cluster_annotation_term_name", aggfunc="first",
        )
        _cell_meta = _cell_meta.join(_taxonomy, on="cluster_alias")
        _cell_meta = _cell_meta[_cell_meta["feature_matrix_label"] == "WMB-10Xv3-HY"]
        _merged = _usages_norm[_programs].join(_cell_meta[["cluster"]], how="inner")
        _cluster_usage = _merged.groupby("cluster")[_programs].mean()

    # Map cluster-level usage onto all MERFISH cells in bbox
    cnmf_cell_usage = cells_df[["cell_id", "cluster", "x", "y", "z_slice"]].merge(
        _cluster_usage, left_on="cluster", right_index=True, how="inner",
    )

    # Load annotations
    _annot_path = _data / "cnmf" / f"program_annotations_{_run_name}_k{_k}.csv"
    if not _annot_path.exists():
        _annot_path = _data / "cnmf" / f"program_annotations_{_run_name}.csv"
    if not _annot_path.exists() and _run_name == "hypo_cnmf":
        _annot_path = _data / "cnmf" / "program_annotations.csv"

    cnmf_annotations = {}
    cnmf_program_labels = {p: p for p in _programs}
    if _annot_path.exists():
        _annot_df = pd.read_csv(_annot_path).set_index("program")
        for p in _programs:
            if p in _annot_df.index:
                cnmf_program_labels[p] = f"{p}: {_annot_df.loc[p, 'name']}"
                cnmf_annotations[p] = _annot_df.loc[p].to_dict()

    cnmf_programs = _programs
    mo.md(f"**cNMF:** {len(cnmf_cell_usage):,} MERFISH cells mapped, {len(_programs)} programs (K={_k})")
    return (
        cnmf_annotations,
        cnmf_cell_usage,
        cnmf_program_labels,
        cnmf_programs,
    )


@app.cell
def _(cnmf_program_labels, cnmf_programs, mo):
    cnmf_checkboxes = mo.ui.dictionary({
        p: mo.ui.checkbox(value=i < 5, label=cnmf_program_labels[p])
        for i, p in enumerate(cnmf_programs)
    })
    cnmf_thresh_slider = mo.ui.slider(
        start=50, stop=99, step=1, value=85,
        label="Contour threshold (percentile)",
    )
    cnmf_smooth_slider = mo.ui.slider(
        start=0.02, stop=0.3, step=0.01, value=0.1,
        label="Smoothing σ (mm)",
    )
    mo.vstack([mo.hstack([cnmf_thresh_slider, cnmf_smooth_slider]), cnmf_checkboxes])
    return cnmf_checkboxes, cnmf_smooth_slider, cnmf_thresh_slider


@app.cell
def _(cnmf_checkboxes):
    cnmf_selected_programs = [p for p, checked in cnmf_checkboxes.value.items() if checked]
    return (cnmf_selected_programs,)


@app.cell(hide_code=True)
def _(
    bbox_slices,
    boundaries,
    centroids,
    cnmf_annotations,
    cnmf_cell_usage,
    cnmf_program_labels,
    cnmf_selected_programs,
    cnmf_smooth_slider,
    cnmf_thresh_slider,
    draw_region_overlays,
    gaussian_filter,
    make_territory_axes,
    mo,
    np,
    plt,
    smooth_sigma_slider,
    x_centers,
    x_edges,
    x_range_slider,
    y_centers,
    y_edges,
    y_range_slider,
):
    _x_range = tuple(x_range_slider.value)
    _y_range = tuple(y_range_slider.value)
    _selected = cnmf_selected_programs
    _sigma = smooth_sigma_slider.value

    _output = mo.md("*Select at least one program above.*")
    if _selected:
        _cmap = plt.cm.tab20
        _prog_colors = {p: _cmap(i % 20) for i, p in enumerate(_selected)}

        _fig, _axes, _nrows, _ncols = make_territory_axes(
            len(bbox_slices), _x_range, _y_range,
        )

        for _i, _z in enumerate(bbox_slices):
            _ax = _axes[_i // _ncols][_i % _ncols]
            _ax.set_facecolor("#f0f0f0")
            _ax.set_title(f"z={_z:.2f}", fontsize=8)
            _ax.tick_params(labelsize=5)

            _slice = cnmf_cell_usage[cnmf_cell_usage["z_slice"] == _z]
            if len(_slice) == 0:
                continue

            # Build padded grid for smoothing (avoids edge artifacts)
            _bin_size = x_edges[1] - x_edges[0]
            _sigma_bins = cnmf_smooth_slider.value / _bin_size
            _pad = int(np.ceil(_sigma_bins * 3))
            _x_edges_pad = np.arange(
                x_edges[0] - _pad * _bin_size,
                x_edges[-1] + (_pad + 1) * _bin_size,
                _bin_size,
            )
            _y_edges_pad = np.arange(
                y_edges[0] - _pad * _bin_size,
                y_edges[-1] + (_pad + 1) * _bin_size,
                _bin_size,
            )

            # Use ALL cells in this slice (not just bbox) for padded histogram
            _slice_all = cnmf_cell_usage[cnmf_cell_usage["z_slice"] == _z]

            for _p in _selected:
                if _p not in _slice_all.columns:
                    continue
                _vals = _slice_all[_p].values
                _sx = _slice_all["x"].values
                _sy = _slice_all["y"].values

                _sum_h, _, _ = np.histogram2d(_sx, _sy, bins=[_x_edges_pad, _y_edges_pad], weights=_vals)
                _cnt_h, _, _ = np.histogram2d(_sx, _sy, bins=[_x_edges_pad, _y_edges_pad])
                with np.errstate(invalid="ignore"):
                    _mean_h = np.where(_cnt_h > 0, _sum_h / _cnt_h, 0)
                _mean_h = gaussian_filter(_mean_h, sigma=_sigma_bins)

                # Crop back to original bbox grid
                _mean_h = _mean_h[_pad:_pad + len(x_centers), _pad:_pad + len(y_centers)]

                _nonzero = _mean_h[_mean_h > 0]
                if len(_nonzero) == 0:
                    continue
                _thresh = np.percentile(_nonzero, cnmf_thresh_slider.value)
                _color = _prog_colors[_p]
                _ax.contourf(
                    x_centers, y_centers, _mean_h.T,
                    levels=[_thresh, _mean_h.max() * 10],
                    colors=[_color], alpha=0.2,
                )
                _ax.contour(
                    x_centers, y_centers, _mean_h.T,
                    levels=[_thresh], colors=[_color], linewidths=0.8,
                )

            draw_region_overlays(_ax, _z, boundaries, centroids, _x_range, _y_range)

        _fig.suptitle("cNMF program territories", fontsize=11)

        # Remove tick labels from all subplots for consistency
        for _row in _axes:
            for _ax2 in _row:
                _ax2.tick_params(labelbottom=False, labelleft=False)

        _fig.tight_layout(rect=[0, 0, 1.0, 0.96])

        # Find rightmost subplot edge in figure coords to place labels tight
        _right_edge = max(
            _axes[_i // _ncols][_i % _ncols].get_position().x1
            for _i in range(len(bbox_slices))
        )

        # Build description sidebar text just past the subplots
        _desc_lines = []
        for _p in _selected:
            _color = _prog_colors[_p]
            _label = cnmf_program_labels.get(_p, _p)
            _desc = ""
            if _p in cnmf_annotations:
                _desc = cnmf_annotations[_p].get("description", "")
                if len(_desc) > 80:
                    _desc = _desc[:77] + "..."
            _desc_lines.append((_label, _desc, _color))

        _x_text = _right_edge + 0.02
        _y_start = 0.88
        for _j, (_label, _desc, _color) in enumerate(_desc_lines):
            _y_pos = _y_start - _j * 0.045
            _fig.text(
                _x_text, _y_pos, _label,
                fontsize=6, fontweight="bold", color=_color,
                transform=_fig.transFigure, va="top",
            )
            if _desc:
                _fig.text(
                    _x_text, _y_pos - 0.018, _desc,
                    fontsize=5, color="#444444",
                    transform=_fig.transFigure, va="top",
                )

        _output = plt.gca()
    _output
    return


if __name__ == "__main__":
    app.run()
