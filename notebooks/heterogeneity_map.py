# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "hypomap",
#     "auto-mix-prep",
#     "marimo>=0.19.11",
#     "matplotlib",
#     "numpy",
#     "pandas",
#     "plotly",
#     "pyarrow",
#     "scipy",
# ]
#
# [tool.uv.sources]
# hypomap = { path = "../", editable = true }
#
# [tool.hatch.build.targets.wheel]
# packages = ["hypomap"]
# ///

import marimo

__generated_with = "0.19.4"
app = marimo.App(width="full")


@app.cell
def _():
    import math
    from pathlib import Path

    import marimo as mo
    import msgpack
    import numpy as np
    import pandas as pd
    import plotly.graph_objects as go
    from hypomap.diversity import (
        build_diversity_sats,
        build_local_diversity_grid,
        find_max_mean_diversity_box,
        mask_excluded_regions,
    )
    from plotly.subplots import make_subplots
    from scipy.spatial import cKDTree
    return cKDTree, go, make_subplots, math, mo, msgpack, np, pd


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Subcortical Heterogeneity Map

    Interactive visualization of **cell-type diversity** across subcortical structures
    (HY, TH, STR, PAL) in the mouse ABC atlas (~746k cells).

    For each coronal slice, diversity is computed within an adjustable integration radius
    at a chosen taxonomic granularity (class / subclass / supertype / cluster / dominant neuropeptide).

    ## Metrics

    | Metric | Formula | Interpretation |
    |--------|---------|---------------|
    | Shannon diversity | H = -sum(p_i * log2(p_i)) | Sensitive to rare types; higher = more diverse |
    | Simpson diversity | D = 1 - sum(p_i^2) | Probability two random cells are different types |
    | Type count | len(unique types in radius) | Simple count of distinct types at chosen granularity |
    """)
    return


@app.cell
def _(mo, msgpack, np, pd):
    # Load cell data
    _data_dir = mo.notebook_dir() / ".." / "data" / "processed" / "mouse_abc_whole_brain"
    cells_df = pd.read_parquet(_data_dir / "cells_with_coords.parquet")

    # Filter unassigned regions
    cells_df = cells_df[
        ~cells_df["region"].str.contains("unassigned", case=False)
    ].copy()
    cells_df["z_slice"] = cells_df["z"].round(2)

    # Compute dominant neuropeptide per cell
    np_expr = pd.read_parquet(_data_dir / "neuropeptide_expression.parquet")
    _np_vals = np_expr.loc[cells_df["cell_id"]].values.astype(np.float32)
    _np_genes = np.array(np_expr.columns.tolist())

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

    del _np_vals, _zscored  # free memory

    # Load neuropeptide-receptor pair map
    np_map = pd.read_csv(
        mo.notebook_dir() / ".." / "data" / "processed" / "mouse_common" / "np_map.csv"
    )

    # Load region boundaries (compact msgpack with float32 binary arrays)
    with open(_data_dir / "coronal_atlas_regions.msgpack", "rb") as _f:
        _region_data = msgpack.unpackb(_f.read(), raw=False)

    all_slices = sorted(cells_df["z_slice"].unique())
    boundaries = {
        float(k): {
            r: np.frombuffer(buf, dtype=np.float32).reshape(-1, 2).tolist()
            for r, buf in regions.items()
        }
        for k, regions in _region_data["boundaries"].items()
    }
    centroids = {
        float(k): {
            r: tuple(np.frombuffer(buf, dtype=np.float32).tolist())
            for r, buf in regions.items()
        }
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
    return all_slices, boundaries, cells_df, centroids, np_expr, np_map


@app.cell(hide_code=True)
def _(cKDTree, np):
    def compute_diversity(labels, metric):
        """Compute a diversity metric from an array of integer labels.

        Args:
            labels: 1D array of integer category labels
            metric: one of 'Shannon', 'Simpson', 'Type count'

        Returns:
            float diversity value, or NaN if fewer than 2 labels
        """
        if len(labels) < 2:
            return float("nan")

        _counts = np.bincount(labels)
        _counts = _counts[_counts > 0]
        _n = _counts.sum()

        if metric == "Type count":
            return float(len(_counts))

        _p = _counts / _n
        if metric == "Shannon":
            return float(-np.sum(_p * np.log2(_p)))
        elif metric == "Simpson":
            return float(1.0 - np.sum(_p**2))
        return float("nan")

    def compute_heterogeneity_grid(xy, int_labels, radius, metric, grid_res):
        """Compute diversity on a regular grid over cell positions.

        Args:
            xy: (N, 2) array of cell positions
            int_labels: (N,) integer-encoded category labels
            radius: integration radius in mm
            metric: diversity metric name
            grid_res: number of grid points per axis

        Returns:
            (grid_x, grid_y, grid_values) arrays for plotting
        """
        _xmin, _xmax = xy[:, 0].min(), xy[:, 0].max()
        _ymin, _ymax = xy[:, 1].min(), xy[:, 1].max()
        _gx = np.linspace(_xmin, _xmax, grid_res)
        _gy = np.linspace(_ymin, _ymax, grid_res)
        _gxx, _gyy = np.meshgrid(_gx, _gy)
        _grid_points = np.column_stack([_gxx.ravel(), _gyy.ravel()])

        _tree = cKDTree(xy)
        _neighbors = _tree.query_ball_point(_grid_points, radius)

        _values = np.full(len(_grid_points), float("nan"))
        for _i, _nbrs in enumerate(_neighbors):
            if len(_nbrs) >= 2:
                _values[_i] = compute_diversity(int_labels[_nbrs], metric)

        return _grid_points[:, 0], _grid_points[:, 1], _values

    def compute_heterogeneity_percell(xy, int_labels, radius, metric):
        """Compute diversity at each cell position.

        Args:
            xy: (N, 2) array of cell positions
            int_labels: (N,) integer-encoded category labels
            radius: integration radius in mm
            metric: diversity metric name

        Returns:
            (cell_x, cell_y, cell_values) arrays for plotting
        """
        _tree = cKDTree(xy)
        _neighbors = _tree.query_ball_point(xy, radius)

        _values = np.full(len(xy), float("nan"))
        for _i, _nbrs in enumerate(_neighbors):
            if len(_nbrs) >= 2:
                _values[_i] = compute_diversity(int_labels[_nbrs], metric)

        return xy[:, 0], xy[:, 1], _values
    return (
        compute_diversity,
        compute_heterogeneity_grid,
        compute_heterogeneity_percell,
    )


@app.cell(hide_code=True)
def _(
    all_slices,
    cells_df,
    compute_heterogeneity_grid,
    compute_heterogeneity_percell,
    granularity_dropdown,
    grid_res_slider,
    metric_dropdown,
    mode_dropdown,
    np,
    pd,
    radius_slider,
):
    _radius = radius_slider.value
    _metric = metric_dropdown.value
    _granularity = granularity_dropdown.value
    _mode = mode_dropdown.value
    _grid_res = grid_res_slider.value

    # Integer-encode labels once for the entire dataset
    _codes, _uniques = pd.factorize(cells_df[_granularity])
    _int_labels = _codes.astype(np.int32)

    _raw_xy = cells_df[["x", "y"]].values
    _raw_z = cells_df["z_slice"].values

    slice_results = {}
    slice_images = {}
    _all_vals = []

    for _z in all_slices:
        _mask = _raw_z == _z
        _xy = _raw_xy[_mask]
        _labels = _int_labels[_mask]

        if len(_xy) < 10:
            continue

        if _mode == "Grid":
            _sx, _sy, _sv = compute_heterogeneity_grid(
                _xy, _labels, _radius, _metric, _grid_res
            )
            # Store 2D image before NaN filtering (meshgrid shape is grid_res×grid_res)
            _xmin_s, _xmax_s = float(_xy[:, 0].min()), float(_xy[:, 0].max())
            _ymin_s, _ymax_s = float(_xy[:, 1].min()), float(_xy[:, 1].max())
            slice_images[_z] = (
                _sv.reshape(_grid_res, _grid_res),
                [_xmin_s, _xmax_s, _ymin_s, _ymax_s],
            )
        else:
            _sx, _sy, _sv = compute_heterogeneity_percell(
                _xy, _labels, _radius, _metric
            )

        # Filter NaN
        _valid = ~np.isnan(_sv)
        _sx = _sx[_valid]
        _sy = _sy[_valid]
        _sv = _sv[_valid]

        if len(_sx) > 0:
            slice_results[_z] = (_sx, _sy, _sv)
            _all_vals.append(_sv)

    # Global color scale from 2nd to 98th percentile
    if _all_vals:
        _combined = np.concatenate(_all_vals)
        color_min = float(np.percentile(_combined, 2))
        color_max = float(np.percentile(_combined, 98))
    else:
        color_min = 0.0
        color_max = 1.0

    print(
        f"Computed {_metric} ({_granularity}) for {len(slice_results)} slices, "
        f"mode={_mode}, radius={_radius}mm"
    )
    print(f"Color range: [{color_min:.2f}, {color_max:.2f}]")
    return color_max, color_min, slice_images, slice_results


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    ## Multi-Slice Overview

    Each panel shows one coronal slice. Color = diversity metric value (Viridis scale).
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    radius_slider = mo.ui.slider(
        start=0.05, stop=0.5, step=0.05, value=0.15, label="Integration radius (mm)"
    )
    metric_dropdown = mo.ui.dropdown(
        options=["Shannon", "Simpson", "Type count"],
        value="Shannon",
        label="Diversity metric",
    )
    granularity_dropdown = mo.ui.dropdown(
        options=["class", "subclass", "supertype", "cluster", "dominant_neuropeptide"],
        value="subclass",
        label="Taxonomic granularity",
    )
    mode_dropdown = mo.ui.dropdown(
        options=["Grid", "Per-cell"], value="Grid", label="Computation mode"
    )
    grid_res_slider = mo.ui.slider(
        start=20, stop=200, step=10, value=100, label="Grid resolution"
    )

    mo.md(
        f"""
    ## Controls

    {mo.hstack([radius_slider, metric_dropdown, granularity_dropdown], justify="start")}

    {mo.hstack([mode_dropdown, grid_res_slider], justify="start")}

    **Grid mode** evaluates diversity on a regular grid (fast). **Per-cell mode** evaluates at each cell position (slower but gives per-cell values).
    """
    )
    return (
        granularity_dropdown,
        grid_res_slider,
        metric_dropdown,
        mode_dropdown,
        radius_slider,
    )


@app.cell(hide_code=True)
def _(
    color_max,
    color_min,
    go,
    granularity_dropdown,
    make_subplots,
    math,
    metric_dropdown,
    np,
    radius_slider,
    slice_results,
):
    # Multi-slice overview figure
    _n_slices = len(slice_results)
    _ncols = 4
    _nrows = math.ceil(_n_slices / _ncols)

    _sorted_slices = sorted(slice_results.keys())

    _fig = make_subplots(
        rows=_nrows,
        cols=_ncols,
        subplot_titles=[f"y={z:.2f}" for z in _sorted_slices],
        horizontal_spacing=0.03,
        vertical_spacing=0.06,
    )

    # Compute aspect ratio from data to set figure height
    _all_x = np.concatenate([slice_results[z][0] for z in _sorted_slices])
    _all_y = np.concatenate([slice_results[z][1] for z in _sorted_slices])
    _x_range = _all_x.max() - _all_x.min()
    _y_range = _all_y.max() - _all_y.min()
    _aspect = _y_range / _x_range if _x_range > 0 else 1.0

    # Build traces as plain dicts for performance
    _all_traces = []
    for _i, _z in enumerate(_sorted_slices):
        _sx, _sy, _sv = slice_results[_z]
        _row = _i // _ncols + 1
        _col = _i % _ncols + 1
        _xaxis = f"x{_i + 1}" if _i > 0 else "x"
        _yaxis = f"y{_i + 1}" if _i > 0 else "y"

        # Subsample if too many points for overview
        if len(_sx) > 5000:
            _idx = np.random.default_rng(42).choice(len(_sx), 5000, replace=False)
            _sx = _sx[_idx]
            _sy = _sy[_idx]
            _sv = _sv[_idx]

        _all_traces.append(
            {
                "type": "scattergl",
                "x": np.round(_sx, 3).tolist(),
                "y": np.round(_sy, 3).tolist(),
                "mode": "markers",
                "marker": {
                    "size": 3,
                    "color": np.round(_sv, 3).tolist(),
                    "colorscale": "Viridis",
                    "cmin": color_min,
                    "cmax": color_max,
                    "showscale": (_i == 0),
                    "colorbar": (
                        {
                            "title": metric_dropdown.value,
                            "len": 0.5,
                        }
                        if _i == 0
                        else None
                    ),
                },
                "showlegend": False,
                "hovertemplate": "x=%{x:.2f}, y=%{y:.2f}<br>value=%{marker.color:.3f}<extra></extra>",
                "xaxis": _xaxis,
                "yaxis": _yaxis,
            }
        )

    # Build figure dict (bypass go.Figure overhead)
    _fig_dict = _fig.to_dict()
    _fig_dict["data"] = _all_traces

    _col_width = 230
    _row_height = int(_col_width * _aspect)
    _fig_height = min(16000, _nrows * _row_height + 100)

    _fig_dict["layout"]["height"] = _fig_height
    _fig_dict["layout"]["width"] = _ncols * _col_width + 150
    _fig_dict["layout"]["title"] = {
        "text": f"{metric_dropdown.value} diversity ({granularity_dropdown.value}), "
        f"r={radius_slider.value}mm"
    }
    _fig_dict["layout"]["plot_bgcolor"] = "#f8f8f8"

    # Match y axes to x axes per subplot for 1:1 aspect
    for _i in range(len(_sorted_slices)):
        _xa = f"xaxis{_i + 1}" if _i > 0 else "xaxis"
        _ya = f"yaxis{_i + 1}" if _i > 0 else "yaxis"
        if _ya in _fig_dict["layout"]:
            _fig_dict["layout"][_ya]["scaleanchor"] = f"x{_i + 1}" if _i > 0 else "x"
            _fig_dict["layout"][_ya]["autorange"] = "reversed"

    overview_fig = go.Figure(_fig_dict)
    overview_fig
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    ## Single-Slice Detail

    Select a slice for a detailed view with region boundaries and labels.
    """)
    return


@app.cell
def _():
    return


@app.cell
def _(all_slices, mo):
    detail_slice_slider = mo.ui.slider(
        steps=[round(float(z), 2) for z in all_slices],
        value=round(float(all_slices[len(all_slices) // 2]), 2),
        label="Detail slice (y mm)",
        show_value=True,
        debounce=True,
    )
    mo.md(f"**Select slice** (drag or use arrow keys): {detail_slice_slider}")
    return (detail_slice_slider,)


@app.cell(hide_code=True)
def _(
    boundaries,
    centroids,
    color_max,
    color_min,
    detail_slice_slider,
    go,
    granularity_dropdown,
    metric_dropdown,
    np,
    radius_slider,
    slice_results,
):
    _z_target = float(detail_slice_slider.value)
    _metric_name = metric_dropdown.value
    _gran_name = granularity_dropdown.value

    # Find nearest slice key (avoids float32/float64 mismatch)
    _z_keys = sorted(slice_results.keys())
    _z = min(_z_keys, key=lambda k: abs(float(k) - _z_target)) if _z_keys else None

    detail_fig = None
    if _z is not None and _z in slice_results:
        _sx, _sy, _sv = slice_results[_z]

        _traces = [
            {
                "type": "scattergl",
                "x": np.round(_sx, 3).tolist(),
                "y": np.round(_sy, 3).tolist(),
                "mode": "markers",
                "marker": {
                    "size": 8,
                    "color": np.round(_sv, 3).tolist(),
                    "colorscale": "Viridis",
                    "cmin": color_min,
                    "cmax": color_max,
                    "showscale": True,
                    "colorbar": {"title": _metric_name},
                },
                "showlegend": False,
                "hovertemplate": "x=%{x:.2f}, y=%{y:.2f}<br>value=%{marker.color:.3f}<extra></extra>",
            }
        ]

        # Find nearest boundary/centroid key (avoids float mismatch with msgpack keys)
        _bnd_keys = sorted(boundaries.keys())
        _z_bnd = min(_bnd_keys, key=lambda k: abs(float(k) - float(_z))) if _bnd_keys else None

        _shapes = []
        if _z_bnd is not None and abs(float(_z_bnd) - float(_z)) < 0.01:
            for _region, _coords in boundaries[_z_bnd].items():
                _path = "M " + " L ".join(f"{c[0]},{c[1]}" for c in _coords) + " Z"
                _shapes.append(
                    {
                        "type": "path",
                        "path": _path,
                        "line": {"color": "rgba(100,100,100,0.6)", "width": 1.5},
                        "fillcolor": "rgba(0,0,0,0)",
                    }
                )

        _cen_keys = sorted(centroids.keys())
        _z_cen = min(_cen_keys, key=lambda k: abs(float(k) - float(_z))) if _cen_keys else None

        _annotations = []
        if _z_cen is not None and abs(float(_z_cen) - float(_z)) < 0.01:
            for _region, (_cx, _cy) in centroids[_z_cen].items():
                _annotations.append(
                    {
                        "x": _cx,
                        "y": _cy,
                        "text": _region,
                        "showarrow": False,
                        "font": {"size": 9, "color": "white"},
                        "bgcolor": "rgba(0,0,0,0.5)",
                        "borderpad": 2,
                    }
                )

        _x_range = [float(np.min(_sx)) - 0.2, float(np.max(_sx)) + 0.2]

        detail_fig = go.Figure(
            {
                "data": _traces,
                "layout": {
                    "title": {
                        "text": f"{_metric_name} diversity ({_gran_name}), z={_z}, "
                        f"r={radius_slider.value}mm"
                    },
                    "xaxis": {"title": "x_ras (mm)", "range": _x_range},
                    "yaxis": {
                        "title": "z_ras (mm)",
                        "scaleanchor": "x",
                        "autorange": "reversed",
                    },
                    "shapes": _shapes,
                    "annotations": _annotations,
                    "height": 600,
                    "width": 900,
                    "plot_bgcolor": "#f8f8f8",
                },
            }
        )
    else:
        print(f"No data for slice z={_z}")

    detail_fig
    return


@app.cell
def _(
    boundaries,
    centroids,
    color_max,
    color_min,
    granularity_dropdown,
    math,
    metric_dropdown,
    np,
    radius_slider,
    slice_images,
    slice_results,
):
    import matplotlib

    matplotlib.rcParams["figure.dpi"] = 150
    import matplotlib.colors as mcolors
    import matplotlib.pyplot as plt

    _metric_name = metric_dropdown.value
    _gran_name = granularity_dropdown.value
    _sorted_slices = sorted(slice_results.keys())
    _n = len(_sorted_slices)
    _use_images = len(slice_images) > 0

    # Global extent across all slices for locked axes
    _all_x = np.concatenate([slice_results[z][0] for z in _sorted_slices])
    _all_y = np.concatenate([slice_results[z][1] for z in _sorted_slices])
    _pad = 0.15
    _gx0 = float(_all_x.min()) - _pad
    _gx1 = float(_all_x.max()) + _pad
    _gy0 = float(_all_y.min()) - _pad
    _gy1 = float(_all_y.max()) + _pad

    _ncols = 4
    _nrows = math.ceil(_n / _ncols)
    _data_aspect = (_gy1 - _gy0) / (_gx1 - _gx0)
    _subplot_w = 4
    _subplot_h = _subplot_w * _data_aspect

    _fig, _axes = plt.subplots(
        _nrows,
        _ncols,
        figsize=(_ncols * _subplot_w, _nrows * _subplot_h + 0.8),
        squeeze=False,
        sharex=True,
        sharey=True,
    )

    _norm = mcolors.Normalize(vmin=color_min, vmax=color_max)
    _cmap = plt.cm.viridis.copy()
    _cmap.set_bad(color="#f0f0f0")
    _bnd_keys = sorted(boundaries.keys())
    _cen_keys = sorted(centroids.keys())

    for _i, _z in enumerate(_sorted_slices):
        _ax = _axes[_i // _ncols][_i % _ncols]

        if _use_images and _z in slice_images:
            # Grid mode: use pre-built 2D array directly
            _img, (_xmin_s, _xmax_s, _ymin_s, _ymax_s) = slice_images[_z]
            _nr, _nc = _img.shape
            _half_dx = (_xmax_s - _xmin_s) / max(1, _nc - 1) / 2
            _half_dy = (_ymax_s - _ymin_s) / max(1, _nr - 1) / 2
            _ax.imshow(
                _img,
                cmap=_cmap,
                norm=_norm,
                extent=[
                    _xmin_s - _half_dx,
                    _xmax_s + _half_dx,
                    _ymax_s + _half_dy,
                    _ymin_s - _half_dy,
                ],
                aspect="equal",
                interpolation="nearest",
                origin="upper",
            )
        else:
            # Per-cell mode: scatter fallback
            _sx, _sy, _sv = slice_results[_z]
            _ax.scatter(
                _sx,
                _sy,
                c=_sv,
                s=4,
                cmap="viridis",
                vmin=color_min,
                vmax=color_max,
                edgecolors="none",
                rasterized=True,
            )

        # Region boundaries (nearest-key match to handle float mismatch)
        _z_bnd = min(_bnd_keys, key=lambda k: abs(float(k) - float(_z))) if _bnd_keys else None
        if _z_bnd is not None and abs(float(_z_bnd) - float(_z)) < 0.01:
            for _region, _coords in boundaries[_z_bnd].items():
                _xs = [c[0] for c in _coords] + [_coords[0][0]]
                _ys = [c[1] for c in _coords] + [_coords[0][1]]
                _ax.plot(_xs, _ys, color="gray", linewidth=0.5, alpha=0.6)

        # Region labels
        _z_cen = min(_cen_keys, key=lambda k: abs(float(k) - float(_z))) if _cen_keys else None
        if _z_cen is not None and abs(float(_z_cen) - float(_z)) < 0.01:
            for _region, (_cx, _cy) in centroids[_z_cen].items():
                _ax.text(
                    _cx,
                    _cy,
                    _region,
                    fontsize=5,
                    color="white",
                    ha="center",
                    va="center",
                    bbox=dict(boxstyle="round,pad=0.15", fc="black", alpha=0.5, lw=0),
                )

        _ax.set_title(f"y={_z:.2f}", fontsize=9)
        _ax.tick_params(labelsize=6)
        for _spine in _ax.spines.values():
            _spine.set_visible(False)

    # Set shared limits (y inverted) on the shared axes
    _axes[0][0].set_xlim(_gx0, _gx1)
    _axes[0][0].set_ylim(_gy1, _gy0)
    _axes[0][0].set_aspect("equal")

    # Hide unused axes
    for _i in range(_n, _nrows * _ncols):
        _axes[_i // _ncols][_i % _ncols].set_visible(False)

    _fig.suptitle(
        f"{_metric_name} diversity ({_gran_name}), r={radius_slider.value}mm",
        fontsize=12,
    )
    _fig.tight_layout(rect=[0, 0, 1.0, 0.96])

    # Floating colorbar in the bottom-right corner
    _cbar_ax = _fig.add_axes([0.77, 0.01, 0.2, 0.015])
    _fig.colorbar(
        plt.cm.ScalarMappable(norm=_norm, cmap=_cmap),
        cax=_cbar_ax,
        orientation="horizontal",
        label=_metric_name,
    )
    _cbar_ax.tick_params(labelsize=6)

    plt.gca()
    return


@app.cell
def _(mo):
    mo.md(r"""
    ---
    ## Region Diversity Table

    Cell-type diversity computed over **all cells** within each anatomical unit,
    at every level of the Allen CCF hierarchy:

    - **Substructure**: finest parcellation (e.g. CEAm, MGv)
    - **Region**: parent grouping (e.g. CEA, MG)
    - **Division**: broadest level (HY, TH, STR, PAL)

    Metrics are computed at the currently selected taxonomic granularity.
    Sort any column by clicking its header.
    """)
    return


@app.cell
def _(cells_df, compute_diversity, granularity_dropdown, np, pd):
    _granularity = granularity_dropdown.value

    # Build a mapping from each cell to its hierarchy levels
    _hier_df = cells_df[
        [_granularity, "parcellation_division", "region", "parcellation_substructure"]
    ].copy()

    def _diversity_for_group(group_labels):
        """Compute all three metrics for a group of cell-type labels."""
        _codes, _ = pd.factorize(group_labels)
        _int = _codes.astype(np.int32)
        return {
            "n_cells": len(_int),
            "Shannon": compute_diversity(_int, "Shannon"),
            "Simpson": compute_diversity(_int, "Simpson"),
            "Type count": compute_diversity(_int, "Type count"),
        }

    _rows = []

    # Level 1: substructure (finest)
    for (_div, _reg, _sub), _grp in _hier_df.groupby(
        ["parcellation_division", "region", "parcellation_substructure"]
    ):
        _m = _diversity_for_group(_grp[_granularity])
        _rows.append(
            {
                "level": "substructure",
                "division": _div,
                "region": _reg,
                "substructure": _sub,
                "name": _sub,
                **_m,
            }
        )

    # Level 2: region
    for (_div, _reg), _grp in _hier_df.groupby(["parcellation_division", "region"]):
        _m = _diversity_for_group(_grp[_granularity])
        _rows.append(
            {
                "level": "region",
                "division": _div,
                "region": _reg,
                "substructure": "",
                "name": _reg,
                **_m,
            }
        )

    # Level 3: division
    for _div, _grp in _hier_df.groupby("parcellation_division"):
        _m = _diversity_for_group(_grp[_granularity])
        _rows.append(
            {
                "level": "division",
                "division": _div,
                "region": "",
                "substructure": "",
                "name": _div,
                **_m,
            }
        )

    diversity_table = pd.DataFrame(_rows)

    # Round for display
    for _c in ["Shannon", "Simpson"]:
        diversity_table[_c] = diversity_table[_c].round(3)
    diversity_table["Type count"] = diversity_table["Type count"].astype(int)

    # Sort: division rows first (by name), then regions, then substructures
    _level_order = {"division": 0, "region": 1, "substructure": 2}
    diversity_table["_level_sort"] = diversity_table["level"].map(_level_order)
    diversity_table = (
        diversity_table.sort_values(
            ["division", "_level_sort", "region", "substructure"]
        )
        .drop(columns="_level_sort")
        .reset_index(drop=True)
    )

    print(
        f"Diversity table: {len(diversity_table)} rows "
        f"(granularity: {_granularity})"
    )
    return (diversity_table,)


@app.cell
def _(diversity_table, mo):
    mo.ui.table(
        diversity_table[
            ["level", "division", "name", "n_cells", "Shannon", "Simpson", "Type count"]
        ],
        label=f"Region diversity (sortable)",
        selection=None,
    )
    return


@app.cell
def _(diversity_table, go):
    _df = diversity_table
    _colors = {"division": "#e41a1c", "region": "#377eb8", "substructure": "#4daf4a"}
    _color_list = [_colors[l] for l in _df["level"]]
    _traces = [
        {
            "type": "scattergl",
            "x": _df["n_cells"].tolist(),
            "y": _df["Shannon"].tolist(),
            "mode": "markers+text",
            "marker": {"size": 8, "color": _color_list},
            "text": [
                f"{x} > {y}"
                for x, y in zip(_df["division"].tolist(), _df["name"].tolist())
            ],
            "textposition": "top center",
            "textfont": {"size": 8},
            "showlegend": False,
            "hovertemplate": (
                "<b>%{text}</b><br>"
                "Cells: %{x:,}<br>"
                "Shannon: %{y:.3f}<extra></extra>"
            ),
        }
    ]
    go.Figure(
        {
            "data": _traces,
            "layout": {
                "title": {"text": "Region diversity vs. cell count"},
                "xaxis": {"title": "Number of cells", "type": "log"},
                "yaxis": {"title": "Shannon diversity"},
                "height": 500,
                "width": 800,
                "hovermode": "closest",
            },
        }
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    ## Interpretation Guide

    **High diversity zones** (bright/yellow in Viridis):
    - Multiple cell types intermingle within the integration radius
    - Often found at region borders or in transition zones
    - May indicate areas of functional integration across cell types

    **Low diversity zones** (dark/purple in Viridis):
    - Dominated by one or few cell types
    - Often correspond to well-defined nuclei or homogeneous structures
    - Regions like CP (caudoputamen) tend to have low class-level diversity but higher cluster-level diversity

    **Effect of granularity**:
    - **Class**: Broadest grouping. Low diversity = all cells belong to one major class (e.g., all GABAergic)
    - **Subclass**: Moderate specificity. Reveals diversity within a class
    - **Supertype**: Finer grouping. Most brain regions show moderate-to-high diversity
    - **Cluster**: Finest level. Almost everywhere shows high diversity since clusters are very specific
    - **Dominant neuropeptide**: Each cell is assigned its top neuropeptide (highest z-scored expression after filtering low-variance genes). Diversity of dominant neuropeptides reveals peptidergic heterogeneity independent of cell-type taxonomy

    **Effect of radius**:
    - **Small radius** (0.05-0.1mm): Very local neighborhoods, noisy but captures fine structure
    - **Medium radius** (0.15-0.25mm): Good balance of signal and smoothing
    - **Large radius** (0.3-0.5mm): Smooth map, averages over larger areas, good for region-level patterns

    **Tips**:
    - Start with subclass granularity and 0.15mm radius for a good initial view
    - Compare Shannon vs Simpson: if they disagree, it often means rare types are present (Shannon is more sensitive to rare types)
    - Use the multi-granularity comparison to identify regions where diversity depends strongly on resolution
    """)
    return


@app.cell
def _(cells_df, np, np_expr, np_map):
    _available_genes = set(np_expr.columns)
    _expr_vals = np_expr.loc[cells_df["cell_id"]].values.astype(np.float32)
    _gene_idx = {g: i for i, g in enumerate(np_expr.columns)}
    _threshold = 3.0

    # Group ligand/receptor genes by system
    system_masks = {}
    for _sys, _grp in np_map.groupby("System"):
        # Collect unique ligand genes available in expression data
        _lig_genes = set()
        for _g in _grp["Ligand_Gene"].dropna().unique():
            if _g in _available_genes:
                _lig_genes.add(_g)
        # Build receptor complexes: semicolon = heterodimer (AND subunits),
        # different rows within system are OR'd
        _rec_complexes = []
        for _raw in _grp["Receptor_Gene"].dropna().unique():
            _subunits = [_g.strip() for _g in _raw.split(";")]
            if all(_g in _available_genes for _g in _subunits):
                _rec_complexes.append(_subunits)
        if not _lig_genes or not _rec_complexes:
            continue
        # OR logic: sender if ANY ligand gene > threshold
        _s = np.zeros(len(cells_df), dtype=bool)
        for _g in _lig_genes:
            _s |= _expr_vals[:, _gene_idx[_g]] > _threshold
        # Receptor: AND subunits within each complex, OR across complexes
        _r = np.zeros(len(cells_df), dtype=bool)
        for _subunits in _rec_complexes:
            _complex_mask = np.ones(len(cells_df), dtype=bool)
            for _g in _subunits:
                _complex_mask &= _expr_vals[:, _gene_idx[_g]] > _threshold
            _r |= _complex_mask
        # Cells expressing both ligand and receptor count as senders only
        _r &= ~_s
        system_masks[_sys] = (_s, _r)

    npp_systems = sorted(system_masks.keys())
    print(f"NPP systems with data: {len(npp_systems)}")
    for _sys in npp_systems[:5]:
        _s, _r = system_masks[_sys]
        print(f"  {_sys}: {_s.sum()} senders, {_r.sum()} receivers")
    return npp_systems, system_masks


@app.cell(hide_code=True)
def _(cKDTree, np):
    def compute_bipartite_entropy_grid(xy, sender, receiver, radius, grid_res):
        """Compute interaction index for one neuropeptide system on a grid.

        For each grid point, compute 4 * p_sender * p_receiver where
        p_sender and p_receiver are the fractions of local cells expressing
        any ligand or receptor gene respectively. Peaks at 1.0 when both
        proportions are 0.5, drops to 0 when either population is absent.
        """
        _xmin, _xmax = xy[:, 0].min(), xy[:, 0].max()
        _ymin, _ymax = xy[:, 1].min(), xy[:, 1].max()
        _gx = np.linspace(_xmin, _xmax, grid_res)
        _gy = np.linspace(_ymin, _ymax, grid_res)
        _gxx, _gyy = np.meshgrid(_gx, _gy)
        _grid_points = np.column_stack([_gxx.ravel(), _gyy.ravel()])

        _tree = cKDTree(xy)
        _neighbors = _tree.query_ball_point(_grid_points, radius)

        _values = np.full(len(_grid_points), float("nan"))
        for _i, _nbrs in enumerate(_neighbors):
            if len(_nbrs) < 2:
                continue
            _nbrs = np.array(_nbrs)
            _n = len(_nbrs)
            _p_s = sender[_nbrs].sum() / _n
            _p_r = receiver[_nbrs].sum() / _n
            _values[_i] = 4.0 * _p_s * _p_r

        return _grid_points[:, 0], _grid_points[:, 1], _values
    return (compute_bipartite_entropy_grid,)


@app.cell(hide_code=True)
def _(mo, npp_systems):
    npp_system_dropdown = mo.ui.dropdown(
        options=npp_systems,
        value=npp_systems[0] if npp_systems else None,
        label="Neuropeptide system",
    )
    npp_radius_slider = mo.ui.slider(
        start=0.05, stop=0.5, step=0.05, value=0.15, label="NPP radius (mm)"
    )
    npp_grid_res_slider = mo.ui.slider(
        start=20, stop=200, step=10, value=100, label="NPP grid resolution"
    )
    mo.md(f"""
    ---
    ## NPP Interaction Index

    For the selected neuropeptide system, each cell is classified as a **sender**
    (expresses any ligand gene) or **receiver** (expresses any receptor gene).
    Interaction index = 4 × p_sender × p_receiver, where p_sender and p_receiver
    are local fractions. Peaks at 1.0 when both proportions are 0.5, drops to 0
    when either population is absent.

    {mo.hstack([npp_system_dropdown, npp_radius_slider, npp_grid_res_slider], justify="start")}
    """)
    return npp_grid_res_slider, npp_radius_slider, npp_system_dropdown


@app.cell(hide_code=True)
def _(
    all_slices,
    cells_df,
    compute_bipartite_entropy_grid,
    np,
    npp_grid_res_slider,
    npp_radius_slider,
    npp_system_dropdown,
    system_masks,
):
    _radius = npp_radius_slider.value
    _grid_res = npp_grid_res_slider.value
    _system = npp_system_dropdown.value

    _raw_xy = cells_df[["x", "y"]].values
    _raw_z = cells_df["z_slice"].values
    _sender_all, _receiver_all = system_masks[_system]

    npp_slice_results = {}
    npp_slice_images = {}
    _all_vals = []

    for _z in all_slices:
        _mask = _raw_z == _z
        _xy = _raw_xy[_mask]
        _s = _sender_all[_mask]
        _r = _receiver_all[_mask]

        if len(_xy) < 10:
            continue

        _sx, _sy, _sv = compute_bipartite_entropy_grid(
            _xy, _s, _r, _radius, _grid_res
        )
        _xmin_s, _xmax_s = float(_xy[:, 0].min()), float(_xy[:, 0].max())
        _ymin_s, _ymax_s = float(_xy[:, 1].min()), float(_xy[:, 1].max())
        npp_slice_images[_z] = (
            _sv.reshape(_grid_res, _grid_res),
            [_xmin_s, _xmax_s, _ymin_s, _ymax_s],
        )

        _valid = ~np.isnan(_sv)
        _sx, _sy, _sv = _sx[_valid], _sy[_valid], _sv[_valid]
        if len(_sx) > 0:
            npp_slice_results[_z] = (_sx, _sy, _sv)
            _all_vals.append(_sv)

    if _all_vals:
        _combined = np.concatenate(_all_vals)
        npp_color_min = 0.0
        npp_color_max = float(_combined.max()) if _combined.max() > 0 else 1.0
    else:
        npp_color_min, npp_color_max = 0.0, 1.0

    print(
        f"NPP interaction index ({_system}) for {len(npp_slice_results)} slices, "
        f"radius={_radius}mm"
    )
    print(f"Color range: [{npp_color_min:.3f}, {npp_color_max:.3f}]")
    return npp_color_max, npp_color_min, npp_slice_images, npp_slice_results


@app.cell(hide_code=True)
def _(
    boundaries,
    centroids,
    math,
    np,
    npp_color_max,
    npp_color_min,
    npp_radius_slider,
    npp_slice_images,
    npp_slice_results,
    npp_system_dropdown,
):
    import matplotlib as _matplotlib
    import matplotlib.colors as _mcolors
    import matplotlib.pyplot as _plt

    _matplotlib.rcParams["figure.dpi"] = 150

    _system_name = npp_system_dropdown.value
    _radius = npp_radius_slider.value
    _sorted_slices = sorted(npp_slice_results.keys())
    _n = len(_sorted_slices)
    _use_images = len(npp_slice_images) > 0

    # Global extent across all slices for locked axes
    _all_x = np.concatenate([npp_slice_results[z][0] for z in _sorted_slices])
    _all_y = np.concatenate([npp_slice_results[z][1] for z in _sorted_slices])
    _pad = 0.15
    _gx0 = float(_all_x.min()) - _pad
    _gx1 = float(_all_x.max()) + _pad
    _gy0 = float(_all_y.min()) - _pad
    _gy1 = float(_all_y.max()) + _pad

    _ncols = 4
    _nrows = math.ceil(_n / _ncols)
    _data_aspect = (_gy1 - _gy0) / (_gx1 - _gx0)
    _subplot_w = 4
    _subplot_h = _subplot_w * _data_aspect

    _fig, _axes = _plt.subplots(
        _nrows,
        _ncols,
        figsize=(_ncols * _subplot_w, _nrows * _subplot_h + 0.8),
        squeeze=False,
        sharex=True,
        sharey=True,
    )

    _norm = _mcolors.Normalize(vmin=0.0, vmax=0.1)
    _cmap = _plt.cm.inferno.copy()
    _cmap.set_bad(color="#f0f0f0")
    _bnd_keys = sorted(boundaries.keys())
    _cen_keys = sorted(centroids.keys())

    for _i, _z in enumerate(_sorted_slices):
        _ax = _axes[_i // _ncols][_i % _ncols]

        if _use_images and _z in npp_slice_images:
            _img, (_xmin_s, _xmax_s, _ymin_s, _ymax_s) = npp_slice_images[_z]
            _nr, _nc = _img.shape
            _half_dx = (_xmax_s - _xmin_s) / max(1, _nc - 1) / 2
            _half_dy = (_ymax_s - _ymin_s) / max(1, _nr - 1) / 2
            _ax.imshow(
                _img,
                cmap=_cmap,
                norm=_norm,
                extent=[
                    _xmin_s - _half_dx,
                    _xmax_s + _half_dx,
                    _ymax_s + _half_dy,
                    _ymin_s - _half_dy,
                ],
                aspect="equal",
                interpolation="nearest",
                origin="upper",
            )
        else:
            _sx, _sy, _sv = npp_slice_results[_z]
            _ax.scatter(
                _sx,
                _sy,
                c=_sv,
                s=4,
                cmap="inferno",
                vmin=npp_color_min,
                vmax=npp_color_max,
                edgecolors="none",
                rasterized=True,
            )

        # Region boundaries (nearest-key match to handle float mismatch)
        _z_bnd = min(_bnd_keys, key=lambda k: abs(float(k) - float(_z))) if _bnd_keys else None
        if _z_bnd is not None and abs(float(_z_bnd) - float(_z)) < 0.01:
            for _region, _coords in boundaries[_z_bnd].items():
                _xs = [c[0] for c in _coords] + [_coords[0][0]]
                _ys = [c[1] for c in _coords] + [_coords[0][1]]
                _ax.plot(_xs, _ys, color="gray", linewidth=0.5, alpha=0.6)

        # Region labels
        _z_cen = min(_cen_keys, key=lambda k: abs(float(k) - float(_z))) if _cen_keys else None
        if _z_cen is not None and abs(float(_z_cen) - float(_z)) < 0.01:
            for _region, (_cx, _cy) in centroids[_z_cen].items():
                _ax.text(
                    _cx,
                    _cy,
                    _region,
                    fontsize=5,
                    color="white",
                    ha="center",
                    va="center",
                    bbox=dict(
                        boxstyle="round,pad=0.15", fc="black", alpha=0.5, lw=0
                    ),
                )

        _ax.set_title(f"y={_z:.2f}", fontsize=9)
        _ax.tick_params(labelsize=6)
        for _spine in _ax.spines.values():
            _spine.set_visible(False)

    _axes[0][0].set_xlim(_gx0, _gx1)
    _axes[0][0].set_ylim(_gy1, _gy0)
    _axes[0][0].set_aspect("equal")

    for _i in range(_n, _nrows * _ncols):
        _axes[_i // _ncols][_i % _ncols].set_visible(False)

    _fig.suptitle(
        f"NPP interaction index: {_system_name}, r={_radius}mm",
        fontsize=12,
    )
    _fig.tight_layout(rect=[0, 0, 1.0, 0.96])

    _cbar_ax = _fig.add_axes([0.77, 0.01, 0.2, 0.015])
    _fig.colorbar(
        _plt.cm.ScalarMappable(norm=_norm, cmap=_cmap),
        cax=_cbar_ax,
        orientation="horizontal",
        label="4·p_s·p_r",
    )
    _cbar_ax.tick_params(labelsize=6)

    _plt.gca()
    return


@app.cell
def _(cKDTree, go, np, npp_slice_results, npp_system_dropdown, slice_results):
    # Collect (shannon_val, npp_val) across all common slices
    _system_name = npp_system_dropdown.value
    _common_slices = sorted(set(slice_results.keys()) & set(npp_slice_results.keys()))
    _all_shannon = []
    _all_npp = []

    for _z in _common_slices:
        _sx, _sy, _sv = slice_results[_z]
        _nx, _ny, _nv = npp_slice_results[_z]

        # Match grid points: find nearest NPP point for each Shannon point
        _tree = cKDTree(np.column_stack([_nx, _ny]))
        _dist, _idx = _tree.query(np.column_stack([_sx, _sy]))
        # Keep only close matches (within 0.05mm)
        _close = _dist < 0.05
        _all_shannon.append(_sv[_close])
        _all_npp.append(_nv[_idx[_close]])

    if _all_shannon:
        _s = np.concatenate(_all_shannon)
        _n = np.concatenate(_all_npp)
        _valid = np.isfinite(_s) & np.isfinite(_n)
        _s, _n = _s[_valid], _n[_valid]
    else:
        _s, _n = np.array([]), np.array([])

    go.Figure(
        {
            "data": [
                {
                    "type": "scattergl",
                    "x": _s.tolist(),
                    "y": _n.tolist(),
                    "mode": "markers",
                    "marker": {"size": 3, "opacity": 0.3, "color": "#2196F3"},
                    "showlegend": False,
                    "hovertemplate": "Shannon: %{x:.3f}<br>Interaction: %{y:.3f}<extra></extra>",
                }
            ],
            "layout": {
                "title": {"text": f"Cluster diversity vs. {_system_name} interaction index"},
                "xaxis": {"title": "Shannon diversity (cluster)"},
                "yaxis": {"title": f"{_system_name} interaction index"},
                "height": 500,
                "width": 600,
            },
        }
    )
    return


if __name__ == "__main__":
    app.run()
