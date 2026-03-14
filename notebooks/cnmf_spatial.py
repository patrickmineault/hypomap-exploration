# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "hypomap",
#     "marimo>=0.19.11",
#     "msgpack",
#     "numpy",
#     "pandas",
#     "plotly",
#     "pyarrow",
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
    from plotly.subplots import make_subplots
    return Path, go, make_subplots, math, mo, msgpack, np, pd


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # cNMF Gene Programs — Spatial Atlas

    Visualize cNMF gene program activity projected onto MERFISH spatial transcriptomics data.

    **Bridge**: cNMF was run on 10X scRNA-seq cells. Mean program usage per ABC taxonomy cluster
    is computed from 10X, then mapped onto MERFISH cells sharing the same cluster identity.
    """)
    return


@app.cell
def _(mo):
    decomposition_selector = mo.ui.dropdown(
        options={
            "All cells (162k)": "hypo_cnmf",
            "Neurons only (96k)": "hypo_cnmf_neurons",
            "ARH/ME/VMH MERFISH (7k)": "nmf_arh_me_vmh",
        },
        value="ARH/ME/VMH MERFISH (7k)",
        label="Decomposition",
    )
    decomposition_selector
    return (decomposition_selector,)


@app.cell
def _(Path, decomposition_selector, mo):
    import re
    _run_name = decomposition_selector.value
    _cnmf_dir = Path(mo.notebook_dir()) / ".." / "data" / "processed" / "mouse_abc" / "cnmf" / _run_name
    _available_ks = sorted({
        int(m.group(1))
        for f in _cnmf_dir.glob(f"{_run_name}.usages.k_*.dt_0_1.consensus.txt")
        if (m := re.search(r"\.k_(\d+)\.", f.name))
    })
    k_selector = mo.ui.dropdown(
        options={str(k): k for k in _available_ks},
        value=str(_available_ks[-1]) if _available_ks else None,
        label="K (number of programs)",
    )
    k_selector
    return (k_selector,)


@app.cell
def _(Path, decomposition_selector, k_selector, mo, pd):
    _root = Path(mo.notebook_dir()) / ".."
    _data = _root / "data" / "processed" / "mouse_abc"
    _run_name = decomposition_selector.value
    _cnmf_dir = _data / "cnmf" / _run_name
    _k = k_selector.value

    # Load cNMF usages
    _usages = pd.read_csv(
        _cnmf_dir / f"{_run_name}.usages.k_{_k}.dt_0_1.consensus.txt",
        sep="\t",
        index_col=0,
    )
    _usages.columns = [f"P{c}" for c in _usages.columns]
    _usages_norm = _usages.div(_usages.sum(axis=1), axis=0)
    _programs = [c for c in _usages_norm.columns if c.startswith("P")]

    # Detect MERFISH vs 10X
    _coords_path = _root / "data" / "processed" / "mouse_abc_extended" / "cells_with_coords.parquet"
    _sample_ids = set(_usages_norm.index[:100])
    _coord_ids = set(pd.read_parquet(_coords_path, columns=["cell_id"])["cell_id"]) if _coords_path.exists() else set()
    is_merfish = len(_sample_ids & _coord_ids) > len(_sample_ids) * 0.5

    if is_merfish:
        # MERFISH run: cell IDs in usages ARE MERFISH cells — use cluster from coords
        _coords = pd.read_parquet(_coords_path, columns=["cell_id", "cluster"]).set_index("cell_id")
        _merged = _usages_norm[_programs].join(_coords[["cluster"]], how="inner")
        cluster_program_usage = _merged.groupby("cluster")[_programs].mean()
        _status = f"MERFISH direct: {_merged.shape[0]:,} cells, {len(cluster_program_usage)} clusters"
    else:
        # 10X run: bridge via ABC taxonomy cluster
        from abc_atlas_access.abc_atlas_cache.abc_project_cache import AbcProjectCache

        _cache = AbcProjectCache.from_cache_dir(_root / "data" / "raw" / "abc_atlas_cache")
        _cache.load_latest_manifest()
        _cell_meta = _cache.get_metadata_dataframe(
            directory="WMB-10X",
            file_name="cell_metadata",
            dtype={"cell_label": str},
        ).set_index("cell_label")
        _membership = _cache.get_metadata_dataframe(
            directory="WMB-taxonomy",
            file_name="cluster_to_cluster_annotation_membership",
        )
        _taxonomy = _membership.pivot_table(
            index="cluster_alias",
            columns="cluster_annotation_term_set_name",
            values="cluster_annotation_term_name",
            aggfunc="first",
        )
        _cell_meta = _cell_meta.join(_taxonomy, on="cluster_alias")
        _cell_meta = _cell_meta[_cell_meta["feature_matrix_label"] == "WMB-10Xv3-HY"]
        _merged = _usages_norm[_programs].join(_cell_meta[["cluster"]], how="inner")
        cluster_program_usage = _merged.groupby("cluster")[_programs].mean()
        _status = f"10X bridge: {_merged.shape[0]:,} cells, {len(cluster_program_usage)} clusters"

    # Load program annotations (K-specific first, then fallback)
    _annot_path = _data / "cnmf" / f"program_annotations_{_run_name}_k{_k}.csv"
    if not _annot_path.exists():
        _annot_path = _data / "cnmf" / f"program_annotations_{_run_name}.csv"
    if not _annot_path.exists() and _run_name == "hypo_cnmf":
        _annot_path = _data / "cnmf" / "program_annotations.csv"

    program_labels = {p: p for p in _programs}
    if _annot_path.exists():
        program_annotations = pd.read_csv(_annot_path).set_index("program")
        for _p in _programs:
            if _p in program_annotations.index:
                program_labels[_p] = f"{_p} — {program_annotations.loc[_p, 'name']}"
    else:
        program_annotations = pd.DataFrame()

    mo.md(f"**{_status}**")
    return cluster_program_usage, program_annotations, program_labels


@app.cell
def _(Path, mo, msgpack, np, pd):
    _root = Path(mo.notebook_dir()) / ".."
    _data = _root / "data" / "processed" / "mouse_abc"

    # Load MERFISH spatial cells
    merfish = pd.read_parquet(_data / "cells_with_coords.parquet")

    # Load region boundaries
    with open(_data / "coronal_atlas_regions.msgpack", "rb") as f:
        _atlas = msgpack.unpack(f, raw=False)
    atlas_slices = _atlas["slices"]
    region_boundaries = {}
    for z_str, regions in _atlas["boundaries"].items():
        z_val = float(z_str)
        region_boundaries[z_val] = {
            r: np.frombuffer(b, dtype=np.float32).reshape(-1, 2).tolist()
            for r, b in regions.items()
        }
    region_centroids = {}
    for z_str, regions in _atlas["centroids"].items():
        z_val = float(z_str)
        region_centroids[z_val] = {
            r: np.frombuffer(b, dtype=np.float32).tolist()
            for r, b in regions.items()
        }

    mo.md(f"**MERFISH:** {len(merfish):,} cells, {len(atlas_slices)} slices")
    return atlas_slices, merfish, region_boundaries, region_centroids


@app.cell
def _(cluster_program_usage, merfish, mo, np):
    # Join mean program usage onto MERFISH cells via cluster
    _programs = cluster_program_usage.columns.tolist()
    _mapped = merfish[["cluster"]].join(cluster_program_usage, on="cluster", how="left")
    _n_mapped = _mapped[_programs[0]].notna().sum()
    _n_total = len(_mapped)

    # Store as numpy arrays for fast rendering
    program_values = {}
    for p in _programs:
        program_values[p] = _mapped[p].fillna(0).values

    # Also compute dominant program per MERFISH cell
    _usage_matrix = np.column_stack([program_values[p] for p in _programs])
    dominant_program_idx = np.argmax(_usage_matrix, axis=1)
    dominant_program_names = np.array(_programs)[dominant_program_idx]

    mo.md(f"**Mapped:** {_n_mapped:,}/{_n_total:,} MERFISH cells have cNMF usage ({_n_total - _n_mapped:,} unmapped)")
    return (program_values,)


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ## Spatial program expression
    """)
    return


@app.cell
def _(mo, program_labels):
    _programs = sorted(program_labels.keys(), key=lambda x: int(x[1:]))
    _options = {program_labels[p]: p for p in _programs}
    spatial_program = mo.ui.dropdown(
        options=_options,
        value=program_labels[_programs[0]],
        label="Program",
    )
    spatial_program
    return (spatial_program,)


@app.cell(hide_code=True)
def _(mo, pd, program_annotations, spatial_program):
    _prog = spatial_program.value
    _text = ""
    if len(program_annotations) and _prog in program_annotations.index:
        _a = program_annotations.loc[_prog]
        _certainty_icon = {"high": "●", "medium": "◐", "low": "○"}.get(str(_a.get("certainty", "")), "")
        _parts = [
            f"**{_a.get('name', _prog)}** — _{_a.get('category', '')}_  {_certainty_icon} {_a.get('certainty', '')}",
            f"> {_a.get('description', '')}",
        ]
        if "detailed_information" in _a and pd.notna(_a["detailed_information"]):
            _parts.append("")
            _parts.append(f"<details><summary>Detailed information</summary>\n\n{_a['detailed_information']}\n\n</details>")
        _text = "\n\n".join(_parts)
    mo.md(_text) if _text else None
    return


@app.cell(hide_code=True)
def _(
    atlas_slices,
    go,
    make_subplots,
    math,
    merfish,
    np,
    program_labels,
    program_values,
    region_boundaries,
    region_centroids,
    spatial_program,
):
    _prog = spatial_program.value
    _vals = program_values[_prog]
    _prog_label = program_labels.get(_prog, _prog)

    _n_slices = len(atlas_slices)
    _n_cols = 3
    _n_rows = math.ceil(_n_slices / _n_cols)

    _raw_x = merfish["x"].values
    _raw_y = merfish["y"].values
    _raw_z = merfish["z_slice"].values

    # Shared bin edges across all slices (fixed grid for consistent appearance)
    _bin_size = 0.1  # mm per bin
    _x_min, _x_max = float(_raw_x.min()), float(_raw_x.max())
    _y_min, _y_max = float(_raw_y.min()), float(_raw_y.max())
    _x_edges = np.arange(_x_min, _x_max + _bin_size, _bin_size)
    _y_edges = np.arange(_y_min, _y_max + _bin_size, _bin_size)
    _x_centers = (_x_edges[:-1] + _x_edges[1:]) / 2
    _y_centers = (_y_edges[:-1] + _y_edges[1:]) / 2

    # Color scale: clip at 98th percentile for contrast
    _vmax = float(np.percentile(_vals[_vals > 0], 98)) if (_vals > 0).any() else 0.1

    _fig = make_subplots(
        rows=_n_rows, cols=_n_cols,
        horizontal_spacing=0.02, vertical_spacing=0.01,
    )

    _all_traces = []
    _all_shapes = []
    _all_annotations = []

    for _si, _z in enumerate(atlas_slices):
        _row = _si // _n_cols + 1
        _col = _si % _n_cols + 1
        _xaxis = f"x{_si + 1}" if _si > 0 else "x"
        _yaxis = f"y{_si + 1}" if _si > 0 else "y"

        # Region boundaries (white lines, visible over heatmap)
        if _z in region_boundaries:
            for _region, _hull in region_boundaries[_z].items():
                _path_parts = [f"M {_hull[0][0]},{_hull[0][1]}"]
                for _px, _py in _hull[1:]:
                    _path_parts.append(f"L {_px},{_py}")
                _path_parts.append("Z")
                _all_shapes.append({
                    "type": "path", "path": " ".join(_path_parts),
                    "fillcolor": "rgba(0,0,0,0)",
                    "line": {"color": "rgba(200,200,200,0.7)", "width": 1.5},
                    "layer": "above", "xref": _xaxis, "yref": _yaxis,
                })

        # Region labels at centroids
        if _z in region_centroids:
            for _region, (_cx, _cy) in region_centroids[_z].items():
                _all_annotations.append({
                    "x": _cx, "y": _cy,
                    "xref": _xaxis, "yref": _yaxis,
                    "text": _region,
                    "showarrow": False,
                    "font": {"size": 8, "color": "white", "family": "Arial Black"},
                    "bgcolor": "rgba(0,0,0,0.5)",
                    "opacity": 0.9,
                })

        # Bin cells for this slice
        _mask = _raw_z == _z
        if _mask.sum() == 0:
            continue
        _sx = _raw_x[_mask]
        _sy = _raw_y[_mask]
        _sv = _vals[_mask]

        # Weighted histogram: sum of usage values per bin
        _sum_vals, _, _ = np.histogram2d(_sx, _sy, bins=[_x_edges, _y_edges], weights=_sv)
        _counts, _, _ = np.histogram2d(_sx, _sy, bins=[_x_edges, _y_edges])
        # Mean usage per bin (avoid division by zero)
        with np.errstate(invalid="ignore"):
            _mean_vals = np.where(_counts > 0, _sum_vals / _counts, np.nan)

        # Transpose: histogram2d returns (x_bins, y_bins), heatmap expects (y, x)
        _z_data = _mean_vals.T

        _all_traces.append({
            "type": "heatmap",
            "z": np.where(np.isnan(_z_data), None, np.round(_z_data, 4)).tolist(),
            "x": np.round(_x_centers, 3).tolist(),
            "y": np.round(_y_centers, 3).tolist(),
            "colorscale": [[0, "#000000"], [0.33, "#800000"], [0.67, "#FF4500"], [1, "#FFFF00"]],
            "zmin": 0,
            "zmax": _vmax,
            "showscale": _si == 0,
            "colorbar": {"title": "Usage", "len": 0.3, "y": 0.85} if _si == 0 else None,
            "hovertemplate": f"{_prog}: " + "%{z:.3f}<extra></extra>",
            "xaxis": _xaxis,
            "yaxis": _yaxis,
        })

        _all_annotations.append({
            "xref": f"{_xaxis} domain", "yref": f"{_yaxis} domain",
            "x": 0.98, "y": 0.02, "text": f"Z={_z}",
            "showarrow": False, "font": {"size": 9, "color": "#666"},
            "xanchor": "right", "yanchor": "bottom",
        })

    # Equalize axes
    _span = max(_x_max - _x_min, _y_max - _y_min) * 1.05
    _xc, _yc = (_x_min + _x_max) / 2, (_y_min + _y_max) / 2
    _half = _span / 2
    _ax_x = {"showticklabels": False, "showgrid": False, "zeroline": False, "matches": "x", "range": [_xc - _half, _xc + _half]}
    _ax_y = {"showticklabels": False, "showgrid": False, "zeroline": False, "matches": "y", "range": [_yc + _half, _yc - _half]}

    _fig_dict = _fig.to_dict()
    _layout = _fig_dict["layout"]
    for _i in range(1, _n_rows * _n_cols + 1):
        _layout.setdefault(f"xaxis{_i}" if _i > 1 else "xaxis", {}).update(_ax_x)
        _layout.setdefault(f"yaxis{_i}" if _i > 1 else "yaxis", {}).update(_ax_y)

    _fig_height = max(400, _n_rows * 350)
    _layout.update(
        height=min(_fig_height, 16000),
        showlegend=False,
        margin={"l": 10, "r": 10, "t": 30, "b": 10},
        title={"text": _prog_label, "font": {"size": 14}},
        paper_bgcolor="white", plot_bgcolor="white",
        annotations=_all_annotations, shapes=_all_shapes,
    )
    _fig_dict["data"] = _all_traces
    go.Figure(_fig_dict)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ## Program usage by subclass (heatmap)
    """)
    return


@app.cell
def _(cluster_program_usage, go, merfish, program_labels):
    # Compute mean program usage per subclass (via MERFISH cells for spatial relevance)
    _programs = sorted(program_labels.keys(), key=lambda x: int(x[1:]))
    _prog_display = [program_labels[p] for p in _programs]

    # Map cluster-level usage to subclass by averaging across clusters within each subclass
    _cluster_sub = merfish[["cluster", "subclass"]].drop_duplicates()
    _usage_with_sub = cluster_program_usage[_programs].join(
        _cluster_sub.set_index("cluster"), how="inner"
    )
    _mean_by_sub = _usage_with_sub.groupby("subclass")[_programs].mean()

    # Z-score across subclasses
    _zscore = (_mean_by_sub - _mean_by_sub.mean()) / _mean_by_sub.std()

    # Sort subclasses by peak program
    _peak = _zscore.idxmax(axis=1)
    _sub_order = _peak.sort_values(key=lambda s: s.map(lambda x: int(x[1:]))).index.tolist()

    _z = _zscore.loc[_sub_order, _programs].values
    _fig = go.Figure(data=go.Heatmap(
        z=_z, x=_prog_display, y=_sub_order,
        colorscale="RdBu_r", zmin=-3, zmax=3,
        colorbar={"title": "Z-score"},
        hovertemplate="%{y}<br>%{x}<br>z=%{z:.2f}<extra></extra>",
    ))
    _fig.update_layout(
        title="Program specificity by subclass (z-scored mean usage)",
        xaxis_title="Gene Program",
        yaxis_title="Subclass",
        height=max(len(_sub_order) * 16, 400),
        yaxis={"autorange": "reversed"},
    )
    _fig
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ## Program usage by hypothalamic region
    """)
    return


@app.cell
def _(cluster_program_usage, go, merfish, program_labels):
    # Mean program usage per anatomical region
    _programs = sorted(program_labels.keys(), key=lambda x: int(x[1:]))
    _prog_display = [program_labels[p] for p in _programs]

    # Join usage to MERFISH cells, then average by region
    _cells = merfish[["cluster", "region"]].copy()
    _usage_map = cluster_program_usage[_programs]
    _cells = _cells.join(_usage_map, on="cluster", how="inner")
    _mean_by_region = _cells.groupby("region")[_programs].mean()

    # Filter regions with at least 100 cells
    _region_counts = merfish["region"].value_counts()
    _keep = _region_counts[_region_counts >= 100].index
    _mean_by_region = _mean_by_region.loc[_mean_by_region.index.isin(_keep)]

    # Z-score
    _zscore = (_mean_by_region - _mean_by_region.mean()) / _mean_by_region.std()

    # Sort by peak program
    _peak = _zscore.idxmax(axis=1)
    _region_order = _peak.sort_values(key=lambda s: s.map(lambda x: int(x[1:]))).index.tolist()

    _z = _zscore.loc[_region_order, _programs].values
    _fig = go.Figure(data=go.Heatmap(
        z=_z, x=_prog_display, y=_region_order,
        colorscale="RdBu_r", zmin=-3, zmax=3,
        colorbar={"title": "Z-score"},
        hovertemplate="%{y}<br>%{x}<br>z=%{z:.2f}<extra></extra>",
    ))
    _fig.update_layout(
        title="Program specificity by hypothalamic region (z-scored)",
        xaxis_title="Gene Program",
        yaxis_title="Region",
        height=max(len(_region_order) * 18, 400),
        yaxis={"autorange": "reversed"},
    )
    _fig
    return


if __name__ == "__main__":
    app.run()
