# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "hypomap",
#     "marimo>=0.19.11",
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
    from pathlib import Path

    import marimo as mo
    import numpy as np
    import pandas as pd
    import plotly.express as px
    import plotly.graph_objects as go
    from scipy.cluster.hierarchy import leaves_list, linkage, optimal_leaf_ordering
    from scipy.spatial.distance import squareform
    return (
        Path,
        go,
        leaves_list,
        linkage,
        mo,
        np,
        optimal_leaf_ordering,
        pd,
        px,
        squareform,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # cNMF Gene Programs — Mouse Hypothalamus

    Explore the 50 consensus gene programs discovered by cNMF on hypothalamus
    scRNA-seq cells (10Xv3, Allen ABC Atlas).

    - **Usages**: how much each cell uses each program (cells x 50)
    - **Gene spectra**: which genes define each program (50 x ~32k genes)
    - **Cell metadata**: cell-type annotations from the ABC taxonomy
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
    if not _available_ks:
        _available_ks = [50]

    k_selector = mo.ui.dropdown(
        options={str(k): k for k in _available_ks},
        value=str(_available_ks[-1]),
        label="K (number of programs)",
    )
    k_selector
    return (k_selector,)


@app.cell
def _(Path, decomposition_selector, k_selector, mo, pd):
    # Load cNMF results
    _run_name = decomposition_selector.value
    _k = k_selector.value
    _cnmf_dir = Path(mo.notebook_dir()) / ".." / "data" / "processed" / "mouse_abc" / "cnmf" / _run_name

    usages = pd.read_csv(
        _cnmf_dir / f"{_run_name}.usages.k_{_k}.dt_0_1.consensus.txt",
        sep="\t",
        index_col=0,
    )
    usages.columns = [f"P{c}" for c in usages.columns]

    gene_spectra_score = pd.read_csv(
        _cnmf_dir / f"{_run_name}.gene_spectra_score.k_{_k}.dt_0_1.txt",
        sep="\t",
        index_col=0,
    )
    gene_spectra_score.index = [f"P{i}" for i in gene_spectra_score.index]

    gene_spectra_tpm = pd.read_csv(
        _cnmf_dir / f"{_run_name}.gene_spectra_tpm.k_{_k}.dt_0_1.txt",
        sep="\t",
        index_col=0,
    )
    gene_spectra_tpm.index = [f"P{i}" for i in gene_spectra_tpm.index]

    # Load program annotations if available
    # Check K-specific path first, then legacy paths
    _cnmf_base = Path(mo.notebook_dir()) / ".." / "data" / "processed" / "mouse_abc" / "cnmf"
    _annot_path = _cnmf_base / f"program_annotations_{_run_name}_k{_k}.csv"
    if not _annot_path.exists():
        _annot_path = _cnmf_base / f"program_annotations_{_run_name}.csv"
    if not _annot_path.exists() and _run_name == "hypo_cnmf":
        _annot_path = _cnmf_base / "program_annotations.csv"

    if _annot_path.exists():
        _annot_df = pd.read_csv(_annot_path).set_index("program")
        # Only use annotations if programs align with current K
        if all(p in usages.columns for p in _annot_df.index):
            program_annotations = _annot_df
            program_labels = {
                p: f"{p} — {program_annotations.loc[p, 'name']}"
                for p in program_annotations.index
                if p in usages.columns
            }
        else:
            program_annotations = pd.DataFrame()
            program_labels = {p: p for p in usages.columns}
    else:
        program_annotations = pd.DataFrame()
        program_labels = {p: p for p in usages.columns}

    # Normalize usages to fractions per cell
    usages_norm = usages.div(usages.sum(axis=1), axis=0)

    # Dominant program per cell
    dominant_program = usages_norm.idxmax(axis=1).rename("dominant_program")

    mo.md(f"""
    **Loaded ({_run_name}, K={_k}):**
    - Usages: {usages.shape[0]:,} cells x {usages.shape[1]} programs
    - Gene spectra: {gene_spectra_score.shape[0]} programs x {gene_spectra_score.shape[1]:,} genes
    - Annotations: {"✓ " + str(len(program_annotations)) + " programs" if len(program_annotations) else "none (K mismatch or not generated)"}
    """)
    return (
        dominant_program,
        gene_spectra_score,
        program_annotations,
        program_labels,
        usages_norm,
    )


@app.cell
def _(Path, decomposition_selector, mo, pd, usages_norm):
    _run_name = decomposition_selector.value
    _data_dir = Path(mo.notebook_dir()) / ".." / "data"

    # Detect MERFISH vs 10X by checking cell ID overlap with cells_with_coords
    _coords_path = _data_dir / "processed" / "mouse_abc_extended" / "cells_with_coords.parquet"
    _sample_ids = set(usages_norm.index[:100])
    _is_merfish = False
    if _coords_path.exists():
        _coord_ids = set(pd.read_parquet(_coords_path, columns=["cell_id"])["cell_id"])
        _is_merfish = len(_sample_ids & _coord_ids) > len(_sample_ids) * 0.5

    if _is_merfish:
        # MERFISH: metadata comes directly from cells_with_coords
        cell_meta = pd.read_parquet(
            _coords_path,
            columns=["cell_id", "class", "subclass", "supertype", "cluster", "region",
                     "class_color", "subclass_color", "supertype_color", "cluster_color"],
        )
        cell_meta = cell_meta[~cell_meta["region"].str.contains("unassigned", case=False)]
        cell_meta = cell_meta.set_index("cell_id")
    else:
        # 10X: load from ABC Atlas cache
        from abc_atlas_access.abc_atlas_cache.abc_project_cache import AbcProjectCache

        _cache_dir = _data_dir / "raw" / "abc_atlas_cache"
        _cache = AbcProjectCache.from_cache_dir(_cache_dir)
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

        _colors = _membership.pivot_table(
            index="cluster_alias",
            columns="cluster_annotation_term_set_name",
            values="color_hex_triplet",
            aggfunc="first",
        )
        _colors.columns = [f"{c}_color" for c in _colors.columns]

        cell_meta = _cell_meta.join(_taxonomy, on="cluster_alias").join(
            _colors, on="cluster_alias"
        )
        cell_meta = cell_meta[cell_meta["feature_matrix_label"] == "WMB-10Xv3-HY"]

    mo.md(f"""
    **Cell metadata** ({'MERFISH' if _is_merfish else '10X'}): {len(cell_meta):,} cells

    | Level | Unique |
    |-------|--------|
    | class | {cell_meta['class'].nunique()} |
    | subclass | {cell_meta['subclass'].nunique()} |
    | supertype | {cell_meta['supertype'].nunique()} |
    | cluster | {cell_meta['cluster'].nunique()} |
    """)
    return (cell_meta,)


@app.cell
def _(cell_meta, dominant_program, mo, pd, usages_norm):
    # Merge usages with cell metadata
    cells = pd.concat(
        [usages_norm, dominant_program, cell_meta.reindex(usages_norm.index)],
        axis=1,
    )

    # Drop cells without annotations (should be very few)
    _n_before = len(cells)
    cells = cells.dropna(subset=["class"])
    _n_after = len(cells)

    mo.md(f"**Merged:** {_n_after:,} cells with both usages and annotations ({_n_before - _n_after} dropped)")
    return (cells,)


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ## Program-Cell Type Associations
    """)
    return


@app.cell
def _(mo):
    annotation_level = mo.ui.dropdown(
        options=["class", "subclass", "supertype", "cluster"],
        value="subclass",
        label="Annotation level",
    )
    annotation_level
    return (annotation_level,)


@app.cell
def _(annotation_level, cells, go, program_labels):
    # Mean usage of each program per cell type
    _level = annotation_level.value
    _programs = [c for c in cells.columns if c.startswith("P")]

    _mean_usage = cells.groupby(_level)[_programs].mean()

    # Z-score across cell types for each program (highlights specificity)
    _zscore = (_mean_usage - _mean_usage.mean()) / _mean_usage.std()

    # Sort programs and cell types
    _prog_order = sorted(_programs, key=lambda x: int(x[1:]))
    _prog_display = [program_labels.get(p, p) for p in _prog_order]
    _peak_prog = _zscore.idxmax(axis=1)
    _type_order = _peak_prog.sort_values(key=lambda s: s.map(lambda x: int(x[1:]))).index.tolist()

    _z = _zscore.loc[_type_order, _prog_order].values
    _fig = go.Figure(data=go.Heatmap(
        z=_z,
        x=_prog_display,
        y=_type_order,
        colorscale="RdBu_r",
        zmin=-3,
        zmax=3,
        colorbar=dict(title="Z-score"),
        hovertemplate="%{y}<br>%{x}<br>z=%{z:.2f}<extra></extra>",
    ))
    _fig.update_layout(
        title=f"Program specificity by {_level} (z-scored mean usage)",
        xaxis_title="Gene Program",
        yaxis_title=_level.capitalize(),
        height=max(len(_type_order) * 18, 400),
        yaxis=dict(autorange="reversed"),
    )
    _fig
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ## Top Genes per Program
    """)
    return


@app.cell
def _(mo, program_labels, usages_norm):
    _programs = sorted(usages_norm.columns, key=lambda x: int(x[1:]))
    _options = {program_labels.get(p, p): p for p in _programs}
    program_selector = mo.ui.dropdown(
        options=_options,
        value=program_labels.get("P1", "P1"),
        label="Select program",
    )
    n_top_genes = mo.ui.slider(
        start=10, stop=50, value=20, step=5, label="Top N genes"
    )
    mo.hstack([program_selector, n_top_genes])
    return n_top_genes, program_selector


@app.cell
def _(mo, pd, program_annotations, program_selector):
    _prog = program_selector.value
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


@app.cell
def _(gene_spectra_score, go, n_top_genes, program_labels, program_selector):
    _prog = program_selector.value
    _n = n_top_genes.value
    _prog_name = program_labels.get(_prog, _prog)

    # Top genes by z-score
    _scores = gene_spectra_score.loc[_prog].sort_values(ascending=False)
    _top = _scores.head(_n)

    _fig = go.Figure(data=go.Bar(
        x=_top.values[::-1],
        y=_top.index[::-1],
        orientation="h",
        marker_color=_top.values[::-1],
        marker_colorscale="Reds",
        hovertemplate="%{y}: %{x:.2f}<extra></extra>",
    ))
    _fig.update_layout(
        title=f"Top {_n} genes for {_prog_name}",
        xaxis_title="Gene Spectra Score (z)",
        yaxis_title="Gene",
        height=max(_n * 25, 300),
    )
    _fig
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ## Program Usage Distribution by Cell Type
    """)
    return


@app.cell
def _(annotation_level, cells, program_labels, program_selector, px):
    # Distribution of selected program's usage across cell types
    _prog = program_selector.value
    _level = annotation_level.value
    _prog_name = program_labels.get(_prog, _prog)

    # Get top cell types by mean usage of this program
    _mean = cells.groupby(_level)[_prog].mean().sort_values(ascending=False)
    _top_types = _mean.head(15).index.tolist()

    _plot_df = cells[cells[_level].isin(_top_types)][[_prog, _level]].copy()
    _plot_df = _plot_df.rename(columns={_prog: "usage"})

    _fig = px.box(
        _plot_df,
        x="usage",
        y=_level,
        category_orders={_level: _top_types},
        title=f"{_prog_name} usage across top {_level}s",
        labels={"usage": f"{_prog} usage (fraction)", _level: _level.capitalize()},
    )
    _fig.update_layout(
        height=max(len(_top_types) * 35, 300),
        showlegend=False,
    )
    _fig
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ## Program Correlations
    """)
    return


@app.cell
def _(go, program_labels, usages_norm, v_order_programs):
    # U matrix: program-program correlation from cell usages, ordered by V matrix clustering
    _corr_u = usages_norm.corr()

    _prog_display = [program_labels.get(p, p) for p in v_order_programs]
    _z_u = _corr_u.loc[v_order_programs, v_order_programs].values

    _fig_u = go.Figure(data=go.Heatmap(
        z=_z_u,
        x=_prog_display,
        y=_prog_display,
        colorscale="RdBu_r",
        zmin=-1,
        zmax=1,
        colorbar=dict(title="Correlation"),
        hovertemplate="%{x} vs %{y}<br>r=%{z:.3f}<extra></extra>",
    ))
    _fig_u.update_layout(
        title="U matrix — program correlation by cell usages (V-matrix order)",
        width=800,
        height=800,
    )
    _fig_u
    return


@app.cell
def _(
    gene_spectra_score,
    go,
    leaves_list,
    linkage,
    np,
    optimal_leaf_ordering,
    program_labels,
    squareform,
):
    # V matrix: program-program correlation from gene spectra (rows of V)
    _corr_v = gene_spectra_score.T.corr()

    # Hierarchical clustering with optimal leaf ordering
    _dist_v = 1 - _corr_v.values
    np.fill_diagonal(_dist_v, 0)
    _dist_v = (_dist_v + _dist_v.T) / 2
    _linkage_v = linkage(squareform(_dist_v, checks=False), method="average")
    _linkage_v = optimal_leaf_ordering(_linkage_v, squareform(_dist_v, checks=False))
    _order_v = leaves_list(_linkage_v)

    # Export the V-matrix ordering for use by U matrix and description similarity
    v_order_programs = [_corr_v.columns[i] for i in _order_v]
    _prog_display = [program_labels.get(p, p) for p in v_order_programs]
    _z_v = _corr_v.loc[v_order_programs, v_order_programs].values

    _fig_v = go.Figure(data=go.Heatmap(
        z=_z_v,
        x=_prog_display,
        y=_prog_display,
        colorscale="RdBu_r",
        zmin=-1,
        zmax=1,
        colorbar=dict(title="Correlation"),
        hovertemplate="%{x} vs %{y}<br>r=%{z:.3f}<extra></extra>",
    ))
    _fig_v.update_layout(
        title="V matrix — program correlation by gene spectra (clustered)",
        width=800,
        height=800,
    )
    _fig_v
    return (v_order_programs,)


@app.cell
def _(Path, decomposition_selector, go, k_selector, mo, np, v_order_programs):
    # Load precomputed sentence-BERT embeddings for the current K
    # Generated by: uv run python scripts/embed_annotations.py --run-name <run_name>
    _run_name = decomposition_selector.value
    _k = k_selector.value
    _cnmf_base = Path(mo.notebook_dir()) / ".." / "data" / "processed" / "mouse_abc" / "cnmf"
    _emb_path = _cnmf_base / f"annotation_embeddings_{_run_name}.npz"

    _output = mo.md(f"**No precomputed embeddings found.** Run `uv run python scripts/embed_annotations.py --run-name {_run_name}`")

    if _emb_path.exists():
        _data = np.load(_emb_path, allow_pickle=True)
        _key = f"embeddings_k{_k}"
        if _key not in _data:
            _output = mo.md(f"**No embeddings for K={_k}.** Re-run `scripts/embed_annotations.py` to include this K.")
        else:
            _embeddings = _data[_key]
            _labels = list(_data[f"labels_k{_k}"])

            # Cosine similarity (embeddings are already L2-normalized)
            _sim = _embeddings @ _embeddings.T

            # Reorder to match V-matrix clustering
            # Labels are like "P0: name", extract program id to map to v_order
            _label_by_prog = {lbl.split(":")[0]: i for i, lbl in enumerate(_labels)}
            _order = [_label_by_prog[p] for p in v_order_programs if p in _label_by_prog]

            _labels_ordered = [_labels[i] for i in _order]
            _sim_ordered = _sim[np.ix_(_order, _order)]

            _output = go.Figure(data=go.Heatmap(
                z=_sim_ordered,
                x=_labels_ordered,
                y=_labels_ordered,
                colorscale="RdBu_r",
                zmin=-0.5,
                zmax=1,
                colorbar=dict(title="Cosine sim"),
                hovertemplate="%{y}<br>vs<br>%{x}<br>sim=%{z:.3f}<extra></extra>",
            ))
            _output.update_layout(
                title=f"Description similarity (sentence-BERT) — K={_k} (V-matrix order)",
                width=800,
                height=800,
            )
    _output
    return


@app.cell
def _(
    Path,
    decomposition_selector,
    go,
    mo,
    np,
    pd,
    program_labels,
    usages_norm,
    v_order_programs,
):
    # Spatial footprint similarity: correlate programs by their binned spatial expression
    _run_name = decomposition_selector.value
    _data_dir = Path(mo.notebook_dir()) / ".." / "data"
    _coords_path = _data_dir / "processed" / "mouse_abc_extended" / "cells_with_coords.parquet"

    # Check if this is a MERFISH run (has spatial coords)
    _sample_ids = set(usages_norm.index[:100])
    _is_merfish = False
    if _coords_path.exists():
        _coord_ids = set(pd.read_parquet(_coords_path, columns=["cell_id"])["cell_id"])
        _is_merfish = len(_sample_ids & _coord_ids) > len(_sample_ids) * 0.5

    if not _is_merfish:
        _output = mo.md("*Spatial footprint similarity requires MERFISH data.*")
    else:
        _coords = pd.read_parquet(_coords_path, columns=["cell_id", "x", "y", "z_slice"]).set_index("cell_id")
        _merged = usages_norm.join(_coords, how="inner")
        _programs = [p for p in v_order_programs if p in usages_norm.columns]

        # Build spatial bins
        _bin_size = 0.1
        _x_edges = np.arange(_merged["x"].min(), _merged["x"].max() + _bin_size, _bin_size)
        _y_edges = np.arange(_merged["y"].min(), _merged["y"].max() + _bin_size, _bin_size)
        _slices = sorted(_merged["z_slice"].unique())

        # Compute spatial footprint vector per program: concatenated mean-usage histograms across slices
        _footprints = {}
        for _p in _programs:
            _bins_all = []
            for _z in _slices:
                _mask = _merged["z_slice"] == _z
                _sx = _merged.loc[_mask, "x"].values
                _sy = _merged.loc[_mask, "y"].values
                _sv = _merged.loc[_mask, _p].values
                _sum_vals, _, _ = np.histogram2d(_sx, _sy, bins=[_x_edges, _y_edges], weights=_sv)
                _counts, _, _ = np.histogram2d(_sx, _sy, bins=[_x_edges, _y_edges])
                with np.errstate(invalid="ignore"):
                    _mean_vals = np.where(_counts > 0, _sum_vals / _counts, 0)
                _bins_all.append(_mean_vals.ravel())
            _footprints[_p] = np.concatenate(_bins_all)

        # Correlation matrix of spatial footprints
        _fp_matrix = np.column_stack([_footprints[p] for p in _programs])
        _corr = np.corrcoef(_fp_matrix.T)

        _prog_display = [program_labels.get(p, p) for p in _programs]
        _output = go.Figure(data=go.Heatmap(
            z=_corr,
            x=_prog_display,
            y=_prog_display,
            colorscale="RdBu_r",
            zmin=-1,
            zmax=1,
            colorbar=dict(title="Correlation"),
            hovertemplate="%{x} vs %{y}<br>r=%{z:.3f}<extra></extra>",
        ))
        _output.update_layout(
            title="Spatial footprint correlation (V-matrix order)",
            width=800,
            height=800,
        )
    _output
    return


@app.cell
def _(mo):
    merge_threshold = mo.ui.slider(
        start=0.3, stop=0.95, step=0.05, value=0.7,
        label="Merge threshold (min similarity across metrics)",
    )
    merge_threshold
    return (merge_threshold,)


@app.cell(hide_code=True)
def _(
    Path,
    decomposition_selector,
    gene_spectra_score,
    k_selector,
    merge_threshold,
    mo,
    np,
    pd,
    program_labels,
    usages_norm,
    v_order_programs,
):
    _threshold = merge_threshold.value
    _programs = [p for p in v_order_programs if p in usages_norm.columns]
    _k = k_selector.value
    _n = len(_programs)

    # 1. U matrix correlation
    _corr_u = usages_norm[_programs].corr().values

    # 2. V matrix correlation
    _corr_v = gene_spectra_score.T[_programs].corr().values

    # 3. Description similarity (sentence-BERT)
    _run_name = decomposition_selector.value
    _cnmf_base = Path(mo.notebook_dir()) / ".." / "data" / "processed" / "mouse_abc" / "cnmf"
    _emb_path = _cnmf_base / f"annotation_embeddings_{_run_name}.npz"
    _sim_desc = np.full((_n, _n), np.nan)
    if _emb_path.exists():
        _data = np.load(_emb_path, allow_pickle=True)
        _key = f"embeddings_k{_k}"
        if _key in _data:
            _embeddings = _data[_key]
            _labels = list(_data[f"labels_k{_k}"])
            _label_by_prog = {lbl.split(":")[0]: i for i, lbl in enumerate(_labels)}
            _idx = [_label_by_prog.get(p) for p in _programs]
            if all(i is not None for i in _idx):
                _emb_ordered = _embeddings[_idx]
                _sim_desc = _emb_ordered @ _emb_ordered.T

    # 4. Spatial footprint correlation
    _data_dir = Path(mo.notebook_dir()) / ".." / "data"
    _coords_path = _data_dir / "processed" / "mouse_abc_extended" / "cells_with_coords.parquet"
    _sim_spatial = np.full((_n, _n), np.nan)
    _sample_ids = set(usages_norm.index[:100])
    _has_spatial = False
    if _coords_path.exists():
        _coord_ids = set(pd.read_parquet(_coords_path, columns=["cell_id"])["cell_id"])
        _has_spatial = len(_sample_ids & _coord_ids) > len(_sample_ids) * 0.5
    if _has_spatial:
        _coords = pd.read_parquet(_coords_path, columns=["cell_id", "x", "y", "z_slice"]).set_index("cell_id")
        _merged = usages_norm.join(_coords, how="inner")
        _bin_size = 0.1
        _x_edges = np.arange(_merged["x"].min(), _merged["x"].max() + _bin_size, _bin_size)
        _y_edges = np.arange(_merged["y"].min(), _merged["y"].max() + _bin_size, _bin_size)
        _slices = sorted(_merged["z_slice"].unique())
        _footprints = {}
        for _p in _programs:
            _bins_all = []
            for _z in _slices:
                _mask = _merged["z_slice"] == _z
                _sx = _merged.loc[_mask, "x"].values
                _sy = _merged.loc[_mask, "y"].values
                _sv = _merged.loc[_mask, _p].values
                _sum_vals, _, _ = np.histogram2d(_sx, _sy, bins=[_x_edges, _y_edges], weights=_sv)
                _counts, _, _ = np.histogram2d(_sx, _sy, bins=[_x_edges, _y_edges])
                with np.errstate(invalid="ignore"):
                    _mean_vals = np.where(_counts > 0, _sum_vals / _counts, 0)
                _bins_all.append(_mean_vals.ravel())
            _footprints[_p] = np.concatenate(_bins_all)
        _fp_matrix = np.column_stack([_footprints[p] for p in _programs])
        _sim_spatial = np.corrcoef(_fp_matrix.T)

    # Find merge candidates: pairs where ALL available metrics exceed threshold
    _candidates = []
    for _i in range(_n):
        for _j in range(_i + 1, _n):
            _scores = {
                "U corr": _corr_u[_i, _j],
                "V corr": _corr_v[_i, _j],
            }
            if not np.isnan(_sim_desc[_i, _j]):
                _scores["Description sim"] = _sim_desc[_i, _j]
            if not np.isnan(_sim_spatial[_i, _j]):
                _scores["Spatial corr"] = _sim_spatial[_i, _j]
            _mean = float(np.mean(list(_scores.values())))
            if _mean >= _threshold:
                _candidates.append({
                    "Program A": program_labels.get(_programs[_i], _programs[_i]),
                    "Program B": program_labels.get(_programs[_j], _programs[_j]),
                    **{k: round(v, 3) for k, v in _scores.items()},
                    "Mean": round(_mean, 3),
                })

    _candidates.sort(key=lambda x: -x["Mean"])

    if _candidates:
        _df = pd.DataFrame(_candidates)
        _output = mo.vstack([
            mo.md(f"### Suggested merges (mean similarity ≥ {_threshold})"),
            mo.md(f"**{len(_candidates)} pairs** with mean similarity ≥ {_threshold}:"),
            _df,
        ])
    else:
        _output = mo.md(f"### No merge candidates at mean similarity ≥ {_threshold}")
    _output
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ## Dominant Program per Cell Type
    """)
    return


@app.cell
def _(annotation_level, cells, pd):
    _level = annotation_level.value

    _ct = pd.crosstab(cells[_level], cells["dominant_program"], normalize="index")

    #_ct.sort_index()
    #_ct.sort_index(in_)
    return


@app.cell
def _(annotation_level, cells, go, pd, program_labels):
    # What fraction of each cell type is dominated by each program?
    _level = annotation_level.value

    _ct = pd.crosstab(cells[_level], cells["dominant_program"], normalize="index")
    _ct.sort_index(inplace=True)

    # Sort cell types by their most common dominant program
    _dominant = _ct.idxmax(axis=1)
    #_type_order = _dominant.sort_values(key=lambda s: s.map(lambda x: int(x[1:]))).index.tolist()
    _prog_order = sorted(_ct.columns, key=lambda x: int(x[1:]))
    _prog_display = [program_labels.get(p, p) for p in _prog_order]

    #_z = _ct.loc[:, _prog_order].values
    _z = _ct.loc[:, _prog_order].values
    _fig = go.Figure(data=go.Heatmap(
        z=_z,
        x=_prog_display,
        y=_ct.index,
        colorscale="Viridis",
        colorbar=dict(title="Fraction"),
        hovertemplate="%{y}<br>%{x}<br>fraction=%{z:.2%}<extra></extra>",
    ))
    _fig.update_layout(
        title=f"Fraction of cells in each {_level} dominated by each program",
        xaxis_title="Dominant Program",
        yaxis_title=_level.capitalize(),
        height=max(len(_ct.index) * 18, 400),
        yaxis=dict(autorange="reversed"),
    )
    _fig
    return


if __name__ == "__main__":
    app.run()
