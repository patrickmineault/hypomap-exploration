# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "hypomap",
#     "marimo>=0.19.11",
#     "matplotlib",
#     "msgpack-python",
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
    import marimo as mo
    import matplotlib
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    from matplotlib.patches import Polygon
    from matplotlib.collections import PatchCollection
    import msgpack
    import numpy as np
    import pandas as pd
    import plotly.graph_objects as go

    matplotlib.rcParams.update({
        "figure.facecolor": "white",
        "axes.facecolor": "#f5f5f5",
        "font.size": 8,
    })
    return Polygon, go, mo, msgpack, np, pd, plt


@app.cell
def _(mo, msgpack, np, pd):
    _data_dir = mo.notebook_dir() / ".." / "data" / "processed" / "mouse_abc_extended"

    _all_cells = pd.read_parquet(_data_dir / "cells_with_coords.parquet")

    FOCUS_REGIONS = ["PAG", "SCm", "SCs"]
    brainstem_df = _all_cells[_all_cells["region"].isin(FOCUS_REGIONS)].copy()
    brainstem_df = brainstem_df.reset_index(drop=True)

    if "neurotransmitter" not in brainstem_df.columns:
        _meta = pd.read_parquet(
            _data_dir / "cell_metadata.parquet",
            columns=["cell_id", "neurotransmitter"],
        )
        brainstem_df = brainstem_df.merge(_meta, on="cell_id", how="left")

    _np_expr_full = pd.read_parquet(_data_dir / "neuropeptide_expression.parquet")
    bs_expr = _np_expr_full.loc[brainstem_df["cell_id"]].reset_index(drop=True)

    del _all_cells, _np_expr_full

    with open(_data_dir / "coronal_atlas_regions.msgpack", "rb") as _f:
        _region_data = msgpack.unpackb(_f.read(), raw=False)

    boundaries = {
        float(k): {
            r: np.frombuffer(buf, dtype=np.float32).reshape(-1, 2)
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

    print(f"Loaded {len(brainstem_df):,} brainstem cells "
          f"(PAG={sum(brainstem_df['region']=='PAG'):,}, "
          f"SCm={sum(brainstem_df['region']=='SCm'):,}, "
          f"SCs={sum(brainstem_df['region']=='SCs'):,})")
    print(f"Expression matrix: {bs_expr.shape[0]:,} cells x {bs_expr.shape[1]} genes")
    return FOCUS_REGIONS, boundaries, brainstem_df, bs_expr


@app.cell
def _(brainstem_df, mo):
    _all_z = sorted(brainstem_df["z_slice"].unique())
    _three_region_slices = []
    for _z in _all_z:
        _regions = brainstem_df[brainstem_df["z_slice"] == _z]["region"].unique()
        if all(r in _regions for r in ["PAG", "SCm", "SCs"]):
            _three_region_slices.append(_z)

    DEFAULT_SLICES = [-4.67, -4.31, -3.94, -3.58, -3.22, -3.04]
    HEATMAP_N_BINS = 50
    HEATMAP_SIGMA = 0.8
    HEATMAP_THRESHOLD = 3.0
    threshold = 3.0

    slice_selector = mo.ui.multiselect(
        options={f"z={z:.2f}": z for z in _three_region_slices},
        value=[f"z={z:.2f}" for z in DEFAULT_SLICES],
        label="Coronal slices to display",
    )
    mo.md(f"""
    ## Slice Selection
    {slice_selector}
    """)
    return (
        DEFAULT_SLICES,
        HEATMAP_N_BINS,
        HEATMAP_SIGMA,
        HEATMAP_THRESHOLD,
        slice_selector,
        threshold,
    )


@app.cell
def _(DEFAULT_SLICES, slice_selector):
    selected_slices = list(slice_selector.value) if slice_selector.value else DEFAULT_SLICES
    selected_slices = sorted(selected_slices)
    return (selected_slices,)


@app.cell(hide_code=True)
def _(
    FOCUS_REGIONS,
    HEATMAP_N_BINS,
    HEATMAP_SIGMA,
    HEATMAP_THRESHOLD,
    Polygon,
    np,
    plt,
):
    def _draw_boundaries(ax, z, boundaries):
        """Draw PAG/SCm/SCs region boundaries on an axis."""
        bnd_keys = sorted(boundaries.keys())
        z_bnd = min(bnd_keys, key=lambda k: abs(k - z))
        if abs(z_bnd - z) > 0.05:
            return
        for region, coords in boundaries[z_bnd].items():
            base = region.split("-")[0] if "-" in region else region
            if base not in FOCUS_REGIONS:
                continue
            poly = Polygon(coords, closed=True, fill=False,
                           edgecolor="k", linewidth=0.8, alpha=0.6)
            ax.add_patch(poly)
            cx, cy = coords.mean(axis=0)
            ax.text(cx, cy, region, fontsize=5, ha="center", va="center",
                    color="white", fontweight="bold",
                    bbox=dict(boxstyle="round,pad=0.15", fc="#333333", alpha=0.85, lw=0))

    def _strip_spines(ax):
        """Remove all four spines from an axis."""
        for sp in ax.spines.values():
            sp.set_visible(False)

    def make_spatial_figure(
        cells_df, color_values, selected_slices, boundaries,
        title, colorscale="magma", cmin=None, cmax=None,
        marker_size=1, n_cols=3, colorbar_title="Expression",
    ):
        """Multi-slice spatial expression map (matplotlib)."""
        n_slices = len(selected_slices)
        n_rows = int(np.ceil(n_slices / n_cols))

        arr_x = cells_df["x"].values
        arr_y = cells_df["y"].values
        arr_z = cells_df["z_slice"].values
        arr_c = np.asarray(color_values, dtype=np.float32)

        if cmin is None or cmax is None:
            pos = arr_c[arr_c > 0]
            auto_min = float(np.percentile(pos, 2)) if len(pos) > 0 else 0.0
            auto_max = float(np.percentile(pos, 98)) if len(pos) > 0 else 1.0
            if cmin is None:
                cmin = auto_min
            if cmax is None:
                cmax = auto_max

        fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols * 3, n_rows * 3),
                                 sharex=True, sharey=True, squeeze=False,
                                 layout="constrained")

        im = None
        for i, z in enumerate(selected_slices):
            ax = axes[i // n_cols, i % n_cols]
            mask = arr_z == z
            sx, sy, sc = arr_x[mask], arr_y[mask], arr_c[mask]

            bg = (sc <= 0) | np.isnan(sc)
            fg = ~bg

            if bg.sum() > 0:
                ax.scatter(sx[bg], sy[bg], s=0.3, c="#ddd", alpha=0.3, rasterized=True)
            if fg.sum() > 0:
                im = ax.scatter(sx[fg], sy[fg], s=marker_size, c=sc[fg],
                                cmap=colorscale, vmin=cmin, vmax=cmax,
                                alpha=0.8, rasterized=True)
            _draw_boundaries(ax, z, boundaries)
            ax.set_title(f"z = {z:.2f}", fontsize=7)

            ax.set_aspect("equal")
            ax.tick_params(labelsize=5)
            _strip_spines(ax)

        axes[0, 0].invert_yaxis()

        for i in range(n_slices, n_rows * n_cols):
            axes[i // n_cols, i % n_cols].set_visible(False)

        if im is not None:
            fig.colorbar(im, ax=axes.ravel().tolist(), label=colorbar_title,
                         shrink=0.6)
        fig.suptitle(title, fontsize=10)
        return plt.gca()

    def make_categorical_figure(
        cells_df, category_values, category_colors, selected_slices, boundaries,
        title, marker_size=1, n_cols=3,
    ):
        """Multi-slice categorical spatial map (matplotlib)."""
        n_slices = len(selected_slices)
        n_rows = int(np.ceil(n_slices / n_cols))

        arr_x = cells_df["x"].values
        arr_y = cells_df["y"].values
        arr_z = cells_df["z_slice"].values
        cats = np.asarray(category_values)

        # Plot order: background categories last so foreground is on top
        ordered_cats = list(category_colors.keys())

        fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols * 3, n_rows * 3),
                                 sharex=True, sharey=True, squeeze=False,
                                 layout="constrained")

        for i, z in enumerate(selected_slices):
            ax = axes[i // n_cols, i % n_cols]
            mask = arr_z == z

            for cat in ordered_cats:
                cat_mask = mask & (cats == cat)
                if cat_mask.sum() == 0:
                    continue
                color = category_colors[cat]
                alpha = 0.15 if color in ("#e0e0e0", "#d0d0d0", "#ddd") else 0.7
                s = 0.3 if color in ("#e0e0e0", "#d0d0d0", "#ddd") else marker_size
                ax.scatter(arr_x[cat_mask], arr_y[cat_mask], s=s, c=color,
                           alpha=alpha, label=cat if i == 0 else None, rasterized=True)

            _draw_boundaries(ax, z, boundaries)
            ax.set_title(f"z = {z:.2f}", fontsize=7)

            ax.set_aspect("equal")
            ax.tick_params(labelsize=5)
            _strip_spines(ax)

        axes[0, 0].invert_yaxis()

        for i in range(n_slices, n_rows * n_cols):
            axes[i // n_cols, i % n_cols].set_visible(False)

        handles, labels = axes[0, 0].get_legend_handles_labels()
        if handles:
            fig.legend(handles, labels, loc="outside lower center",
                       ncol=min(6, len(handles)), fontsize=6, markerscale=4, frameon=False)
        fig.suptitle(title, fontsize=10)
        return plt.gca()

    def _bin_proportion(arr_x, arr_y, arr_expressing, x_edges, y_edges, sigma):
        """Bin cells into a 2D grid and return smoothed p(expressing) per bin.

        Applies gaussian smoothing (sigma in bin units) to both numerator
        and denominator before dividing, giving a smooth local average.
        """
        from scipy.ndimage import gaussian_filter
        counts_all, _, _ = np.histogram2d(arr_x, arr_y, bins=[x_edges, y_edges])
        counts_expr, _, _ = np.histogram2d(
            arr_x[arr_expressing], arr_y[arr_expressing], bins=[x_edges, y_edges]
        )
        smooth_all = gaussian_filter(counts_all.astype(np.float64), sigma=sigma)
        smooth_expr = gaussian_filter(counts_expr.astype(np.float64), sigma=sigma)
        with np.errstate(divide="ignore", invalid="ignore"):
            grid = np.where(smooth_all > 0.5, smooth_expr / smooth_all, np.nan)
        return grid

    def make_heatmap_figure(
        cells_df, expression_values, selected_slices, boundaries,
        title, n_cols=3, colorbar_title="P(expr > threshold)",
        threshold=HEATMAP_THRESHOLD, n_bins=HEATMAP_N_BINS, sigma=HEATMAP_SIGMA,
    ):
        """Multi-slice binned proportion heatmap.

        For each spatial bin, shows the fraction of cells with
        expression > threshold.
        """
        n_slices = len(selected_slices)
        n_rows = int(np.ceil(n_slices / n_cols))

        arr_x = cells_df["x"].values
        arr_y = cells_df["y"].values
        arr_z = cells_df["z_slice"].values
        arr_e = np.asarray(expression_values, dtype=np.float32) > threshold

        # Global bin edges from all displayed slices
        disp = np.isin(arr_z, selected_slices)
        x_edges = np.linspace(arr_x[disp].min(), arr_x[disp].max(), n_bins + 1)
        y_edges = np.linspace(arr_y[disp].min(), arr_y[disp].max(), n_bins + 1)

        fig, axes = plt.subplots(n_rows, n_cols, figsize=(n_cols * 3, n_rows * 3),
                                 sharex=True, sharey=True, squeeze=False,
                                 layout="constrained")
        im = None
        for i, z in enumerate(selected_slices):
            ax = axes[i // n_cols, i % n_cols]
            mask = arr_z == z
            grid = _bin_proportion(arr_x[mask], arr_y[mask], arr_e[mask], x_edges, y_edges, sigma=sigma)

            im = ax.pcolormesh(
                x_edges, y_edges, grid.T,
                cmap="magma", vmin=0, vmax=1, rasterized=True,
            )
            _draw_boundaries(ax, z, boundaries)
            ax.set_title(f"z = {z:.2f}", fontsize=7)
            ax.set_aspect("equal")
            ax.tick_params(labelsize=5)
            _strip_spines(ax)

        axes[0, 0].invert_yaxis()

        for i in range(n_slices, n_rows * n_cols):
            axes[i // n_cols, i % n_cols].set_visible(False)

        if im is not None:
            fig.colorbar(im, ax=axes.ravel().tolist(), label=colorbar_title,
                         shrink=0.6)
        fig.suptitle(title, fontsize=10)
        return plt.gca()

    def make_ligand_receptor_heatmap(
        cells_df, bs_expr, lig_gene, rec_gene, lig_label, rec_label,
        selected_slices, boundaries,
        threshold=HEATMAP_THRESHOLD, n_bins=HEATMAP_N_BINS, sigma=HEATMAP_SIGMA,
    ):
        """2-row heatmap: ligand on top, receptor on bottom."""
        n_slices = len(selected_slices)
        arr_x = cells_df["x"].values
        arr_y = cells_df["y"].values
        arr_z = cells_df["z_slice"].values

        lig_vals = bs_expr[lig_gene].values.astype(np.float32) if lig_gene in bs_expr.columns else np.zeros(len(cells_df))
        rec_vals = bs_expr[rec_gene].values.astype(np.float32) if rec_gene in bs_expr.columns else np.zeros(len(cells_df))

        # Global bin edges
        disp = np.isin(arr_z, selected_slices)
        x_edges = np.linspace(arr_x[disp].min(), arr_x[disp].max(), n_bins + 1)
        y_edges = np.linspace(arr_y[disp].min(), arr_y[disp].max(), n_bins + 1)

        fig, axes = plt.subplots(2, n_slices, figsize=(n_slices * 2.5, 5),
                                 sharex=True, sharey=True, squeeze=False,
                                 layout="constrained")

        for row_idx, (vals, gene, label) in enumerate([
            (lig_vals, lig_gene, lig_label),
            (rec_vals, rec_gene, rec_label),
        ]):
            expressing = vals > threshold
            for col_idx, z in enumerate(selected_slices):
                ax = axes[row_idx, col_idx]
                mask = arr_z == z
                grid = _bin_proportion(arr_x[mask], arr_y[mask], expressing[mask], x_edges, y_edges, sigma=sigma)

                ax.pcolormesh(
                    x_edges, y_edges, grid.T,
                    cmap="magma", vmin=0, vmax=1, rasterized=True,
                )
                _draw_boundaries(ax, z, boundaries)
                ax.set_aspect("equal")
                ax.tick_params(labelsize=4)
                _strip_spines(ax)
                if row_idx == 0:
                    ax.set_title(f"z={z:.2f}", fontsize=6)
                if col_idx == 0:
                    ax.set_ylabel(f"{gene}\n({label})", fontsize=6)

        axes[0, 0].invert_yaxis()
        fig.suptitle(f"{lig_gene} / {rec_gene}  (P(log2 expr > {threshold}))", fontsize=10)
        return plt.gca()
    return (
        make_categorical_figure,
        make_heatmap_figure,
        make_ligand_receptor_heatmap,
        make_spatial_figure,
    )


@app.cell(hide_code=True)
def _(bs_expr, np):
    def get_expr(gene):
        if gene in bs_expr.columns:
            return bs_expr[gene].values.astype(np.float32)
        return np.zeros(bs_expr.shape[0], dtype=np.float32)

    def gene_available(gene):
        return gene in bs_expr.columns
    return gene_available, get_expr


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    # Section 1: Anatomical Orientation -- Reference Maps for PAG, SCm, and SCs

    - **PAG** (periaqueductal gray) -- core affective pain and aversion circuitry
    - **SCm** (superior colliculus, motor/deep layers) -- defensive orienting, approach/avoidance
    - **SCs** (superior colliculus, sensory/superficial layers) -- primarily visual/sensory, likely **safe to preserve**

    **Caveat:** PAG is NOT segmented into functional columns (dlPAG, lPAG, vlPAG).
    If aversive markers cluster spatially within PAG, a follow-up notebook should validate sub-regional boundaries.
    """)
    return


@app.cell
def _(boundaries, brainstem_df, make_categorical_figure, selected_slices):
    _cats = brainstem_df["is_neuron"].map({True: "Neuron", False: "Non-neuron"}).values
    make_categorical_figure(
        brainstem_df, _cats,
        {"Non-neuron": "#e0e0e0", "Neuron": "#2b5ea7"},
        selected_slices, boundaries,
        title="Plot 1: Pan-neuronal Reference Map",
    )
    return


@app.cell
def _(boundaries, brainstem_df, make_categorical_figure, np, selected_slices):
    _nt_type = np.full(len(brainstem_df), "Non-neuron", dtype=object)
    _cls = brainstem_df["class"].values
    _is_n = brainstem_df["is_neuron"].values
    for _i in range(len(brainstem_df)):
        if not _is_n[_i]:
            continue
        _c = str(_cls[_i])
        if "Glut" in _c:
            _nt_type[_i] = "Glutamatergic"
        elif "GABA" in _c:
            _nt_type[_i] = "GABAergic"
        elif "Sero" in _c:
            _nt_type[_i] = "Serotonergic"
        elif "Dopa" in _c:
            _nt_type[_i] = "Dopaminergic"
        else:
            _nt_type[_i] = "Other neuron"

    make_categorical_figure(
        brainstem_df, _nt_type,
        {
            "Non-neuron": "#e0e0e0",
            "Other neuron": "#999999",
            "Glutamatergic": "#d62728",
            "GABAergic": "#1f77b4",
            "Serotonergic": "#ff7f0e",
            "Dopaminergic": "#2ca02c",
        },
        selected_slices, boundaries,
        title="Plot 2: Excitatory (Glut) vs Inhibitory (GABA) Landscape",
    )
    return


@app.cell
def _(brainstem_df, go, mo):
    _neurons = brainstem_df[brainstem_df["is_neuron"]].copy()
    _cls = _neurons["class"].astype(str)
    _neurons["nt_type"] = "Other"
    _neurons.loc[_cls.str.contains("Glut"), "nt_type"] = "Glutamatergic"
    _neurons.loc[_cls.str.contains("GABA"), "nt_type"] = "GABAergic"
    _neurons.loc[_cls.str.contains("Sero"), "nt_type"] = "Serotonergic"
    _neurons.loc[_cls.str.contains("Dopa"), "nt_type"] = "Dopaminergic"

    _counts = _neurons.groupby(["region", "nt_type"]).size().unstack(fill_value=0)
    _colors = {
        "Glutamatergic": "#d62728", "GABAergic": "#1f77b4",
        "Serotonergic": "#ff7f0e", "Dopaminergic": "#2ca02c", "Other": "#999999",
    }
    _fig = go.Figure()
    for _nt in ["Glutamatergic", "GABAergic", "Serotonergic", "Dopaminergic", "Other"]:
        if _nt in _counts.columns:
            _fig.add_trace(go.Bar(x=_counts.index.tolist(), y=_counts[_nt].tolist(),
                                  name=_nt, marker_color=_colors.get(_nt, "#999")))
    _fig.update_layout(barmode="stack", title="Plot 3: Neuron Composition by Region",
                       xaxis_title="Region", yaxis_title="Count", height=350, width=500)
    mo.hstack([_fig])
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    # Section 2: The Opioid System -- The Core of Pain Affect Modulation

    - **Enkephalin / mu-opioid (Penk -> Oprm1):** Analgesia in vlPAG. The system morphine hijacks.
    - **Enkephalin / delta-opioid (Penk -> Oprd1):** Emotional regulation, anxiolysis.
    - **Dynorphin / kappa-opioid (Pdyn -> Oprk1):** **Aversion/dysphoria pathway.**
      The single highest-priority target to remove for eliminating suffering.
    """)
    return


@app.cell
def _(
    boundaries,
    brainstem_df,
    bs_expr,
    make_ligand_receptor_heatmap,
    selected_slices,
):
    make_ligand_receptor_heatmap(
        brainstem_df, bs_expr, "Penk", "Oprm1",
        "proenkephalin", "mu-opioid receptor",
        selected_slices, boundaries,
    )
    return


@app.cell(hide_code=True)
def _(
    boundaries,
    brainstem_df,
    bs_expr,
    make_ligand_receptor_heatmap,
    mo,
    selected_slices,
):
    mo.md(r"""
    **HIGH-PRIORITY EXCLUSION TARGET:** Wherever Pdyn and Oprk1 co-localize in PAG,
    that is a candidate suffering node.
    """)
    make_ligand_receptor_heatmap(
        brainstem_df, bs_expr, "Pdyn", "Oprk1",
        "prodynorphin -- aversion", "kappa-opioid -- suffering",
        selected_slices, boundaries,
    )
    return


@app.cell
def _(
    FOCUS_REGIONS,
    HEATMAP_N_BINS,
    HEATMAP_SIGMA,
    HEATMAP_THRESHOLD,
    Polygon,
    boundaries,
    brainstem_df,
    bs_expr,
    np,
    plt,
):
    # Plot 6: Opioid receptor comparison on single slice — binned heatmap
    _z_ref = -3.94
    _receptors = [("Oprm1", "mu -- analgesia"), ("Oprd1", "delta -- emotional"),
                  ("Oprk1", "kappa -- suffering"), ("Oprl1", "nociceptin")]
    _arr_x = brainstem_df["x"].values
    _arr_y = brainstem_df["y"].values
    _arr_z = brainstem_df["z_slice"].values
    _mask = _arr_z == _z_ref

    _x_edges = np.linspace(_arr_x[_mask].min(), _arr_x[_mask].max(), HEATMAP_N_BINS + 1)
    _y_edges = np.linspace(_arr_y[_mask].min(), _arr_y[_mask].max(), HEATMAP_N_BINS + 1)

    _fig, _axes = plt.subplots(1, 4, figsize=(14, 3.5), sharex=True, sharey=True)
    for _i, (_gene, _label) in enumerate(_receptors):
        _ax = _axes[_i]
        _vals = bs_expr[_gene].values.astype(np.float32) if _gene in bs_expr.columns else np.zeros(len(brainstem_df))
        _expressing = (_vals > HEATMAP_THRESHOLD) & _mask
        from scipy.ndimage import gaussian_filter as _gf
        _counts_all, _, _ = np.histogram2d(_arr_x[_mask], _arr_y[_mask], bins=[_x_edges, _y_edges])
        _counts_expr, _, _ = np.histogram2d(_arr_x[_expressing], _arr_y[_expressing], bins=[_x_edges, _y_edges])
        _smooth_all = _gf(_counts_all.astype(np.float64), sigma=HEATMAP_SIGMA)
        _smooth_expr = _gf(_counts_expr.astype(np.float64), sigma=HEATMAP_SIGMA)
        with np.errstate(divide="ignore", invalid="ignore"):
            _grid = np.where(_smooth_all > 0.5, _smooth_expr / _smooth_all, np.nan)

        _ax.pcolormesh(_x_edges, _y_edges, _grid.T, cmap="magma", vmin=0, vmax=1, rasterized=True)

        # Boundaries
        _bnd_keys = sorted(boundaries.keys())
        _z_bnd = min(_bnd_keys, key=lambda k: abs(k - _z_ref))
        if abs(_z_bnd - _z_ref) < 0.05:
            for _region, _coords in boundaries[_z_bnd].items():
                _base = _region.split("-")[0] if "-" in _region else _region
                if _base in FOCUS_REGIONS:
                    _ax.add_patch(Polygon(_coords, closed=True, fill=False, ec="k", lw=0.8, alpha=0.6))

        _ax.set_title(f"{_gene} ({_label})", fontsize=7)
        _ax.set_aspect("equal")
        _ax.tick_params(labelsize=4)
        for _sp in _ax.spines.values():
            _sp.set_visible(False)

    _axes[0].invert_yaxis()
    _fig.suptitle(f"Plot 6: Opioid Receptor Comparison (z={_z_ref}, P(log2 > {HEATMAP_THRESHOLD}))", fontsize=10)
    _fig.tight_layout()
    plt.gca()
    return


@app.cell
def _(mo):
    mo.md(r"""
    ---
    # Section 3: Neuropeptide Systems Driving Aversion, Panic, and Fear

    - **Substance P (Tac1 -> Tacr1):** Nociceptive-affective processing. Strong exclusion candidate.
    - **CRH (Crh -> Crhr1):** Stress-induced aversion, anxiety, learned helplessness via vlPAG.
    - **CCK (Cck -> Cckbr):** One of the most potent panicogenic substances. Panic attacks from dlPAG.
    - **CGRP (Calca):** Aversive visceral signals (nausea, visceral pain, malaise).
    """)
    return


@app.cell
def _(
    boundaries,
    brainstem_df,
    bs_expr,
    make_ligand_receptor_heatmap,
    selected_slices,
):
    make_ligand_receptor_heatmap(
        brainstem_df, bs_expr, "Tac1", "Tacr1",
        "substance P -- nociceptive affect", "NK1 receptor",
        selected_slices, boundaries,
    )
    return


@app.cell
def _(
    boundaries,
    brainstem_df,
    bs_expr,
    make_ligand_receptor_heatmap,
    selected_slices,
):
    make_ligand_receptor_heatmap(
        brainstem_df, bs_expr, "Crh", "Crhr1",
        "CRH -- stress/aversion", "CRF receptor 1",
        selected_slices, boundaries,
    )
    return


@app.cell
def _(
    boundaries,
    brainstem_df,
    bs_expr,
    make_ligand_receptor_heatmap,
    selected_slices,
):
    make_ligand_receptor_heatmap(
        brainstem_df, bs_expr, "Cck", "Cckbr",
        "CCK -- panic peptide", "CCK-B receptor",
        selected_slices, boundaries,
    )
    return


@app.cell
def _(
    boundaries,
    brainstem_df,
    get_expr,
    make_heatmap_figure,
    mo,
    selected_slices,
):
    mo.vstack([
        mo.md("### CGRP System"),
        make_heatmap_figure(brainstem_df, get_expr("Calca"), selected_slices, boundaries,
                            "Calca (CGRP -- visceral aversion)"),
        make_heatmap_figure(brainstem_df, get_expr("Calcrl"), selected_slices, boundaries,
                            "Calcrl (CLR -- CGRP receptor component)"),
        make_heatmap_figure(brainstem_df, get_expr("Ramp1"), selected_slices, boundaries,
                            "Ramp1 (CGRP receptor component)"),
    ])
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    # Section 4: Monoamine Systems -- Emotional Tone and Arousal

    - **Dorsal raphe (serotonin):** Within/ventral to PAG. Complex: some pro-aversive, some anti-aversive.
    - **VTA (dopamine):** Small population in vlPAG. Positive affect -- might want to KEEP.
    - **Locus coeruleus (norepinephrine):** Adjacent to PAG. Pain suppression + arousal/fear.

    Cell-type markers (Tph2, Th, Chat) are not in the MERFISH panel;
    we use ABC atlas class/neurotransmitter labels as proxies.
    """)
    return


@app.cell
def _(boundaries, brainstem_df, make_categorical_figure, np, selected_slices):
    _nt = brainstem_df["neurotransmitter"].fillna("").values
    _cls = brainstem_df["class"].fillna("").values

    _amino_type = np.full(len(brainstem_df), "Other", dtype=object)
    for _i in range(len(brainstem_df)):
        if "Sero" in str(_cls[_i]) or _nt[_i] == "Sero":
            _amino_type[_i] = "Serotonergic"
        elif "Dopa" in str(_cls[_i]) or _nt[_i] == "Dopa":
            _amino_type[_i] = "Dopaminergic"
        elif _nt[_i] == "Chol":
            _amino_type[_i] = "Cholinergic"
        elif _nt[_i] == "Nora":
            _amino_type[_i] = "Noradrenergic"

    for _t in ["Serotonergic", "Dopaminergic", "Cholinergic", "Noradrenergic"]:
        _c = np.sum(_amino_type == _t)
        _p = np.sum((_amino_type == _t) & (brainstem_df["region"].values == "PAG"))
        print(f"  {_t}: {_c} total ({_p} in PAG)")

    make_categorical_figure(
        brainstem_df, _amino_type,
        {"Other": "#e0e0e0", "Serotonergic": "#ff7f0e", "Dopaminergic": "#2ca02c",
         "Cholinergic": "#9467bd", "Noradrenergic": "#17becf"},
        selected_slices, boundaries,
        title="Plot 14: Combined Aminergic Overview",
        marker_size=3,
    )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    # Section 5: Cholinergic and Arousal Systems -- What to Preserve

    Cholinergic brainstem nuclei are critical for arousal, REM sleep, motor function -- likely "keep" candidates.
    Orexin (Hcrt) and MCH (Pmch) are hypothalamic, but their **receptors** may be expressed here.
    """)
    return


@app.cell
def _(
    boundaries,
    brainstem_df,
    gene_available,
    get_expr,
    make_heatmap_figure,
    mo,
    selected_slices,
):
    _figs = []
    for _gene, _label in [
        ("Hcrtr1", "Orexin receptor 1 -- arousal"),
        ("Hcrtr2", "Orexin receptor 2 -- arousal"),
        ("Pmch", "MCH -- sleep/energy"),
    ]:
        if gene_available(_gene):
            _figs.append(make_heatmap_figure(
                brainstem_df, get_expr(_gene), selected_slices, boundaries,
                f"{_gene} ({_label})",
            ))
        else:
            _figs.append(mo.md(f"*{_gene} not in MERFISH panel*"))
    mo.vstack([mo.md("### Arousal Receptors")] + _figs)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    # Section 6: Co-expression Analysis -- Identifying Core Affective Pain Nodes

    Neurons co-expressing multiple aversive markers are stronger candidates for
    "affective suffering nodes" than any single-marker population.
    """)
    return


@app.cell
def _(
    boundaries,
    brainstem_df,
    get_expr,
    make_categorical_figure,
    np,
    selected_slices,
    threshold,
):
    _tac1 = get_expr("Tac1")
    _oprk1 = get_expr("Oprk1")

    _coexpr = np.full(len(brainstem_df), "Neither", dtype=object)
    _coexpr[(_tac1 > threshold) & (_oprk1 <= threshold)] = "Tac1 only"
    _coexpr[(_tac1 <= threshold) & (_oprk1 > threshold)] = "Oprk1 only"
    _coexpr[(_tac1 > threshold) & (_oprk1 > threshold)] = "Both (suffering node)"

    _n = np.sum(_coexpr == "Both (suffering node)")
    _np = np.sum((_coexpr == "Both (suffering node)") & (brainstem_df["region"].values == "PAG"))
    print(f"Tac1+Oprk1 co-expressing: {_n} total, {_np} in PAG")

    make_categorical_figure(
        brainstem_df, _coexpr,
        {"Neither": "#e0e0e0", "Tac1 only": "#2ca02c",
         "Oprk1 only": "#1f77b4", "Both (suffering node)": "#ff0000"},
        selected_slices, boundaries,
        title="Plot 17: Pain Affect Node -- Tac1 + Oprk1 Co-expression",
        marker_size=2,
    )
    return


@app.cell
def _(
    boundaries,
    brainstem_df,
    get_expr,
    make_categorical_figure,
    np,
    selected_slices,
    threshold,
):
    _cck = get_expr("Cck")
    _crh = get_expr("Crh")

    _coexpr = np.full(len(brainstem_df), "Neither", dtype=object)
    _coexpr[(_cck > threshold) & (_crh <= threshold)] = "Cck only"
    _coexpr[(_cck <= threshold) & (_crh > threshold)] = "Crh only"
    _coexpr[(_cck > threshold) & (_crh > threshold)] = "Both (panic node)"

    _n = np.sum(_coexpr == "Both (panic node)")
    _np = np.sum((_coexpr == "Both (panic node)") & (brainstem_df["region"].values == "PAG"))
    print(f"Cck+Crh co-expressing: {_n} total, {_np} in PAG")

    make_categorical_figure(
        brainstem_df, _coexpr,
        {"Neither": "#e0e0e0", "Cck only": "#ff7f0e",
         "Crh only": "#9467bd", "Both (panic node)": "#ff0000"},
        selected_slices, boundaries,
        title="Plot 18: Panic Node -- Cck + Crh Co-expression",
        marker_size=2,
    )
    return


@app.cell
def _(
    boundaries,
    brainstem_df,
    get_expr,
    make_spatial_figure,
    np,
    selected_slices,
    threshold,
):
    _aversive = ["Pdyn", "Oprk1", "Tac1", "Tacr1", "Crh", "Crhr1", "Cck", "Cckbr", "Calca"]
    _burden = np.zeros(len(brainstem_df), dtype=np.float32)
    for _g in _aversive:
        _vals = get_expr(_g)
        _burden += (_vals > threshold).astype(np.float32)

    _regions = brainstem_df["region"].values
    _neurons = brainstem_df["is_neuron"].values
    for _r in ["PAG", "SCm", "SCs"]:
        _m = (_regions == _r) & _neurons
        print(f"  {_r}: mean burden = {_burden[_m].mean():.2f}")

    make_spatial_figure(
        brainstem_df, _burden, selected_slices, boundaries,
        "Plot 19: Aversive Gene Burden (# expressed per cell)",
        colorscale="YlOrRd", cmin=0, cmax=5,
        colorbar_title="# genes",
    )
    return


@app.cell
def _(brainstem_df, get_expr, go, mo, np):
    _aversive = ["Pdyn", "Oprk1", "Tac1", "Tacr1", "Crh", "Crhr1", "Cck", "Cckbr", "Calca"]
    _regions = ["PAG", "SCm", "SCs"]
    _neurons = brainstem_df["is_neuron"].values
    _region_col = brainstem_df["region"].values

    _pct = np.zeros((len(_aversive), len(_regions)), dtype=np.float32)
    for _gi, _gene in enumerate(_aversive):
        _vals = get_expr(_gene)
        for _ri, _reg in enumerate(_regions):
            _mask = _neurons & (_region_col == _reg)
            _n = _mask.sum()
            _pct[_gi, _ri] = 100.0 * ((_vals > 0) & _mask).sum() / _n if _n > 0 else 0.0

    _text = [[f"{_pct[_gi, _ri]:.1f}%" for _ri in range(3)] for _gi in range(len(_aversive))]

    _fig = go.Figure(data=go.Heatmap(
        z=_pct.tolist(), x=_regions, y=_aversive,
        colorscale="YlOrRd", text=_text, texttemplate="%{text}",
        colorbar={"title": "% neurons"},
    ))
    _fig.update_layout(
        title="Molecular Exclusion Targets for Affective-Pain-Free Bodyoid Design",
        height=450, width=450, yaxis={"autorange": "reversed"},
    )
    mo.vstack([
        mo.md("### Plot 20: Exclusion Summary -- the hit list"),
        _fig,
    ])
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    # Section 7: Social and Affiliative Affect -- A Gray Area

    Oxytocin and vasopressin modulate social bonding and autonomic function.
    Removing them eliminates social affect but also beneficial neuromodulation.
    A philosophical design choice, not a clear-cut exclusion.
    """)
    return


@app.cell
def _(
    boundaries,
    brainstem_df,
    get_expr,
    make_heatmap_figure,
    mo,
    selected_slices,
):
    mo.vstack([
        make_heatmap_figure(brainstem_df, get_expr("Oxtr"), selected_slices, boundaries,
                            "Oxtr (Oxytocin receptor -- social bonding)"),
        make_heatmap_figure(brainstem_df, get_expr("Avpr1a"), selected_slices, boundaries,
                            "Avpr1a (V1a receptor -- vasopressin)"),
        make_heatmap_figure(brainstem_df, get_expr("Avp"), selected_slices, boundaries,
                            "Avp (Vasopressin)"),
    ])
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    # Section 8: Summary and Recommendations for Jarod

    ## Highest-priority exclusion targets
    - **Pdyn / Oprk1** (dynorphin/kappa -- dysphoria and aversion)
    - **Tac1 / Tacr1** (substance P -- nociceptive affect)
    - **Crh / Crhr1** (CRH -- stress-induced aversion)
    - **Cck / Cckbr** (CCK -- panic)

    ## Second-tier targets
    - **Penk / Oprm1** (enkephalin/mu -- pain suppression, but implies pain signal exists)
    - **Serotonin neurons** (mixed anxiogenic/anxiolytic)
    - **Calca / CGRP** (visceral aversion relay)

    ## Likely safe to preserve
    - **SCs neurons** (sensory SC -- visual processing)
    - **Cholinergic neurons** (arousal, motor)
    - **Dopaminergic neurons in ventral PAG** (positive affect)

    ## Gray areas (ethical/philosophical)
    - **Oxtr / Avp / Avpr1a** (social affect)
    - **Opioid analgesia circuitry** (prevent pain entirely, or keep the brake?)
    """)
    return


if __name__ == "__main__":
    app.run()
