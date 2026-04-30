import marimo

__generated_with = "0.19.4"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Cortex Neurotransmitter / Neuromodulator / Neuropeptide Expression

    For each NT/NM/NP system, computes the **proportion of cortical neurons** that:

    - match the **ligand** (NT identity for small-molecule systems, peptide
      gene expression for neuropeptides/hormones), and
    - express **a receptor** (any receptor gene for that system)

    Granularity (class / subclass / supertype / cluster) is selectable.

    ## Data sources

    - **Cells & expression**: existing whole-brain MERFISH bundle at
      `data/processed/mouse_abc_whole_brain/` (`cell_metadata.parquet` +
      `neuropeptide_expression.parquet`), filtered in-notebook to
      `parcellation_division == 'Isocortex'` and to neurons (class number < 30).
    - **Ligand-receptor catalog**: `np_map.csv` (neuropeptides), `hormone_map.csv`
      (hormones), and `nt_map.csv` (small-molecule NTs/NMs).

    ## Definitions

    - **Ligand match for NT systems** uses the canonical ABC `neurotransmitter`
      taxonomy column. A cell with `neurotransmitter` containing the system's
      NT_Class token matches — so `Glut-GABA` co-release cells count as both
      Glut and GABA.
    - **Ligand match for neuropeptides/hormones** uses gene expression: a cell
      counts if it expresses **any** of the system's ligand genes above threshold
      (no canonical metadata exists for peptide release).
    - **Receptor match** always uses gene expression: a cell counts if **any**
      receptor gene for the system exceeds threshold.

    > Note: `neuropeptide_expression.parquet` only contains the ~130 ligand/receptor
    > genes from `np_map.csv` + `hormone_map.csv`. Small-molecule NT receptor
    > genes (Gria*, Gabra*, Drd*, Htr*, etc.) are not in this matrix, so those
    > systems will show only the ligand side here.
    """)
    return


@app.cell
def _():
    EXPRESSION_THRESHOLD = 2.0
    return (EXPRESSION_THRESHOLD,)


@app.cell
def _():
    from pathlib import Path

    import numpy as np
    import pandas as pd
    import plotly.graph_objects as go

    import marimo as mo

    REPO_ROOT = Path(mo.notebook_dir()) / ".."
    WB_DIR = REPO_ROOT / "data" / "processed" / "mouse_abc_whole_brain"
    LR_MAPS_DIR = REPO_ROOT / "data" / "generated" / "mouse_common"
    return LR_MAPS_DIR, WB_DIR, go, mo, np, pd


@app.cell
def _(WB_DIR, pd):
    # Load the whole-brain ABC MERFISH parquets and filter to cortex.
    # `neuropeptide_expression.parquet` is per-cell × gene log2(CPM+1) for the
    # ligand/receptor genes already extracted by build_cluster_ligand_receptor_map.
    # `cell_metadata.parquet` carries class/subclass/supertype/cluster, the
    # canonical `neurotransmitter` taxonomy token, and `parcellation_division`.
    cell_meta_all = pd.read_parquet(WB_DIR / "cell_metadata.parquet")
    if "cell_id" in cell_meta_all.columns:
        cell_meta_all = cell_meta_all.set_index("cell_id")

    expr_all = pd.read_parquet(WB_DIR / "neuropeptide_expression.parquet")

    # Cortex filter: parcellation_division == 'Isocortex'. Adjust if you also
    # want HPF / OLF / CTXsp.
    CORTEX_DIVISIONS = ["Isocortex"]
    cortex_mask = cell_meta_all["parcellation_division"].isin(CORTEX_DIVISIONS)
    cortex_meta = cell_meta_all[cortex_mask]

    # Drop non-neuronal classes (leading number >= 30: Astro-Epen, OPC-Oligo, Vascular, Immune)
    _class_num = cortex_meta["class"].str.extract(r"^(\d+)", expand=False).astype(float)
    cell_meta = cortex_meta[_class_num < 30].copy()

    # Align expression to the filtered cells
    common = expr_all.index.intersection(cell_meta.index)
    expr_df = expr_all.loc[common]
    cell_meta = cell_meta.loc[common]

    print(f"Whole-brain cells: {len(cell_meta_all):,}")
    print(f"Cortex cells (divisions={CORTEX_DIVISIONS}): {len(cortex_meta):,}")
    print(f"Cortex neurons (after non-neuronal filter): {len(cell_meta):,}")
    print(f"Genes in expression matrix: {expr_df.shape[1]}")
    print(f"Class values: {cell_meta['class'].nunique()}")
    print(f"Subclass values: {cell_meta['subclass'].nunique()}")
    print(f"Supertype values: {cell_meta['supertype'].nunique()}")
    print(f"Cluster values: {cell_meta['cluster'].nunique()}")
    return cell_meta, expr_df


@app.cell
def _(LR_MAPS_DIR, pd):
    # Combine all three ligand-receptor maps into one long table.
    # nt_map.csv carries an NT_Class column with the canonical ABC `neurotransmitter`
    # taxonomy token (e.g. "Glut", "GABA"). For those systems the ligand side is
    # measured from the metadata column directly, NOT from gene expression.
    # np_map.csv and hormone_map.csv have no NT_Class — they fall back to gene expression.
    def _load(name, source_label):
        path = LR_MAPS_DIR / name
        if not path.exists():
            return pd.DataFrame()
        df = pd.read_csv(path)
        keep = ["System", "Ligand_Gene", "Receptor_Gene", "Ligand_Class", "NT_Class"]
        keep = [c for c in keep if c in df.columns]
        df = df[keep].copy()
        if "NT_Class" not in df.columns:
            df["NT_Class"] = ""
        df["source"] = source_label
        return df

    parts = [
        _load("nt_map.csv", "neurotransmitter"),
        _load("np_map.csv", "neuropeptide"),
        _load("hormone_map.csv", "hormone"),
    ]
    lr_map = pd.concat([p for p in parts if not p.empty], ignore_index=True)

    # Explode semicolon-separated genes (e.g. heterodimer receptors "Calcrl;Ramp2")
    def _explode(col):
        return (
            lr_map[col]
            .fillna("")
            .astype(str)
            .str.split(";")
            .apply(lambda lst: [g.strip() for g in lst if g.strip()])
        )

    lr_map["Ligand_Gene"] = _explode("Ligand_Gene")
    lr_map["Receptor_Gene"] = _explode("Receptor_Gene")
    lr_map = lr_map.explode("Ligand_Gene").explode("Receptor_Gene").reset_index(drop=True)

    print(f"Total ligand-receptor pairs (after explode): {len(lr_map)}")
    print(f"Systems: {lr_map['System'].nunique()}")
    print(f"By source:\n{lr_map.groupby('source')['System'].nunique()}")
    return (lr_map,)


@app.cell
def _(expr_df, lr_map):
    # Per-system bookkeeping. Each system gets:
    #   - source: 'neurotransmitter' / 'neuropeptide' / 'hormone'
    #   - nt_class: canonical ABC taxonomy token (str), or "" if ligand uses gene expr
    #   - ligands: list of marker genes (only used when nt_class is empty)
    #   - receptors: list of receptor genes
    _available = set(expr_df.columns)

    system_genes = {}
    for _system, _grp in lr_map.groupby("System"):
        _ligands = sorted({g for g in _grp["Ligand_Gene"].unique() if g and g in _available})
        _receptors = sorted({g for g in _grp["Receptor_Gene"].unique() if g and g in _available})
        _nt_class = ""
        for _val in _grp["NT_Class"].dropna().unique():
            if str(_val).strip():
                _nt_class = str(_val).strip()
                break
        if not _ligands and not _receptors and not _nt_class:
            continue
        system_genes[_system] = {
            "source": _grp["source"].iloc[0],
            "nt_class": _nt_class,
            "ligands": _ligands,
            "receptors": _receptors,
        }

    # Sort systems: NT first, then NP, then hormone — alphabetical within
    _SOURCE_ORDER = {"neurotransmitter": 0, "neuropeptide": 1, "hormone": 2}
    sorted_systems = sorted(
        system_genes.keys(),
        key=lambda s: (_SOURCE_ORDER.get(system_genes[s]["source"], 99), s),
    )
    print(f"Systems with at least one available ligand source or receptor gene: {len(sorted_systems)}")
    print(
        f"  Sources: ",
        {k: sum(1 for s in sorted_systems if system_genes[s]['source'] == k) for k in ['neurotransmitter', 'neuropeptide', 'hormone']},
    )
    return sorted_systems, system_genes


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    ## Controls
    """)
    return


@app.cell
def _(mo):
    GRANULARITIES = {
        "All cortex (1 group)": "_all",
        "Class": "class",
        "Subclass": "subclass",
        "Supertype": "supertype",
        "Cluster": "cluster",
    }
    granularity = mo.ui.dropdown(
        options=GRANULARITIES,
        value="Subclass",
        label="Granularity",
    )
    return (granularity,)


@app.cell
def _(EXPRESSION_THRESHOLD, mo):
    threshold_slider = mo.ui.slider(
        start=0.5,
        stop=5.0,
        step=0.25,
        value=EXPRESSION_THRESHOLD,
        label="Expression threshold (log2(CPM+1))",
        show_value=True,
    )
    return (threshold_slider,)


@app.cell
def _(mo):
    SOURCE_FILTERS = {
        "All (NT + NM + NP + hormone)": "all",
        "Small-molecule NT/NM only": "nt_only",
        "Neuropeptides + hormones only": "np_only",
    }
    source_filter = mo.ui.dropdown(
        options=SOURCE_FILTERS,
        value="All (NT + NM + NP + hormone)",
        label="Systems to show",
    )
    return (source_filter,)


@app.cell
def _(granularity, mo, source_filter, threshold_slider):
    mo.vstack([
        granularity,
        threshold_slider,
        source_filter,
    ])
    return


@app.cell
def _(cell_meta, expr_df, granularity, pd, system_genes, threshold_slider):
    # Compute, for each (system, group), the proportion of cells matching the
    # ligand and the proportion expressing any receptor.
    #
    # Ligand match:
    #   - If the system has an `nt_class` (Glut, GABA, ...), use the canonical
    #     ABC `neurotransmitter` taxonomy column. Substring match handles
    #     co-release classes (e.g. nt_class='Glut' matches 'Glut' AND 'Glut-GABA').
    #   - Else (neuropeptide / hormone), a cell counts as expressing the ligand
    #     if ANY ligand gene exceeds threshold.
    # Receptor match: ANY receptor gene exceeds threshold.
    threshold = threshold_slider.value
    grouping = granularity.value

    if grouping == "_all":
        _groups = pd.Series("All cortex neurons", index=expr_df.index, name="group")
    else:
        _groups = cell_meta[grouping].rename("group")

    _expressed = expr_df > threshold  # cells x genes
    _group_n = _groups.value_counts().rename("n_cells")
    _nt_col = cell_meta["neurotransmitter"].fillna("") if "neurotransmitter" in cell_meta.columns else None

    _rows = []
    for _sys, _info in system_genes.items():
        # Ligand mask
        _ligand_mask = None
        if _info["nt_class"]:
            if _nt_col is None:
                continue  # Cannot compute without taxonomy column
            _ligand_mask = _nt_col.str.contains(_info["nt_class"], na=False, regex=False)
        elif _info["ligands"]:
            _ligand_mask = _expressed[_info["ligands"]].any(axis=1)

        # Receptor mask
        _receptor_mask = (
            _expressed[_info["receptors"]].any(axis=1) if _info["receptors"] else None
        )

        for _kind, _mask, _basis in (
            ("ligand", _ligand_mask, "nt_class" if _info["nt_class"] else "gene_expr"),
            ("receptor", _receptor_mask, "gene_expr"),
        ):
            if _mask is None:
                continue
            _counts = _mask.groupby(_groups).sum()
            for _grp_name, _n_match in _counts.items():
                _n = _group_n.get(_grp_name, 0)
                if _n == 0:
                    continue
                _rows.append({
                    "system": _sys,
                    "source": _info["source"],
                    "kind": _kind,
                    "basis": _basis,
                    "group": _grp_name,
                    "n_cells": int(_n),
                    "n_matching": int(_n_match),
                    "pct": 100 * float(_n_match) / float(_n),
                })

    prop_df = pd.DataFrame(_rows)
    print(f"Computed {len(prop_df):,} (system, kind, group) entries at granularity={grouping}, threshold={threshold}")
    if not prop_df.empty:
        print(f"Ligand basis breakdown:\n{prop_df[prop_df['kind']=='ligand'].drop_duplicates('system')['basis'].value_counts().to_string()}")
    return grouping, prop_df, threshold


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    ## Summary table (sorted by ligand %)

    For each system, the proportion of neurons (across all groups within the
    selected granularity) expressing the ligand or any receptor.
    """)
    return


@app.cell
def _(cell_meta, expr_df, np, pd, source_filter, system_genes, threshold_slider):
    # Whole-cortex proportions per system (independent of granularity).
    _thr = threshold_slider.value
    _expressed = expr_df > _thr
    _n_total = len(expr_df)
    _nt_col = cell_meta["neurotransmitter"].fillna("") if "neurotransmitter" in cell_meta.columns else None

    _rows = []
    for _sys, _info in system_genes.items():
        if source_filter.value == "nt_only" and _info["source"] != "neurotransmitter":
            continue
        if source_filter.value == "np_only" and _info["source"] not in ("neuropeptide", "hormone"):
            continue
        if _info["nt_class"] and _nt_col is not None:
            _ligand_pct = 100 * _nt_col.str.contains(_info["nt_class"], na=False, regex=False).sum() / _n_total
            _ligand_basis = f"NT_class={_info['nt_class']}"
        elif _info["ligands"]:
            _ligand_pct = 100 * _expressed[_info["ligands"]].any(axis=1).sum() / _n_total
            _ligand_basis = "gene_expr"
        else:
            _ligand_pct = np.nan
            _ligand_basis = "(none)"
        _receptor_pct = 100 * _expressed[_info["receptors"]].any(axis=1).sum() / _n_total if _info["receptors"] else np.nan
        _rows.append({
            "system": _sys,
            "source": _info["source"],
            "ligand_basis": _ligand_basis,
            "n_receptor_genes": len(_info["receptors"]),
            "pct_cells_expressing_ligand": _ligand_pct,
            "pct_cells_expressing_any_receptor": _receptor_pct,
        })

    summary_df = pd.DataFrame(_rows).sort_values("pct_cells_expressing_ligand", ascending=False).reset_index(drop=True)
    summary_df
    return (summary_df,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    ## Heatmap: % expressing ligand, by group

    Rows = systems. Columns = groups at the selected granularity. Color = % of
    cells in that group expressing **any** ligand gene for the system.
    """)
    return


@app.cell
def _(go, np, pd, prop_df, sorted_systems, source_filter, system_genes):
    def _filter_systems(systems):
        if source_filter.value == "nt_only":
            return [s for s in systems if system_genes[s]["source"] == "neurotransmitter"]
        if source_filter.value == "np_only":
            return [s for s in systems if system_genes[s]["source"] in ("neuropeptide", "hormone")]
        return systems

    def _heatmap(kind, colorscale, title):
        df = prop_df[prop_df["kind"] == kind]
        if df.empty:
            return None
        sys_list = [s for s in _filter_systems(sorted_systems) if s in set(df["system"])]
        if not sys_list:
            return None
        pivot = df.pivot_table(index="system", columns="group", values="pct", aggfunc="first").reindex(sys_list)
        # Sort columns: most-expressing first
        col_order = pivot.mean(axis=0).sort_values(ascending=False).index.tolist()
        pivot = pivot[col_order]

        # Hover text
        hover = [[f"{s}<br>{c}<br>{pivot.loc[s, c]:.1f}% expressing" if pd.notna(pivot.loc[s, c]) else f"{s}<br>{c}<br>(n/a)" for c in pivot.columns] for s in pivot.index]

        zmax = max(5.0, float(np.nanpercentile(pivot.values, 95)))
        n_rows, n_cols = pivot.shape
        fig = go.Figure(data=go.Heatmap(
            z=pivot.values,
            x=pivot.columns.tolist(),
            y=pivot.index.tolist(),
            hovertext=hover,
            hoverinfo="text",
            colorscale=colorscale,
            zmin=0,
            zmax=zmax,
            colorbar=dict(title="% expressing"),
        ))
        fig.update_layout(
            title=title,
            xaxis_title="Group",
            yaxis_title="System",
            height=max(450, 18 * n_rows),
            width=max(700, 12 * n_cols + 250),
            yaxis=dict(autorange="reversed", tickfont=dict(size=10)),
            xaxis=dict(tickangle=45, tickfont=dict(size=9)),
            margin=dict(l=200, b=180),
        )
        return fig

    fig_ligand = _heatmap("ligand", "YlOrRd", "% cortical neurons expressing system LIGAND")
    fig_ligand
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Heatmap: % expressing any receptor, by group
    """)
    return


@app.cell
def _(go, np, pd, prop_df, sorted_systems, source_filter, system_genes):
    def _filter_systems(systems):
        if source_filter.value == "nt_only":
            return [s for s in systems if system_genes[s]["source"] == "neurotransmitter"]
        if source_filter.value == "np_only":
            return [s for s in systems if system_genes[s]["source"] in ("neuropeptide", "hormone")]
        return systems

    _df = prop_df[prop_df["kind"] == "receptor"]
    _sys_list = [s for s in _filter_systems(sorted_systems) if s in set(_df["system"])]
    _pivot = _df.pivot_table(index="system", columns="group", values="pct", aggfunc="first").reindex(_sys_list)
    _col_order = _pivot.mean(axis=0).sort_values(ascending=False).index.tolist()
    _pivot = _pivot[_col_order]

    _hover = [
        [f"{s}<br>{c}<br>{_pivot.loc[s, c]:.1f}% expressing" if pd.notna(_pivot.loc[s, c]) else f"{s}<br>{c}<br>(n/a)" for c in _pivot.columns]
        for s in _pivot.index
    ]
    _zmax = max(5.0, float(np.nanpercentile(_pivot.values, 95)))
    _n_rows, _n_cols = _pivot.shape
    fig_receptor = go.Figure(data=go.Heatmap(
        z=_pivot.values,
        x=_pivot.columns.tolist(),
        y=_pivot.index.tolist(),
        hovertext=_hover,
        hoverinfo="text",
        colorscale="YlGnBu",
        zmin=0,
        zmax=_zmax,
        colorbar=dict(title="% expressing"),
    ))
    fig_receptor.update_layout(
        title="% cortical neurons expressing ANY receptor for system",
        xaxis_title="Group",
        yaxis_title="System",
        height=max(450, 18 * _n_rows),
        width=max(700, 12 * _n_cols + 250),
        yaxis=dict(autorange="reversed", tickfont=dict(size=10)),
        xaxis=dict(tickangle=45, tickfont=dict(size=9)),
        margin=dict(l=200, b=180),
    )
    fig_receptor
    return (fig_receptor,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    ## Per-gene breakdown (debugging)

    Pick a system to see, gene-by-gene, the fraction of cells in each group
    expressing that gene. Useful for sanity-checking which marker is driving the
    "any-of" signal above.
    """)
    return


@app.cell
def _(mo, sorted_systems):
    system_picker = mo.ui.dropdown(
        options=sorted_systems,
        value=sorted_systems[0] if sorted_systems else None,
        label="System",
    )
    system_picker
    return (system_picker,)


@app.cell
def _(
    cell_meta,
    expr_df,
    grouping,
    pd,
    system_genes,
    system_picker,
    threshold,
):
    if system_picker.value is None or system_picker.value not in system_genes:
        per_gene = pd.DataFrame()
    else:
        _info = system_genes[system_picker.value]
        # For NT systems with an nt_class, the "ligand" side is metadata-driven,
        # so the gene list here is just the receptors. Neuropeptides/hormones
        # show both ligand and receptor genes.
        if _info["nt_class"]:
            _gene_list = list(_info["receptors"])
            _gene_kind = {g: "receptor" for g in _info["receptors"]}
        else:
            _gene_list = list(_info["ligands"]) + list(_info["receptors"])
            _gene_kind = {g: "ligand" for g in _info["ligands"]}
            _gene_kind.update({g: "receptor" for g in _info["receptors"]})

        if grouping == "_all":
            _groups = pd.Series("All cortex neurons", index=expr_df.index, name="group")
        else:
            _groups = cell_meta[grouping].rename("group")

        _rows = []
        for _g in _gene_list:
            _expr_bool = expr_df[_g] > threshold
            _counts = _expr_bool.groupby(_groups).agg(["sum", "count"])
            _counts["pct"] = 100 * _counts["sum"] / _counts["count"]
            for _grp, _row in _counts.iterrows():
                _rows.append({
                    "gene": _g,
                    "kind": _gene_kind[_g],
                    "group": _grp,
                    "n_cells": int(_row["count"]),
                    "n_expressing": int(_row["sum"]),
                    "pct": float(_row["pct"]),
                })
        per_gene = pd.DataFrame(_rows).sort_values(["kind", "gene", "pct"], ascending=[True, True, False])
    per_gene
    return (per_gene,)


if __name__ == "__main__":
    app.run()
