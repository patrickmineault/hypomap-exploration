import marimo

__generated_with = "0.19.4"
app = marimo.App()


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Distribution of Neuropeptide & Neuromodulator Receptors Across Neurons

    For each neuron, count how many NP **systems** it expresses a ligand for
    (and separately a receptor for), using the same definitions and thresholds
    as the HypoMap app. Also includes classical **neuromodulator** receptors
    (dopamine, serotonin, norepinephrine, acetylcholine, histamine).

    - **NP systems** defined in `np_map.csv` (75 systems, 72 ligand genes, 89 receptor genes)
    - **Neuromodulator receptors** from the MERFISH 500-gene panel
    - **Expression threshold**: 3.0 log₂(CPM) (same as preprocessing pipeline)
    - **Receptor heterodimer logic**: AND — all subunit genes must exceed threshold
    - **Scope**: neurons only (`is_neuron == True`)

    Five regional slices:
    1. **Whole brain** — all ~2.1M neurons in the ABC Atlas MERFISH dataset
    2. **Isocortex** — ~658k cortical neurons
    3. **Hippocampus** — ~177k HPF neurons
    4. **Hypothalamus** — ~74k HY neurons
    5. **ARH + ME + VMH** — ~8k neurons in arcuate, median eminence, and ventromedial hypothalamus
    """)
    return


@app.cell
def _():
    import marimo as mo
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    from pathlib import Path

    DATA_DIR = Path(mo.notebook_dir()) / ".." / "data"
    EXPRESSION_THRESHOLD = 3.0

    # Genes to exclude: ubiquitously expressed, uninformative for NP signaling
    BLACKLISTED_GENES = {
        "Pcsk1n",    # ProSAAS — 100% of neurons, proprotein convertase not a signaling NP
        "Adcyap1r1", # PACAP receptor — 96% of neurons
        "Pdyn", "Adcyap1" # For argument's sake
    }
    return BLACKLISTED_GENES, DATA_DIR, EXPRESSION_THRESHOLD, mo, np, pd, plt


@app.cell
def _(DATA_DIR, pd):
    # Load NP system definitions (same parsing as app.py load_np_systems)
    np_map = pd.read_csv(DATA_DIR / "processed" / "mouse_common" / "np_map.csv")

    systems = {}
    for _system in np_map["System"].unique():
        _system_df = np_map[np_map["System"] == _system]

        _ligand_genes = set()
        for _lg in _system_df["Ligand_Gene"].dropna().unique():
            for _g in _lg.split(";"):
                _ligand_genes.add(_g.strip())

        _receptor_complexes = []
        for _rg in _system_df["Receptor_Gene"].dropna().unique():
            _genes = tuple(_g.strip() for _g in _rg.split(";"))
            if _genes not in _receptor_complexes:
                _receptor_complexes.append(_genes)

        systems[_system] = {
            "ligands": list(_ligand_genes),
            "receptors": _receptor_complexes,
        }

    # Add classical neuromodulator systems
    # Receptors available in the MERFISH 500-gene panel
    neuromod_systems = {
        "Dopamine": {
            "ligands": [],  # Th/Ddc/Slc6a3 not in panel
            "receptors": [("Drd1",), ("Drd3",), ("Drd5",)],
        },
        "Serotonin": {
            "ligands": [],  # Tph2/Slc6a4 not in panel
            "receptors": [("Htr1a",), ("Htr1b",), ("Htr1d",), ("Htr1f",),
                          ("Htr2a",), ("Htr2c",), ("Htr4",), ("Htr5b",), ("Htr7",)],
        },
        "Norepinephrine": {
            "ligands": [],  # Dbh/Slc6a2 not in panel
            "receptors": [("Adra1a",), ("Adra1b",), ("Adra2a",), ("Adra2b",),
                          ("Adra2c",), ("Adrb1",)],
        },
        "Acetylcholine (muscarinic)": {
            "ligands": [],  # Chat/Slc18a3 not in panel
            "receptors": [("Chrm1",), ("Chrm2",), ("Chrm3",), ("Chrm5",)],
        },
        "Histamine": {
            "ligands": [],  # Hdc not in panel
            "receptors": [("Hrh1",), ("Hrh2",), ("Hrh3",)],
        },
    }
    for _k, _v in neuromod_systems.items():
        systems[_k] = _v

    print(f"Loaded {len(systems)} systems ({len(systems) - len(neuromod_systems)} NP + {len(neuromod_systems)} neuromodulator)")
    return (systems,)


@app.cell
def _(DATA_DIR, pd):
    # Load whole-brain cell metadata and expression
    cells = pd.read_parquet(
        DATA_DIR / "processed" / "mouse_abc_whole_brain" / "cells_with_coords.parquet"
    )
    expr = pd.read_parquet(
        DATA_DIR / "processed" / "mouse_abc_whole_brain" / "neuropeptide_expression.parquet"
    )

    # Filter to neurons
    neuron_ids = cells.loc[cells["is_neuron"], "cell_id"]
    expr_neurons = expr.loc[expr.index.isin(neuron_ids)]
    cells_neurons = cells[cells["is_neuron"]].set_index("cell_id")

    print(f"Total cells: {len(cells):,}")
    print(f"Neurons: {len(expr_neurons):,}")
    print(f"Expression matrix: {expr_neurons.shape[0]:,} neurons × {expr_neurons.shape[1]} genes")
    return cells_neurons, expr_neurons


@app.cell
def _(BLACKLISTED_GENES, EXPRESSION_THRESHOLD, expr_neurons, np, systems):
    # For each neuron, count how many NP systems it expresses a ligand / receptor for
    # Blacklisted genes are excluded before counting
    _available_genes = set(expr_neurons.columns) - BLACKLISTED_GENES

    n_ligand_systems = np.zeros(len(expr_neurons), dtype=int)
    n_receptor_systems = np.zeros(len(expr_neurons), dtype=int)

    _expr_vals = expr_neurons.values
    _col_idx = {g: i for i, g in enumerate(expr_neurons.columns)}

    for _sys_name, _sys_info in systems.items():
        # Ligand: any non-blacklisted ligand gene above threshold
        _lig_cols = [_col_idx[g] for g in _sys_info["ligands"] if g in _available_genes]
        if _lig_cols:
            _max_lig = _expr_vals[:, _lig_cols].max(axis=1)
            n_ligand_systems += (_max_lig >= EXPRESSION_THRESHOLD).astype(int)

        # Receptor: skip any complex containing a blacklisted gene
        _any_complex_above = np.zeros(len(expr_neurons), dtype=bool)
        for _receptor_complex in _sys_info["receptors"]:
            if any(g in BLACKLISTED_GENES for g in _receptor_complex):
                continue
            _complex_gene_cols = [_col_idx[g] for g in _receptor_complex if g in _available_genes]
            if len(_complex_gene_cols) != len(_receptor_complex):
                continue
            _min_across_complex = _expr_vals[:, _complex_gene_cols].min(axis=1)
            _any_complex_above |= (_min_across_complex >= EXPRESSION_THRESHOLD)
        n_receptor_systems += _any_complex_above.astype(int)

    print(f"Blacklisted genes: {sorted(BLACKLISTED_GENES)}")
    print(f"Computed per-neuron NP system counts for {len(expr_neurons):,} neurons")
    print(f"Ligand systems per neuron: mean={n_ligand_systems.mean():.1f}, max={n_ligand_systems.max()}")
    print(f"Receptor systems per neuron: mean={n_receptor_systems.mean():.1f}, max={n_receptor_systems.max()}")
    return n_ligand_systems, n_receptor_systems


@app.cell
def _(cells_neurons, expr_neurons, n_ligand_systems, n_receptor_systems, pd):
    # Build analysis dataframe
    neuron_df = pd.DataFrame({
        "cell_id": expr_neurons.index,
        "n_ligand_systems": n_ligand_systems,
        "n_receptor_systems": n_receptor_systems,
    }).set_index("cell_id")

    neuron_df = neuron_df.join(cells_neurons[["region", "class", "subclass", "parcellation_division"]])

    # Define region subsets
    subsets = {
        "Whole brain": neuron_df,
        "Isocortex": neuron_df[neuron_df["parcellation_division"] == "Isocortex"],
        "Hippocampus": neuron_df[neuron_df["parcellation_division"] == "HPF"],
        "Hypothalamus": neuron_df[neuron_df["parcellation_division"] == "HY"],
        "ARH + ME + VMH": neuron_df[neuron_df["region"].isin({"ARH", "ME", "VMH"})],
    }

    for _name, _df in subsets.items():
        print(f"{_name}: {len(_df):,} neurons")
    return (subsets,)


@app.cell(hide_code=True)
def _(np, plt, subsets):
    _subset_keys = ["Whole brain", "Isocortex", "Hippocampus", "Hypothalamus", "ARH + ME + VMH"]
    _metrics = [("n_ligand_systems", "Ligand systems"), ("n_receptor_systems", "Receptor systems")]
    _colors = {
        "Whole brain": "#636EFA", "Isocortex": "#00CC96", "Hippocampus": "#FFA15A",
        "Hypothalamus": "#AB63FA", "ARH + ME + VMH": "#EF553B",
    }

    # Shared x-axis limit per row (same metric across all regions)
    _xlims = {}
    for _metric, _ in _metrics:
        _xmax = max(int(subsets[s][_metric].max()) for s in _subset_keys)
        _xlims[_metric] = (-0.8, _xmax + 0.8)

    _fig, _axes = plt.subplots(2, 5, figsize=(18, 5.5))

    for _col_i, _subset_name in enumerate(_subset_keys):
        _df = subsets[_subset_name]
        _color = _colors[_subset_name]

        for _row_i, (_metric, _label) in enumerate(_metrics):
            _ax = _axes[_row_i, _col_i]
            _vals = _df[_metric].values
            _max_val = int(_vals.max())
            _bins = np.arange(0, 21) - 0.5
            _counts, _ = np.histogram(_vals, bins=_bins)
            _pcts = 100 * _counts / _counts.sum()
            _x_vals = np.arange(0, int(_bins[-1] + .5))

            _ax.bar(_x_vals, _pcts, color=_color, width=0.8, edgecolor="white", linewidth=0.3)

            _ax.set_xlim(_xlims[_metric])
            _ax.set_xticks(np.arange(0, int(_xlims[_metric][1]) + 1, 2))

            if _col_i == 0:
                _ax.set_ylabel("% of neurons")
            if _row_i == 1:
                _ax.set_xlabel("# NP systems expressed")
            if _row_i == 0:
                _ax.set_title(f"{_subset_name}\n(n={len(_df):,})", fontsize=10)

            _ax.spines[["top", "right"]].set_visible(False)

    # Second pass: draw median/mean markers after y-limits settle
    for _col_i, _subset_name in enumerate(_subset_keys):
        _df = subsets[_subset_name]
        for _row_i, (_metric, _) in enumerate(_metrics):
            _ax = _axes[_row_i, _col_i]
            _vals = _df[_metric].values
            _ymax = _ax.get_ylim()[1]
            _median = int(np.median(_vals))
            _mean = _vals.mean()

            # Median: black dashed line + black triangle
            _ax.axvline(_median, color="k", linewidth=1, linestyle="--", zorder=3)
            _ax.plot(_median, _ymax * 0.95, marker="v", color="k", markersize=7,
                     zorder=4, clip_on=False)

            # Mean: gray solid line + gray triangle
            _ax.axvline(_mean, color="gray", linewidth=1, linestyle="-", zorder=3)
            _ax.plot(_mean, _ymax * 0.95, marker="v", color="gray", markersize=7,
                     zorder=4, clip_on=False)

    # Row labels on right side
    for _row_i, (_, _label) in enumerate(_metrics):
        _axes[_row_i, 4].annotate(_label, xy=(1.02, 0.5), xycoords="axes fraction",
                                   fontsize=11, fontweight="bold", rotation=-90,
                                   va="center", ha="left")

    _fig.suptitle(r"NP systems per neuron (threshold = 3.0 log$_2$CPM)"
                  "\nblack dashed = median, gray solid = mean", fontsize=12, y=1.02)
    _fig.tight_layout()
    plt.gca()
    return


@app.cell(hide_code=True)
def _(mo, np, pd, subsets):
    # Summary statistics table
    _rows = []
    for _subset_name, _df in subsets.items():
        for _metric, _label in [
            ("n_ligand_systems", "Ligand"),
            ("n_receptor_systems", "Receptor"),
        ]:
            _vals = _df[_metric].values
            _rows.append({
                "Region": _subset_name,
                "Type": _label,
                "N neurons": f"{len(_df):,}",
                "Mean": f"{_vals.mean():.2f}",
                "Median": f"{int(np.median(_vals))}",
                "Max": f"{int(_vals.max())}",
                "% with 0": f"{100 * (_vals == 0).mean():.1f}%",
                "% with ≥1": f"{100 * (_vals >= 1).mean():.1f}%",
                "% with ≥3": f"{100 * (_vals >= 3).mean():.1f}%",
                "% with ≥5": f"{100 * (_vals >= 5).mean():.1f}%",
                "% with ≥10": f"{100 * (_vals >= 10).mean():.1f}%",
            })

    summary_df = pd.DataFrame(_rows)
    mo.md("## Summary statistics")
    return (summary_df,)


@app.cell(hide_code=True)
def _(mo, summary_df):
    mo.ui.table(summary_df, selection=None)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ## Cumulative distributions (CCDFs)
    """)
    return


@app.cell(hide_code=True)
def _(np, plt, subsets):
    _colors2 = {
        "Whole brain": "#636EFA", "Isocortex": "#00CC96", "Hippocampus": "#FFA15A",
        "Hypothalamus": "#AB63FA", "ARH + ME + VMH": "#EF553B",
    }

    _fig2, _axes2 = plt.subplots(1, 2, figsize=(10, 4))

    for _col_i, (_metric, _label) in enumerate([
        ("n_ligand_systems", "NP Ligand Systems"),
        ("n_receptor_systems", "NP Receptor Systems"),
    ]):
        _ax = _axes2[_col_i]
        for _subset_name, _df in subsets.items():
            _vals = _df[_metric].values
            _max_val = int(_vals.max())
            _x = np.arange(0, _max_val + 1)
            _ccdf = np.array([100 * (_vals >= k).mean() for k in _x])

            _ax.plot(_x, _ccdf, "-o", markersize=3, color=_colors2[_subset_name],
                     label=_subset_name)

        _ax.set_xlabel("k (number of NP systems)")
        _ax.set_ylabel("% neurons with ≥ k systems")
        _ax.set_title(_label, fontsize=11)
        _ax.set_xticks(np.arange(0, _ax.get_xlim()[1] + 1, 2))
        _ax.spines[["top", "right"]].set_visible(False)
        if _col_i == 0:
            _ax.legend(fontsize=8)

    _fig2.suptitle("Complementary CDF: NP systems per neuron", fontsize=12, y=1.01)
    _fig2.tight_layout()
    plt.gca()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ## Per-gene prevalence: % of neurons expressing each ligand / receptor gene

    Sorted by whole-sample prevalence. Use this to identify suspiciously ubiquitous genes
    that may warrant blacklisting.
    """)
    return


@app.cell
def _(EXPRESSION_THRESHOLD, cells_neurons, expr_neurons, np, pd, systems):
    # Build per-gene prevalence table
    _available = set(expr_neurons.columns)
    _expr_vals = expr_neurons.values
    _col_idx = {g: i for i, g in enumerate(expr_neurons.columns)}

    _div_col = cells_neurons.loc[expr_neurons.index, "parcellation_division"]
    _region_col = cells_neurons.loc[expr_neurons.index, "region"]
    _masks = {
        "% Whole brain": np.ones(len(expr_neurons), dtype=bool),
        "% Isocortex": (_div_col == "Isocortex").values,
        "% Hippocampus": (_div_col == "HPF").values,
        "% Hypothalamus": (_div_col == "HY").values,
        "% ARH+ME+VMH": _region_col.isin({"ARH", "ME", "VMH"}).values,
    }

    # Collect all ligand and receptor genes with their system membership
    _ligand_gene_systems = {}
    _receptor_gene_systems = {}
    for _sys_name, _sys_info in systems.items():
        for _g in _sys_info["ligands"]:
            if _g in _available:
                _ligand_gene_systems.setdefault(_g, []).append(_sys_name)
        for _complex in _sys_info["receptors"]:
            for _g in _complex:
                if _g in _available:
                    _receptor_gene_systems.setdefault(_g, []).append(_sys_name)

    _gene_rows = []
    for _gene, _sys_list in sorted(_ligand_gene_systems.items()):
        _ci = _col_idx[_gene]
        _expr_col = _expr_vals[:, _ci]
        _row = {
            "Gene": _gene,
            "Type": "Ligand",
            "Systems": ", ".join(sorted(set(_sys_list))),
        }
        for _mname, _mask in _masks.items():
            _row[_mname] = round(100 * ((_expr_col >= EXPRESSION_THRESHOLD) & _mask).sum() / _mask.sum(), 1)
        _gene_rows.append(_row)
    for _gene, _sys_list in sorted(_receptor_gene_systems.items()):
        _ci = _col_idx[_gene]
        _expr_col = _expr_vals[:, _ci]
        _row = {
            "Gene": _gene,
            "Type": "Receptor",
            "Systems": ", ".join(sorted(set(_sys_list))),
        }
        for _mname, _mask in _masks.items():
            _row[_mname] = round(100 * ((_expr_col >= EXPRESSION_THRESHOLD) & _mask).sum() / _mask.sum(), 1)
        _gene_rows.append(_row)

    gene_prevalence_df = pd.DataFrame(_gene_rows)
    gene_prevalence_df
    return (gene_prevalence_df,)


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ### Ligand genes sorted by % neurons expressing (whole brain)
    """)
    return


@app.cell(hide_code=True)
def _(gene_prevalence_df, mo):
    _lig_table = (
        gene_prevalence_df[gene_prevalence_df["Type"] == "Ligand"]
        .sort_values("% Whole brain", ascending=False)
        .reset_index(drop=True)
    )
    mo.ui.table(_lig_table, selection=None)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ### Receptor genes sorted by % neurons expressing (whole brain)
    """)
    return


@app.cell(hide_code=True)
def _(gene_prevalence_df, mo):
    _rec_table = (
        gene_prevalence_df[gene_prevalence_df["Type"] == "Receptor"]
        .sort_values("% Whole brain", ascending=False)
        .reset_index(drop=True)
    )
    mo.ui.table(_rec_table, selection=None)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ## Top NP systems by prevalence
    """)
    return


@app.cell
def _(EXPRESSION_THRESHOLD, cells_neurons, expr_neurons, np, pd, systems):
    # For each NP system, what fraction of neurons express its ligand / receptor?
    _available = set(expr_neurons.columns)
    _expr_vals = expr_neurons.values
    _col_idx = {g: i for i, g in enumerate(expr_neurons.columns)}

    # Region masks
    _div_col = cells_neurons.loc[expr_neurons.index, "parcellation_division"]
    _region_col = cells_neurons.loc[expr_neurons.index, "region"]
    _region_masks = [
        (np.ones(len(expr_neurons), dtype=bool), "Whole brain"),
        ((_div_col == "Isocortex").values, "Isocortex"),
        ((_div_col == "HPF").values, "Hippocampus"),
        ((_div_col == "HY").values, "Hypothalamus"),
        (_region_col.isin({"ARH", "ME", "VMH"}).values, "ARH+ME+VMH"),
    ]

    _system_prevalence = []
    for _sys_name, _sys_info in systems.items():
        # Ligand expressing mask
        _lig_cols = [_col_idx[g] for g in _sys_info["ligands"] if g in _available]
        if _lig_cols:
            _lig_mask = _expr_vals[:, _lig_cols].max(axis=1) >= EXPRESSION_THRESHOLD
        else:
            _lig_mask = np.zeros(len(expr_neurons), dtype=bool)

        # Receptor expressing mask
        _rec_mask = np.zeros(len(expr_neurons), dtype=bool)
        for _receptor_complex in _sys_info["receptors"]:
            _cg = [_col_idx[g] for g in _receptor_complex if g in _available]
            if len(_cg) != len(_receptor_complex):
                continue
            _rec_mask |= (_expr_vals[:, _cg].min(axis=1) >= EXPRESSION_THRESHOLD)

        for _mask, _region_name in _region_masks:
            _n = _mask.sum()
            _system_prevalence.append({
                "System": _sys_name,
                "Region": _region_name,
                "Ligand genes": ", ".join(sorted(set(_sys_info["ligands"]) & _available)),
                "% expressing ligand": 100 * (_lig_mask & _mask).sum() / _n if _n > 0 else 0,
                "% expressing receptor": 100 * (_rec_mask & _mask).sum() / _n if _n > 0 else 0,
            })

    prevalence_df = pd.DataFrame(_system_prevalence)
    prevalence_df
    return (prevalence_df,)


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ### Top 15 ligand systems by % neurons expressing
    """)
    return


@app.cell(hide_code=True)
def _(mo, prevalence_df):
    _top_lig = (
        prevalence_df
        .pivot(index="System", columns="Region", values="% expressing ligand")
        .sort_values("Whole brain", ascending=False)
        .head(15)
        .round(1)
        .reset_index()
    )
    mo.ui.table(_top_lig, selection=None)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md("""
    ### Top 15 receptor systems by % neurons expressing
    """)
    return


@app.cell(hide_code=True)
def _(mo, prevalence_df):
    _top_rec = (
        prevalence_df
        .pivot(index="System", columns="Region", values="% expressing receptor")
        .sort_values("Whole brain", ascending=False)
        .head(15)
        .round(1)
        .reset_index()
    )
    mo.ui.table(_top_rec, selection=None)
    return


@app.cell
def _(mo):
    mo.md("""
    ## Neuromodulator receptor prevalence by region

    Classical neuromodulator receptor genes (dopamine, serotonin, norepinephrine,
    acetylcholine, histamine) — % of neurons expressing each gene above threshold.
    """)
    return


@app.cell
def _(EXPRESSION_THRESHOLD, cells_neurons, expr_neurons, np, pd):
    _neuromod_genes = {
        "Dopamine": ["Drd1", "Drd3", "Drd5"],
        "Serotonin": ["Htr1a", "Htr1b", "Htr1d", "Htr1f", "Htr2a", "Htr2c", "Htr4", "Htr5b", "Htr7"],
        "Norepinephrine": ["Adra1a", "Adra1b", "Adra2a", "Adra2b", "Adra2c", "Adrb1"],
        "Acetylcholine": ["Chrm1", "Chrm2", "Chrm3", "Chrm5"],
        "Histamine": ["Hrh1", "Hrh2", "Hrh3"],
    }

    _expr_vals = expr_neurons.values
    _col_idx = {g: i for i, g in enumerate(expr_neurons.columns)}
    _div_col = cells_neurons.loc[expr_neurons.index, "parcellation_division"]
    _region_col = cells_neurons.loc[expr_neurons.index, "region"]
    _masks = {
        "% Whole brain": np.ones(len(expr_neurons), dtype=bool),
        "% Isocortex": (_div_col == "Isocortex").values,
        "% Hippocampus": (_div_col == "HPF").values,
        "% Hypothalamus": (_div_col == "HY").values,
        "% ARH+ME+VMH": _region_col.isin({"ARH", "ME", "VMH"}).values,
    }

    _rows = []
    for _system, _genes in _neuromod_genes.items():
        for _gene in _genes:
            if _gene not in _col_idx:
                continue
            _ci = _col_idx[_gene]
            _expr_col = _expr_vals[:, _ci]
            _row = {"System": _system, "Gene": _gene}
            for _mname, _mask in _masks.items():
                _row[_mname] = round(100 * ((_expr_col >= EXPRESSION_THRESHOLD) & _mask).sum() / _mask.sum(), 1)
            _rows.append(_row)

    neuromod_prevalence_df = pd.DataFrame(_rows)
    neuromod_prevalence_df
    return (neuromod_prevalence_df,)


@app.cell
def _(mo, neuromod_prevalence_df):
    mo.ui.table(neuromod_prevalence_df, selection=None)
    return


if __name__ == "__main__":
    app.run()
