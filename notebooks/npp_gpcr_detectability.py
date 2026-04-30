import marimo

__generated_with = "0.19.4"
app = marimo.App()


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Neuropeptide Precursor & GPCR Transcript Detectability

    How detectable are neuropeptide precursor (NPP) and GPCR transcripts in MERFISH data?
    We compare at three spatial scales: **whole brain**, **hypothalamus**, and **ARH millimeter cube**.

    **Prior expectation**: GPCRs should be harder to detect than precursors because
    receptors operate at high sensitivity (few copies needed) while precursors are produced
    in bulk for secretion.

    **Actual finding — a double dissociation**:

    - **Receptors are more broadly expressed**: detected in 2-3x more cells than ligands
      (median % expressing ~1-3% vs ~0.3-0.7% across scales).
    - **Ligands are expressed at higher levels when present**: median expression among
      expressing cells is ~5.4 log2 for ligands vs ~4.3 log2 for receptors.

    This makes biological sense: neuropeptide precursors are produced in bulk by
    *specialized* cell populations (few cells, high expression), while receptors are
    expressed at modest levels across *many* target cells to enable sensitivity.
    The constraint for spatial mapping is therefore different for each side of the
    circuit: for **ligands**, the challenge is spatial sparsity; for **receptors**,
    the challenge is low per-cell expression levels.
    """)
    return


@app.cell
def _(Path, mo):
    _img_path = Path(mo.notebook_dir()) / "assets" / "langlieb2023_np_receptor_histogram.png"
    mo.vstack([
        mo.md(r"""
        ### Prior work: Langlieb et al. (2023)

        Langlieb et al. curated **65 NP precursor genes** with known downstream GPCRs
        (Supplementary Table 8). Their Fig. 3c shows the number of NP precursors and
        NP receptors expressed per cluster across ~5,000 transcriptomic clusters (snRNA-seq).
        The receptor distribution is shifted right — most clusters express more receptor
        types than precursor types — consistent with our MERFISH finding that receptors
        are more broadly detected.

        The Allen ABC Atlas defines a broader functional gene set (`abc_functional_genes.csv`)
        with **49 NP precursors**, **22 additional secreted peptides** (cerebellins,
        neurexophilins, natriuretic peptides, chromogranins, etc.), and **209 GPCRs**
        (59 NP receptors + 40 NT GPCRs + 110 other GPCRs). Our expression file contains
        a subset of these — filtered to genes in `np_map.csv` + `hormone_map.csv`.
        """),
        mo.image(_img_path, width=300),
        mo.md("*Source: Langlieb et al. (2023), Fig. 3c.*"),
    ])
    return


@app.cell
def _():
    import marimo as mo
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    from pathlib import Path
    return Path, mo, np, pd, plt


@app.cell
def _(Path, mo, np, pd):
    _root = Path(mo.notebook_dir()) / ".."

    # --- Gene classification from ABC functional gene list ---
    # Source: data/raw/mouse_abc/abc_functional_genes.csv
    # Sections: ion channels | NP receptors | NT GPCRs | other GPCRs | NP precursors | extra secreted
    _func_lines = (_root / "data/raw/mouse_abc/abc_functional_genes.csv").read_text().splitlines()
    _func_genes = [l.split(",")[0] for l in _func_lines if l.strip() and l.split(",")[0]]

    _np_rec_start = _func_genes.index("Adcyap1r1")
    _nt_gpcr_start = _func_genes.index("Adora1")
    _other_gpcr_start = _func_genes.index("Ackr1")
    _np_prec_start = _func_genes.index("Adcyap1")
    _extra_start = _func_genes.index("Retn")

    _abc_np_precursors = set(_func_genes[_np_prec_start:_extra_start])
    _abc_extra_secreted = set(_func_genes[_extra_start:])
    _abc_np_receptors = set(_func_genes[_np_rec_start:_nt_gpcr_start])
    _abc_nt_gpcrs = set(_func_genes[_nt_gpcr_start:_other_gpcr_start])
    _abc_other_gpcrs = set(_func_genes[_other_gpcr_start:_np_prec_start])
    _abc_all_ligands = _abc_np_precursors | _abc_extra_secreted
    _abc_all_receptors = _abc_np_receptors | _abc_nt_gpcrs | _abc_other_gpcrs

    # --- MERFISH 500-gene panel (directly measured) ---
    # Source: data/raw/mouse_abc/500_gene_panel.csv
    merfish_500 = set(pd.read_csv(_root / "data/raw/mouse_abc/500_gene_panel.csv")["vizgen_gene"])

    # --- Load expression data and cell metadata ---
    # neuropeptide_expression.parquet contains all genes from np_map + hormone_map + abc_functional_genes
    # Regenerate via: uv run snakemake build_cluster_ligand_receptor_map_extended --cores all
    expr = pd.read_parquet(_root / "data/processed/mouse_abc_extended/neuropeptide_expression.parquet")
    cells = pd.read_parquet(_root / "data/processed/mouse_abc_extended/cells_with_coords.parquet")

    # Cast float16 → float32 to avoid overflow in mean/sum on 1.3M rows
    expr = expr.astype(np.float32)

    # Align indices: expression is indexed by cell_id, cells has cell_id as column
    cells = cells.set_index("cell_id")

    # --- Classify genes present in expression data ---
    _expr_genes = set(expr.columns) - {"__index_level_0__"}
    ligand_genes = _abc_all_ligands & _expr_genes
    receptor_genes = _abc_all_receptors & _expr_genes
    both_genes = ligand_genes & receptor_genes

    gene_category = {}
    for g in ligand_genes - both_genes:
        gene_category[g] = "ligand"
    for g in receptor_genes - both_genes:
        gene_category[g] = "receptor"
    for g in both_genes:
        gene_category[g] = "both"

    # Track which genes are directly measured vs imputed
    gene_source = {g: ("measured" if g in merfish_500 else "imputed") for g in gene_category}

    _lig_measured = sum(1 for g in ligand_genes if g in merfish_500)
    _rec_measured = sum(1 for g in receptor_genes if g in merfish_500)

    mo.md(f"""
    **Data loaded**: {expr.shape[0]:,} cells × {expr.shape[1]} genes in expression file

    Gene classification from `abc_functional_genes.csv`:

    | Category | In ABC list | In expression file | Measured (MERFISH 500) | Imputed |
    |----------|------------:|-------------------:|-----------------------:|--------:|
    | NP precursors | {len(_abc_np_precursors)} | {len(_abc_np_precursors & _expr_genes)} | {len(_abc_np_precursors & merfish_500 & _expr_genes)} | {len((_abc_np_precursors & _expr_genes) - merfish_500)} |
    | Extra secreted | {len(_abc_extra_secreted)} | {len(_abc_extra_secreted & _expr_genes)} | {len(_abc_extra_secreted & merfish_500 & _expr_genes)} | {len((_abc_extra_secreted & _expr_genes) - merfish_500)} |
    | NP receptors | {len(_abc_np_receptors)} | {len(_abc_np_receptors & _expr_genes)} | {len(_abc_np_receptors & merfish_500 & _expr_genes)} | {len((_abc_np_receptors & _expr_genes) - merfish_500)} |
    | NT GPCRs | {len(_abc_nt_gpcrs)} | {len(_abc_nt_gpcrs & _expr_genes)} | — | — |
    | Other GPCRs | {len(_abc_other_gpcrs)} | {len(_abc_other_gpcrs & _expr_genes)} | — | — |

    - **Ligands** (NP precursors + extra secreted): **{len(ligand_genes)}** in expression file ({_lig_measured} measured, {len(ligand_genes) - _lig_measured} imputed)
    - **Receptors** (NP + NT + other GPCRs): **{len(receptor_genes)}** in expression file ({_rec_measured} measured, {len(receptor_genes) - _rec_measured} imputed)

    **Key caveat**: Only {_lig_measured + _rec_measured} of {len(gene_category)} genes are directly
    measured by MERFISH. The rest come from imputation (MapMyCells). Key hypothalamic
    NPs like Npy, Pomc, Agrp, Oxt, Avp, Cck, Sst, Hcrt, and Tac1 are all **imputed only**.
    """)
    return cells, expr, gene_category, gene_source


@app.cell
def _(cells, mo):
    # Define spatial subsets
    _whole_brain_idx = cells.index
    _hypo_idx = cells[cells["parcellation_division"] == "HY"].index
    _arh_idx = cells[
        (cells["x"] >= -0.1) & (cells["x"] <= 0.8) &
        (cells["y"] >= 1.2) & (cells["y"] <= 2.05) &
        (cells["z"] >= -2.5) & (cells["z"] <= -1.55)
    ].index

    spatial_subsets = {
        "Whole brain": _whole_brain_idx,
        "Hypothalamus": _hypo_idx,
        "ARH bbox": _arh_idx,
    }

    mo.md(f"""
    **Spatial subsets**:

    | Scale | Cells |
    |-------|------:|
    | Whole brain | {len(_whole_brain_idx):,} |
    | Hypothalamus (HY division) | {len(_hypo_idx):,} |
    | ARH bbox | {len(_arh_idx):,} |
    """)
    return (spatial_subsets,)


@app.cell
def _(expr, gene_category, gene_source, np, pd, spatial_subsets):
    # Compute per-gene stats for each spatial scale
    THRESHOLD = 3.0
    _rows = []

    for scale_name, idx in spatial_subsets.items():
        _sub = expr.loc[idx]
        _n_cells = len(idx)
        for gene, cat in gene_category.items():
            _vals = _sub[gene].dropna()
            _expressing = _vals[_vals > THRESHOLD]
            _n_expressing = len(_expressing)
            _pct = 100.0 * _n_expressing / _n_cells if _n_cells > 0 else 0.0
            _mean_all = _vals.mean() if len(_vals) > 0 else np.nan
            _mean_expr = _expressing.mean() if _n_expressing > 0 else np.nan

            _rows.append({
                "gene": gene,
                "scale": scale_name,
                "category": cat,
                "source": gene_source[gene],
                "n_cells": _n_cells,
                "n_expressing": _n_expressing,
                "pct_expressing": _pct,
                "mean_expr_all": _mean_all,
                "mean_expr_expressing": _mean_expr,
            })

    stats_df = pd.DataFrame(_rows)
    stats_df
    return (stats_df,)


@app.cell
def _(mo, np, stats_df):
    # Summary table — broken out by category AND measured vs imputed
    _summary = (
        stats_df[stats_df["category"].isin(["ligand", "receptor"])]
        .groupby(["category", "source", "scale"])
        .agg(
            n_genes=("gene", "nunique"),
            median_pct_expressing=("pct_expressing", "median"),
            mean_pct_expressing=("pct_expressing", "mean"),
            median_mean_expr_expressing=("mean_expr_expressing", lambda x: np.nanmedian(x)),
        )
        .reset_index()
    )

    _lines = [
        "| Category | Source | Scale | Genes | Median % expressing | Mean % expressing | Median expr (expressing) |",
        "|----------|--------|-------|------:|--------------------:|------------------:|-------------------------:|",
    ]
    for _, _r in _summary.iterrows():
        _lines.append(
            f"| {_r['category']} | {_r['source']} | {_r['scale']} | {_r['n_genes']} | "
            f"{_r['median_pct_expressing']:.2f}% | {_r['mean_pct_expressing']:.2f}% | "
            f"{_r['median_mean_expr_expressing']:.2f} |"
        )

    mo.md("## Summary: Ligand vs Receptor Detectability\n\n" + "\n".join(_lines))
    return


@app.cell
def _(mo, np, stats_df):
    # Interpretation of the double dissociation
    _lig_wb = stats_df[(stats_df["scale"] == "Whole brain") & (stats_df["category"] == "ligand")]
    _rec_wb = stats_df[(stats_df["scale"] == "Whole brain") & (stats_df["category"] == "receptor")]
    _lig_arh = stats_df[(stats_df["scale"] == "ARH bbox") & (stats_df["category"] == "ligand")]
    _rec_arh = stats_df[(stats_df["scale"] == "ARH bbox") & (stats_df["category"] == "receptor")]

    _lig_arh_detected = _lig_arh["mean_expr_expressing"].dropna()
    _rec_arh_detected = _rec_arh["mean_expr_expressing"].dropna()

    mo.md(f"""
    ### Interpretation: Double Dissociation

    The data reveal a clear **double dissociation** between ligands and receptors:

    | Metric | Ligands | Receptors | Ratio |
    |--------|--------:|----------:|------:|
    | Median % expressing (whole brain) | {_lig_wb['pct_expressing'].median():.2f}% | {_rec_wb['pct_expressing'].median():.2f}% | {_rec_wb['pct_expressing'].median() / _lig_wb['pct_expressing'].median():.1f}x more receptor cells |
    | Median % expressing (ARH bbox) | {_lig_arh['pct_expressing'].median():.2f}% | {_rec_arh['pct_expressing'].median():.2f}% | {_rec_arh['pct_expressing'].median() / _lig_arh['pct_expressing'].median():.1f}x more receptor cells |
    | Median expr when expressing (whole brain) | {np.nanmedian(_lig_wb['mean_expr_expressing']):.2f} | {np.nanmedian(_rec_wb['mean_expr_expressing']):.2f} | {np.nanmedian(_lig_wb['mean_expr_expressing']) / np.nanmedian(_rec_wb['mean_expr_expressing']):.1f}x higher ligand expression |
    | Genes with any expressing cells (ARH) | {len(_lig_arh_detected)}/{len(_lig_arh)} | {len(_rec_arh_detected)}/{len(_rec_arh)} | |

    **Why this matters for spatial mapping**:

    - **Ligand source mapping** is limited by *spatial sparsity* — precursors are produced by
      a small number of specialized cells, so you need sufficient spatial sampling to capture them.
      In the ARH bbox, {len(_lig_arh) - len(_lig_arh_detected)} of {len(_lig_arh)} ligand genes
      have zero expressing cells above threshold.
    - **Receptor target mapping** is limited by *low per-cell abundance* — receptors are present
      in many cells but at modest levels, making them harder to distinguish from imputation noise.
    - The implication is that **ligand-receptor circuit mapping from MERFISH is asymmetric**:
      we can confidently locate receptor-expressing target populations (broad expression),
      but ligand sources require either larger spatial volumes or lower detection thresholds
      (at the cost of more false positives from imputation artifacts).
    """)
    return


@app.cell
def _(plt, stats_df):
    # Histogram: pct_expressing by category
    _cats = ["ligand", "receptor"]
    _scales = ["Whole brain", "Hypothalamus", "ARH bbox"]

    fig1, axes1 = plt.subplots(1, 3, figsize=(14, 4), sharey=True)
    for _i, _scale in enumerate(_scales):
        _ax = axes1[_i]
        for _cat, _color in zip(_cats, ["#2196F3", "#FF5722"]):
            _vals = stats_df[
                (stats_df["scale"] == _scale) & (stats_df["category"] == _cat)
            ]["pct_expressing"]
            _ax.hist(_vals, bins=20, alpha=0.6, label=_cat.capitalize(), color=_color)
        _ax.set_title(_scale)
        _ax.set_xlabel("% cells expressing")
        if _i == 0:
            _ax.set_ylabel("Number of genes")
        _ax.legend()

    fig1.suptitle("Detectability: % cells expressing (threshold > 3.0 log2)", y=1.02)
    fig1.tight_layout()
    fig1
    return


@app.cell
def _(np, plt, stats_df):
    # Histogram: mean expression among expressing cells
    _cats = ["ligand", "receptor"]
    _scales = ["Whole brain", "Hypothalamus", "ARH bbox"]

    fig2, axes2 = plt.subplots(1, 3, figsize=(14, 4), sharey=True)
    for _i, _scale in enumerate(_scales):
        _ax = axes2[_i]
        for _cat, _color in zip(_cats, ["#2196F3", "#FF5722"]):
            _vals = stats_df[
                (stats_df["scale"] == _scale) & (stats_df["category"] == _cat)
            ]["mean_expr_expressing"]
            _vals = _vals[np.isfinite(_vals)]
            if len(_vals) > 0:
                _ax.hist(_vals, bins=20, alpha=0.6, label=_cat.capitalize(), color=_color)
        _ax.set_title(_scale)
        _ax.set_xlabel("Mean expression (expressing cells, log2)")
        if _i == 0:
            _ax.set_ylabel("Number of genes")
        _ax.legend()

    fig2.suptitle("Expression Level Among Expressing Cells", y=1.02)
    fig2.tight_layout()
    fig2
    return


@app.cell
def _(plt, stats_df):
    # Scatter: pct_expressing vs mean expression, with measured vs imputed distinction
    _cats = ["ligand", "receptor"]
    _scales = ["Whole brain", "Hypothalamus", "ARH bbox"]

    fig3, axes3 = plt.subplots(1, 3, figsize=(14, 4))
    for _i, _scale in enumerate(_scales):
        _ax = axes3[_i]
        for _cat, _color in zip(_cats, ["#2196F3", "#FF5722"]):
            for _src, _marker, _alpha in [("measured", "o", 0.9), ("imputed", "x", 0.4)]:
                _sub = stats_df[
                    (stats_df["scale"] == _scale)
                    & (stats_df["category"] == _cat)
                    & (stats_df["source"] == _src)
                ]
                if len(_sub) == 0:
                    continue
                _ax.scatter(
                    _sub["pct_expressing"], _sub["mean_expr_expressing"],
                    alpha=_alpha, label=f"{_cat} ({_src})", color=_color,
                    marker=_marker, s=40 if _src == "measured" else 25,
                )
        _ax.set_title(_scale)
        _ax.set_xlabel("% cells expressing")
        _ax.set_ylabel("Mean expr (expressing)")
        _ax.legend()

    fig3.suptitle("Detectability vs Expression Level", y=1.02)
    fig3.tight_layout()
    fig3
    return


@app.cell
def _(mo):
    # Per-gene bar chart for ARH bbox
    pct_threshold = mo.ui.slider(
        0, 50, value=0, step=1, label="Min % expressing to show"
    )
    pct_threshold
    return (pct_threshold,)


@app.cell
def _(pct_threshold, plt, stats_df):
    _arh = stats_df[stats_df["scale"] == "ARH bbox"].copy()
    _arh = _arh[_arh["category"].isin(["ligand", "receptor"])]
    _arh = _arh[_arh["pct_expressing"] >= pct_threshold.value]
    _arh = _arh.sort_values("pct_expressing", ascending=True)

    # Color by category, hatch pattern for imputed
    _colors = _arh["category"].map({"ligand": "#2196F3", "receptor": "#FF5722"}).values
    _hatches = _arh["source"].map({"measured": "", "imputed": "///"}).values
    _height = max(6, len(_arh) * 0.22)

    fig4, ax4 = plt.subplots(figsize=(8, _height))
    _bars = ax4.barh(range(len(_arh)), _arh["pct_expressing"], color=_colors)
    for _bar, _hatch in zip(_bars, _hatches):
        _bar.set_hatch(_hatch)
    ax4.set_yticks(range(len(_arh)))
    ax4.set_yticklabels(
        [f"{g} {'*' if s == 'imputed' else ''}"
         for g, s in zip(_arh["gene"], _arh["source"])],
        fontsize=8,
    )
    ax4.set_xlabel("% cells expressing (ARH bbox)")
    ax4.set_title("Per-Gene Detectability in ARH Millimeter Cube (* = imputed)")

    # Legend
    from matplotlib.patches import Patch
    ax4.legend(
        handles=[
            Patch(color="#2196F3", label="Ligand (measured)"),
            Patch(color="#2196F3", hatch="///", label="Ligand (imputed)"),
            Patch(color="#FF5722", label="Receptor (measured)"),
            Patch(color="#FF5722", hatch="///", label="Receptor (imputed)"),
        ],
        loc="lower right",
        fontsize=7,
    )
    fig4.tight_layout()
    fig4
    return


@app.cell
def _(mo, stats_df):
    # Interactive gene table for ARH bbox
    _arh_table = (
        stats_df[stats_df["scale"] == "ARH bbox"]
        .sort_values("pct_expressing", ascending=False)
        [["gene", "category", "source", "n_expressing", "pct_expressing", "mean_expr_all", "mean_expr_expressing"]]
        .round(3)
    )
    mo.ui.table(_arh_table, label="ARH bbox gene stats")
    return


@app.cell
def _(Path, expr, mo, pd):
    # Missing genes note
    _root = Path(mo.notebook_dir()) / ".."
    _missing_path = _root / "data/processed/mouse_abc/neuropeptide_genes_missing.txt"
    _missing_lines = _missing_path.read_text().strip().split("\n")
    _missing_genes = [l.strip() for l in _missing_lines if l.strip() and not l.startswith("#")]

    # Re-parse ABC functional gene list for the coverage table
    _func_lines = (_root / "data/raw/mouse_abc/abc_functional_genes.csv").read_text().splitlines()
    _func_genes = [l.split(",")[0] for l in _func_lines if l.strip() and l.split(",")[0]]
    _abc_np_precursors = set(_func_genes[_func_genes.index("Adcyap1"):_func_genes.index("Retn")])
    _abc_extra_secreted = set(_func_genes[_func_genes.index("Retn"):])
    _abc_np_receptors = set(_func_genes[_func_genes.index("Adcyap1r1"):_func_genes.index("Adora1")])
    _abc_nt_gpcrs = set(_func_genes[_func_genes.index("Adora1"):_func_genes.index("Ackr1")])
    _abc_other_gpcrs = set(_func_genes[_func_genes.index("Ackr1"):_func_genes.index("Adcyap1")])
    _expr_genes = set(expr.columns) - {"__index_level_0__"}
    _merfish_500 = set(pd.read_csv(_root / "data/raw/mouse_abc/500_gene_panel.csv")["vizgen_gene"])

    mo.md(f"""
    ## Missing Genes (not in MERFISH panel)

    **{len(_missing_genes)} genes** from the neuropeptide map are not present in the MERFISH 500-gene panel
    and would require imputed (MapMyCells) data for spatial expression analysis:

    {', '.join(f'`{g}`' for g in sorted(_missing_genes))}

    These are predominantly receptors and less-common ligands. Their absence further
    limits the ability to map complete ligand-receptor circuits from MERFISH alone.

    ## ABC functional gene panel — what's NOT in our expression file

    The Allen ABC Atlas defines a comprehensive functional gene set
    (`abc_functional_genes.csv`) that is much broader than what our expression
    file covers:

    | Category | In ABC list | In our expression file | Not included |
    |----------|------------:|-----------------------:|-------------:|
    | NP precursors | 49 | {len(_abc_np_precursors & _expr_genes)} | {len(_abc_np_precursors - _expr_genes)} |
    | Extra secreted (Cbln, Nxph, Npp, Scg, etc.) | 22 | {len(_abc_extra_secreted & _expr_genes)} | {len(_abc_extra_secreted - _expr_genes)} |
    | NP receptors | 59 | {len(_abc_np_receptors & _expr_genes)} | {len(_abc_np_receptors - _expr_genes)} |
    | NT GPCRs (Adora, Chrm, Drd, Grm, Htr, etc.) | 40 | {len(_abc_nt_gpcrs & _expr_genes)} | {len(_abc_nt_gpcrs - _expr_genes)} |
    | Other GPCRs (orphan, adhesion, Fzd, etc.) | 110 | {len(_abc_other_gpcrs & _expr_genes)} | {len(_abc_other_gpcrs - _expr_genes)} |
    | Ion channels | 261 | 0 | 261 |

    The extra secreted peptides not in our file include: {', '.join(f'`{g}`' for g in sorted(_abc_extra_secreted - _expr_genes))}.

    Many of these (especially cerebellins, neurexophilins, and orphan GPCRs) are
    **directly measured** in the MERFISH 500-gene panel and could extend the analysis
    of paracrine signaling beyond classical neuropeptide-GPCR pairs.
    """)
    return


if __name__ == "__main__":
    app.run()
