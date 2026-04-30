import marimo

__generated_with = "0.19.4"
app = marimo.App(width="full")


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Validation: Langlieb Ligand/Receptor Coverage

    Investigates why neuropeptide coverage in the Langlieb whole-brain snRNA-seq data
    is hit-and-miss compared to the Allen Brain Cell (ABC) MERFISH data.

    **Key differences between datasets:**
    - **ABC**: ~385k cells, 550-gene MERFISH panel + imputed full genome, hypothalamus-focused
    - **Langlieb**: ~4.4M nuclei snRNA-seq, cluster-level averages only (no per-cell data), whole-brain

    **Known issues from Langlieb et al. 2023:**
    - OXT, AVP, PMCH, AGRP show contamination across clusters (ambient RNA)
    - Paper uses dual threshold: ≥30% expressing AND avg ≥0.5 (≥80% + ≥5 for contaminated NPs)
    """)
    return


@app.cell
def _():
    import pandas as pd
    import numpy as np
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    import marimo as mo
    return go, make_subplots, mo, np, pd


@app.cell
def _(pd):
    # Load data
    lang_profile = pd.read_parquet(
        "../data/processed/mouse_langlieb/cluster_ligand_receptor_profile.parquet"
    )
    abc_profile = pd.read_parquet(
        "../data/processed/mouse_abc/cluster_ligand_receptor_profile.parquet"
    )
    lang_np_expr = pd.read_parquet(
        "../data/processed/mouse_langlieb/cluster_np_expression.parquet"
    )
    abc_np_expr = pd.read_parquet(
        "../data/processed/mouse_abc/cluster_np_expression.parquet"
    )
    np_map = pd.read_csv("../data/processed/mouse_common/np_map.csv")

    print(f"Langlieb profile: {lang_profile.shape[0]:,} cluster-gene pairs, {lang_profile['cluster'].nunique()} clusters")
    print(f"ABC profile: {abc_profile.shape[0]:,} cluster-gene pairs, {abc_profile['cluster'].nunique()} clusters")
    print(f"NP systems in np_map: {np_map['System'].nunique()}")
    return abc_np_expr, abc_profile, lang_np_expr, lang_profile, np_map


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    ## 1. Gene Availability Check

    Are all expected genes present in the Langlieb snRNA-seq data?
    """)
    return


@app.cell
def _(lang_profile, np_map, pd):
    # Check gene availability
    _all_ligands = set()
    _all_receptors = set()
    for _g in np_map["Ligand_Gene"].dropna().unique():
        for _x in _g.split(";"):
            _all_ligands.add(_x.strip())
    for _g in np_map["Receptor_Gene"].dropna().unique():
        for _x in _g.split(";"):
            _all_receptors.add(_x.strip())

    _all_expected = _all_ligands | _all_receptors
    _found = set(lang_profile["gene"].unique())
    _missing = _all_expected - _found

    availability = pd.DataFrame(
        {
            "Category": ["Ligand genes", "Receptor genes", "Total unique", "Found in Langlieb", "Missing"],
            "Count": [
                len(_all_ligands),
                len(_all_receptors),
                len(_all_expected),
                len(_all_expected & _found),
                len(_missing),
            ],
        }
    )
    print(availability.to_string(index=False))
    if _missing:
        print(f"\nMissing genes: {sorted(_missing)}")
    else:
        print("\nAll genes found - naming is not the issue.")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    ## 2. Expression Level Comparison: Langlieb vs ABC

    For each NP system, compare the maximum ligand and receptor expression values
    between the two datasets. This reveals the fundamental dynamic range difference.
    """)
    return


@app.cell
def _(abc_np_expr, lang_np_expr, np, pd):
    # Build per-system comparison
    _systems = sorted(lang_np_expr["system"].unique())
    _rows = []
    for _sys in _systems:
        _lang = lang_np_expr[lang_np_expr["system"] == _sys]
        _abc = abc_np_expr[abc_np_expr["system"] == _sys]

        _rows.append(
            {
                "system": _sys,
                "lang_lig_max": _lang["max_ligand_expr"].max(),
                "lang_lig_mean": _lang["max_ligand_expr"].mean(),
                "lang_lig_gt01": (_lang["max_ligand_expr"] > 0.1).sum(),
                "lang_rec_max": _lang["max_receptor_expr"].max(),
                "lang_rec_mean": _lang["max_receptor_expr"].mean(),
                "lang_rec_gt01": (_lang["max_receptor_expr"] > 0.1).sum(),
                "abc_lig_max": _abc["max_ligand_expr"].max() if len(_abc) > 0 else np.nan,
                "abc_lig_gt01": (_abc["max_ligand_expr"] > 0.1).sum() if len(_abc) > 0 else 0,
                "abc_rec_max": _abc["max_receptor_expr"].max() if len(_abc) > 0 else np.nan,
                "abc_rec_gt01": (_abc["max_receptor_expr"] > 0.1).sum() if len(_abc) > 0 else 0,
                "lang_n_clusters": len(_lang),
                "abc_n_clusters": len(_abc) if len(_abc) > 0 else 0,
            }
        )

    system_comparison = pd.DataFrame(_rows)
    return (system_comparison,)


@app.cell
def _(go, make_subplots, system_comparison):
    # Scatter: Langlieb max vs ABC max for ligands and receptors
    _df = system_comparison.dropna(subset=["abc_lig_max", "abc_rec_max"])

    _fig = make_subplots(
        rows=1,
        cols=2,
        subplot_titles=["Ligand max expression", "Receptor max expression"],
    )

    # Ligands
    _fig.add_trace(
        go.Scatter(
            x=_df["abc_lig_max"],
            y=_df["lang_lig_max"],
            mode="markers+text",
            text=_df["system"],
            textposition="top center",
            textfont=dict(size=8),
            marker=dict(size=8, color="crimson"),
            name="Ligands",
            hovertext=[
                f"{r['system']}<br>ABC: {r['abc_lig_max']:.2f}<br>Langlieb: {r['lang_lig_max']:.2f}"
                for _, r in _df.iterrows()
            ],
            hoverinfo="text",
        ),
        row=1,
        col=1,
    )

    # Receptors
    _fig.add_trace(
        go.Scatter(
            x=_df["abc_rec_max"],
            y=_df["lang_rec_max"],
            mode="markers+text",
            text=_df["system"],
            textposition="top center",
            textfont=dict(size=8),
            marker=dict(size=8, color="steelblue"),
            name="Receptors",
            hovertext=[
                f"{r['system']}<br>ABC: {r['abc_rec_max']:.2f}<br>Langlieb: {r['lang_rec_max']:.2f}"
                for _, r in _df.iterrows()
            ],
            hoverinfo="text",
        ),
        row=1,
        col=2,
    )

    # Add diagonal
    _max_val = max(
        _df["abc_lig_max"].max(),
        _df["lang_lig_max"].max(),
        _df["abc_rec_max"].max(),
        _df["lang_rec_max"].max(),
    )
    for _col in [1, 2]:
        _fig.add_trace(
            go.Scatter(
                x=[0, _max_val],
                y=[0, _max_val],
                mode="lines",
                line=dict(dash="dash", color="gray"),
                showlegend=False,
            ),
            row=1,
            col=_col,
        )

    _fig.update_xaxes(title_text="ABC max expression", row=1, col=1)
    _fig.update_xaxes(title_text="ABC max expression", row=1, col=2)
    _fig.update_yaxes(title_text="Langlieb max expression", row=1, col=1)
    _fig.update_yaxes(title_text="Langlieb max expression", row=1, col=2)
    _fig.update_layout(
        title="Max Expression Comparison: Langlieb vs ABC<br><sup>Points below diagonal = lower Langlieb expression</sup>",
        height=600,
        width=1200,
        showlegend=False,
    )
    _fig
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    ## 3. Per-System Coverage Grade

    Grade each NP system on whether it has usable signal in the Langlieb data.

    **Criteria:**
    - Ligand "present": ≥5 clusters with max_ligand_expr > 0.1
    - Receptor "present": ≥5 clusters with max_receptor_expr > 0.1
    - **A**: Both ligand and receptor well-represented (≥50 clusters each >0.1)
    - **B**: Both present but sparse (≥5 clusters each >0.1)
    - **C**: Only one of ligand/receptor present
    - **F**: Neither usable
    """)
    return


@app.cell
def _(lang_np_expr, pd):
    # Grade each system
    _systems = sorted(lang_np_expr["system"].unique())
    _grade_rows = []

    for _sys in _systems:
        _data = lang_np_expr[lang_np_expr["system"] == _sys]

        _lig_gt01 = (_data["max_ligand_expr"] > 0.1).sum()
        _rec_gt01 = (_data["max_receptor_expr"] > 0.1).sum()
        _lig_gt05 = (_data["max_ligand_expr"] > 0.5).sum()
        _rec_gt05 = (_data["max_receptor_expr"] > 0.5).sum()
        _lig_max = _data["max_ligand_expr"].max()
        _rec_max = _data["max_receptor_expr"].max()

        _lig_present = _lig_gt01 >= 5
        _rec_present = _rec_gt01 >= 5

        if _lig_gt01 >= 50 and _rec_gt01 >= 50:
            _grade = "A"
        elif _lig_present and _rec_present:
            _grade = "B"
        elif _lig_present or _rec_present:
            _grade = "C"
        else:
            _grade = "F"

        _which = []
        if _lig_present:
            _which.append("L")
        if _rec_present:
            _which.append("R")

        _grade_rows.append(
            {
                "System": _sys,
                "Grade": _grade,
                "Has": "+".join(_which) if _which else "none",
                "Lig>0.1": _lig_gt01,
                "Lig>0.5": _lig_gt05,
                "Lig_max": round(_lig_max, 2),
                "Rec>0.1": _rec_gt01,
                "Rec>0.5": _rec_gt05,
                "Rec_max": round(_rec_max, 2),
            }
        )

    grade_df = pd.DataFrame(_grade_rows).sort_values(
        ["Grade", "System"]
    )

    # Print summary
    _grade_counts = grade_df["Grade"].value_counts().sort_index()
    print("=== GRADE DISTRIBUTION ===")
    for _g, _c in _grade_counts.items():
        print(f"  Grade {_g}: {_c} systems")

    print(f"\n=== SYSTEMS WITH ISSUES (Grade C or F) ===")
    _issues = grade_df[grade_df["Grade"].isin(["C", "F"])]
    for _, _r in _issues.iterrows():
        print(
            f"  {_r['System']}: Grade {_r['Grade']} ({_r['Has']}) "
            f"- Lig max={_r['Lig_max']}, Rec max={_r['Rec_max']}"
        )
    return (grade_df,)


@app.cell
def _(grade_df, mo):
    mo.ui.table(grade_df, label="Langlieb NP System Coverage Grades")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    ## 4. Contamination Analysis

    The Langlieb paper explicitly flags OXT, AVP, PMCH, and AGRP as contamination-prone.
    Here we quantify the "specificity" of each ligand: what fraction of clusters show
    expression above noise? A highly specific ligand is expressed in few clusters
    (e.g. Hcrt in orexin neurons), while a contaminated one appears everywhere.
    """)
    return


@app.cell
def _(lang_profile, np_map, pd):
    # Compute specificity metrics for each ligand gene
    _ligand_genes = set()
    for _g in np_map["Ligand_Gene"].dropna().unique():
        for _x in _g.split(";"):
            _ligand_genes.add(_x.strip())

    _spec_rows = []
    for _gene in sorted(_ligand_genes):
        _data = lang_profile[lang_profile["gene"] == _gene]
        if len(_data) == 0:
            continue

        _total = len(_data)
        _gt0 = (_data["mean_expr"] > 0).sum()
        _gt01 = (_data["mean_expr"] > 0.1).sum()
        _gt05 = (_data["mean_expr"] > 0.5).sum()
        _gt1 = (_data["mean_expr"] > 1.0).sum()
        _max_expr = _data["mean_expr"].max()
        _mean_expr = _data["mean_expr"].mean()
        _max_pct = _data["pct_expressing"].max()

        # Specificity: fraction of clusters NOT expressing (higher = more specific)
        _specificity = 1.0 - _gt01 / _total if _total > 0 else 0

        # Get system name
        _sys = np_map[np_map["Ligand_Gene"].str.contains(_gene, na=False)][
            "System"
        ].iloc[0] if len(np_map[np_map["Ligand_Gene"].str.contains(_gene, na=False)]) > 0 else "?"

        # Flag contaminated genes (from Langlieb paper)
        _contaminated = _gene in ["Oxt", "Avp", "Pmch", "Agrp"]

        _spec_rows.append(
            {
                "gene": _gene,
                "system": _sys,
                "specificity": round(_specificity, 3),
                "clusters>0.1": _gt01,
                "clusters>0.5": _gt05,
                "clusters>1.0": _gt1,
                "max_expr": round(_max_expr, 3),
                "mean_expr": round(_mean_expr, 4),
                "max_pct_expressing": round(_max_pct, 1),
                "contaminated": _contaminated,
            }
        )

    specificity_df = pd.DataFrame(_spec_rows).sort_values("specificity")
    return (specificity_df,)


@app.cell
def _(go, specificity_df):
    # Bar chart of specificity
    _df = specificity_df.sort_values("specificity")
    _colors = [
        "red" if c else "steelblue" for c in _df["contaminated"]
    ]

    _fig = go.Figure(
        data=go.Bar(
            x=_df["gene"],
            y=_df["specificity"],
            marker_color=_colors,
            hovertext=[
                f"{r['gene']} ({r['system']})<br>"
                f"Specificity: {r['specificity']:.3f}<br>"
                f"Clusters>0.1: {r['clusters>0.1']}<br>"
                f"Max expr: {r['max_expr']:.3f}<br>"
                f"{'CONTAMINATED (flagged by Langlieb)' if r['contaminated'] else ''}"
                for _, r in _df.iterrows()
            ],
            hoverinfo="text",
        )
    )
    _fig.update_layout(
        title="Ligand Gene Specificity in Langlieb snRNA-seq<br>"
        "<sup>Red = flagged as contaminated by Langlieb paper. Higher = more specific expression.</sup>",
        xaxis_title="Ligand Gene",
        yaxis_title="Specificity (1 - fraction of clusters with expr > 0.1)",
        height=500,
        width=1200,
        xaxis=dict(tickangle=45, tickfont=dict(size=8)),
    )
    _fig
    return


@app.cell
def _(specificity_df):
    # Show the worst offenders
    _low_spec = specificity_df[specificity_df["specificity"] < 0.5].sort_values(
        "specificity"
    )
    print("=== LOW SPECIFICITY LIGANDS (expressed in >50% of clusters) ===")
    print("These are likely contaminated or ubiquitously expressed:\n")
    for _, _r in _low_spec.iterrows():
        _flag = " *** KNOWN CONTAMINATED ***" if _r["contaminated"] else ""
        print(
            f"  {_r['gene']} ({_r['system']}): specificity={_r['specificity']:.3f}, "
            f"clusters>0.1={_r['clusters>0.1']}, mean_expr={_r['mean_expr']:.4f}{_flag}"
        )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    ## 5. Receptor Coverage Deep Dive

    For NP systems where the receptor signal is weak, diagnose whether it's a
    heterodimer AND-logic issue (components present individually but not together)
    or genuinely low expression.
    """)
    return


@app.cell
def _(lang_profile, np_map, pd):
    # Check heterodimer receptor complexes
    _rows = []
    for _, _map_row in np_map.drop_duplicates("Receptor_Gene").iterrows():
        _receptor_str = _map_row["Receptor_Gene"]
        _system = _map_row["System"]
        if pd.isna(_receptor_str):
            continue
        _genes = [g.strip() for g in _receptor_str.split(";")]
        _is_heterodimer = len(_genes) > 1

        _gene_exprs = {}
        for _g in _genes:
            _data = lang_profile[lang_profile["gene"] == _g]
            if len(_data) > 0:
                _gene_exprs[_g] = {
                    "max_expr": _data["mean_expr"].max(),
                    "clusters_gt01": (_data["mean_expr"] > 0.1).sum(),
                    "mean_expr": _data["mean_expr"].mean(),
                }
            else:
                _gene_exprs[_g] = {"max_expr": 0, "clusters_gt01": 0, "mean_expr": 0}

        # For heterodimers, effective expression = min of components
        if _is_heterodimer:
            _eff_max = min(v["max_expr"] for v in _gene_exprs.values())
            _component_detail = " & ".join(
                f"{g}(max={v['max_expr']:.2f}, n={v['clusters_gt01']})"
                for g, v in _gene_exprs.items()
            )
        else:
            _g = _genes[0]
            _eff_max = _gene_exprs[_g]["max_expr"]
            _component_detail = f"{_g}(max={_gene_exprs[_g]['max_expr']:.2f}, n={_gene_exprs[_g]['clusters_gt01']})"

        _rows.append(
            {
                "system": _system,
                "receptor": _receptor_str,
                "is_heterodimer": _is_heterodimer,
                "effective_max": round(_eff_max, 3),
                "components": _component_detail,
            }
        )

    receptor_detail_df = pd.DataFrame(_rows).sort_values(
        ["is_heterodimer", "effective_max"], ascending=[False, True]
    )
    return (receptor_detail_df,)


@app.cell
def _(receptor_detail_df):
    # Show heterodimer receptors and their component expression
    _hetero = receptor_detail_df[receptor_detail_df["is_heterodimer"]]
    print(f"=== HETERODIMER RECEPTORS ({len(_hetero)} receptor complexes) ===\n")
    for _, _r in _hetero.iterrows():
        print(f"  {_r['system']}: {_r['receptor']}")
        print(f"    Effective max: {_r['effective_max']:.3f}")
        print(f"    Components: {_r['components']}")
        print()

    # Show weakest non-heterodimer receptors
    _mono = receptor_detail_df[
        (~receptor_detail_df["is_heterodimer"])
        & (receptor_detail_df["effective_max"] < 0.5)
    ]
    print(f"=== WEAK SINGLE-GENE RECEPTORS (max < 0.5, {len(_mono)} receptors) ===\n")
    for _, _r in _mono.iterrows():
        print(f"  {_r['system']}: {_r['receptor']} - max={_r['effective_max']:.3f}")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    ## 6. Side-by-Side System Profiles

    Compare expression distributions for selected systems between ABC and Langlieb.
    """)
    return


@app.cell
def _(mo):
    FOCUS_SYSTEMS = [
        "Orexin",
        "AgRP",
        "Adrenomedullin",
        "Oxytocin",
        "Pomc",
        "Neuropeptide Y",
        "CRH",
        "Vasopressin",
        "Kisspeptin",
        "MCH",
        "TRH",
        "Somatostatin",
    ]
    system_selector = mo.ui.dropdown(
        options={s: s for s in FOCUS_SYSTEMS},
        value="Orexin",
        label="Select NP system",
    )
    mo.md(f"**Select system for detailed comparison:** {system_selector}")
    return (system_selector,)


@app.cell(hide_code=True)
def _(abc_profile, go, lang_profile, make_subplots, np_map, system_selector):
    _sys = system_selector.value
    _sys_map = np_map[np_map["System"] == _sys]

    # Get ligand and receptor genes
    _lig_genes = set()
    _rec_genes = set()
    for _g in _sys_map["Ligand_Gene"].dropna().unique():
        for _x in _g.split(";"):
            _lig_genes.add(_x.strip())
    for _g in _sys_map["Receptor_Gene"].dropna().unique():
        for _x in _g.split(";"):
            _rec_genes.add(_x.strip())

    _all_genes = sorted(_lig_genes | _rec_genes)

    # Build comparison data
    _fig = make_subplots(
        rows=len(_all_genes),
        cols=2,
        subplot_titles=[
            f"{g} - Langlieb" if i % 2 == 0 else f"{g} - ABC"
            for g in _all_genes
            for i in range(2)
        ],
        vertical_spacing=0.05,
    )

    for _i, _gene in enumerate(_all_genes):
        _row = _i + 1

        # Langlieb
        _lang_data = lang_profile[lang_profile["gene"] == _gene]
        if len(_lang_data) > 0:
            _lang_vals = _lang_data["mean_expr"].values
            _lang_vals = _lang_vals[_lang_vals > 0]  # Remove zeros for histogram
            if len(_lang_vals) > 0:
                _fig.add_trace(
                    go.Histogram(
                        x=_lang_vals,
                        nbinsx=50,
                        marker_color="crimson" if _gene in _lig_genes else "steelblue",
                        name=f"{_gene} (Lang)",
                        showlegend=False,
                    ),
                    row=_row,
                    col=1,
                )

        # ABC
        _abc_data = abc_profile[abc_profile["gene"] == _gene]
        if len(_abc_data) > 0:
            _abc_vals = _abc_data["mean_expr"].values
            _abc_vals = _abc_vals[_abc_vals > 0]
            if len(_abc_vals) > 0:
                _fig.add_trace(
                    go.Histogram(
                        x=_abc_vals,
                        nbinsx=50,
                        marker_color="crimson" if _gene in _lig_genes else "steelblue",
                        name=f"{_gene} (ABC)",
                        showlegend=False,
                    ),
                    row=_row,
                    col=2,
                )

    _fig.update_layout(
        title=f"{_sys} System: Gene Expression Distributions<br>"
        "<sup>Red = ligand, Blue = receptor. Histogram of mean_expr across clusters (zeros excluded).</sup>",
        height=250 * len(_all_genes),
        width=1000,
    )
    _fig
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    ## 7. Langlieb Paper Threshold Analysis

    The paper uses dual thresholds for NP assignment:
    - **Standard NPs**: ≥30% expressing AND avg ≥0.5
    - **Contaminated NPs** (OXT, AVP, PMCH, AGRP): ≥80% expressing AND avg ≥5.0

    How many clusters pass these thresholds for each ligand?
    """)
    return


@app.cell
def _(lang_profile, np_map, pd):
    # Apply Langlieb paper thresholds
    _contaminated_genes = {"Oxt", "Avp", "Pmch", "Agrp"}

    _ligand_genes_map = {}
    for _, _r in np_map.drop_duplicates("Ligand_Gene").iterrows():
        _lg = _r["Ligand_Gene"]
        if pd.isna(_lg):
            continue
        for _g in _lg.split(";"):
            _g = _g.strip()
            _ligand_genes_map[_g] = _r["System"]

    _thresh_rows = []
    for _gene, _sys in sorted(_ligand_genes_map.items()):
        _data = lang_profile[lang_profile["gene"] == _gene]
        if len(_data) == 0:
            continue

        _is_contam = _gene in _contaminated_genes
        _pct_thresh = 80 if _is_contam else 30
        _expr_thresh = 5.0 if _is_contam else 0.5

        # Count clusters passing each threshold individually and both
        _pass_pct = (_data["pct_expressing"] >= _pct_thresh).sum()
        _pass_expr = (_data["mean_expr"] >= _expr_thresh).sum()
        _pass_both = (
            (_data["pct_expressing"] >= _pct_thresh)
            & (_data["mean_expr"] >= _expr_thresh)
        ).sum()

        _thresh_rows.append(
            {
                "gene": _gene,
                "system": _sys,
                "contaminated": _is_contam,
                "pct_threshold": _pct_thresh,
                "expr_threshold": _expr_thresh,
                "pass_pct": _pass_pct,
                "pass_expr": _pass_expr,
                "pass_both": _pass_both,
                "total_clusters": len(_data),
            }
        )

    threshold_df = pd.DataFrame(_thresh_rows).sort_values("pass_both", ascending=False)
    return (threshold_df,)


@app.cell
def _(go, threshold_df):
    # Bar chart of clusters passing Langlieb thresholds
    _df = threshold_df.sort_values("pass_both", ascending=True)
    _colors = [
        "red" if c else "steelblue" for c in _df["contaminated"]
    ]

    _fig = go.Figure(
        data=go.Bar(
            y=_df["gene"],
            x=_df["pass_both"],
            orientation="h",
            marker_color=_colors,
            hovertext=[
                f"{r['gene']} ({r['system']})<br>"
                f"Pass both: {r['pass_both']}<br>"
                f"Pass pct≥{r['pct_threshold']}%: {r['pass_pct']}<br>"
                f"Pass expr≥{r['expr_threshold']}: {r['pass_expr']}<br>"
                f"{'CONTAMINATED' if r['contaminated'] else ''}"
                for _, r in _df.iterrows()
            ],
            hoverinfo="text",
        )
    )
    _fig.update_layout(
        title="Clusters Passing Langlieb Paper Thresholds<br>"
        "<sup>Red = contaminated NPs (higher thresholds). Shows clusters where ligand is reliably detected.</sup>",
        xaxis_title="Number of clusters passing both thresholds",
        yaxis_title="Ligand gene",
        height=max(500, len(_df) * 18),
        width=900,
    )
    _fig
    return


@app.cell
def _(threshold_df):
    # Show genes with 0 clusters passing
    _zero = threshold_df[threshold_df["pass_both"] == 0]
    _few = threshold_df[
        (threshold_df["pass_both"] > 0) & (threshold_df["pass_both"] <= 5)
    ]

    print("=== LIGANDS WITH NO CLUSTERS PASSING LANGLIEB THRESHOLDS ===\n")
    for _, _r in _zero.iterrows():
        print(
            f"  {_r['gene']} ({_r['system']}): "
            f"pct≥{_r['pct_threshold']}%: {_r['pass_pct']}, "
            f"expr≥{_r['expr_threshold']}: {_r['pass_expr']}"
        )

    print(f"\n=== LIGANDS WITH ≤5 CLUSTERS PASSING ===\n")
    for _, _r in _few.iterrows():
        print(
            f"  {_r['gene']} ({_r['system']}): {_r['pass_both']} clusters"
        )
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    ## 8. Pct-Expressing vs Mean Expression (Diagnostic Scatter)

    For each cluster-gene pair, plot pct_expressing vs mean_expr. This helps identify
    genes with ambient RNA contamination (high mean but low pct_expressing) vs
    genuine expression (high both).
    """)
    return


@app.cell
def _(mo):
    DIAG_GENES = [
        "Agrp",
        "Hcrt",
        "Oxt",
        "Avp",
        "Pomc",
        "Npy",
        "Crh",
        "Adm",
        "Pmch",
        "Kiss1",
    ]
    diag_gene_selector = mo.ui.dropdown(
        options={g: g for g in DIAG_GENES},
        value="Agrp",
        label="Gene",
    )
    mo.md(f"**Select gene for diagnostic plot:** {diag_gene_selector}")
    return (diag_gene_selector,)


@app.cell
def _(diag_gene_selector, go, lang_profile):
    _gene = diag_gene_selector.value
    _data = lang_profile[lang_profile["gene"] == _gene].copy()

    _is_contam = _gene in ["Oxt", "Avp", "Pmch", "Agrp"]
    _pct_thresh = 80 if _is_contam else 30
    _expr_thresh = 5.0 if _is_contam else 0.5

    _fig = go.Figure()

    # All clusters
    _fig.add_trace(
        go.Scatter(
            x=_data["mean_expr"],
            y=_data["pct_expressing"],
            mode="markers",
            marker=dict(size=3, color="lightgray", opacity=0.5),
            name="All clusters",
            hovertext=[
                f"{r['cluster']}<br>mean_expr={r['mean_expr']:.3f}<br>pct={r['pct_expressing']:.1f}%"
                for _, r in _data.iterrows()
            ],
            hoverinfo="text",
        )
    )

    # Passing clusters
    _passing = _data[
        (_data["pct_expressing"] >= _pct_thresh)
        & (_data["mean_expr"] >= _expr_thresh)
    ]
    if len(_passing) > 0:
        _fig.add_trace(
            go.Scatter(
                x=_passing["mean_expr"],
                y=_passing["pct_expressing"],
                mode="markers",
                marker=dict(size=5, color="crimson"),
                name=f"Pass thresholds ({len(_passing)})",
                hovertext=[
                    f"{r['cluster']}<br>mean_expr={r['mean_expr']:.3f}<br>pct={r['pct_expressing']:.1f}%"
                    for _, r in _passing.iterrows()
                ],
                hoverinfo="text",
            )
        )

    # Threshold lines
    _fig.add_hline(y=_pct_thresh, line_dash="dash", line_color="gray",
                   annotation_text=f"pct≥{_pct_thresh}%")
    _fig.add_vline(x=_expr_thresh, line_dash="dash", line_color="gray",
                   annotation_text=f"expr≥{_expr_thresh}")

    _fig.update_layout(
        title=f"{_gene}: pct_expressing vs mean_expr across {len(_data)} clusters<br>"
        f"<sup>{'CONTAMINATED - higher thresholds' if _is_contam else 'Standard thresholds'}</sup>",
        xaxis_title="Mean expression (log2)",
        yaxis_title="% cells expressing",
        height=500,
        width=800,
    )
    _fig
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    ## 9. Summary & Recommendations

    Based on the analysis above, here are the key findings and recommended actions.
    """)
    return


@app.cell
def _(grade_df, mo, specificity_df, threshold_df):
    # Merge all metrics for final summary
    _n_grade_a = (grade_df["Grade"] == "A").sum()
    _n_grade_b = (grade_df["Grade"] == "B").sum()
    _n_grade_c = (grade_df["Grade"] == "C").sum()
    _n_grade_f = (grade_df["Grade"] == "F").sum()
    _n_total = len(grade_df)

    _low_spec = specificity_df[specificity_df["specificity"] < 0.5]
    _no_pass = threshold_df[threshold_df["pass_both"] == 0]

    mo.md(f"""
    ### Key Findings

    **1. Gene naming is NOT the issue.** All {len(specificity_df)} ligand genes from np_map.csv are
    found in the Langlieb snRNA-seq data with correct capitalization (mouse gene symbols).

    **2. Expression dynamic range is much lower.** Langlieb snRNA-seq cluster averages have
    ~3-5x lower max expression than ABC imputed MERFISH for most genes. This is inherent to
    the data type (snRNA-seq averages vs spatially-imputed full-genome expression).

    **3. System coverage grades:**
    - Grade A (both L+R well-represented): **{_n_grade_a}/{_n_total}** systems
    - Grade B (both present, sparse): **{_n_grade_b}/{_n_total}** systems
    - Grade C (only L or R): **{_n_grade_c}/{_n_total}** systems
    - Grade F (neither usable): **{_n_grade_f}/{_n_total}** systems

    **4. Contamination:** {len(_low_spec)} ligand genes have specificity < 0.5 (expressed
    in >50% of clusters). The Langlieb paper flags OXT, AVP, PMCH, AGRP explicitly.

    **5. {len(_no_pass)} ligand genes** have zero clusters passing the paper's own thresholds.

    ### Recommendations

    1. **For the app**: Consider using `pct_expressing` alongside `mean_expr` as a dual
       threshold (matching the paper's approach) rather than just mean expression.

    2. **Contaminated NPs** (Agrp, Oxt, Avp, Pmch): Apply the paper's stricter thresholds
       (≥80% expressing AND avg ≥5.0) to filter noise.

    3. **Heterodimer receptors** (Calcrl;Ramp2 etc.): The AND-logic correctly handles these
       but reduces effective signal. Individual components often have decent expression -
       the low combined score is expected biology (co-expression requirement).

    4. **Low-expression ligands** (Hcrt, Adm2, many others): These are genuinely rare cell
       populations. The snRNA-seq cluster averages dilute their signal. For these systems,
       only the top 1-5 clusters will show meaningful expression.

    5. **Whole-brain perspective**: Unlike ABC (hypothalamus-focused), Langlieb covers the
       entire brain. Many NP systems that are hypothalamus-specific will appear as a tiny
       fraction of the 5,030 whole-brain clusters.
    """)
    return


@app.cell
def _(mo):
    mo.md(r"""
    ---
    ## 10. Per-Gene Scale Factor: Q-Q Analysis

    The Langlieb and ABC expression profiles look similar but scaled down. Here we
    quantify the per-gene relationship: for each common gene, we compute a Q-Q plot
    of cluster-level `mean_expr` distributions and fit `langlieb = scale * abc`.

    After normalizing each gene's distribution by its 95th percentile, the shape
    match should be nearly perfect if the difference is purely multiplicative.
    """)
    return


@app.cell
def _(abc_profile, lang_profile, np, np_map, pd):
    # Identify all ligand/receptor genes from np_map
    _all_genes = set()
    _ligand_genes_set = set()
    _receptor_genes_set = set()
    for _g in np_map["Ligand_Gene"].dropna().unique():
        for _x in _g.split(";"):
            _x = _x.strip()
            _all_genes.add(_x)
            _ligand_genes_set.add(_x)
    for _g in np_map["Receptor_Gene"].dropna().unique():
        for _x in _g.split(";"):
            _x = _x.strip()
            _all_genes.add(_x)
            _receptor_genes_set.add(_x)

    # Common genes present in both datasets
    _common = sorted(
        _all_genes
        & set(lang_profile["gene"].unique())
        & set(abc_profile["gene"].unique())
    )

    _quantile_points = np.linspace(0.01, 0.99, 50)
    _rows = []

    for _gene in _common:
        # Extract mean_expr values, filter inf/nan, cast to float64
        _lang_vals = lang_profile.loc[
            lang_profile["gene"] == _gene, "mean_expr"
        ].values.astype(np.float64)
        _lang_vals = _lang_vals[np.isfinite(_lang_vals)]

        _abc_vals = abc_profile.loc[
            abc_profile["gene"] == _gene, "mean_expr"
        ].values.astype(np.float64)
        _abc_vals = _abc_vals[np.isfinite(_abc_vals)]

        if len(_lang_vals) < 10 or len(_abc_vals) < 10:
            continue

        # 95th percentiles
        _lang_95 = np.percentile(_lang_vals, 95)
        _abc_95 = np.percentile(_abc_vals, 95)

        # Skip genes with no signal
        if _lang_95 < 0.01 or _abc_95 < 0.01:
            continue

        _scale_factor = _lang_95 / _abc_95

        # Compute quantiles
        _lang_q = np.quantile(_lang_vals, _quantile_points)
        _abc_q = np.quantile(_abc_vals, _quantile_points)

        # Normalize by 95th percentile
        _lang_q_norm = _lang_q / _lang_95
        _abc_q_norm = _abc_q / _abc_95

        # Fit Q-Q slope on normalized data (intercept forced to 0)
        # Mask points where both values are > 0.01 (skip near-zero noise)
        _mask = (_lang_q_norm > 0.01) & (_abc_q_norm > 0.01)
        if _mask.sum() < 5:
            continue

        _x = _abc_q_norm[_mask]
        _y = _lang_q_norm[_mask]

        # Slope = sum(x*y) / sum(x*x) for zero-intercept regression
        _norm_slope = np.sum(_x * _y) / np.sum(_x * _x)

        # R² for zero-intercept model
        _y_pred = _norm_slope * _x
        _ss_res = np.sum((_y - _y_pred) ** 2)
        _ss_tot = np.sum(_y ** 2)  # for zero-intercept, SS_tot = sum(y^2)
        _norm_r2 = 1.0 - _ss_res / _ss_tot if _ss_tot > 0 else 0.0

        # Zero fractions
        _lang_zero_frac = np.mean(_lang_vals == 0)
        _abc_zero_frac = np.mean(_abc_vals == 0)

        _is_ligand = _gene in _ligand_genes_set
        _is_receptor = _gene in _receptor_genes_set

        _rows.append({
            "gene": _gene,
            "scale_factor": round(_scale_factor, 4),
            "norm_slope": round(_norm_slope, 4),
            "norm_r2": round(_norm_r2, 4),
            "abc_95": round(_abc_95, 4),
            "lang_95": round(_lang_95, 4),
            "lang_zero_frac": round(_lang_zero_frac, 3),
            "abc_zero_frac": round(_abc_zero_frac, 3),
            "is_ligand": _is_ligand,
            "is_receptor": _is_receptor,
        })

    scale_df = pd.DataFrame(_rows).sort_values("scale_factor")

    print(f"Analyzed {len(scale_df)} common genes")
    print(f"Median scale factor (lang_95 / abc_95): {scale_df['scale_factor'].median():.3f}")
    print(f"Scale factor IQR: {scale_df['scale_factor'].quantile(0.25):.3f} - {scale_df['scale_factor'].quantile(0.75):.3f}")
    print(f"Median normalized Q-Q slope: {scale_df['norm_slope'].median():.3f}")
    print(f"Median normalized R²: {scale_df['norm_r2'].median():.3f}")
    return (scale_df,)


@app.cell
def _(mo):
    mo.md(r"""
    ---
    ## 11. Scale Factor Distribution

    Horizontal bar chart of all genes sorted by their scale factor (`langlieb_95th / abc_95th`).
    Color indicates how well the shape matches after rescaling (normalized R²):
    - **Green**: R² > 0.9 (good shape match)
    - **Yellow**: R² 0.7–0.9 (moderate)
    - **Red**: R² < 0.7 (shape mismatch — contamination or dropout)
    """)
    return


@app.cell
def _(go, scale_df):
    _df = scale_df.sort_values("scale_factor", ascending=True)

    # Color by norm_r2 quality
    _colors = []
    for _r2 in _df["norm_r2"]:
        if _r2 > 0.9:
            _colors.append("green")
        elif _r2 > 0.7:
            _colors.append("goldenrod")
        else:
            _colors.append("red")

    # Annotate ligand/receptor
    _labels = []
    for _, _r in _df.iterrows():
        _tag = ""
        if _r["is_ligand"] and _r["is_receptor"]:
            _tag = " [L+R]"
        elif _r["is_ligand"]:
            _tag = " [L]"
        elif _r["is_receptor"]:
            _tag = " [R]"
        _labels.append(f"{_r['gene']}{_tag}")

    _fig = go.Figure(
        data=go.Bar(
            y=_labels,
            x=_df["scale_factor"],
            orientation="h",
            marker_color=_colors,
            hovertext=[
                f"{r['gene']}<br>"
                f"Scale factor: {r['scale_factor']:.3f}<br>"
                f"Norm slope: {r['norm_slope']:.3f}<br>"
                f"Norm R²: {r['norm_r2']:.3f}<br>"
                f"ABC 95th: {r['abc_95']:.3f}<br>"
                f"Lang 95th: {r['lang_95']:.3f}<br>"
                f"{'Ligand' if r['is_ligand'] else ''}{'/' if r['is_ligand'] and r['is_receptor'] else ''}{'Receptor' if r['is_receptor'] else ''}"
                for _, r in _df.iterrows()
            ],
            hoverinfo="text",
        )
    )
    _fig.update_layout(
        title="Per-Gene Scale Factor (Langlieb 95th / ABC 95th)<br>"
        "<sup>Green = good shape match (R²>0.9), Yellow = moderate (0.7-0.9), Red = poor (<0.7). [L]=ligand, [R]=receptor</sup>",
        xaxis_title="Scale factor (langlieb_95th / abc_95th)",
        yaxis_title="Gene",
        height=max(600, len(_df) * 18),
        width=900,
    )
    _fig
    return


@app.cell
def _(mo):
    mo.md(r"""
    ---
    ## 12. Per-Gene Q-Q Explorer

    Select a gene to see:
    - **Left**: Raw Q-Q plot (ABC quantiles vs Langlieb quantiles) with fitted `y = scale * x`
    - **Right**: Normalized Q-Q (after dividing by 95th percentile) with identity line
    """)
    return


@app.cell
def _(mo, scale_df):
    _gene_options = {g: g for g in sorted(scale_df["gene"].unique())}
    qq_gene_selector = mo.ui.dropdown(
        options=_gene_options,
        value=scale_df.iloc[len(scale_df) // 2]["gene"],  # start at median
        label="Gene",
    )
    mo.md(f"**Select gene for Q-Q comparison:** {qq_gene_selector}")
    return (qq_gene_selector,)


@app.cell
def _(
    abc_profile,
    go,
    lang_profile,
    make_subplots,
    mo,
    np,
    qq_gene_selector,
    scale_df,
):
    _gene = qq_gene_selector.value
    _row = scale_df[scale_df["gene"] == _gene].iloc[0]

    # Re-extract values for this gene
    _lang_vals = lang_profile.loc[
        lang_profile["gene"] == _gene, "mean_expr"
    ].values.astype(np.float64)
    _lang_vals = _lang_vals[np.isfinite(_lang_vals)]
    _abc_vals = abc_profile.loc[
        abc_profile["gene"] == _gene, "mean_expr"
    ].values.astype(np.float64)
    _abc_vals = _abc_vals[np.isfinite(_abc_vals)]

    _quantile_points = np.linspace(0.01, 0.99, 50)
    _lang_q = np.quantile(_lang_vals, _quantile_points)
    _abc_q = np.quantile(_abc_vals, _quantile_points)

    _lang_95 = _row["lang_95"]
    _abc_95 = _row["abc_95"]
    _scale = _row["scale_factor"]

    _lang_q_norm = _lang_q / _lang_95
    _abc_q_norm = _abc_q / _abc_95

    _fig = make_subplots(
        rows=1, cols=2,
        subplot_titles=[
            f"Raw Q-Q (scale={_scale:.3f})",
            f"Normalized Q-Q (slope={_row['norm_slope']:.3f}, R²={_row['norm_r2']:.3f})",
        ],
    )

    # Left: Raw Q-Q
    _fig.add_trace(
        go.Scatter(
            x=_abc_q, y=_lang_q,
            mode="markers",
            marker=dict(size=6, color="steelblue"),
            name="Quantiles",
            hovertext=[
                f"p={p:.2f}<br>ABC={a:.3f}<br>Lang={l:.3f}"
                for p, a, l in zip(_quantile_points, _abc_q, _lang_q)
            ],
            hoverinfo="text",
        ),
        row=1, col=1,
    )
    # Fitted line y = scale * x
    _x_max_raw = max(_abc_q.max(), 0.01)
    _fig.add_trace(
        go.Scatter(
            x=[0, _x_max_raw], y=[0, _scale * _x_max_raw],
            mode="lines",
            line=dict(dash="dash", color="crimson"),
            name=f"y = {_scale:.3f}x",
        ),
        row=1, col=1,
    )

    # Right: Normalized Q-Q
    _fig.add_trace(
        go.Scatter(
            x=_abc_q_norm, y=_lang_q_norm,
            mode="markers",
            marker=dict(size=6, color="green"),
            name="Norm quantiles",
            hovertext=[
                f"p={p:.2f}<br>ABC_norm={a:.3f}<br>Lang_norm={l:.3f}"
                for p, a, l in zip(_quantile_points, _abc_q_norm, _lang_q_norm)
            ],
            hoverinfo="text",
        ),
        row=1, col=2,
    )
    # Identity line
    _x_max_norm = max(_abc_q_norm.max(), _lang_q_norm.max(), 1.1)
    _fig.add_trace(
        go.Scatter(
            x=[0, _x_max_norm], y=[0, _x_max_norm],
            mode="lines",
            line=dict(dash="dash", color="gray"),
            name="y = x",
        ),
        row=1, col=2,
    )

    _fig.update_xaxes(title_text="ABC quantile", row=1, col=1)
    _fig.update_yaxes(title_text="Langlieb quantile", row=1, col=1)
    _fig.update_xaxes(title_text="ABC quantile (normalized)", row=1, col=2)
    _fig.update_yaxes(title_text="Langlieb quantile (normalized)", row=1, col=2)
    _fig.update_layout(height=500, width=1100, showlegend=True)

    # Stats text
    _stats = mo.md(f"""
    **{_gene}** | Scale factor: **{_scale:.3f}** | Norm slope: **{_row['norm_slope']:.3f}** | Norm R²: **{_row['norm_r2']:.3f}**

    ABC 95th: {_abc_95:.3f} | Lang 95th: {_lang_95:.3f} | ABC zero fraction: {_row['abc_zero_frac']:.1%} | Lang zero fraction: {_row['lang_zero_frac']:.1%}
    """)

    mo.vstack([_fig, _stats])
    return


@app.cell
def _(mo):
    mo.md(r"""
    ---
    ## 13. Outlier Genes

    Genes where the scale factor or shape match is anomalous:
    - **High scale factor (>0.8)**: Langlieb expression nearly matches ABC — possible ambient RNA contamination
    - **Low scale factor (<0.08)**: Extreme compression — dropout-prone genes
    - **Low R² (<0.7)**: Shape doesn't match even after rescaling — different biology or technical artifact

    Known contaminated genes (Agrp, Oxt, Avp, Pmch) are flagged.
    """)
    return


@app.cell
def _(mo, scale_df):
    _contaminated = {"Agrp", "Oxt", "Avp", "Pmch"}

    _outliers = scale_df[
        (scale_df["scale_factor"] > 0.8)
        | (scale_df["scale_factor"] < 0.08)
        | (scale_df["norm_r2"] < 0.7)
    ].copy()

    # Add reason and contamination flag
    _reasons = []
    for _, _r in _outliers.iterrows():
        _parts = []
        if _r["scale_factor"] > 0.8:
            _parts.append("high scale (>0.8)")
        if _r["scale_factor"] < 0.08:
            _parts.append("low scale (<0.08)")
        if _r["norm_r2"] < 0.7:
            _parts.append("low R² (<0.7)")
        _reasons.append("; ".join(_parts))

    _outliers = _outliers.assign(
        reason=_reasons,
        known_contaminated=_outliers["gene"].isin(_contaminated),
    )

    # Sort: contaminated first, then by scale_factor descending
    _outliers = _outliers.sort_values(
        ["known_contaminated", "scale_factor"], ascending=[False, False]
    )

    _display_cols = [
        "gene", "scale_factor", "norm_slope", "norm_r2",
        "abc_95", "lang_95", "reason", "known_contaminated",
    ]
    outlier_df = _outliers[_display_cols].reset_index(drop=True)

    # Summary
    _n_high = (scale_df["scale_factor"] > 0.8).sum()
    _n_low = (scale_df["scale_factor"] < 0.08).sum()
    _n_bad_r2 = (scale_df["norm_r2"] < 0.7).sum()
    _summary = mo.md(f"""
    **Outlier summary**: {_n_high} genes with high scale factor (>0.8),
    {_n_low} with low scale factor (<0.08), {_n_bad_r2} with poor shape match (R²<0.7).
    Total unique outliers: {len(outlier_df)}.
    """)

    mo.vstack([_summary, mo.ui.table(outlier_df, label="Outlier Genes")])
    return


if __name__ == "__main__":
    app.run()
