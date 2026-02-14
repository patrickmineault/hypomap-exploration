import marimo

__generated_with = "0.19.4"
app = marimo.App()


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Somatodendritic Release Index (SDI) Validation

    Can we build a per-cell transcriptomic index that distinguishes neurons biased toward
    **somatodendritic neuropeptide release** from those biased toward **axonal terminal release**,
    using only somatic mRNA (MERFISH)?

    ## The Contrast

    | Group | Population | Release Mode | Marker | Region |
    |-------|-----------|-------------|--------|--------|
    | **Somatodendritic** | OXT/AVP magnocellular | Dendritic release | Oxt (+ Avp) | PVH, SO |
    | **Axonal** | HCRT (orexin) | Long-range axonal | Hcrt | LHA |
    | **Axonal** | CRH parvocellular | Median eminence axonal | Crh | PVH |

    Note: In SO, Oxt and Avp are 100% co-expressed in imputed data, so we combine
    OXT/AVP magnocellular as a single somatodendritic population.

    ## The Index

    Based on synaptotagmin compartmentalization (Bhatt et al. 2022):
    - **SYT4**: enriched in somatodendritic compartment
    - **SYT1**: canonical fast Ca2+ sensor for axonal terminal release
    - **SDI = log2((SYT4 + e) / (SYT1 + e))**

    **Caveat**: The SYT4/SYT1 story comes from DA neurons. If genes are imputed
    rather than measured, spatial smoothing could create or destroy the signal.
    """)
    return


@app.cell
def _():
    import marimo as mo
    import pandas as pd
    import numpy as np
    import plotly.graph_objects as go
    import plotly.express as px
    from plotly.subplots import make_subplots
    from scipy import stats
    from pathlib import Path
    return Path, go, make_subplots, mo, np, pd, stats


@app.cell
def _(Path):
    EXPRESSION_THRESHOLD = 3.0  # log2 scale

    # All release machinery genes we want (from imputed data)
    RELEASE_GENES_ALL = [
        'Syt1', 'Syt4', 'Syt7',
        'Syt2', 'Syt6', 'Syt10', 'Syt11', 'Syt17',
        'Scg2', 'Scg3', 'Scg5', 'Chga', 'Chgb',
        'Pcsk1', 'Pcsk2',
        'Snap25', 'Snap23', 'Snap29',
        'Vamp2', 'Vamp7',
        'Stx1a', 'Stx1b',
        'Cadps', 'Cadps2',
        'Sv2a', 'Sv2b', 'Sv2c',
        'Cplx1', 'Cplx2',
        'Unc13a', 'Unc13b',
        'Doc2a', 'Doc2b',
        'Rph3a', 'Stxbp1',
    ]

    # Subset available in the measured MERFISH 550-gene panel
    RELEASE_GENES_MEASURED = [
        'Syt2', 'Syt6', 'Syt10', 'Syt17',
        'Sv2b', 'Sv2c',
        'Rab3b',
        'Unc13c',
        'Synpr',
        'Baiap3',
    ]

    CACHE_PARQUET = Path('data/processed/mouse_abc/release_machinery_expression.parquet')
    MEASURED_H5AD = Path(
        'data/raw/abc_atlas_cache/expression_matrices/'
        'MERFISH-C57BL6J-638850/20230830/C57BL6J-638850-log2.h5ad'
    )
    return (
        CACHE_PARQUET,
        EXPRESSION_THRESHOLD,
        MEASURED_H5AD,
        RELEASE_GENES_ALL,
        RELEASE_GENES_MEASURED,
    )


@app.cell
def _(pd):
    cell_metadata = pd.read_parquet('data/processed/mouse_abc/cell_metadata.parquet')
    expr_df = pd.read_parquet('data/processed/mouse_abc/neuropeptide_expression.parquet')

    # Join expression with metadata, filter to neurons
    neurons_with_expr = expr_df.join(
        cell_metadata.set_index('cell_id')[
            ['region', 'class', 'subclass', 'supertype', 'cluster', 'x', 'y', 'z']
        ],
        how='inner'
    )
    NON_NEURONAL = ['30 Astro-Epen', '31 OPC-Oligo', '33 Vascular', '34 Immune']
    neurons_with_expr = neurons_with_expr[
        ~neurons_with_expr['class'].isin(NON_NEURONAL)
    ].copy()

    print(f"Neurons: {len(neurons_with_expr):,}")
    print(f"Regions: {neurons_with_expr['region'].nunique()}")
    return cell_metadata, neurons_with_expr


@app.cell
def _(
    CACHE_PARQUET,
    MEASURED_H5AD,
    RELEASE_GENES_ALL,
    RELEASE_GENES_MEASURED,
    cell_metadata,
    np,
    pd,
):
    # Load release machinery expression.
    # Priority: 1) cached parquet from imputed extraction, 2) measured MERFISH h5ad
    _release_expr = None
    gene_source = {}
    data_mode = 'none'

    # Try 1: Cached parquet (from previous imputed extraction)
    if CACHE_PARQUET.exists():
        _release_expr = pd.read_parquet(CACHE_PARQUET)
        data_mode = 'imputed_cached'
        for _g in _release_expr.columns:
            gene_source[_g] = 'imputed'
        print(f"Loaded cached imputed data: {len(_release_expr.columns)} genes, "
              f"{len(_release_expr):,} cells")

    # Try 2: Extract from measured MERFISH h5ad
    if _release_expr is None and MEASURED_H5AD.exists():
        print("No imputed cache found. Extracting from measured MERFISH h5ad...")
        try:
            import anndata
            _adata = anndata.read_h5ad(str(MEASURED_H5AD), backed='r')
            _panel_genes = set(_adata.var['gene_symbol'].tolist())
            _to_extract = [g for g in RELEASE_GENES_MEASURED if g in _panel_genes]

            if _to_extract:
                _hy_ids = set(cell_metadata['cell_id'].values)
                _gmask = _adata.var['gene_symbol'].isin(_to_extract)
                _gidx = np.where(_gmask)[0]
                _cmask = pd.Series(_adata.obs.index).isin(_hy_ids)
                _cidx = np.where(_cmask)[0]

                print(f"  Extracting {len(_to_extract)} genes for {len(_cidx)} cells...")

                _chunks = []
                _ids = []
                for _s in range(0, len(_cidx), 50000):
                    _e = min(_s + 50000, len(_cidx))
                    _ci = _cidx[_s:_e]
                    _ch = _adata.X[_ci, :][:, _gidx]
                    if hasattr(_ch, 'toarray'):
                        _ch = _ch.toarray()
                    _chunks.append(_ch)
                    _ids.extend(_adata.obs.index[_ci].tolist())

                _gnames = _adata.var.loc[_gmask, 'gene_symbol'].tolist()
                _adata.file.close()

                _release_expr = pd.DataFrame(
                    np.vstack(_chunks), index=_ids, columns=_gnames
                )
                data_mode = 'measured'
                for _g in _to_extract:
                    gene_source[_g] = 'measured'
                print(f"  Extracted: {sorted(_to_extract)}")
            else:
                _adata.file.close()
                print("  No target genes in measured panel.")
        except Exception as _ex:
            print(f"  Extraction failed: {_ex}")

    # Try 3: Download imputed data (large, ~50GB)
    if _release_expr is None:
        print("\nNo local data available. Attempting imputed data download...")
        print("(This requires ~50GB download from Allen Brain Cell Atlas)")
        try:
            from abc_atlas_access.abc_atlas_cache.abc_project_cache import (
                AbcProjectCache,
            )

            _cache = AbcProjectCache.from_cache_dir(
                CACHE_PARQUET.parent.parent.parent / "raw" / "abc_atlas_cache"
            )
            _cache.load_latest_manifest()
            _expr_path = _cache.get_file_path(
                directory="MERFISH-C57BL6J-638850-imputed",
                file_name="C57BL6J-638850-imputed/log2"
            )

            import anndata
            _adata = anndata.read_h5ad(str(_expr_path), backed='r')
            _imp_genes = set(_adata.var['gene_symbol'].tolist())
            _to_extract = [g for g in RELEASE_GENES_ALL if g in _imp_genes]

            _hy_ids = set(cell_metadata['cell_id'].values)
            _gmask = _adata.var['gene_symbol'].isin(_to_extract)
            _gidx = np.where(_gmask)[0]
            _cmask = pd.Series(_adata.obs.index).isin(_hy_ids)
            _cidx = np.where(_cmask)[0]

            print(f"  Extracting {len(_to_extract)} genes for {len(_cidx)} cells...")
            _chunks = []
            _ids = []
            for _s in range(0, len(_cidx), 50000):
                _e = min(_s + 50000, len(_cidx))
                _ci = _cidx[_s:_e]
                _ch = _adata.X[_ci, :][:, _gidx]
                if hasattr(_ch, 'toarray'):
                    _ch = _ch.toarray()
                _chunks.append(_ch)
                _ids.extend(_adata.obs.index[_ci].tolist())

            _gnames = _adata.var.loc[_gmask, 'gene_symbol'].tolist()
            _adata.file.close()

            _release_expr = pd.DataFrame(
                np.vstack(_chunks), index=_ids, columns=_gnames
            )
            data_mode = 'imputed_downloaded'
            for _g in _to_extract:
                gene_source[_g] = 'imputed'

            # Cache for next time
            CACHE_PARQUET.parent.mkdir(parents=True, exist_ok=True)
            _release_expr.to_parquet(CACHE_PARQUET)
            print(f"  Cached to {CACHE_PARQUET}")
        except Exception as _ex:
            print(f"  Download/extraction failed: {_ex}")

    if _release_expr is None:
        _release_expr = pd.DataFrame()

    release_expr = _release_expr

    print(f"\nData mode: {data_mode}")
    if len(release_expr.columns) > 0:
        print(f"Available genes: {sorted(release_expr.columns.tolist())}")
    else:
        print("WARNING: No release machinery expression data available.")
    return data_mode, gene_source, release_expr


@app.cell(hide_code=True)
def _(data_mode, gene_source, mo, release_expr):
    _measured = [g for g, s in gene_source.items() if s == 'measured']
    _imputed = [g for g, s in gene_source.items() if s == 'imputed']
    _has_core = all(g in release_expr.columns for g in ['Syt1', 'Syt4'])

    mo.md(f"""
    ---
    ## Data Availability

    **Mode**: `{data_mode}` | **Genes loaded**: {len(release_expr.columns)}

    {"- **Measured** (MERFISH panel): " + ", ".join(sorted(_measured)) if _measured else ""}
    {"- **Imputed**: " + ", ".join(sorted(_imputed)) if _imputed else ""}

    {"**Core SDI genes (Syt1, Syt4) available.**" if _has_core else
     "**Core SDI genes (Syt1, Syt4) NOT available.** The primary SDI cannot be computed. "
     "Analysis will use whatever genes are available. To get full results, "
     "run with imputed MERFISH data (~50GB download from Allen Brain Cell Atlas)."}
    """)
    return


@app.cell
def _(EXPRESSION_THRESHOLD, neurons_with_expr, pd):
    # Identify three populations (OXT/AVP co-express in SO, so no separate AVP_magno)
    _n = neurons_with_expr

    pop_labels = pd.Series('other', index=_n.index, name='population')

    # OXT/AVP magnocellular: Oxt+ in PVH/SO/PVa
    # Note: In SO, 100% of Oxt+ neurons also express Avp (and vice versa) in imputed data.
    # These are the textbook somatodendritic releasers.
    _oxt_mask = (
        (_n['Oxt'] > EXPRESSION_THRESHOLD) &
        (_n['region'].isin(['PVH', 'SO', 'PVa']))
    )
    pop_labels[_oxt_mask] = 'OXT_magno'

    # HCRT neurons: Hcrt+ in LHA/DMH/PeF
    _hcrt_mask = (
        (_n['Hcrt'] > EXPRESSION_THRESHOLD) &
        (_n['region'].isin(['LHA', 'DMH', 'PeF']))
    )
    pop_labels[_hcrt_mask] = 'HCRT'

    # CRH parvocellular: Crh+ in PVH, excluding Oxt+ (likely magnocellular)
    # Avp-only PVH neurons (supertype 0587) also express Crh and release at
    # median eminence (axonal), so they belong in the axonal group.
    _crh_parvo_mask = (
        (_n['Crh'] > EXPRESSION_THRESHOLD) &
        (_n['region'] == 'PVH') &
        (~_oxt_mask) &
        (_n['Oxt'] <= EXPRESSION_THRESHOLD)
    )
    pop_labels[_crh_parvo_mask] = 'CRH_parvo'

    # Release mode
    release_mode = pd.Series('other', index=_n.index, name='release_mode')
    release_mode[pop_labels == 'OXT_magno'] = 'somatodendritic'
    release_mode[pop_labels.isin(['HCRT', 'CRH_parvo'])] = 'axonal'

    # Build population dataframe
    pop_df = pd.DataFrame({
        'population': pop_labels,
        'release_mode': release_mode,
    })

    # Add spatial coords, region, marker genes
    for _col in ['Oxt', 'Avp', 'Hcrt', 'Crh', 'x', 'y', 'z', 'region', 'supertype']:
        if _col in _n.columns:
            pop_df[_col] = _n[_col]

    # Summary
    print("=== Population Identification ===\n")
    for _pop in ['OXT_magno', 'HCRT', 'CRH_parvo']:
        _mask = pop_df['population'] == _pop
        _n_cells = _mask.sum()
        _regions = pop_df.loc[_mask, 'region'].value_counts().head(3)
        _rstr = ', '.join([f"{r}({c})" for r, c in _regions.items()])
        print(f"  {_pop}: {_n_cells} cells  [{_rstr}]")

    # Show SO co-expression note
    _so_oxt = pop_df[(pop_df['region'] == 'SO') & (pop_df['population'] == 'OXT_magno')]
    print(f"\n  Note: {len(_so_oxt)} SO cells in OXT_magno (all co-express Avp)")

    _n_sd = (release_mode == 'somatodendritic').sum()
    _n_ax = (release_mode == 'axonal').sum()
    print(f"\n  Somatodendritic group: {_n_sd} cells")
    print(f"  Axonal group: {_n_ax} cells")
    return pop_df, release_mode


@app.cell(hide_code=True)
def _(mo, release_mode):
    _sd = (release_mode == 'somatodendritic').sum()
    _ax = (release_mode == 'axonal').sum()
    mo.md(f"""
    ---
    ## Population Summary

    **Somatodendritic** ({_sd} cells): OXT/AVP magnocellular neurons in PVH/SO/PVa.
    In SO, all Oxt+ neurons co-express Avp (100% overlap in imputed data), consistent
    with the known co-expression in magnocellular neurons.

    **Axonal** ({_ax} cells): HCRT neurons (long-range LHA projectors) and CRH
    parvocellular neurons (axonal release at median eminence). The CRH_parvo group
    includes Avp-only PVH neurons (supertype 0587) that co-express Crh and release
    CRH/AVP at median eminence terminals.
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    ## Release Machinery Gene Expression by Population
    """)
    return


@app.cell
def _(go, np, pop_df, release_expr):
    # Expression heatmap
    if len(release_expr.columns) == 0:
        print("No release machinery expression data. Skipping heatmap.")
    else:
        _pops = ['OXT_magno', 'HCRT', 'CRH_parvo']
        _genes = sorted(release_expr.columns.tolist())

        _mean_matrix = []
        _pct_matrix = []
        for _pop in _pops:
            _cells = pop_df.index[pop_df['population'] == _pop]
            _data = release_expr.loc[release_expr.index.isin(_cells), _genes]
            _mean_matrix.append(_data.mean().values)
            _pct_matrix.append((_data > 3.0).mean().values * 100)

        _mean_arr = np.array(_mean_matrix)
        _pct_arr = np.array(_pct_matrix)

        _hover = [
            [f'{_genes[j]}<br>{_pops[i]}<br>Mean: {_mean_arr[i,j]:.2f}<br>'
             f'{_pct_arr[i,j]:.1f}% expressing'
             for j in range(len(_genes))]
            for i in range(len(_pops))
        ]

        fig_heatmap = go.Figure(data=go.Heatmap(
            z=_mean_arr,
            x=_genes,
            y=_pops,
            hovertext=_hover,
            hoverinfo='text',
            colorscale='Viridis',
            colorbar=dict(title='Mean log2'),
        ))
        fig_heatmap.update_layout(
            title='Release Machinery: Mean Expression by Population',
            xaxis_title='Gene',
            yaxis_title='Population',
            height=350,
            width=max(800, len(_genes) * 35),
            xaxis=dict(tickangle=45, tickfont=dict(size=10)),
        )
        fig_heatmap.show()

        # Also show % expressing heatmap
        fig_pct = go.Figure(data=go.Heatmap(
            z=_pct_arr,
            x=_genes,
            y=_pops,
            hovertext=_hover,
            hoverinfo='text',
            colorscale='YlOrRd',
            colorbar=dict(title='% expressing'),
        ))
        fig_pct.update_layout(
            title='Release Machinery: % Expressing by Population',
            xaxis_title='Gene',
            yaxis_title='Population',
            height=350,
            width=max(800, len(_genes) * 35),
            xaxis=dict(tickangle=45, tickfont=dict(size=10)),
        )
        fig_pct.show()
    return


@app.cell(hide_code=True)
def _(mo, release_expr):
    _has_core = all(g in release_expr.columns for g in ['Syt1', 'Syt4'])
    mo.md(f"""
    ---
    ## Somatodendritic Release Index (SDI)

    {"Computing **SDI = log2((Syt4 + e) / (Syt1 + e))** for each neuron." if _has_core else
     "**Syt1/Syt4 not available.** Computing alternative indices from available genes."}

    Also exploring alternative ratios using other gene pairs.
    """)
    return


@app.cell
def _(np, pop_df, release_expr):
    # Compute SDI indices
    sdi_df = pop_df[['population', 'release_mode', 'region', 'x', 'y', 'z']].copy()

    epsilon = 0.5
    indices_computed = []

    _has = lambda g: len(release_expr.columns) > 0 and g in release_expr.columns

    def _get(gene):
        return release_expr.reindex(sdi_df.index)[gene].fillna(0)

    if _has('Syt1') and _has('Syt4'):
        sdi_df['SDI_Syt4_Syt1'] = np.log2((_get('Syt4') + epsilon) / (_get('Syt1') + epsilon))
        indices_computed.append('SDI_Syt4_Syt1')
        print("Computed SDI_Syt4_Syt1 = log2((Syt4 + e) / (Syt1 + e))")

    if _has('Syt7') and _has('Syt1'):
        sdi_df['SDI_Syt7_Syt1'] = np.log2((_get('Syt7') + epsilon) / (_get('Syt1') + epsilon))
        indices_computed.append('SDI_Syt7_Syt1')
        print("Computed SDI_Syt7_Syt1 = log2((Syt7 + e) / (Syt1 + e))")

    if _has('Syt4') and _has('Syt7') and _has('Syt1'):
        sdi_df['SDI_combined'] = np.log2(
            (_get('Syt4') + _get('Syt7') + epsilon) / (_get('Syt1') + epsilon)
        )
        indices_computed.append('SDI_combined')
        print("Computed SDI_combined = log2((Syt4 + Syt7 + e) / (Syt1 + e))")

    # Fallback: Syt2 as Syt1 proxy (if only measured data)
    if _has('Syt2') and not _has('Syt1'):
        if _has('Syt17'):
            sdi_df['SDI_Syt17_Syt2'] = np.log2(
                (_get('Syt17') + epsilon) / (_get('Syt2') + epsilon)
            )
            indices_computed.append('SDI_Syt17_Syt2')
            print("Computed SDI_Syt17_Syt2 (measured-only fallback)")

    # Dense-core vesicle indices
    if _has('Cadps') and _has('Snap25'):
        sdi_df['SDI_Cadps_Snap25'] = np.log2(
            (_get('Cadps') + epsilon) / (_get('Snap25') + epsilon)
        )
        indices_computed.append('SDI_Cadps_Snap25')
        print("Computed SDI_Cadps_Snap25")

    if _has('Pcsk2') and _has('Pcsk1'):
        sdi_df['SDI_Pcsk2_Pcsk1'] = np.log2(
            (_get('Pcsk2') + epsilon) / (_get('Pcsk1') + epsilon)
        )
        indices_computed.append('SDI_Pcsk2_Pcsk1')
        print("Computed SDI_Pcsk2_Pcsk1")

    # Add all release genes to sdi_df for downstream use
    if len(release_expr.columns) > 0:
        for _g in release_expr.columns:
            if _g not in sdi_df.columns:
                sdi_df[_g] = release_expr.reindex(sdi_df.index)[_g].fillna(0)

    if not indices_computed:
        print("WARNING: No SDI indices could be computed.")
        print("Core genes (Syt1, Syt4) are required for the primary SDI.")
    else:
        print(f"\nIndices computed: {indices_computed}")
    return indices_computed, sdi_df


@app.cell
def _(go, indices_computed, make_subplots, sdi_df, stats):
    # SDI box plots by population
    if not indices_computed:
        print("No SDI indices to plot.")
    else:
        _n_idx = len(indices_computed)
        fig_sdi = make_subplots(
            rows=1, cols=_n_idx,
            subplot_titles=indices_computed,
        )

        _colors = {
            'OXT_magno': '#e41a1c',
            'HCRT': '#4daf4a',
            'CRH_parvo': '#984ea3',
        }
        _pops = ['OXT_magno', 'HCRT', 'CRH_parvo']

        for _ci, _idx in enumerate(indices_computed):
            for _pop in _pops:
                _vals = sdi_df.loc[sdi_df['population'] == _pop, _idx].dropna()
                if len(_vals) > 0:
                    fig_sdi.add_trace(
                        go.Box(
                            y=_vals,
                            name=_pop,
                            marker_color=_colors[_pop],
                            showlegend=(_ci == 0),
                            boxmean=True,
                        ),
                        row=1, col=_ci + 1
                    )

        fig_sdi.update_layout(
            title='SDI by Population (red/blue = SD, green/purple = axonal)',
            height=500,
            width=max(600, 350 * _n_idx),
        )
        fig_sdi.show()

        # Statistical tests
        print("\n=== Statistical Tests: Somatodendritic vs Axonal ===\n")
        for _idx in indices_computed:
            _sd = sdi_df.loc[sdi_df['release_mode'] == 'somatodendritic', _idx].dropna()
            _ax = sdi_df.loc[sdi_df['release_mode'] == 'axonal', _idx].dropna()

            if len(_sd) > 0 and len(_ax) > 0:
                _pooled = ((_sd.var() + _ax.var()) / 2) ** 0.5
                _d = (_sd.mean() - _ax.mean()) / _pooled if _pooled > 0 else 0
                _u, _p_mw = stats.mannwhitneyu(_sd, _ax, alternative='two-sided')
                _auroc = _u / (len(_sd) * len(_ax))

                print(f"{_idx}:")
                print(f"  SD:     mean={_sd.mean():.3f} sd={_sd.std():.3f} n={len(_sd)}")
                print(f"  Axonal: mean={_ax.mean():.3f} sd={_ax.std():.3f} n={len(_ax)}")
                print(f"  Cohen's d = {_d:+.3f}")
                print(f"  AUROC = {_auroc:.3f}")
                print(f"  Mann-Whitney p = {_p_mw:.2e}")
                print()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    ## Effect Sizes: Individual Release Machinery Genes

    Which individual genes best separate somatodendritic from axonal populations?
    Positive Cohen's d = higher in somatodendritic group.
    """)
    return


@app.cell
def _(pd, pop_df, release_expr, stats):
    # Effect sizes for individual genes
    if len(release_expr.columns) == 0:
        print("No release machinery data.")
        effect_sizes = pd.DataFrame()
    else:
        _results = []
        for _gene in sorted(release_expr.columns):
            _sd_cells = pop_df.index[pop_df['release_mode'] == 'somatodendritic']
            _ax_cells = pop_df.index[pop_df['release_mode'] == 'axonal']

            _sd_vals = release_expr.loc[
                release_expr.index.isin(_sd_cells), _gene
            ].dropna()
            _ax_vals = release_expr.loc[
                release_expr.index.isin(_ax_cells), _gene
            ].dropna()

            if len(_sd_vals) > 10 and len(_ax_vals) > 10:
                _pooled = ((_sd_vals.var() + _ax_vals.var()) / 2) ** 0.5
                _d = ((_sd_vals.mean() - _ax_vals.mean()) / _pooled) if _pooled > 0 else 0
                _u, _p = stats.mannwhitneyu(_sd_vals, _ax_vals, alternative='two-sided')
                _auroc = _u / (len(_sd_vals) * len(_ax_vals))

                _results.append({
                    'gene': _gene,
                    'mean_SD': _sd_vals.mean(),
                    'mean_axonal': _ax_vals.mean(),
                    'diff': _sd_vals.mean() - _ax_vals.mean(),
                    'cohens_d': _d,
                    'auroc': _auroc,
                    'p_value': _p,
                })

        effect_sizes = pd.DataFrame(_results).sort_values('cohens_d', ascending=False)

        print("=== Effect Sizes: Somatodendritic vs Axonal ===\n")
        print(f"{'Gene':12s} {'d':>8s} {'AUROC':>7s} {'mean_SD':>9s} {'mean_ax':>9s}  sig")
        print("-" * 55)
        for _, _r in effect_sizes.iterrows():
            _star = ('***' if _r['p_value'] < 0.001 else
                     '**' if _r['p_value'] < 0.01 else
                     '*' if _r['p_value'] < 0.05 else '')
            print(f"{_r['gene']:12s} {_r['cohens_d']:+8.3f} {_r['auroc']:7.3f} "
                  f"{_r['mean_SD']:9.2f} {_r['mean_axonal']:9.2f}  {_star}")
    return (effect_sizes,)


@app.cell
def _(effect_sizes, go):
    # Effect size bar chart
    if len(effect_sizes) == 0:
        print("No effect size data.")
    else:
        _sorted = effect_sizes.sort_values('cohens_d')
        _colors = ['#e41a1c' if d > 0 else '#377eb8' for d in _sorted['cohens_d']]

        fig_effect = go.Figure(data=go.Bar(
            x=_sorted['cohens_d'],
            y=_sorted['gene'],
            orientation='h',
            marker_color=_colors,
            hovertext=[
                f"{r['gene']}: d={r['cohens_d']:.3f}, p={r['p_value']:.2e}"
                for _, r in _sorted.iterrows()
            ],
            hoverinfo='text',
        ))
        fig_effect.update_layout(
            title="Cohen's d: Somatodendritic vs Axonal<br>"
                  "<sup>Red = higher in SD, Blue = higher in axonal</sup>",
            xaxis_title="Cohen's d",
            yaxis_title='Gene',
            height=max(400, len(_sorted) * 22),
            width=700,
        )
        fig_effect.show()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    ## Sanity Check: Is SDI Just Measuring "PVH-ness"?

    If SDI merely recapitulates a PVH vs LHA transcriptomic signature (rather than
    release biology), it would separate ALL PVH neurons from ALL LHA neurons equally well.

    **Test 1**: Compare SDI across ALL PVH vs ALL LHA neurons (regardless of peptide identity).
    The region effect should be smaller than the population-specific effect.

    **Test 2**: Compare OXT_magno vs CRH_parvo WITHIN PVH. If SDI reflects release mode,
    it should separate these two populations even though they share the same region.
    """)
    return


@app.cell
def _(indices_computed, np, sdi_df, stats):
    # Sanity checks
    if not indices_computed:
        print("No SDI indices for sanity check.")
    else:
        print("=== Test 1: ALL PVH vs ALL LHA ===\n")
        for _idx in indices_computed:
            _pvh = sdi_df.loc[sdi_df['region'] == 'PVH', _idx].dropna()
            _lha = sdi_df.loc[sdi_df['region'] == 'LHA', _idx].dropna()

            if len(_pvh) > 10 and len(_lha) > 10:
                _d_region = (
                    (_pvh.mean() - _lha.mean()) /
                    ((_pvh.var() + _lha.var()) / 2) ** 0.5
                )
                _, _p_region = stats.mannwhitneyu(
                    _pvh, _lha, alternative='two-sided'
                )

                _sd = sdi_df.loc[
                    sdi_df['release_mode'] == 'somatodendritic', _idx
                ].dropna()
                _ax = sdi_df.loc[
                    sdi_df['release_mode'] == 'axonal', _idx
                ].dropna()
                _pooled = ((_sd.var() + _ax.var()) / 2) ** 0.5
                _d_pop = (_sd.mean() - _ax.mean()) / _pooled if _pooled > 0 else 0

                print(f"{_idx}:")
                print(f"  ALL PVH vs ALL LHA:    d={_d_region:.3f} "
                      f"(n={len(_pvh)},{len(_lha)}) p={_p_region:.2e}")
                print(f"  Specific SD vs Axonal: d={_d_pop:.3f}")

                _ratio = abs(_d_region) / abs(_d_pop) if abs(_d_pop) > 0 else np.inf
                if _ratio < 0.5:
                    print(f"  --> GOOD: Region effect much smaller "
                          f"({abs(_d_region):.3f} vs {abs(_d_pop):.3f})")
                elif _ratio > 0.8:
                    print(f"  --> WARNING: Region effect similar to population effect. "
                          f"SDI may reflect region identity.")
                else:
                    print(f"  --> MODERATE: Region effect present but smaller.")
                print()

        print("=== Test 2: Within-PVH (OXT_magno vs CRH_parvo) ===\n")
        for _idx in indices_computed:
            _oxt = sdi_df.loc[sdi_df['population'] == 'OXT_magno', _idx].dropna()
            _crh = sdi_df.loc[sdi_df['population'] == 'CRH_parvo', _idx].dropna()

            if len(_oxt) > 10 and len(_crh) > 10:
                _pooled = ((_oxt.var() + _crh.var()) / 2) ** 0.5
                _d = (_oxt.mean() - _crh.mean()) / _pooled if _pooled > 0 else 0
                _, _p = stats.mannwhitneyu(_oxt, _crh, alternative='two-sided')

                print(f"{_idx}:")
                print(f"  OXT_magno vs CRH_parvo (both PVH): "
                      f"d={_d:.3f}, p={_p:.2e}")
                if abs(_d) > 0.3:
                    print(f"  --> SDI separates populations WITHIN PVH")
                else:
                    print(f"  --> SDI does NOT separate well within PVH")
                print()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    ## Spatial Distribution of SDI
    """)
    return


@app.cell
def _(go, indices_computed, np, sdi_df):
    # Spatial plot
    if not indices_computed:
        print("No SDI indices for spatial plot.")
    else:
        _idx = indices_computed[0]
        _focus_regions = ['PVH', 'SO', 'PVa', 'LHA', 'DMH', 'PeF']
        _focus = sdi_df[sdi_df['region'].isin(_focus_regions)].dropna(
            subset=[_idx, 'x', 'y']
        ).copy()

        if len(_focus) > 0:
            _vmin = np.percentile(_focus[_idx], 5)
            _vmax = np.percentile(_focus[_idx], 95)

            fig_spatial = go.Figure()

            # Background: all neurons in focus regions
            _other = _focus[_focus['population'] == 'other']
            fig_spatial.add_trace(go.Scattergl(
                x=_other['x'], y=_other['y'],
                mode='markers',
                marker=dict(
                    size=2,
                    color=_other[_idx],
                    colorscale='RdBu_r',
                    cmin=_vmin, cmax=_vmax,
                    opacity=0.3,
                    colorbar=dict(title=_idx),
                ),
                name='Other neurons',
            ))

            # Labeled populations
            _pop_colors = {
                'OXT_magno': '#e41a1c',
                'HCRT': '#4daf4a',
                'CRH_parvo': '#984ea3',
            }
            for _pop, _color in _pop_colors.items():
                _pd = _focus[_focus['population'] == _pop]
                if len(_pd) > 0:
                    _mean = _pd[_idx].mean()
                    fig_spatial.add_trace(go.Scattergl(
                        x=_pd['x'], y=_pd['y'],
                        mode='markers',
                        marker=dict(
                            size=5, color=_color, opacity=0.7,
                            line=dict(width=0.5, color='black'),
                        ),
                        name=f'{_pop} (mean={_mean:.2f})',
                    ))

            fig_spatial.update_layout(
                title=f'Spatial Distribution of {_idx}',
                xaxis_title='X (mm)', yaxis_title='Y (mm)',
                height=600, width=800,
                yaxis=dict(scaleanchor='x'),
            )
            fig_spatial.show()
    return


@app.cell(hide_code=True)
def _(indices_computed, mo, sdi_df):
    if not indices_computed:
        _text = """
    ---
    ## Conclusion

    **No SDI indices could be computed.** The analysis requires Syt1 and Syt4 expression
    from the imputed MERFISH dataset (~50GB from Allen Brain Cell Atlas).

    The measured MERFISH panel (550 genes) does not include Syt1 or Syt4.
    """
    else:
        _idx = indices_computed[0]
        _sd = sdi_df.loc[sdi_df['release_mode'] == 'somatodendritic', _idx].dropna()
        _ax = sdi_df.loc[sdi_df['release_mode'] == 'axonal', _idx].dropna()
        _pooled = ((_sd.var() + _ax.var()) / 2) ** 0.5
        _d = (_sd.mean() - _ax.mean()) / _pooled if len(_sd) > 0 and _pooled > 0 else 0
        _dir = "higher" if _sd.mean() > _ax.mean() else "lower"

        if abs(_d) > 0.5:
            _verdict = "**SDI discriminates the two release mode groups.**"
        elif abs(_d) > 0.2:
            _verdict = ("**SDI effect is modest** - may not reliably discriminate "
                        "at single-cell level.")
        else:
            _verdict = "**SDI does NOT meaningfully discriminate the two groups.**"

        _text = f"""
    ---
    ## Conclusion

    **Primary SDI ({_idx})**: SD mean = {_sd.mean():.3f}, Axonal mean = {_ax.mean():.3f},
    Cohen's d = {_d:.3f}. Somatodendritic group has **{_dir}** SDI.

    {_verdict}

    ### Caveats
    - Release machinery genes are **imputed**, not directly measured by MERFISH
    - Imputation spatial smoothing could create or destroy signal
    - The SYT4/SYT1 compartmental story comes from DA neurons, not peptidergic hypothalamic neurons
    - CRH parvocellular identification is approximate (Crh+ PVH, excluding Oxt co-expressors)
    """
    mo.md(_text)
    return


if __name__ == "__main__":
    app.run()
