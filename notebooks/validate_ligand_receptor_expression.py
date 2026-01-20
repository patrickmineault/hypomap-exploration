import marimo

__generated_with = "0.19.4"
app = marimo.App()


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Validation: Neuropeptide Ligand/Receptor Expression by Region

    This notebook validates expression profiles against known neuroscience facts.

    **Two validation sources:**
    1. **np_map.csv** - Curated neuropeptide systems with expected hypothalamic nuclei
    2. **hypothalamic_nuclei_genes.csv** - Manual gene-region annotations

    **Scoring**: 1 point for each gene where the expected region ranks in the **top 5** for expression.
    """)
    return


@app.cell
def _():
    EXPRESSION_THRESHOLD = 2.0
    return (EXPRESSION_THRESHOLD,)


@app.cell
def _():
    import pandas as pd
    import numpy as np
    import plotly.graph_objects as go

    # Load data
    expr_df = pd.read_parquet('../data/processed/mouse_abc/neuropeptide_expression.parquet')
    cell_metadata = pd.read_parquet('../data/processed/mouse_abc/cell_metadata.parquet')

    # Join expression with region and cell type hierarchy, filter to neurons only
    expr_with_region = expr_df.join(
        cell_metadata.set_index('cell_id')[['region', 'class', 'subclass', 'supertype', 'cluster']],
        how='inner'
    )
    NON_NEURONAL = ['30 Astro-Epen', '31 OPC-Oligo', '33 Vascular', '34 Immune']
    expr_with_region = expr_with_region[~expr_with_region['class'].isin(NON_NEURONAL)]

    # Load expected mappings
    expected_df = pd.read_csv('../data/generated/mouse_common/hypothalamic_nuclei_genes.csv')

    # Load ligand-receptor list to know what's tracked
    np_map = pd.read_csv('../data/generated/mouse_common/np_map.csv')
    tracked_genes = set(np_map['Ligand_Gene'].unique()) | set(np_map['Receptor_Gene'].unique())

    # Available genes in expression data
    available_genes = set(expr_df.columns)

    print(f"Neurons with expression data: {len(expr_with_region):,}")
    print(f"Tracked genes in LR list: {len(tracked_genes)}")
    print(f"Available in expression data: {len(available_genes)}")
    return (
        available_genes,
        expected_df,
        expr_with_region,
        go,
        np,
        pd,
        tracked_genes,
    )


@app.cell
def _(EXPRESSION_THRESHOLD, expr_with_region):
    def get_region_ranks(_gene: str) -> dict:
        """Get rank of each region for a gene (1 = highest expression)."""
        if _gene not in expr_with_region.columns:
            return {}
        stats = expr_with_region.groupby('region').agg(total_cells=('region', 'count'), expressing_cells=(_gene, lambda x: (x > EXPRESSION_THRESHOLD).sum()))
        stats['pct_expressing'] = 100 * stats['expressing_cells'] / stats['total_cells']
        stats = stats.sort_values('pct_expressing', ascending=False)
        stats['rank'] = range(1, len(stats) + 1)
        return stats['rank'].to_dict()
    return (get_region_ranks,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    ## Validation 1: np_map.csv (Curated Neuropeptide Systems)

    Validates ligand and receptor expression against expected hypothalamic nuclei from the curated neuropeptide map.
    """)
    return


@app.cell
def _(available_genes, get_region_ranks, pd):
    # Load np_map.csv and validate LIGANDS
    np_map_df = pd.read_csv('../data/generated/mouse_common/np_map.csv')

    # Extract unique (gene, region) pairs from Hypothalamic_Nuclei column for LIGANDS
    np_map_ligand_results = []
    for _, _row in np_map_df.iterrows():
        _gene = _row['Ligand_Gene']
        _system = _row['System']
        _regions_str = _row['Hypothalamic_Nuclei']

        if pd.isna(_regions_str) or _gene not in available_genes:
            continue

        # Parse semicolon-separated regions
        _expected_regions = [r.strip() for r in _regions_str.split(';')]
        _ranks = get_region_ranks(_gene)

        if not _ranks:
            continue

        for _region in _expected_regions:
            _rank = _ranks.get(_region, 999)
            _in_top5 = _rank <= 5
            np_map_ligand_results.append({
                'system': _system,
                'gene': _gene,
                'expected_region': _region,
                'rank': _rank,
                'in_top5': _in_top5,
                'score': 1 if _in_top5 else 0
            })

    np_map_ligand_df = pd.DataFrame(np_map_ligand_results)

    # Deduplicate (same gene-region pair may appear multiple times)
    np_map_ligand_df = np_map_ligand_df.drop_duplicates(subset=['gene', 'expected_region'])

    print(f"=== np_map.csv LIGAND Validation ===")
    print(f"Testable gene-region pairs: {len(np_map_ligand_df)}")
    print(f"Unique ligand genes tested: {np_map_ligand_df['gene'].nunique()}")
    print(f"Unique systems: {np_map_ligand_df['system'].nunique()}")
    return np_map_df, np_map_ligand_df


@app.cell
def _(available_genes, get_region_ranks, np_map_df, pd):
    # Validate RECEPTORS from np_map.csv
    np_map_receptor_results = []
    for _, _row in np_map_df.iterrows():
        _gene = _row['Receptor_Gene']
        _system = _row['System']
        _regions_str = _row['Hypothalamic_Nuclei']

        if pd.isna(_regions_str) or pd.isna(_gene) or _gene not in available_genes:
            continue

        # Parse semicolon-separated regions
        _expected_regions = [r.strip() for r in _regions_str.split(';')]
        _ranks = get_region_ranks(_gene)

        if not _ranks:
            continue

        for _region in _expected_regions:
            _rank = _ranks.get(_region, 999)
            _in_top5 = _rank <= 5
            np_map_receptor_results.append({
                'system': _system,
                'gene': _gene,
                'expected_region': _region,
                'rank': _rank,
                'in_top5': _in_top5,
                'score': 1 if _in_top5 else 0
            })

    np_map_receptor_df = pd.DataFrame(np_map_receptor_results)

    # Deduplicate (same gene-region pair may appear multiple times)
    np_map_receptor_df = np_map_receptor_df.drop_duplicates(subset=['gene', 'expected_region'])

    print(f"=== np_map.csv RECEPTOR Validation ===")
    print(f"Testable gene-region pairs: {len(np_map_receptor_df)}")
    print(f"Unique receptor genes tested: {np_map_receptor_df['gene'].nunique()}")
    print(f"Unique systems: {np_map_receptor_df['system'].nunique()}")
    return (np_map_receptor_df,)


@app.cell
def _(np_map_ligand_df):
    # Score LIGANDS by system
    ligand_system_scores = np_map_ligand_df.groupby('system').agg(
        total=('score', 'count'),
        passed=('score', 'sum')
    ).reset_index()
    ligand_system_scores['pct'] = 100 * ligand_system_scores['passed'] / ligand_system_scores['total']
    ligand_system_scores['grade'] = ligand_system_scores['pct'].apply(
        lambda x: 'A' if x >= 80 else 'B' if x >= 60 else 'C' if x >= 40 else 'D' if x >= 20 else 'F'
    )
    ligand_system_scores = ligand_system_scores.sort_values('pct', ascending=False)

    print("=== LIGAND GRADES BY NEUROPEPTIDE SYSTEM ===\n")
    for _, _row in ligand_system_scores.iterrows():
        print(f"{_row['system']}: {_row['passed']}/{_row['total']} ({_row['pct']:.0f}%) - Grade {_row['grade']}")

    # Overall
    _total = np_map_ligand_df['score'].sum()
    _possible = len(np_map_ligand_df)
    print(f"\n=== OVERALL LIGANDS (np_map.csv) ===")
    print(f"Total: {_total}/{_possible} ({100*_total/_possible:.0f}%)")
    return


@app.cell
def _(np_map_receptor_df):
    # Score RECEPTORS by system
    receptor_system_scores = np_map_receptor_df.groupby('system').agg(
        total=('score', 'count'),
        passed=('score', 'sum')
    ).reset_index()
    receptor_system_scores['pct'] = 100 * receptor_system_scores['passed'] / receptor_system_scores['total']
    receptor_system_scores['grade'] = receptor_system_scores['pct'].apply(
        lambda x: 'A' if x >= 80 else 'B' if x >= 60 else 'C' if x >= 40 else 'D' if x >= 20 else 'F'
    )
    receptor_system_scores = receptor_system_scores.sort_values('pct', ascending=False)

    print("=== RECEPTOR GRADES BY NEUROPEPTIDE SYSTEM ===\n")
    for _, _row in receptor_system_scores.iterrows():
        print(f"{_row['system']}: {_row['passed']}/{_row['total']} ({_row['pct']:.0f}%) - Grade {_row['grade']}")

    # Overall
    _total = np_map_receptor_df['score'].sum()
    _possible = len(np_map_receptor_df)
    print(f"\n=== OVERALL RECEPTORS (np_map.csv) ===")
    print(f"Total: {_total}/{_possible} ({100*_total/_possible:.0f}%)")
    return


@app.cell
def _(np_map_ligand_df):
    # Score LIGANDS by gene
    ligand_gene_scores = np_map_ligand_df.groupby('gene').agg(
        total=('score', 'count'),
        passed=('score', 'sum')
    ).reset_index()
    ligand_gene_scores['pct'] = 100 * ligand_gene_scores['passed'] / ligand_gene_scores['total']
    ligand_gene_scores = ligand_gene_scores.sort_values('pct', ascending=False)

    print("=== TOP PERFORMING LIGAND GENES (all expected regions in top 5) ===\n")
    _perfect = ligand_gene_scores[ligand_gene_scores['pct'] == 100]
    print(f"Perfect scores ({len(_perfect)} ligands): {', '.join(_perfect['gene'].tolist())}")

    print("\n=== WORST PERFORMING LIGAND GENES ===\n")
    _worst = ligand_gene_scores[ligand_gene_scores['pct'] < 50].sort_values('pct')
    for _, _row in _worst.head(10).iterrows():
        print(f"{_row['gene']}: {_row['passed']}/{_row['total']} ({_row['pct']:.0f}%)")
    return


@app.cell
def _(np_map_receptor_df):
    # Score RECEPTORS by gene
    receptor_gene_scores = np_map_receptor_df.groupby('gene').agg(
        total=('score', 'count'),
        passed=('score', 'sum')
    ).reset_index()
    receptor_gene_scores['pct'] = 100 * receptor_gene_scores['passed'] / receptor_gene_scores['total']
    receptor_gene_scores = receptor_gene_scores.sort_values('pct', ascending=False)

    print("=== TOP PERFORMING RECEPTOR GENES (all expected regions in top 5) ===\n")
    _perfect = receptor_gene_scores[receptor_gene_scores['pct'] == 100]
    print(f"Perfect scores ({len(_perfect)} receptors): {', '.join(_perfect['gene'].tolist())}")

    print("\n=== WORST PERFORMING RECEPTOR GENES ===\n")
    _worst = receptor_gene_scores[receptor_gene_scores['pct'] < 50].sort_values('pct')
    for _, _row in _worst.head(10).iterrows():
        print(f"{_row['gene']}: {_row['passed']}/{_row['total']} ({_row['pct']:.0f}%)")
    return


@app.cell
def _(np_map_ligand_df):
    # Show failed LIGAND pairs
    _failed = np_map_ligand_df[~np_map_ligand_df['in_top5']].sort_values(['system', 'gene', 'rank'])

    print("=== FAILED LIGAND GENE-REGION PAIRS (not in top 5) ===\n")
    for _system in _failed['system'].unique():
        _sys_failed = _failed[_failed['system'] == _system]
        _pairs = [f"{r['gene']}@{r['expected_region']}(#{r['rank']})" for _, r in _sys_failed.iterrows()]
        print(f"{_system}: {', '.join(_pairs)}")
    return


@app.cell
def _(np_map_receptor_df):
    # Show failed RECEPTOR pairs
    _failed = np_map_receptor_df[~np_map_receptor_df['in_top5']].sort_values(['system', 'gene', 'rank'])

    print("=== FAILED RECEPTOR GENE-REGION PAIRS (not in top 5) ===\n")
    for _system in _failed['system'].unique():
        _sys_failed = _failed[_failed['system'] == _system]
        _pairs = [f"{r['gene']}@{r['expected_region']}(#{r['rank']})" for _, r in _sys_failed.iterrows()]
        print(f"{_system}: {', '.join(_pairs)}")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    ## Validation 2: hypothalamic_nuclei_genes.csv (Manual Annotations)

    Original validation against manually curated gene-region pairs.
    """)
    return


@app.cell
def _(available_genes, expected_df, get_region_ranks, pd, tracked_genes):
    # Score each expected gene-region pair
    results = []
    for _, _row in expected_df.iterrows():
        _area = _row['area']
        gene_field = _row['gene']
        gene_type = _row['type']
        _genes = [g.strip() for g in gene_field.split('/')]
        for _gene in _genes:
            if _gene.startswith('Gabr') or _gene.startswith('Grin'):  # Handle multi-gene entries like "Npy1r/Npy5r" or "Gabra1/Gabrb2/Gabrg2"
                continue
            if _gene not in tracked_genes:
                continue
            if _gene not in available_genes:  # Skip GABA receptors
                continue
            _ranks = get_region_ranks(_gene)
            if not _ranks:
                continue  # Skip if not in tracked genes
            _rank = _ranks.get(_area, 999)
            _in_top5 = _rank <= 5
            _score = 1 if _in_top5 else 0
            results.append({'area': _area, 'gene': _gene, 'type': gene_type, 'expected_rank': _rank, 'in_top5': _in_top5, 'score': _score})  # Skip if not in available expression data
    results_df = pd.DataFrame(results)
    print(f'Testable gene-region pairs: {len(results_df)}')
    print(f'\nBreakdown by area:')
    print(results_df.groupby('area').size())  # Get ranks
    return (results_df,)


@app.cell
def _(results_df):
    # Compute grade per area
    area_grades = results_df.groupby('area').agg(total_tests=('score', 'count'), passed=('score', 'sum')).reset_index()
    area_grades['grade_pct'] = 100 * area_grades['passed'] / area_grades['total_tests']
    area_grades['grade'] = area_grades['grade_pct'].apply(lambda x: 'A' if x >= 80 else 'B' if x >= 60 else 'C' if x >= 40 else 'D' if x >= 20 else 'F')
    print('=== VALIDATION GRADES BY AREA ===\n')
    for _, _row in area_grades.sort_values('grade_pct', ascending=False).iterrows():
        print(f"{_row['area']}: {_row['passed']}/{_row['total_tests']} ({_row['grade_pct']:.0f}%) - Grade {_row['grade']}")
    print(f'\n=== OVERALL ===')
    _total = results_df['score'].sum()
    _possible = len(results_df)
    _overall_pct = 100 * _total / _possible
    print(f'Total: {_total}/{_possible} ({_overall_pct:.0f}%)')
    return


@app.cell
def _(results_df):
    # Show details: what passed and what failed
    print('=== PASSED (in top 5) ===\n')
    passed = results_df[results_df['in_top5']].sort_values(['area', 'expected_rank'])
    for _area in passed['area'].unique():
        area_passed = passed[passed['area'] == _area]
        _genes = [f"{r['gene']}(#{r['expected_rank']})" for _, r in area_passed.iterrows()]
        print(f"{_area}: {', '.join(_genes)}")
    print('\n=== FAILED (not in top 5) ===\n')
    _failed = results_df[~results_df['in_top5']].sort_values(['area', 'expected_rank'])
    for _area in _failed['area'].unique():
        _area_failed = _failed[_failed['area'] == _area]
        _genes = [f"{r['gene']}(#{r['expected_rank']})" for _, r in _area_failed.iterrows()]
        print(f"{_area}: {', '.join(_genes)}")
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    ## Interactive Ranked Heatmaps

    Visualize expression rankings across all regions.
    """)
    return


@app.cell
def _(available_genes, expr_with_region, pd):
    # Get all regions and prepare data for heatmaps
    region_counts = expr_with_region['region'].value_counts()
    top_regions = region_counts.index.tolist()

    # Get cluster profiles for gene classification
    cluster_profiles = pd.read_parquet('../data/processed/mouse_abc/cluster_ligand_receptor_profile.parquet')
    gene_info = cluster_profiles.groupby('gene')[['is_ligand', 'is_receptor']].first()

    # Get ligand and receptor genes
    ligand_genes = sorted([g for g in available_genes if g in gene_info.index and gene_info.loc[g, 'is_ligand']])
    receptor_genes = sorted([g for g in available_genes if g in gene_info.index and gene_info.loc[g, 'is_receptor']])

    print(f"Ligands: {len(ligand_genes)}")
    print(f"Receptors: {len(receptor_genes)}")
    return ligand_genes, receptor_genes, top_regions


@app.cell
def _(get_region_ranks, go, ligand_genes, np, pd, top_regions):
    # Ranked LIGAND heatmap
    rank_data_lig = []
    for _gene in ligand_genes:
        _ranks = get_region_ranks(_gene)
        for _region in top_regions:
            _rank = _ranks.get(_region, len(top_regions) + 1)
            rank_data_lig.append({'gene': _gene, 'region': _region, 'rank': _rank})
    rank_df_lig = pd.DataFrame(rank_data_lig)
    rank_pivot_lig = rank_df_lig.pivot(index='gene', columns='region', values='rank').loc[ligand_genes]
    text_lig = [[f'#{int(rank_pivot_lig.loc[g, r])}' if rank_pivot_lig.loc[g, r] <= 5 else '' for r in rank_pivot_lig.columns] for g in ligand_genes]
    hover_lig = [[f'Ligand: {g}<br>Region: {r}<br>Rank: #{int(rank_pivot_lig.loc[g, r])}' for r in rank_pivot_lig.columns] for g in ligand_genes]
    # Text annotations (top 5 only)
    color_lig = 21 - np.clip(rank_pivot_lig.values, 1, 20)
    fig_lig = go.Figure(data=go.Heatmap(z=color_lig, x=rank_pivot_lig.columns.tolist(), y=ligand_genes, text=text_lig, texttemplate='%{text}', textfont={'size': 9}, hovertext=hover_lig, hoverinfo='text', colorscale='YlOrRd', colorbar=dict(title='Rank', tickvals=[1, 5, 10, 15, 20], ticktext=['#20', '#15', '#10', '#5', '#1']), zmin=1, zmax=20))
    fig_lig.update_layout(title='LIGAND Expression by Region (RANKED)<br><sup>Brighter = higher rank</sup>', xaxis_title='Region', yaxis_title='Ligand', height=max(500, len(ligand_genes) * 20), width=1100, yaxis=dict(tickfont=dict(size=10), autorange='reversed'), xaxis=dict(tickangle=45, tickfont=dict(size=9)))
    # Hover text
    # Color (inverted rank)
    fig_lig.show()
    return


@app.cell
def _(get_region_ranks, go, np, pd, receptor_genes, top_regions):
    # Ranked RECEPTOR heatmap
    rank_data_rec = []
    for _gene in receptor_genes:
        _ranks = get_region_ranks(_gene)
        for _region in top_regions:
            _rank = _ranks.get(_region, len(top_regions) + 1)
            rank_data_rec.append({'gene': _gene, 'region': _region, 'rank': _rank})
    rank_df_rec = pd.DataFrame(rank_data_rec)
    rank_pivot_rec = rank_df_rec.pivot(index='gene', columns='region', values='rank').loc[receptor_genes]
    text_rec = [[f'#{int(rank_pivot_rec.loc[g, r])}' if rank_pivot_rec.loc[g, r] <= 5 else '' for r in rank_pivot_rec.columns] for g in receptor_genes]
    hover_rec = [[f'Receptor: {g}<br>Region: {r}<br>Rank: #{int(rank_pivot_rec.loc[g, r])}' for r in rank_pivot_rec.columns] for g in receptor_genes]
    # Text annotations (top 5 only)
    color_rec = 21 - np.clip(rank_pivot_rec.values, 1, 20)
    fig_rec = go.Figure(data=go.Heatmap(z=color_rec, x=rank_pivot_rec.columns.tolist(), y=receptor_genes, text=text_rec, texttemplate='%{text}', textfont={'size': 9}, hovertext=hover_rec, hoverinfo='text', colorscale='YlGnBu', colorbar=dict(title='Rank', tickvals=[1, 5, 10, 15, 20], ticktext=['#20', '#15', '#10', '#5', '#1']), zmin=1, zmax=20))
    fig_rec.update_layout(title='RECEPTOR Expression by Region (RANKED)<br><sup>Brighter = higher rank</sup>', xaxis_title='Region', yaxis_title='Receptor', height=max(500, len(receptor_genes) * 18), width=1100, yaxis=dict(tickfont=dict(size=10), autorange='reversed'), xaxis=dict(tickangle=45, tickfont=dict(size=9)))
    # Hover text
    # Color (inverted rank)
    fig_rec.show()
    return


@app.cell
def _(mo):
    mo.md(r"""
    ---
    ## Gene Co-expression Analysis

    Analyze correlations between ligand and receptor expression patterns across different granularity levels.

    For each level, we compute the **mean expression** per gene per group, then calculate **Pearson correlations** between gene pairs across groups.
    """)
    return


@app.cell
def _(expr_with_region, np, pd):
    # Define granularity levels for co-expression analysis
    GRANULARITY_LEVELS = {
        'Cell Class (27 groups)': 'class',
        'Cell Subclass (164 groups)': 'subclass',
        'Cell Supertype (483 groups)': 'supertype',
        'Cell Cluster (1,745 groups)': 'cluster',
        'Brain Region (45 groups)': 'region',
    }

    def compute_mean_expression_matrix(genes, groupby_col):
        """Compute mean expression per gene per group."""
        gene_cols = [g for g in genes if g in expr_with_region.columns]
        grouped = expr_with_region.groupby(groupby_col)[gene_cols].mean()
        return grouped  # rows = groups, columns = genes

    def compute_coexpression_matrix(genes_row, genes_col, groupby_col):
        """Compute Pearson correlation matrix between two sets of genes."""
        all_genes = list(set(genes_row) | set(genes_col))
        mean_expr = compute_mean_expression_matrix(all_genes, groupby_col)

        # Filter to genes that exist
        genes_row = [g for g in genes_row if g in mean_expr.columns]
        genes_col = [g for g in genes_col if g in mean_expr.columns]

        # Compute correlation matrix
        corr_matrix = np.zeros((len(genes_row), len(genes_col)))
        for i, g1 in enumerate(genes_row):
            for j, g2 in enumerate(genes_col):
                v1 = mean_expr[g1].values
                v2 = mean_expr[g2].values
                # Handle constant vectors
                if np.std(v1) == 0 or np.std(v2) == 0:
                    corr_matrix[i, j] = 0
                else:
                    corr_matrix[i, j] = np.corrcoef(v1, v2)[0, 1]

        return pd.DataFrame(corr_matrix, index=genes_row, columns=genes_col)
    return GRANULARITY_LEVELS, compute_coexpression_matrix


@app.cell
def _(GRANULARITY_LEVELS, mo):
    # UI control for granularity level
    granularity_dropdown = mo.ui.dropdown(
        options=GRANULARITY_LEVELS,
        value='Brain Region (45 groups)',
        label='Granularity Level'
    )
    mo.md(f"**Select grouping level for co-expression analysis:** {granularity_dropdown}")
    return (granularity_dropdown,)


@app.cell(hide_code=True)
def _(compute_coexpression_matrix, go, granularity_dropdown, ligand_genes):
    # Ligand-Ligand co-expression matrix
    _level = granularity_dropdown.value
    _corr_ll = compute_coexpression_matrix(ligand_genes, ligand_genes, _level)

    _hover_ll = [[f'{g1} vs {g2}<br>r = {_corr_ll.loc[g1, g2]:.3f}'
                  for g2 in _corr_ll.columns] for g1 in _corr_ll.index]

    fig_ll = go.Figure(data=go.Heatmap(
        z=_corr_ll.values,
        x=_corr_ll.columns.tolist(),
        y=_corr_ll.index.tolist(),
        hovertext=_hover_ll,
        hoverinfo='text',
        colorscale='RdBu_r',
        zmid=0,
        zmin=-1,
        zmax=1,
        colorbar=dict(title='Pearson r')
    ))
    fig_ll.update_layout(
        title=f'Ligand-Ligand Co-expression<br><sup>Grouped by {_level}</sup>',
        xaxis_title='Ligand',
        yaxis_title='Ligand',
        height=max(500, len(ligand_genes) * 15),
        width=max(600, len(ligand_genes) * 15),
        yaxis=dict(tickfont=dict(size=9), autorange='reversed'),
        xaxis=dict(tickangle=45, tickfont=dict(size=9))
    )
    fig_ll.show()
    return


@app.cell(hide_code=True)
def _(compute_coexpression_matrix, go, granularity_dropdown, receptor_genes):
    # Receptor-Receptor co-expression matrix
    _level = granularity_dropdown.value
    _corr_rr = compute_coexpression_matrix(receptor_genes, receptor_genes, _level)

    _hover_rr = [[f'{g1} vs {g2}<br>r = {_corr_rr.loc[g1, g2]:.3f}'
                  for g2 in _corr_rr.columns] for g1 in _corr_rr.index]

    fig_rr = go.Figure(data=go.Heatmap(
        z=_corr_rr.values,
        x=_corr_rr.columns.tolist(),
        y=_corr_rr.index.tolist(),
        hovertext=_hover_rr,
        hoverinfo='text',
        colorscale='RdBu_r',
        zmid=0,
        zmin=-1,
        zmax=1,
        colorbar=dict(title='Pearson r')
    ))
    fig_rr.update_layout(
        title=f'Receptor-Receptor Co-expression<br><sup>Grouped by {_level}</sup>',
        xaxis_title='Receptor',
        yaxis_title='Receptor',
        height=max(500, len(receptor_genes) * 12),
        width=max(600, len(receptor_genes) * 12),
        yaxis=dict(tickfont=dict(size=8), autorange='reversed'),
        xaxis=dict(tickangle=45, tickfont=dict(size=8))
    )
    fig_rr.show()
    return


@app.cell(hide_code=True)
def _(
    compute_coexpression_matrix,
    go,
    granularity_dropdown,
    ligand_genes,
    receptor_genes,
):
    # Ligand-Receptor co-expression matrix
    _level = granularity_dropdown.value
    _corr_lr = compute_coexpression_matrix(ligand_genes, receptor_genes, _level)

    _hover_lr = [[f'{g1} (L) vs {g2} (R)<br>r = {_corr_lr.loc[g1, g2]:.3f}'
                  for g2 in _corr_lr.columns] for g1 in _corr_lr.index]

    fig_lr = go.Figure(data=go.Heatmap(
        z=_corr_lr.values,
        x=_corr_lr.columns.tolist(),
        y=_corr_lr.index.tolist(),
        hovertext=_hover_lr,
        hoverinfo='text',
        colorscale='RdBu_r',
        zmid=0,
        zmin=-1,
        zmax=1,
        colorbar=dict(title='Pearson r')
    ))
    fig_lr.update_layout(
        title=f'Ligand-Receptor Co-expression<br><sup>Grouped by {_level}</sup>',
        xaxis_title='Receptor',
        yaxis_title='Ligand',
        height=max(500, len(ligand_genes) * 15),
        width=max(600, len(receptor_genes) * 12),
        yaxis=dict(tickfont=dict(size=9), autorange='reversed'),
        xaxis=dict(tickangle=45, tickfont=dict(size=8))
    )
    fig_lr.show()
    return


@app.cell
def _(mo):
    mo.md(r"""
    ---
    ## Expression Prevalence Heatmaps

    Visualize the **percentage of cells expressing** each gene above threshold, grouped by brain region, cell class, or cell subclass.
    """)
    return


@app.cell
def _(mo):
    # UI control for expression heatmap grouping
    EXPRESSION_GROUPINGS = {
        'Brain Region (45 groups)': 'region',
        'Cell Class (27 groups)': 'class',
        'Cell Subclass (164 groups)': 'subclass',
    }
    expression_grouping_dropdown = mo.ui.dropdown(
        options=EXPRESSION_GROUPINGS,
        value='Brain Region (45 groups)',
        label='Group by'
    )
    mo.md(f"**Select grouping for expression prevalence:** {expression_grouping_dropdown}")
    return (expression_grouping_dropdown,)


@app.cell(hide_code=True)
def _(
    EXPRESSION_THRESHOLD,
    expr_with_region,
    expression_grouping_dropdown,
    go,
    ligand_genes,
    np,
    pd,
):
    # Ligand expression prevalence heatmap
    _groupby = expression_grouping_dropdown.value
    _genes = [g for g in ligand_genes if g in expr_with_region.columns]

    # Compute % expressing per group
    _pct_data = []
    for _gene in _genes:
        _stats = expr_with_region.groupby(_groupby).agg(
            total=(_groupby, 'count'),
            expressing=(_gene, lambda x: (x > EXPRESSION_THRESHOLD).sum())
        )
        _stats['pct'] = 100 * _stats['expressing'] / _stats['total']
        for _group in _stats.index:
            _pct_data.append({'gene': _gene, 'group': _group, 'pct': _stats.loc[_group, 'pct']})

    _pct_df = pd.DataFrame(_pct_data)
    _pct_pivot = _pct_df.pivot(index='gene', columns='group', values='pct').loc[_genes]

    # Sort columns by total expression
    _col_order = _pct_pivot.mean(axis=0).sort_values(ascending=False).index.tolist()
    _pct_pivot = _pct_pivot[_col_order]

    _hover = [[f'{g}<br>{c}<br>{_pct_pivot.loc[g, c]:.1f}% expressing'
               for c in _pct_pivot.columns] for g in _pct_pivot.index]

    fig_lig_pct = go.Figure(data=go.Heatmap(
        z=_pct_pivot.values,
        x=_pct_pivot.columns.tolist(),
        y=_pct_pivot.index.tolist(),
        hovertext=_hover,
        hoverinfo='text',
        colorscale='YlOrRd',
        zmin=0,
        zmax=np.percentile(_pct_pivot.values, 95),
        colorbar=dict(title='% Expressing')
    ))
    fig_lig_pct.update_layout(
        title=f'LIGAND Expression Prevalence<br><sup>% cells above threshold, grouped by {_groupby}</sup>',
        xaxis_title=_groupby.capitalize(),
        yaxis_title='Ligand',
        height=max(500, len(_genes) * 18),
        width=max(800, len(_pct_pivot.columns) * 20),
        yaxis=dict(tickfont=dict(size=9), autorange='reversed'),
        xaxis=dict(tickangle=45, tickfont=dict(size=8))
    )
    fig_lig_pct.show()
    pct_df = _pct_df
    return (pct_df,)


@app.cell
def _(pct_df):
    pct_df[['gene', 'pct']].groupby('gene').mean()
    return


@app.cell(hide_code=True)
def _(
    EXPRESSION_THRESHOLD,
    expr_with_region,
    expression_grouping_dropdown,
    go,
    np,
    pd,
    receptor_genes,
):
    # Receptor expression prevalence heatmap
    _groupby = expression_grouping_dropdown.value
    _genes = [g for g in receptor_genes if g in expr_with_region.columns]

    # Compute % expressing per group
    _pct_data = []
    for _gene in _genes:
        _stats = expr_with_region.groupby(_groupby).agg(
            total=(_groupby, 'count'),
            expressing=(_gene, lambda x: (x > EXPRESSION_THRESHOLD).sum())
        )
        _stats['pct'] = 100 * _stats['expressing'] / _stats['total']
        for _group in _stats.index:
            _pct_data.append({'gene': _gene, 'group': _group, 'pct': _stats.loc[_group, 'pct']})

    _pct_df = pd.DataFrame(_pct_data)
    _pct_pivot = _pct_df.pivot(index='gene', columns='group', values='pct').loc[_genes]

    # Sort columns by total expression
    _col_order = _pct_pivot.mean(axis=0).sort_values(ascending=False).index.tolist()
    _pct_pivot = _pct_pivot[_col_order]

    _hover = [[f'{g}<br>{c}<br>{_pct_pivot.loc[g, c]:.1f}% expressing'
               for c in _pct_pivot.columns] for g in _pct_pivot.index]

    fig_rec_pct = go.Figure(data=go.Heatmap(
        z=_pct_pivot.values,
        x=_pct_pivot.columns.tolist(),
        y=_pct_pivot.index.tolist(),
        hovertext=_hover,
        hoverinfo='text',
        colorscale='YlGnBu',
        zmin=0,
        zmax=np.percentile(_pct_pivot.values, 95),
        colorbar=dict(title='% Expressing')
    ))
    fig_rec_pct.update_layout(
        title=f'RECEPTOR Expression Prevalence<br><sup>% cells above threshold, grouped by {_groupby}</sup>',
        xaxis_title=_groupby.capitalize(),
        yaxis_title='Receptor',
        height=max(500, len(_genes) * 14),
        width=max(800, len(_pct_pivot.columns) * 20),
        yaxis=dict(tickfont=dict(size=8), autorange='reversed'),
        xaxis=dict(tickangle=45, tickfont=dict(size=8))
    )
    fig_rec_pct.show()

    pct_df_receptor = _pct_df
    return (pct_df_receptor,)


@app.cell
def _(pct_df_receptor):
    pct_df_receptor[['gene', 'pct']].groupby('gene').mean()
    return


@app.cell
def _():
    import marimo as mo
    return (mo,)


if __name__ == "__main__":
    app.run()
