import marimo

__generated_with = "0.19.4"
app = marimo.App()


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Wired vs Wireless Connectome Comparison

    This notebook compares two types of hypothalamic connectivity:

    1. **Wired Connectome**: Axonal projections from the Allen Brain Connectivity Atlas (50 viral tracing experiments)
       - **Important**: Most experiments have injections OUTSIDE the hypothalamus, showing **afferent inputs** to hypothalamic regions
       - Only ~8 experiments have hypothalamic injection sites (true internal wiring)
    2. **Wireless Connectome**: Neuropeptide signaling derived from ligand/receptor co-expression (65 systems, 129 L-R pairs)

    **Key methodological choices:**
    - Connection strength: Geometric mean `sqrt(L × R)` where L=ligand expression, R=receptor expression
    - Cluster→Region aggregation: Cell-weighted mean
    - Region granularity: Aggregated to ABC level (44 common regions)
    """)
    return


@app.cell
def _():
    # Constants
    EXPRESSION_THRESHOLD = 3.0  # log2 scale threshold for "expressing"
    return


@app.cell
def _():
    import pandas as pd
    import numpy as np
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    from pathlib import Path
    from scipy import stats

    # Paths
    DATA_DIR = Path('../data')
    PROCESSED_ABC = DATA_DIR / 'processed' / 'mouse_abc'
    PROCESSED_COMMON = DATA_DIR / 'processed' / 'mouse_common'
    GENERATED = DATA_DIR / 'generated' / 'mouse_common'
    return (
        GENERATED,
        PROCESSED_ABC,
        PROCESSED_COMMON,
        go,
        make_subplots,
        np,
        pd,
        stats,
    )


@app.cell
def _(GENERATED, PROCESSED_ABC, PROCESSED_COMMON, pd):
    # Load all datasets
    wired_df = pd.read_csv(PROCESSED_COMMON / 'hypothalamus_connectivity.csv')
    cluster_profiles = pd.read_parquet(PROCESSED_ABC / 'cluster_ligand_receptor_profile.parquet')
    cell_metadata = pd.read_parquet(PROCESSED_ABC / 'cell_metadata.parquet')
    np_map = pd.read_csv(GENERATED / 'np_map.csv')

    # Load injection site metadata for experiments
    exp_metadata = pd.read_csv(PROCESSED_COMMON / 'connectivity_experiment_metadata.csv')

    # Define hypothalamic regions for filtering
    HYPO_REGIONS = {'ARH', 'VMH', 'DMH', 'PVH', 'PVHd', 'LHA', 'SCH', 'SO', 'AHN', 'PH',
                    'TU', 'ME', 'ZI', 'MPO', 'MPN', 'LPO', 'AVPV', 'MEPO', 'MM', 'LM',
                    'SUM', 'TM', 'PMd', 'PMv', 'PVa', 'PVi', 'PVpo'}
    exp_metadata['is_hypothalamic'] = exp_metadata['injection_acronym'].isin(HYPO_REGIONS)

    print(f"Wired connectome: {wired_df.shape[0]} experiments x {wired_df.shape[1]-1} regions")
    print(f"  - Hypothalamic injections: {exp_metadata['is_hypothalamic'].sum()}")
    print(f"  - Non-hypothalamic injections: {(~exp_metadata['is_hypothalamic']).sum()}")
    print(f"Cluster profiles: {cluster_profiles.shape[0]:,} rows, {cluster_profiles['cluster'].nunique()} clusters")
    print(f"Cell metadata: {cell_metadata.shape[0]:,} cells, {cell_metadata['region'].nunique()} regions")
    print(f"NP systems: {np_map['System'].nunique()} systems, {len(np_map)} L-R pairs")
    return cell_metadata, cluster_profiles, exp_metadata, np_map, wired_df


@app.cell
def _(cell_metadata, wired_df):
    # Build region mappings
    conn_regions = set(wired_df.columns) - {'experiment_id'}
    abc_regions = set(cell_metadata['region'].unique())

    # Map connectivity subregions to ABC parent regions
    # e.g., PVHd -> PVH, VMHc -> VMH, AHNa -> AHN
    SUBREGION_TO_PARENT = {}
    for _r in conn_regions - abc_regions:
        # Try various parent lengths
        for _length in [3, 2, 4]:
            _candidate = _r[:_length]
            if _candidate in abc_regions:
                SUBREGION_TO_PARENT[_r] = _candidate
                break
        # Special cases
        if _r.startswith('MM') and 'MM' in abc_regions and _r not in SUBREGION_TO_PARENT:
            SUBREGION_TO_PARENT[_r] = 'MM'
        elif _r.startswith('SUM') and 'SUM' in abc_regions and _r not in SUBREGION_TO_PARENT:
            SUBREGION_TO_PARENT[_r] = 'SUM'

    # Common regions for comparison (exclude HY-unassigned)
    common_regions = sorted((abc_regions & conn_regions) - {'HY-unassigned'})

    print(f"Common regions: {len(common_regions)}")
    print(f"Subregion mappings: {len(SUBREGION_TO_PARENT)}")
    print(f"\nCommon regions: {', '.join(common_regions[:20])}...")
    return SUBREGION_TO_PARENT, common_regions


@app.cell
def _(cell_metadata, common_regions, np, pd):
    # Compute region centroids from cell coordinates, split by hemisphere
    # Allen CCF: x-axis is left-right, midline approximately at x=5.7
    MIDLINE_X = 5.7

    _region_cells = cell_metadata[cell_metadata['region'].isin(common_regions)].copy()
    _region_cells['hemisphere'] = np.where(_region_cells['x'] < MIDLINE_X, 'L', 'R')

    # Compute centroids per region-hemisphere pair
    region_centroids = (
        _region_cells
        .groupby(['region', 'hemisphere'])
        .agg({
            'x': 'mean',
            'y': 'mean',
            'z': 'mean',
            'cell_id': 'count'  # cell count per hemisphere
        })
        .reset_index()
    )
    region_centroids.columns = ['region', 'hemisphere', 'centroid_x', 'centroid_y', 'centroid_z', 'n_cells']

    # Compute mean centroid per region (for sorting)
    region_mean_centroids = (
        _region_cells
        .groupby('region')
        .agg({'x': 'mean', 'y': 'mean', 'z': 'mean'})
        .reset_index()
    )
    region_mean_centroids.columns = ['region', 'mean_x', 'mean_y', 'mean_z']

    # Add medial distance (distance from midline) for M-L sorting
    region_mean_centroids['medial_dist'] = np.abs(region_mean_centroids['mean_x'] - MIDLINE_X)

    print(f"Region-hemisphere centroids computed: {len(region_centroids)}")
    print(f"  Left hemisphere: {(region_centroids['hemisphere'] == 'L').sum()} regions")
    print(f"  Right hemisphere: {(region_centroids['hemisphere'] == 'R').sum()} regions")

    # Build centroid dict for distance computation
    _centroid_dict = {}
    for _, _row in region_centroids.iterrows():
        _key = (_row['region'], _row['hemisphere'])
        _centroid_dict[_key] = {
            'x': _row['centroid_x'],
            'y': _row['centroid_y'],
            'z': _row['centroid_z']
        }

    def _compute_dist(_c1, _c2):
        return np.sqrt(
            (_c1['x'] - _c2['x'])**2 +
            (_c1['y'] - _c2['y'])**2 +
            (_c1['z'] - _c2['z'])**2
        )

    # Compute pairwise distance matrix (use alphabetical order, will be reordered later)
    _regions_alpha = sorted(common_regions)
    _n = len(_regions_alpha)
    _dist_matrix = np.zeros((_n, _n))

    for _i, _r1 in enumerate(_regions_alpha):
        for _j, _r2 in enumerate(_regions_alpha):
            # Get all hemisphere combinations and take minimum distance
            _distances = []
            for _h1 in ['L', 'R']:
                for _h2 in ['L', 'R']:
                    _k1 = (_r1, _h1)
                    _k2 = (_r2, _h2)
                    if _k1 in _centroid_dict and _k2 in _centroid_dict:
                        _distances.append(_compute_dist(_centroid_dict[_k1], _centroid_dict[_k2]))

            if _distances:
                _dist_matrix[_i, _j] = min(_distances)

    region_distances_unsorted = pd.DataFrame(
        _dist_matrix,
        index=_regions_alpha,
        columns=_regions_alpha
    )

    print(f"Distance matrix shape: {region_distances_unsorted.shape}")
    print(f"Distance range: {region_distances_unsorted.values[region_distances_unsorted.values > 0].min():.2f} - {region_distances_unsorted.values.max():.2f} mm")
    return region_distances_unsorted, region_mean_centroids


@app.cell
def _(mo):
    # UI selector for sorting axis
    sort_axis_selector = mo.ui.dropdown(
        options={
            'Posterior → Anterior (z)': 'z',
            'Dorsal → Ventral (y)': 'y',
            'Medial → Lateral (|x - midline|)': 'medial'
        },
        value='Posterior → Anterior (z)',
        label='Sort regions by'
    )
    mo.md(f"**Region sorting for visualizations:** {sort_axis_selector}")
    return (sort_axis_selector,)


@app.cell
def _(region_distances_unsorted, region_mean_centroids, sort_axis_selector):
    # Apply selected sorting to get ordered regions
    _axis = sort_axis_selector.value

    if _axis == 'z':
        _sorted = region_mean_centroids.sort_values('mean_z', ascending=True)
        sort_axis_label = 'P → A'
    elif _axis == 'y':
        _sorted = region_mean_centroids.sort_values('mean_y', ascending=False)  # Higher y = more dorsal
        sort_axis_label = 'D → V'
    else:  # medial
        _sorted = region_mean_centroids.sort_values('medial_dist', ascending=True)
        sort_axis_label = 'M → L'

    regions_sorted = _sorted['region'].tolist()

    # Reorder distance matrix
    region_distances = region_distances_unsorted.reindex(
        index=regions_sorted,
        columns=regions_sorted
    )

    print(f"Regions sorted by: {sort_axis_label}")
    return region_distances, sort_axis_label


@app.cell
def _(go, region_distances, sort_axis_label):
    # Visualize region distance matrix
    _regions = region_distances.index.tolist()

    fig_distances = go.Figure(data=go.Heatmap(
        z=region_distances.values,
        x=_regions,
        y=_regions,
        colorscale='Viridis',
        colorbar=dict(title='Distance (mm)')
    ))
    fig_distances.update_layout(
        title='Region Distance Matrix<br><sup>Minimum distance across L/R hemispheres</sup>',
        xaxis_title=f'Target Region ({sort_axis_label})',
        yaxis_title=f'Source Region ({sort_axis_label})',
        height=700,
        width=800,
        xaxis=dict(tickangle=45, tickfont=dict(size=8)),
        yaxis=dict(tickfont=dict(size=8))
    )
    fig_distances.show()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    ## Section 2: Wired Connectome Processing

    Aggregate subregions to parent regions and visualize the projection density matrix.
    """)
    return


@app.cell
def _(SUBREGION_TO_PARENT, common_regions, exp_metadata, pd, wired_df):
    # Aggregate wired connectivity: sum subregions into parent regions
    wired_aggregated = pd.DataFrame({'experiment_id': wired_df['experiment_id']})

    for _region in common_regions:
        # Find all columns that map to this region
        _cols = [_region] if _region in wired_df.columns else []
        _cols += [_sub for _sub, _parent in SUBREGION_TO_PARENT.items()
                 if _parent == _region and _sub in wired_df.columns]

        if _cols:
            wired_aggregated[_region] = wired_df[_cols].sum(axis=1)
        else:
            wired_aggregated[_region] = 0.0

    # Add injection site info to the aggregated data
    wired_aggregated = wired_aggregated.merge(
        exp_metadata[['experiment_id', 'injection_acronym', 'is_hypothalamic']],
        on='experiment_id',
        how='left'
    )

    # Average across experiments for canonical connectivity
    wired_avg = wired_aggregated[common_regions].mean()

    print(f"Aggregated wired matrix: {wired_aggregated.shape}")
    print(f"\nTop 10 regions by average incoming projection density:")
    print(wired_avg.sort_values(ascending=False).head(10))
    return wired_aggregated, wired_avg


@app.cell
def _(go, np, region_distances, sort_axis_label, wired_aggregated):
    # Plot: Wired connectivity heatmap (experiments x regions)
    # Get spatially-sorted region order from distance matrix
    _sorted_regions = region_distances.index.tolist()

    # Sort experiments by hypothalamic vs non-hypothalamic, then by injection site
    _sorted = wired_aggregated.sort_values(['is_hypothalamic', 'injection_acronym'], ascending=[False, True])

    # Reorder columns to spatially-sorted regions
    _regions_in_data = [r for r in _sorted_regions if r in _sorted.columns]
    _data = _sorted[_regions_in_data].values
    _data_log = np.log10(_data + 1e-6)  # Log scale for better visibility

    # Create labels showing injection site (mark hypothalamic ones with *)
    _labels = [
        f"{'*' if row['is_hypothalamic'] else ''}{row['injection_acronym']}"
        for _, row in _sorted.iterrows()
    ]

    fig_wired = go.Figure(data=go.Heatmap(
        z=_data_log,
        x=_regions_in_data,
        y=_labels,
        colorscale='Blues',
        colorbar=dict(title='log10(Proj. Density)')
    ))
    fig_wired.update_layout(
        title='Wired Connectome: Projections TO Hypothalamus by Injection Site',
        xaxis_title=f'Target Region ({sort_axis_label})',
        yaxis_title='Injection Site',
        height=800,
        width=1100,
        xaxis=dict(tickangle=45, tickfont=dict(size=9)),
        yaxis=dict(tickfont=dict(size=9))
    )
    fig_wired.show()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    ## Section 3: Wireless Connectome Derivation

    Aggregate cluster-level ligand/receptor expression to region level using cell-weighted means, then compute connection strengths.

    **Distance weighting**: Neuropeptide signals diffuse through extracellular space, so nearby regions receive stronger signals than distant ones. Choose a distance model below.
    """)
    return


@app.cell
def _(mo):
    # UI controls for distance-based weighting
    distance_mode = mo.ui.dropdown(
        options=['None', 'Threshold', 'Exponential decay'],
        value='Threshold',
        label='Distance mode'
    )

    distance_threshold = mo.ui.slider(
        start=0.2,
        stop=2.0,
        step=0.2,
        value=.6,
        label='Distance threshold (mm)'
    )

    decay_constant = mo.ui.slider(
        start=0.1,
        stop=2.0,
        step=0.1,
        value=0.2,
        label='Decay constant λ (mm)'
    )

    mo.md(f"""
    **Distance Weighting Controls**

    {distance_mode}

    - **None**: No distance correction (all connections equally weighted)
    - **Threshold**: Zero out connections beyond threshold distance
    - **Exponential decay**: Weight = exp(-distance / λ)

    {mo.hstack([distance_threshold, decay_constant])}

    *Typical hypothalamic dimensions: ~2-3 mm across*
    """)
    return decay_constant, distance_mode, distance_threshold


@app.cell
def _(cell_metadata):
    # Compute cluster-region cell counts for weighted aggregation
    # Note: use 'n_cells_region' to avoid conflict with 'n_cells' in cluster_profiles
    cluster_region_counts = (
        cell_metadata[cell_metadata['region'] != 'HY-unassigned']
        .groupby(['cluster', 'region'])
        .size()
        .reset_index(name='n_cells_region')
    )

    # Total cells per cluster
    cluster_totals = cluster_region_counts.groupby('cluster')['n_cells_region'].sum().reset_index()
    cluster_totals.columns = ['cluster', 'total_cells']

    print(f"Cluster-region pairs: {len(cluster_region_counts):,}")
    print(f"Clusters with region info: {len(cluster_totals):,}")
    return (cluster_region_counts,)


@app.cell
def _(cluster_profiles, cluster_region_counts, np):
    def compute_region_expression(profiles, cr_counts):
        """
        Aggregate cluster-level expression to region-level using cell-weighted mean.

        For each (gene, region):
            region_expr = sum(cluster_expr * cluster_cells_in_region) / total_cells_in_region
        """
        # Filter out inf values in expression data
        profiles_clean = profiles[~np.isinf(profiles['mean_expr'])].copy()

        # Join cluster profiles with region counts
        joined = profiles_clean.merge(cr_counts, on='cluster', how='inner')

        # Weight expression by cell count in region
        joined['weighted_expr'] = joined['mean_expr'] * joined['n_cells_region']

        # Aggregate to region level
        region_expr = (
            joined
            .groupby(['region', 'gene', 'is_ligand', 'is_receptor'])
            .agg({
                'weighted_expr': 'sum',
                'n_cells_region': 'sum'
            })
            .reset_index()
        )
        region_expr['region_mean_expr'] = region_expr['weighted_expr'] / region_expr['n_cells_region']

        return region_expr

    # Report inf values filtered
    n_inf = np.sum(np.isinf(cluster_profiles['mean_expr']))
    if n_inf > 0:
        print(f"Note: Filtered {n_inf:,} inf values from cluster profiles")

    region_expression = compute_region_expression(cluster_profiles, cluster_region_counts)
    print(f"Region-gene expression pairs: {len(region_expression):,}")
    return (region_expression,)


@app.cell
def _(
    common_regions,
    decay_constant,
    distance_mode,
    distance_threshold,
    np,
    np_map,
    pd,
    region_distances,
    region_expression,
):
    def compute_wireless_connectome(region_expr, lr_map, regions, dist_matrix, mode, threshold, decay):
        """
        Compute wireless connectivity matrix for each L-R pair with distance weighting.

        Connection strength = sqrt(L_source * R_target) * distance_weight

        Distance weighting modes:
        - None: weight = 1 (no distance correction)
        - Threshold: weight = 1 if distance <= threshold, else 0
        - Exponential decay: weight = exp(-distance / decay)

        Returns dict: system -> DataFrame with source, target, strength
        """
        connectomes = {}

        for _system in lr_map['System'].unique():
            _system_pairs = lr_map[lr_map['System'] == _system]
            _system_connections = []

            for _, _row in _system_pairs.iterrows():
                _ligand = _row['Ligand_Gene']
                _receptor = _row['Receptor_Gene']

                # Handle heterodimer receptors (e.g., "Calcrl;Ramp2")
                _receptor_genes = [g.strip() for g in _receptor.split(';')]

                # Get ligand expression per region
                _lig_data = region_expr[
                    (region_expr['gene'] == _ligand) & region_expr['is_ligand']
                ][['region', 'region_mean_expr']].rename(
                    columns={'region_mean_expr': 'ligand_expr'}
                )

                if _lig_data.empty:
                    continue

                # Get receptor expression (min across subunits for heterodimers)
                _rec_dfs = []
                for _rg in _receptor_genes:
                    _rec = region_expr[
                        (region_expr['gene'] == _rg) & region_expr['is_receptor']
                    ][['region', 'region_mean_expr']]
                    if not _rec.empty:
                        _rec_dfs.append(_rec.set_index('region')['region_mean_expr'])

                if not _rec_dfs:
                    continue

                # For heterodimers, take minimum (AND logic)
                if len(_rec_dfs) == 1:
                    _rec_series = _rec_dfs[0]
                else:
                    _rec_combined = pd.concat(_rec_dfs, axis=1)
                    _rec_series = _rec_combined.min(axis=1)

                _rec_data = _rec_series.reset_index()
                _rec_data.columns = ['region', 'receptor_expr']

                # Compute connections for all source-target pairs
                for _src in regions:
                    _src_lig = _lig_data[_lig_data['region'] == _src]
                    if _src_lig.empty:
                        continue
                    _L = max(0, _src_lig['ligand_expr'].values[0])

                    for _tgt in regions:
                        _tgt_rec = _rec_data[_rec_data['region'] == _tgt]
                        if _tgt_rec.empty:
                            continue
                        _R = max(0, _tgt_rec['receptor_expr'].values[0])

                        # Compute base strength (geometric mean)
                        _base_strength = np.sqrt(_L * _R)
                        if _base_strength <= 0:
                            continue

                        # Get distance between regions
                        _dist = 0.0
                        if _src in dist_matrix.index and _tgt in dist_matrix.columns:
                            _dist = dist_matrix.loc[_src, _tgt]

                        # Apply distance weighting
                        if mode == 'None':
                            _weight = 1.0
                        elif mode == 'Threshold':
                            _weight = 1.0 if _dist <= threshold else 0.0
                        else:  # Exponential decay
                            _weight = np.exp(-_dist / decay) if decay > 0 else 0.0

                        _strength = _base_strength * _weight
                        if _strength > 0:
                            _system_connections.append({
                                'source': _src,
                                'target': _tgt,
                                'ligand': _ligand,
                                'receptor': _receptor,
                                'ligand_expr': _L,
                                'receptor_expr': _R,
                                'distance': _dist,
                                'distance_weight': _weight,
                                'strength': _strength
                            })

            if _system_connections:
                connectomes[_system] = pd.DataFrame(_system_connections)

        return connectomes

    # Get UI values
    _mode = distance_mode.value
    _threshold = distance_threshold.value
    _decay = decay_constant.value

    wireless_by_system = compute_wireless_connectome(
        region_expression, np_map, common_regions,
        region_distances, _mode, _threshold, _decay
    )
    print(f"Distance mode: {_mode}" +
          (f" (threshold={_threshold:.1f}mm)" if _mode == 'Threshold' else "") +
          (f" (λ={_decay:.2f}mm)" if _mode == 'Exponential decay' else ""))
    print(f"Systems with wireless connectomes: {len(wireless_by_system)}")
    return (wireless_by_system,)


@app.cell
def _(pd, wireless_by_system):
    # Aggregate wireless connectome across all systems
    _all_connections = []
    for _system, _df in wireless_by_system.items():
        _df_copy = _df.copy()
        _df_copy['system'] = _system
        _all_connections.append(_df_copy)

    wireless_all = pd.concat(_all_connections, ignore_index=True) if _all_connections else pd.DataFrame()

    # Sum strengths per source-target pair
    wireless_total = (
        wireless_all
        .groupby(['source', 'target'])
        .agg({
            'strength': 'sum',
            'system': 'count'
        })
        .reset_index()
        .rename(columns={'strength': 'total_strength', 'system': 'n_pairs'})
    )

    # Pivot to matrix form
    wireless_matrix = wireless_total.pivot(
        index='source', columns='target', values='total_strength'
    ).fillna(0)

    print(f"Total wireless connections: {len(wireless_total):,}")
    print(f"Wireless matrix shape: {wireless_matrix.shape}")
    return (wireless_matrix,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    ## Section 4: Sanity Check Visualizations

    Verify that the cluster-to-region aggregation produces sensible results.
    """)
    return


@app.cell
def _(cluster_profiles, cluster_region_counts, go):
    # Sanity check: Hcrt (orexin) should be highest in LHA
    _test_gene = 'Hcrt'

    # Get cluster-level expression (n_cells is cluster total cells from cluster_profiles)
    _cluster_hcrt = cluster_profiles[cluster_profiles['gene'] == _test_gene][['cluster', 'mean_expr', 'n_cells']]

    # Get primary region for each cluster (based on where most cells are)
    _primary_regions = cluster_region_counts.loc[
        cluster_region_counts.groupby('cluster')['n_cells_region'].idxmax()
    ][['cluster', 'region']].rename(columns={'region': 'primary_region'})

    _cluster_with_region = _cluster_hcrt.merge(_primary_regions, on='cluster', how='inner')

    # Plot
    fig_sanity = go.Figure()
    _key_regions = ['LHA', 'ARH', 'VMH', 'DMH', 'PVH', 'ZI', 'SCH']
    _colors = ['red', 'blue', 'green', 'orange', 'purple', 'brown', 'pink']

    for _region, _color in zip(_key_regions, _colors):
        _subset = _cluster_with_region[_cluster_with_region['primary_region'] == _region]
        if len(_subset) > 0:
            fig_sanity.add_trace(go.Scatter(
                x=_subset['n_cells'],
                y=_subset['mean_expr'],
                mode='markers',
                name=f"{_region} ({len(_subset)} clusters)",
                marker=dict(color=_color, opacity=0.6, size=8)
            ))

    fig_sanity.update_layout(
        title=f'Sanity Check: {_test_gene} (Orexin) Expression by Region<br><sup>Each point is a cluster; Hcrt should be high in LHA clusters</sup>',
        xaxis_title='Cells in Cluster',
        yaxis_title='Mean Expression (log2)',
        xaxis_type='log',
        height=500,
        width=800
    )
    fig_sanity.show()
    return


@app.cell
def _(common_regions, go, region_expression):
    # Region-level expression bar chart for key ligands
    _key_ligands = ['Hcrt', 'Pomc', 'Agrp', 'Oxt', 'Avp', 'Npy', 'Crh']

    fig_ligands = go.Figure()

    for _ligand in _key_ligands:
        _lig_data = region_expression[
            (region_expression['gene'] == _ligand) &
            region_expression['is_ligand'] &
            region_expression['region'].isin(common_regions)
        ]
        if not _lig_data.empty:
            _lig_sorted = _lig_data.sort_values('region_mean_expr', ascending=False).head(10)
            fig_ligands.add_trace(go.Bar(
                x=_lig_sorted['region'],
                y=_lig_sorted['region_mean_expr'],
                name=_ligand
            ))

    fig_ligands.update_layout(
        title='Top 10 Regions by Ligand Expression<br><sup>Validation: Hcrt in LHA, Pomc/Agrp in ARH, Oxt/Avp in PVH</sup>',
        xaxis_title='Region',
        yaxis_title='Mean Expression (log2)',
        barmode='group',
        height=500,
        width=1000
    )
    fig_ligands.show()
    return


@app.cell
def _(mo, wireless_by_system):
    # UI dropdown for NP system selection
    system_options = sorted(wireless_by_system.keys())
    system_selector = mo.ui.dropdown(
        options=system_options,
        value='Orexin' if 'Orexin' in system_options else system_options[0],
        label='Select NP System'
    )
    mo.md(f"**Select neuropeptide system for detailed view:** {system_selector}")
    return (system_selector,)


@app.cell
def _(
    go,
    pd,
    region_distances,
    sort_axis_label,
    system_selector,
    wireless_by_system,
):
    # Per-system wireless heatmap
    _system = system_selector.value
    _df = wireless_by_system.get(_system, pd.DataFrame())

    # Get spatially-sorted region order from distance matrix
    _sorted_regions = region_distances.index.tolist()

    if not _df.empty:
        # Aggregate across L-R pairs within system
        _agg = _df.groupby(['source', 'target'])['strength'].sum().reset_index()
        _matrix = _agg.pivot(index='source', columns='target', values='strength').fillna(0)

        # Reindex to spatially-sorted order (only include regions present in matrix)
        _regions_in_matrix = [r for r in _sorted_regions if r in _matrix.index and r in _matrix.columns]
        _matrix = _matrix.reindex(index=_regions_in_matrix, columns=_regions_in_matrix, fill_value=0)

        _z = _matrix.values
        _sources = _matrix.index.tolist()
        _targets = _matrix.columns.tolist()

        fig_system = go.Figure(data=go.Heatmap(
            z=_z,
            x=_targets,
            y=_sources,
            colorscale='Reds',
            colorbar=dict(title='sqrt(L*R)')
        ))
        fig_system.update_layout(
            title=f'Wireless Connectome: {_system} System',
            xaxis_title=f'Target Region ({sort_axis_label})',
            yaxis_title=f'Source Region ({sort_axis_label})',
            height=600,
            width=700,
            xaxis=dict(tickangle=45, tickfont=dict(size=9)),
            yaxis=dict(tickfont=dict(size=9))
        )
    else:
        # Create empty figure with message
        fig_system = go.Figure()
        fig_system.add_annotation(
            text=f"No data available for {_system}",
            xref="paper", yref="paper",
            x=0.5, y=0.5, showarrow=False,
            font=dict(size=16)
        )
        fig_system.update_layout(height=400, width=600)

    fig_system.show()
    return


@app.cell
def _(
    go,
    make_subplots,
    region_distances,
    sort_axis_label,
    wireless_by_system,
):
    # Multi-system comparison panel
    _key_systems = ['Orexin', 'Neuropeptide Y', 'Oxytocin', 'CRH', 'Pomc', 'GHRH']
    _key_systems = [s for s in _key_systems if s in wireless_by_system]

    # Get spatially-sorted region order
    _sorted_regions = region_distances.index.tolist()

    _n_cols = 3
    _n_rows = (len(_key_systems) + _n_cols - 1) // _n_cols

    fig_multi = make_subplots(
        rows=_n_rows, cols=_n_cols,
        subplot_titles=_key_systems,
        horizontal_spacing=0.08,
        vertical_spacing=0.15
    )

    for _idx, _system in enumerate(_key_systems):
        _row = _idx // _n_cols + 1
        _col = _idx % _n_cols + 1

        _df = wireless_by_system[_system]
        _agg = _df.groupby(['source', 'target'])['strength'].sum().reset_index()
        _matrix = _agg.pivot(index='source', columns='target', values='strength').fillna(0)

        # Reindex to spatially-sorted order
        _regions_in_matrix = [r for r in _sorted_regions if r in _matrix.index and r in _matrix.columns]
        _matrix = _matrix.reindex(index=_regions_in_matrix, columns=_regions_in_matrix, fill_value=0)

        fig_multi.add_trace(
            go.Heatmap(
                z=_matrix.values,
                x=_matrix.columns.tolist(),
                y=_matrix.index.tolist(),
                colorscale='Reds',
                showscale=False
            ),
            row=_row, col=_col
        )

    fig_multi.update_layout(
        title=f'Wireless Connectomes: Key Neuropeptide Systems<br><sup>Regions sorted {sort_axis_label}</sup>',
        height=300 * _n_rows,
        width=1100
    )
    fig_multi.show()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    ## Section 5: Wired vs Wireless Comparison

    Compare the two connectome modalities to identify regions with concordant or discordant connectivity patterns.
    """)
    return


@app.cell
def _(common_regions, np, wired_avg, wireless_matrix):
    # Prepare comparison vectors
    # Wired: average incoming projection density per region
    wired_incoming = wired_avg[common_regions].values

    # Wireless: total incoming signal per region (sum of column)
    _wireless_regions_in_matrix = [r for r in common_regions if r in wireless_matrix.columns]
    wireless_incoming = np.array([
        wireless_matrix[r].sum() if r in wireless_matrix.columns else 0
        for r in common_regions
    ])

    # Wireless: total outgoing signal per region (sum of row)
    wireless_outgoing = np.array([
        wireless_matrix.loc[r].sum() if r in wireless_matrix.index else 0
        for r in common_regions
    ])

    print("Comparison vectors computed")
    print(f"Wired incoming range: {wired_incoming.min():.4f} - {wired_incoming.max():.4f}")
    print(f"Wireless incoming range: {wireless_incoming.min():.2f} - {wireless_incoming.max():.2f}")
    return wired_incoming, wireless_incoming


@app.cell
def _(common_regions, go, np, stats, wired_incoming, wireless_incoming):
    # Correlation analysis
    # Filter out zeros for correlation
    _mask = (wired_incoming > 0) & (wireless_incoming > 0)
    r_value, p_value = stats.pearsonr(wired_incoming[_mask], wireless_incoming[_mask])

    # Spearman (rank) correlation
    r_spearman, p_spearman = stats.spearmanr(wired_incoming, wireless_incoming)

    fig_corr = go.Figure()

    # Main scatter
    fig_corr.add_trace(go.Scatter(
        x=wired_incoming,
        y=wireless_incoming,
        mode='markers+text',
        text=common_regions,
        textposition='top center',
        textfont=dict(size=8),
        marker=dict(size=10, color='steelblue', opacity=0.7),
        name='Regions'
    ))

    # Add trend line
    if np.sum(_mask) > 2:
        _slope, _intercept = np.polyfit(wired_incoming[_mask], wireless_incoming[_mask], 1)
        _x_line = np.linspace(wired_incoming.min(), wired_incoming.max(), 100)
        _y_line = _slope * _x_line + _intercept
        fig_corr.add_trace(go.Scatter(
            x=_x_line,
            y=_y_line,
            mode='lines',
            line=dict(color='red', dash='dash'),
            name=f'Trend (r={r_value:.3f})'
        ))

    fig_corr.update_layout(
        title=f'Wired vs Wireless Connectivity by Region<br><sup>Pearson r={r_value:.3f} (p={p_value:.2e}), Spearman \u03c1={r_spearman:.3f}</sup>',
        xaxis_title='Wired: Avg Incoming Projection Density',
        yaxis_title='Wireless: Total Incoming Signal',
        height=600,
        width=800
    )
    fig_corr.show()
    return r_spearman, r_value


@app.cell
def _(common_regions, go, pd, wired_incoming, wireless_incoming):
    # Create comparison DataFrame with ranks
    comparison_df = pd.DataFrame({
        'region': common_regions,
        'wired': wired_incoming,
        'wireless': wireless_incoming
    })

    # Compute ranks (1 = highest connectivity)
    comparison_df['wired_rank'] = comparison_df['wired'].rank(ascending=False)
    comparison_df['wireless_rank'] = comparison_df['wireless'].rank(ascending=False)
    comparison_df['rank_diff'] = abs(comparison_df['wired_rank'] - comparison_df['wireless_rank'])

    # Bump chart
    fig_bump = go.Figure()

    for _, _row in comparison_df.iterrows():
        _color = 'red' if _row['rank_diff'] >= 10 else 'gray'
        _width = 2 if _row['rank_diff'] >= 10 else 1

        fig_bump.add_trace(go.Scatter(
            x=['Wired Rank', 'Wireless Rank'],
            y=[_row['wired_rank'], _row['wireless_rank']],
            mode='lines+markers+text',
            text=[_row['region'], ''],
            textposition='middle left',
            textfont=dict(size=9),
            line=dict(width=_width, color=_color),
            marker=dict(size=6),
            showlegend=False
        ))

    fig_bump.update_layout(
        title='Region Ranking: Wired vs Wireless Connectivity<br><sup>Red lines = rank difference >= 10</sup>',
        yaxis=dict(autorange='reversed', title='Rank (1 = highest)', dtick=5),
        xaxis=dict(tickfont=dict(size=12)),
        height=700,
        width=600
    )
    fig_bump.show()
    return (comparison_df,)


@app.cell
def _(comparison_df, mo):
    # Display regions with largest rank differences
    discordant = comparison_df.sort_values('rank_diff', ascending=False).head(15)

    mo.md(f"""
    ### Regions with Largest Wired/Wireless Rank Differences

    These regions show the most divergent connectivity patterns between axonal projections and neuropeptide signaling.

    | Region | Wired Rank | Wireless Rank | Rank Diff | Wired Value | Wireless Value |
    |--------|------------|---------------|-----------|-------------|----------------|
    """ + '\n'.join([
        f"| {_r['region']} | {int(_r['wired_rank'])} | {int(_r['wireless_rank'])} | {int(_r['rank_diff'])} | {_r['wired']:.4f} | {_r['wireless']:.2f} |"
        for _, _r in discordant.iterrows()
    ]))
    return


@app.cell
def _(comparison_df, mo):
    # Display concordant regions
    concordant = comparison_df[comparison_df['rank_diff'] <= 3].sort_values('wired_rank')

    mo.md(f"""
    ### Concordant Regions (Rank Difference <= 3)

    These regions have similar connectivity rankings in both wired and wireless modalities.

    | Region | Wired Rank | Wireless Rank | Rank Diff |
    |--------|------------|---------------|-----------|
    """ + '\n'.join([
        f"| {_r['region']} | {int(_r['wired_rank'])} | {int(_r['wireless_rank'])} | {int(_r['rank_diff'])} |"
        for _, _r in concordant.iterrows()
    ]))
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    ## Section 6: Per-System Correlations

    Examine which neuropeptide systems most closely mirror the wired connectome.
    """)
    return


@app.cell
def _(common_regions, go, np, pd, stats, wired_incoming, wireless_by_system):
    # Compute correlation for each NP system
    _system_correlations = []

    for _system, _df in wireless_by_system.items():
        # Aggregate to region-level incoming signal
        _agg = _df.groupby(['source', 'target'])['strength'].sum().reset_index()
        _matrix = _agg.pivot(index='source', columns='target', values='strength').fillna(0)

        # Compute incoming signal per region
        _sys_incoming = np.array([
            _matrix[_r].sum() if _r in _matrix.columns else 0
            for _r in common_regions
        ])

        # Correlation with wired
        if np.sum(_sys_incoming > 0) >= 3:
            _r_val, _p_val = stats.spearmanr(wired_incoming, _sys_incoming)
            _system_correlations.append({
                'system': _system,
                'spearman_r': _r_val,
                'p_value': _p_val,
                'n_connections': len(_df)
            })

    system_corr_df = pd.DataFrame(_system_correlations).sort_values('spearman_r', ascending=False)

    # Plot
    fig_sys_corr = go.Figure(data=go.Bar(
        x=system_corr_df['system'].head(20),
        y=system_corr_df['spearman_r'].head(20),
        marker_color=['green' if _val > 0 else 'red' for _val in system_corr_df['spearman_r'].head(20)]
    ))
    fig_sys_corr.update_layout(
        title='Correlation with Wired Connectome by NP System<br><sup>Spearman \u03c1 between system-specific wireless and wired connectivity</sup>',
        xaxis_title='Neuropeptide System',
        yaxis_title='Spearman \u03c1',
        height=500,
        width=1000,
        xaxis=dict(tickangle=45)
    )
    fig_sys_corr.show()

    print("Top 10 systems most correlated with wired connectome:")
    print(system_corr_df.head(10).to_string(index=False))
    return (system_corr_df,)


@app.cell(hide_code=True)
def _(
    comparison_df,
    mo,
    r_spearman,
    r_value,
    system_corr_df,
    wired_aggregated,
    wireless_by_system,
):
    # Summary statistics
    _n_systems = len(wireless_by_system)
    _n_regions = len(comparison_df)
    _n_concordant = len(comparison_df[comparison_df['rank_diff'] <= 3])
    _n_discordant = len(comparison_df[comparison_df['rank_diff'] >= 10])

    _top_corr_systems = ', '.join(system_corr_df.head(5)['system'].tolist())
    _anti_corr_systems = ', '.join(system_corr_df.tail(3)['system'].tolist())

    mo.md(f"""
    ---
    # Summary: Wired vs Wireless Connectome Comparison

    ## Dataset Overview
    - **Wired experiments**: {len(wired_aggregated)} viral tracing experiments
    - **NP systems analyzed**: {_n_systems}
    - **Regions compared**: {_n_regions}

    ## Overall Correlation
    - **Pearson r**: {r_value:.3f}
    - **Spearman \u03c1**: {r_spearman:.3f}

    ## Region Concordance
    - **Concordant regions** (rank diff <= 3): {_n_concordant} ({100*_n_concordant/_n_regions:.0f}%)
    - **Discordant regions** (rank diff >= 10): {_n_discordant} ({100*_n_discordant/_n_regions:.0f}%)

    ## NP Systems Most Similar to Wired Connectome
    {_top_corr_systems}

    ## NP Systems Least Similar (or Anti-correlated)
    {_anti_corr_systems}

    ## Interpretation

    The correlation between wired and wireless connectomes suggests that **neuropeptide signaling
    patterns partially mirror axonal connectivity patterns**. This could indicate:

    1. Regions with high axonal innervation also tend to receive more neuropeptide signals
    2. There may be common organizational principles governing both communication modes
    3. Discordant regions may rely more heavily on one communication mode over the other

    **Regions with high wired but low wireless connectivity** may rely more on fast synaptic
    transmission, while **regions with high wireless but low wired connectivity** may be
    specialized for volume transmission and neuromodulation.
    """)
    return


@app.cell
def _():
    import marimo as mo
    return (mo,)


if __name__ == "__main__":
    app.run()
