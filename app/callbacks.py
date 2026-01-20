"""Dash callbacks for HypoMap Coronal Atlas Viewer."""

from dash import Input, Output, State, ctx, callback
from dash import html
import dash_bootstrap_components as dbc
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import plotly.express as px
import pandas as pd
import numpy as np
from scipy.spatial import cKDTree


# Extended color palette for many categories
EXTENDED_COLORS = (
    px.colors.qualitative.Dark24 +
    px.colors.qualitative.Light24 +
    px.colors.qualitative.Alphabet
)


def get_color_for_category(categories, category):
    """Get consistent color for a category."""
    if category not in categories:
        return '#808080'
    idx = list(categories).index(category)
    return EXTENDED_COLORS[idx % len(EXTENDED_COLORS)]


def create_slice_figure(
    cells_df,
    slices,
    mode,
    cluster_level,
    nt_system,
    np_system,
    threshold,
    display_options,
    point_size,
    subsample_pct,
    diffusion_enabled,
    diffusion_range,
    nt_mapping,
    np_systems,
    cluster_expression,
    cluster_np_expression,
    region_boundaries,
    region_centroids,
    region_colors,
):
    """Create the 2x9 grid of coronal slices.

    PERFORMANCE NOTE: This function uses several optimizations to achieve <1s render times:
    1. Build traces as plain dicts, not go.Scatter objects (avoids Plotly validation overhead)
    2. Use fig.add_traces() to batch-add all traces at once
    3. Use Scattergl (WebGL) instead of Scatter for cell points
    4. Combine all region boundaries into single trace per slice
    5. Add all annotations via layout dict, not individual add_annotation() calls

    Without these optimizations, render time was 8-10 seconds due to Plotly's slow
    add_trace() calls with subplots.
    """
    # Subsampling is handled per-mode below (NT mode only subsamples background)
    n_slices = len(slices)
    n_cols = 2
    n_rows = 9

    # Create subplot grid (no subplot_titles - we'll add annotations instead)
    fig = make_subplots(
        rows=n_rows,
        cols=n_cols,
        horizontal_spacing=0.02,
        vertical_spacing=0.01,
    )

    show_boundaries = 'show_boundaries' in display_options
    show_labels = 'show_labels' in display_options

    # Get unique categories for consistent coloring
    if mode == 'cluster':
        unique_categories = sorted(cells_df[cluster_level].dropna().unique())
    else:
        unique_categories = []

    # Precompute cell colors based on mode
    if mode == 'cluster':
        # Subsample all cells for cluster mode
        if subsample_pct < 100:
            cells_df = cells_df.sample(frac=subsample_pct / 100, random_state=42)

        # Color by cluster level
        color_map = {cat: get_color_for_category(unique_categories, cat)
                     for cat in unique_categories}
        cells_df = cells_df.copy()
        cells_df['_color'] = cells_df[cluster_level].map(color_map).fillna('#808080')
        cells_df['_legend_group'] = cells_df[cluster_level]

    elif mode == 'nt':
        # Gray all cells, black for selected NT type
        cells_df = cells_df.copy()

        # Map cluster to NT type (direct lookup is fast)
        cells_df['_nt_type'] = cells_df['cluster'].map(nt_mapping).fillna('NA')

        # Only subsample background (non-selected) cells, keep all foreground cells
        foreground_mask = cells_df['_nt_type'] == nt_system
        foreground_df = cells_df[foreground_mask]
        background_df = cells_df[~foreground_mask]

        if subsample_pct < 100:
            background_df = background_df.sample(frac=subsample_pct / 100, random_state=42)

        # Background first, foreground last so active cells render on top
        cells_df = pd.concat([background_df, foreground_df], ignore_index=True)

        cells_df['_color'] = np.where(
            cells_df['_nt_type'] == nt_system,
            '#000000',  # Black for selected
            '#D3D3D3'   # Light gray for others
        )
        cells_df['_legend_group'] = np.where(
            cells_df['_nt_type'] == nt_system,
            nt_system,
            'Other'
        )

    elif mode == 'np':
        # Color by neuropeptide system expression using precomputed lookup
        cells_df = cells_df.copy()

        if np_system and np_system in cluster_np_expression:
            # Use precomputed cluster-system expression lookup
            system_lookup = cluster_np_expression[np_system]

            def get_np_color_group(cluster_name):
                """Get color and legend group for a cluster based on expression."""
                if cluster_name not in system_lookup:
                    return '#D3D3D3', 'Neither'

                ligand_expr, receptor_expr = system_lookup[cluster_name]

                # Apply threshold filter
                ligand_above = ligand_expr >= threshold
                receptor_above = receptor_expr >= threshold

                if ligand_above and receptor_above:
                    # Both - blend to purple, intensity by average
                    intensity = min(1.0, (ligand_expr + receptor_expr) / 10)
                    r = int(128 + 127 * intensity)
                    b = int(128 + 127 * intensity)
                    return f'rgb({r}, 0, {b})', 'Both'
                elif ligand_above:
                    # Red for ligand
                    intensity = min(1.0, ligand_expr / 5)
                    r = int(100 + 155 * intensity)
                    return f'rgb({r}, 50, 50)', 'Ligand'
                elif receptor_above:
                    # Blue for receptor
                    intensity = min(1.0, receptor_expr / 5)
                    b = int(100 + 155 * intensity)
                    return f'rgb(50, 50, {b})', 'Receptor'
                else:
                    return '#D3D3D3', 'Neither'

            # Precompute per-cluster, then map to cells
            unique_clusters = cells_df['cluster'].unique()
            cluster_colors = {c: get_np_color_group(c)[0] for c in unique_clusters}
            cluster_groups = {c: get_np_color_group(c)[1] for c in unique_clusters}

            cells_df['_color'] = cells_df['cluster'].map(cluster_colors)
            cells_df['_legend_group'] = cells_df['cluster'].map(cluster_groups)

            # Only subsample background (Neither) cells, keep all foreground cells
            foreground_mask = cells_df['_legend_group'] != 'Neither'
            foreground_df = cells_df[foreground_mask]
            background_df = cells_df[~foreground_mask]

            if subsample_pct < 100:
                background_df = background_df.sample(frac=subsample_pct / 100, random_state=42)

            cells_df = pd.concat([foreground_df, background_df], ignore_index=True)

            # Apply diffusion range filter if enabled
            if 'enabled' in (diffusion_enabled or []):
                # Z-axis scaling factor for anisotropic distance
                # Gives Z an extra 0.15mm range compared to XY (for slice thickness)
                z_extra = 0.15  # mm
                z_scale = diffusion_range / (diffusion_range + z_extra)

                # Get ligand cell coordinates (include 'Both' as ligand sources)
                ligand_mask = cells_df['_legend_group'].isin(['Ligand', 'Both'])
                ligand_cells = cells_df[ligand_mask]

                if len(ligand_cells) > 0:
                    # Build scaled coordinates for KD-tree
                    ligand_coords = np.column_stack([
                        ligand_cells['x'].values,
                        ligand_cells['y'].values,
                        ligand_cells['z'].values * z_scale
                    ])
                    tree = cKDTree(ligand_coords)

                    # Query receptor cells
                    receptor_mask = cells_df['_legend_group'] == 'Receptor'

                    if receptor_mask.sum() > 0:
                        receptor_indices = cells_df.index[receptor_mask]
                        receptor_coords = np.column_stack([
                            cells_df.loc[receptor_indices, 'x'].values,
                            cells_df.loc[receptor_indices, 'y'].values,
                            cells_df.loc[receptor_indices, 'z'].values * z_scale
                        ])

                        # Find distance to nearest ligand source
                        distances, _ = tree.query(receptor_coords)

                        # Classify receptors as within or outside range
                        within_range = distances <= diffusion_range

                        # Update colors and legend groups
                        cells_df.loc[receptor_indices[within_range], '_color'] = '#4169E1'  # Blue - in range
                        cells_df.loc[receptor_indices[within_range], '_legend_group'] = 'Receptor (in range)'
                        cells_df.loc[receptor_indices[~within_range], '_color'] = '#40E0D0'  # Turquoise - out of range
                        cells_df.loc[receptor_indices[~within_range], '_legend_group'] = 'Receptor (out of range)'

        else:
            if subsample_pct < 100:
                cells_df = cells_df.sample(frac=subsample_pct / 100, random_state=42)
            cells_df['_color'] = '#D3D3D3'
            cells_df['_legend_group'] = 'NA'

    # Build all traces as a list first, then add in batch
    all_traces = []
    all_annotations = []

    for slice_idx, z_slice in enumerate(slices):
        row = slice_idx // n_cols + 1
        col = slice_idx % n_cols + 1
        # Get axis references for this subplot
        xaxis = f'x{slice_idx + 1}' if slice_idx > 0 else 'x'
        yaxis = f'y{slice_idx + 1}' if slice_idx > 0 else 'y'

        slice_df = cells_df[cells_df['z_slice'] == z_slice]

        # Add region boundaries - combine into single trace per slice
        if show_boundaries and z_slice in region_boundaries:
            all_x = []
            all_y = []
            for region, hull_points in region_boundaries[z_slice].items():
                x_coords = [p[0] for p in hull_points]
                y_coords = [p[1] for p in hull_points]
                all_x.extend(x_coords + [None])
                all_y.extend(y_coords + [None])

            if all_x:
                all_traces.append({
                    'type': 'scatter',
                    'x': all_x,
                    'y': all_y,
                    'mode': 'lines',
                    'fill': 'toself',
                    'fillcolor': 'rgba(240, 240, 240, 0.5)',
                    'line': {'color': '#CCCCCC', 'width': 1},
                    'hoverinfo': 'skip',
                    'showlegend': False,
                    'xaxis': xaxis,
                    'yaxis': yaxis,
                })

        # Add cells as single trace with per-point colors
        if len(slice_df) > 0:
            if mode == 'np':
                group_order_map = {'Neither': 0, 'Receptor (out of range)': 1, 'Receptor (in range)': 2,
                                   'Receptor': 3, 'Both': 4, 'Ligand': 5}
                slice_df = slice_df.copy()
                slice_df['_sort'] = slice_df['_legend_group'].map(group_order_map).fillna(0)
                slice_df = slice_df.sort_values('_sort')

            all_traces.append({
                'type': 'scattergl',
                'x': slice_df['x'].tolist(),
                'y': slice_df['y'].tolist(),
                'mode': 'markers',
                'marker': {
                    'size': point_size,
                    'color': slice_df['_color'].tolist(),
                    'opacity': 0.7,
                },
                'showlegend': False,
                'hovertemplate': '<b>%{customdata[0]}</b><br>Region: %{customdata[1]}<br><extra></extra>',
                'customdata': slice_df[['cluster', 'region_display']].values.tolist(),
                'xaxis': xaxis,
                'yaxis': yaxis,
            })

        # Collect annotations
        if show_labels and z_slice in region_centroids:
            for region, (cx, cy) in region_centroids[z_slice].items():
                all_annotations.append({
                    'x': cx,
                    'y': cy,
                    'xref': xaxis,
                    'yref': yaxis,
                    'text': region,
                    'showarrow': False,
                    'font': {'size': 8, 'color': '#333'},
                    'opacity': 0.7,
                })

    # Add all traces at once (much faster than individual add_trace calls)
    fig.add_traces(all_traces)

    # Update layout and add all annotations at once

    # Add Z-value labels to annotations
    for slice_idx, z_slice in enumerate(slices):
        xaxis = f'x{slice_idx + 1}' if slice_idx > 0 else 'x'
        yaxis = f'y{slice_idx + 1}' if slice_idx > 0 else 'y'
        all_annotations.append({
            'xref': f'{xaxis} domain',
            'yref': f'{yaxis} domain',
            'x': 0.98,
            'y': 0.02,
            'text': f'Z={z_slice}',
            'showarrow': False,
            'font': {'size': 9, 'color': '#666'},
            'xanchor': 'right',
            'yanchor': 'bottom',
        })

    fig.update_layout(
        height=2000,
        showlegend=False,
        margin=dict(l=10, r=10, t=10, b=10),
        paper_bgcolor='white',
        plot_bgcolor='white',
        annotations=all_annotations,
    )

    # Compute data range for default zoom (1.2x zoom = tighter range for x)
    x_min, x_max = cells_df['x'].min(), cells_df['x'].max()
    x_center = (x_min + x_max) / 2
    x_range = (x_max - x_min) / 1.2  # Zoom in by 1.2x

    x_range_final = [x_center - x_range / 2, x_center + x_range / 2]
    y_range_final = [5.5, 8.5]  # Fixed y range as requested

    # Update all axes - link them together with matches
    for i in range(1, n_rows * n_cols + 1):
        row_i = (i - 1) // n_cols + 1
        col_i = (i - 1) % n_cols + 1

        fig.update_xaxes(
            showticklabels=False,
            showgrid=False,
            zeroline=False,
            matches='x',  # Link all x-axes together
            range=x_range_final,
            row=row_i,
            col=col_i,
        )
        fig.update_yaxes(
            showticklabels=False,
            showgrid=False,
            zeroline=False,
            matches='y',  # Link all y-axes together
            scaleanchor='x',
            scaleratio=1,
            range=y_range_final,
            row=row_i,
            col=col_i,
        )

    return fig


def create_gene_bar(gene, expr, gene_info, is_ligand=True):
    """Create a horizontal bar for gene expression with hoverable tooltip."""
    # Get system info for this gene
    systems = gene_info.get(gene, [])
    if systems:
        # Build tooltip text
        system_names = list(set(s['system'] for s in systems))
        first_system = systems[0]
        tooltip_lines = [
            f"System: {', '.join(system_names)}",
            f"Ligands: {', '.join(first_system['ligand_genes'])}",
            f"Receptors: {', '.join(first_system['receptor_genes'])}",
        ]
        if first_system['functional_role']:
            tooltip_lines.append(f"Function: {first_system['functional_role']}")
        tooltip_text = '\n'.join(tooltip_lines)
    else:
        tooltip_text = f"{gene}: No system info available"

    # Normalize expression for bar width (assuming max ~10 log2CPM)
    bar_width_pct = min(100, expr * 10)
    color = '#667eea' if is_ligand else '#764ba2'  # Modern gradient colors

    return html.Div([
        # Gene name - fixed width
        html.Abbr(
            gene,
            title=tooltip_text,
            className="gene-tooltip",
            style={
                'fontSize': '0.8rem',
                'width': '70px',
                'flexShrink': '0',
                'textDecoration': 'none',
                'borderBottom': '1px dotted #667eea',
            }
        ),
        # Bar container - flex to fill space
        html.Div(
            style={
                'flex': '1',
                'height': '14px',
                'background': '#e2e8f0',
                'borderRadius': '4px',
                'marginLeft': '8px',
                'marginRight': '8px',
                'overflow': 'hidden',
            },
            children=[
                html.Div(
                    style={
                        'width': f'{bar_width_pct}%',
                        'height': '100%',
                        'background': f'linear-gradient(90deg, {color} 0%, {color}dd 100%)',
                        'borderRadius': '4px',
                    }
                )
            ]
        ),
        # Value - fixed width, right aligned
        html.Span(
            f'{expr:.1f}',
            style={
                'fontSize': '0.75rem',
                'width': '32px',
                'textAlign': 'right',
                'flexShrink': '0',
                'color': '#718096',
                'fontWeight': '500',
            }
        ),
    ], style={'marginBottom': '6px', 'display': 'flex', 'alignItems': 'center'})


def create_region_bar(region_name, count, total, index):
    """Create a horizontal bar for region distribution."""
    pct = (count / total) * 100
    # Gradient colors based on index
    colors = ['#667eea', '#764ba2', '#9f7aea', '#b794f4', '#d6bcfa']
    color = colors[index % len(colors)]

    return html.Div([
        html.Div([
            html.Span(
                region_name,
                style={
                    'fontSize': '0.8rem',
                    'minWidth': '50px',
                    'display': 'inline-block',
                    'fontWeight': '500',
                    'color': '#4a5568',
                }
            ),
            html.Div(
                style={
                    'display': 'inline-block',
                    'flex': '1',
                    'height': '18px',
                    'background': '#e2e8f0',
                    'borderRadius': '4px',
                    'marginLeft': '8px',
                    'marginRight': '8px',
                    'overflow': 'hidden',
                },
                children=[
                    html.Div(
                        style={
                            'width': f'{pct}%',
                            'height': '100%',
                            'background': f'linear-gradient(90deg, {color} 0%, {color}cc 100%)',
                            'borderRadius': '4px',
                        }
                    )
                ]
            ),
            html.Span(
                f'{count}',
                style={
                    'fontSize': '0.75rem',
                    'color': '#718096',
                    'minWidth': '30px',
                    'textAlign': 'right',
                    'fontWeight': '500',
                }
            ),
        ], style={'marginBottom': '6px', 'display': 'flex', 'alignItems': 'center'}),
    ])


def create_cell_details(click_data, cells_df, nt_mapping, cluster_expression, region_descriptions, gene_info):
    """Create cell details panel content."""
    if not click_data or 'points' not in click_data or not click_data['points']:
        return html.P("Click on a cell to see details", className="text-muted")

    point = click_data['points'][0]
    custom_data = point.get('customdata', [])

    if len(custom_data) < 2:
        return html.P("No data available for this cell", className="text-muted")

    cluster_name = custom_data[0]
    region = custom_data[1]

    # Find the cell in the dataframe
    cell_info = cells_df[cells_df['cluster'] == cluster_name].iloc[0] if len(cells_df[cells_df['cluster'] == cluster_name]) > 0 else None

    if cell_info is None:
        return html.P("Cell not found", className="text-muted")

    # Get region info
    region_info = region_descriptions.get(region, {})
    region_full_name = region_info.get('full_name', region)

    # Get NT type (direct lookup)
    nt_type = nt_mapping.get(cluster_name, 'Unknown')

    # Get top ligands and receptors for this cluster
    top_ligands = []
    top_receptors = []
    if cluster_name in cluster_expression:
        cluster_expr = cluster_expression[cluster_name]
        for gene, info in cluster_expr.items():
            if info['mean_expr'] >= 1:  # Expression threshold for display
                if info['is_ligand']:
                    top_ligands.append((gene, info['mean_expr']))
                if info['is_receptor']:
                    top_receptors.append((gene, info['mean_expr']))

    top_ligands = sorted(top_ligands, key=lambda x: -x[1])[:5]
    top_receptors = sorted(top_receptors, key=lambda x: -x[1])[:5]

    # Get region distribution for this cluster
    cluster_regions = cells_df[cells_df['cluster'] == cluster_name]['region'].value_counts().head(5)

    return html.Div([
        # Cluster info
        html.Div([
            html.H6("Cluster", className="text-muted mb-1"),
            html.P(cluster_name, className="fw-bold mb-2", style={'fontSize': '0.9rem'}),
        ]),

        # Hierarchy
        html.Div([
            html.H6("Cell Type Hierarchy", className="text-muted mb-1"),
            html.Ul([
                html.Li(f"Class: {cell_info.get('class', 'NA')}", style={'fontSize': '0.85rem'}),
                html.Li(f"Subclass: {cell_info.get('subclass', 'NA')}", style={'fontSize': '0.85rem'}),
                html.Li(f"Supertype: {cell_info.get('supertype', 'NA')}", style={'fontSize': '0.85rem'}),
            ], className="mb-2", style={'paddingLeft': '1.2rem'}),
        ]),

        html.Hr(),

        # NT system
        html.Div([
            html.H6("Neurotransmitter", className="text-muted mb-1"),
            html.P(nt_type, className="fw-bold mb-2"),
        ]),

        html.Hr(),

        # Top regions bar chart - separate bars for each region
        html.Div([
            html.H6("Top 5 Regions", className="text-muted mb-2", style={'fontSize': '0.85rem', 'fontWeight': '600'}),
            html.Div([
                create_region_bar(region, count, cluster_regions.sum(), i)
                for i, (region, count) in enumerate(cluster_regions.items())
            ]) if len(cluster_regions) > 0 else html.P("No region data", className="text-muted small"),
        ]),

        html.Hr(),

        # Neuropeptide profile - graphical bars with hoverable labels
        html.Div([
            html.H6("Top Ligands [log₂(CPM)]", className="text-muted mb-1"),
            html.Div([
                create_gene_bar(gene, expr, gene_info, is_ligand=True)
                for gene, expr in top_ligands
            ] if top_ligands else [html.P("None expressed", className="text-muted small")]),
        ], className="mb-3"),

        html.Div([
            html.H6("Top Receptors [log₂(CPM)]", className="text-muted mb-1"),
            html.Div([
                create_gene_bar(gene, expr, gene_info, is_ligand=False)
                for gene, expr in top_receptors
            ] if top_receptors else [html.P("None expressed", className="text-muted small")]),
        ]),
    ])


def register_callbacks(app, app_data):
    """Register all callbacks for the application."""

    cells_df = app_data['cells_df']
    slices = app_data['slices']
    nt_mapping = app_data['nt_mapping']
    nt_types = app_data['nt_types']
    np_systems = app_data['np_systems']
    gene_info = app_data['gene_info']
    cluster_expression = app_data['cluster_expression']
    cluster_np_expression = app_data['cluster_np_expression']
    region_boundaries = app_data['region_boundaries']
    region_centroids = app_data['region_centroids']
    region_colors = app_data['region_colors']
    region_descriptions = app_data['region_descriptions']

    @app.callback(
        [
            Output('cluster-controls', 'style'),
            Output('nt-controls', 'style'),
            Output('np-controls', 'style'),
        ],
        Input('viz-mode', 'value'),
    )
    def toggle_controls(mode):
        """Show/hide mode-specific controls."""
        cluster_style = {'display': 'block'} if mode == 'cluster' else {'display': 'none'}
        nt_style = {'display': 'block'} if mode == 'nt' else {'display': 'none'}
        np_style = {'display': 'block'} if mode == 'np' else {'display': 'none'}
        return cluster_style, nt_style, np_style

    @app.callback(
        Output('diffusion-range-container', 'style'),
        Input('diffusion-filter-enabled', 'value'),
    )
    def toggle_diffusion_slider(enabled):
        """Show/hide diffusion range slider based on checkbox."""
        if 'enabled' in (enabled or []):
            return {'display': 'block'}
        return {'display': 'none'}

    @app.callback(
        Output('slice-grid', 'figure'),
        [
            Input('viz-mode', 'value'),
            Input('cluster-level', 'value'),
            Input('nt-system', 'value'),
            Input('np-system', 'value'),
            Input('expression-threshold', 'value'),
            Input('display-options', 'value'),
            Input('point-size', 'value'),
            Input('subsample-pct', 'value'),
            Input('diffusion-filter-enabled', 'value'),
            Input('diffusion-range', 'value'),
        ],
    )
    def update_slices(mode, cluster_level, nt_system, np_system, threshold, display_options, point_size, subsample_pct, diffusion_enabled, diffusion_range):
        """Update the slice grid visualization."""
        return create_slice_figure(
            cells_df=cells_df,
            slices=slices,
            mode=mode,
            cluster_level=cluster_level,
            nt_system=nt_system,
            np_system=np_system,
            threshold=threshold,
            display_options=display_options or [],
            point_size=point_size,
            subsample_pct=subsample_pct or 30,
            diffusion_enabled=diffusion_enabled,
            diffusion_range=diffusion_range or 0.5,
            nt_mapping=nt_mapping,
            np_systems=np_systems,
            cluster_expression=cluster_expression,
            cluster_np_expression=cluster_np_expression,
            region_boundaries=region_boundaries,
            region_centroids=region_centroids,
            region_colors=region_colors,
        )

    @app.callback(
        Output('cell-details', 'children'),
        Input('slice-grid', 'clickData'),
    )
    def update_details(click_data):
        """Update cell details on click."""
        return create_cell_details(
            click_data,
            cells_df,
            nt_mapping,
            cluster_expression,
            region_descriptions,
            gene_info,
        )
