"""Dash callbacks for HypoMap Coronal Atlas Viewer."""

from dash import Input, Output, State, ctx, callback, clientside_callback, ClientsideFunction, ALL, MATCH
from dash import html
from dash.exceptions import PreventUpdate
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


def _region_full_name(region_display, region_descriptions):
    """Get full name for a region acronym, stripping -L/-R suffix for lookup."""
    if not region_descriptions:
        return region_display
    # Strip lateral suffix for lookup
    base = region_display.rsplit('-', 1)[0] if region_display.endswith(('-L', '-R')) else region_display
    info = region_descriptions.get(base, {})
    return info.get('full_name', '') or region_display


def create_slice_figure(
    cells_df,
    slices,
    mode,
    cluster_level,
    nt_system,
    np_system,
    hormone_system,
    threshold,
    display_options,
    point_size,
    subsample_pct,
    diffusion_enabled,
    diffusion_range,
    nt_mapping,
    np_systems,
    hormone_systems,
    cluster_expression,
    cluster_np_expression,
    region_boundaries,
    region_centroids,
    region_colors,
    highlighted_region=None,
    rainbow_mode=False,
    n_cols=2,
    region_descriptions=None,
    fig_height=800,
):
    """Create the grid of coronal slices.

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
    import math
    n_rows = math.ceil(n_slices / n_cols)

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

    # Extract raw arrays from cells_df (zero-copy views — no df.copy() needed).
    # All mode-specific logic operates on numpy arrays instead of adding columns to
    # a copied DataFrame, eliminating the ~0.5s pandas copy + consolidation overhead.
    n_cells = len(cells_df)
    raw_x = cells_df['x'].values
    raw_y = cells_df['y'].values
    raw_z_slice = cells_df['z_slice'].values
    raw_cluster = cells_df['cluster'].values
    raw_region_display = cells_df['region_display'].values

    # Each mode computes:
    #   color_arr  — per-cell color (full size, before subsampling)
    #   keep_idx   — indices to retain after fg/bg subsampling
    #   stroke_arr — per-cell stroke color (rainbow NP mode only)
    #   legend_arr — per-cell legend group (NP mode only, for sort + diffusion)
    stroke_arr = None
    legend_arr = None

    if mode == 'cluster':
        # Subsample all cells uniformly
        if subsample_pct < 100:
            n_keep = max(1, int(n_cells * subsample_pct / 100))
            keep_idx = np.random.RandomState(42).choice(n_cells, size=n_keep, replace=False)
        else:
            keep_idx = np.arange(n_cells)

        raw_ct = cells_df[cluster_level].values
        color_map = {cat: get_color_for_category(unique_categories, cat)
                     for cat in unique_categories}
        color_arr = np.array([color_map.get(c, '#808080') for c in raw_ct])

    elif mode == 'nt':
        # Map cluster to NT type via dict lookup on numpy array
        nt_type_arr = np.array([nt_mapping.get(c, 'NA') for c in raw_cluster])
        fg_mask = nt_type_arr == nt_system
        fg_idx = np.where(fg_mask)[0]
        bg_idx = np.where(~fg_mask)[0]

        if subsample_pct < 100 and len(bg_idx) > 0:
            n_bg = max(1, int(len(bg_idx) * subsample_pct / 100))
            bg_idx = np.random.RandomState(42).choice(bg_idx, size=n_bg, replace=False)

        # Background first, foreground last (fg renders on top)
        keep_idx = np.concatenate([bg_idx, fg_idx])
        color_arr = np.where(fg_mask, '#000000', '#D3D3D3')

    elif mode == 'np':
        if np_system and np_system in cluster_np_expression:
            system_lookup = cluster_np_expression[np_system]

            # Precompute per-cluster color and legend group
            unique_clusters = np.unique(raw_cluster)
            cluster_np_colors = {}
            cluster_groups = {}
            for c in unique_clusters:
                if c not in system_lookup:
                    cluster_np_colors[c] = '#D3D3D3'
                    cluster_groups[c] = 'Neither'
                    continue
                ligand_expr, receptor_expr = system_lookup[c]
                ligand_above = ligand_expr >= threshold
                receptor_above = receptor_expr >= threshold
                if ligand_above and receptor_above:
                    intensity = min(1.0, (ligand_expr + receptor_expr) / 10)
                    r = int(128 + 127 * intensity)
                    b = int(128 + 127 * intensity)
                    cluster_np_colors[c] = f'rgb({r}, 0, {b})'
                    cluster_groups[c] = 'Both'
                elif ligand_above:
                    intensity = min(1.0, ligand_expr / 5)
                    r = int(100 + 155 * intensity)
                    cluster_np_colors[c] = f'rgb({r}, 50, 50)'
                    cluster_groups[c] = 'Ligand'
                elif receptor_above:
                    intensity = min(1.0, receptor_expr / 5)
                    b = int(100 + 155 * intensity)
                    cluster_np_colors[c] = f'rgb(50, 50, {b})'
                    cluster_groups[c] = 'Receptor'
                else:
                    cluster_np_colors[c] = '#D3D3D3'
                    cluster_groups[c] = 'Neither'

            legend_arr = np.array([cluster_groups.get(c, 'Neither') for c in raw_cluster])

            if rainbow_mode:
                unique_clusters_all = sorted(unique_clusters)
                cluster_fill_colors = {cat: get_color_for_category(unique_clusters_all, cat)
                                       for cat in unique_clusters_all}
                stroke_arr = np.array([cluster_np_colors.get(c, '#D3D3D3') for c in raw_cluster])
                color_arr = np.array([
                    cluster_fill_colors.get(c, '#D3D3D3')
                    if cluster_groups.get(c, 'Neither') != 'Neither'
                    else '#D3D3D3'
                    for c in raw_cluster
                ])
            else:
                color_arr = np.array([cluster_np_colors.get(c, '#D3D3D3') for c in raw_cluster])

            # fg/bg split — keep all foreground, subsample background
            fg_mask = legend_arr != 'Neither'
            fg_idx = np.where(fg_mask)[0]
            bg_idx = np.where(~fg_mask)[0]

            if subsample_pct < 100 and len(bg_idx) > 0:
                n_bg = max(1, int(len(bg_idx) * subsample_pct / 100))
                bg_idx = np.random.RandomState(42).choice(bg_idx, size=n_bg, replace=False)

            keep_idx = np.concatenate([fg_idx, bg_idx])
        else:
            if subsample_pct < 100:
                n_keep = max(1, int(n_cells * subsample_pct / 100))
                keep_idx = np.random.RandomState(42).choice(n_cells, size=n_keep, replace=False)
            else:
                keep_idx = np.arange(n_cells)
            color_arr = np.full(n_cells, '#D3D3D3')
            legend_arr = np.full(n_cells, 'NA')

    elif mode == 'hormone':
        if hormone_system and hormone_system in hormone_systems:
            receptor_genes = hormone_systems[hormone_system]['receptors']

            unique_clusters = np.unique(raw_cluster)
            cluster_has_receptor = {}
            for c in unique_clusters:
                if c not in cluster_expression:
                    cluster_has_receptor[c] = False
                    continue
                cluster_genes = cluster_expression[c]
                cluster_has_receptor[c] = any(
                    gene in cluster_genes and cluster_genes[gene]['mean_expr'] >= threshold
                    for gene in receptor_genes
                )

            receptor_mask = np.array([cluster_has_receptor.get(c, False) for c in raw_cluster])
            fg_idx = np.where(receptor_mask)[0]
            bg_idx = np.where(~receptor_mask)[0]

            if subsample_pct < 100 and len(bg_idx) > 0:
                n_bg = max(1, int(len(bg_idx) * subsample_pct / 100))
                bg_idx = np.random.RandomState(42).choice(bg_idx, size=n_bg, replace=False)

            # Background first, foreground last
            keep_idx = np.concatenate([bg_idx, fg_idx])
            color_arr = np.where(receptor_mask, '#000000', '#D3D3D3')
        else:
            if subsample_pct < 100:
                n_keep = max(1, int(n_cells * subsample_pct / 100))
                keep_idx = np.random.RandomState(42).choice(n_cells, size=n_keep, replace=False)
            else:
                keep_idx = np.arange(n_cells)
            color_arr = np.full(n_cells, '#D3D3D3')

    # Apply subsampling to all arrays at once (fancy indexing returns copies)
    arr_x = np.round(raw_x[keep_idx], 3)
    arr_y = np.round(raw_y[keep_idx], 3)
    arr_z = raw_z_slice[keep_idx]
    arr_color = color_arr[keep_idx]
    arr_cluster = raw_cluster[keep_idx]
    arr_region = raw_region_display[keep_idx]
    has_stroke = stroke_arr is not None
    if has_stroke:
        arr_stroke = stroke_arr[keep_idx]
    else:
        arr_stroke = None

    # Diffusion range filter (NP mode only — modifies arr_color/arr_stroke in place)
    if mode == 'np' and legend_arr is not None and 'enabled' in (diffusion_enabled or []):
        arr_legend = legend_arr[keep_idx]
        raw_z = cells_df['z'].values
        arr_z_coord = raw_z[keep_idx]

        z_extra = 0.15  # mm
        z_scale = diffusion_range / (diffusion_range + z_extra)

        ligand_mask = np.isin(arr_legend, ['Ligand', 'Both'])
        if ligand_mask.any():
            ligand_coords = np.column_stack([
                arr_x[ligand_mask], arr_y[ligand_mask],
                arr_z_coord[ligand_mask] * z_scale
            ])
            tree = cKDTree(ligand_coords)

            receptor_mask = arr_legend == 'Receptor'
            if receptor_mask.any():
                receptor_coords = np.column_stack([
                    arr_x[receptor_mask], arr_y[receptor_mask],
                    arr_z_coord[receptor_mask] * z_scale
                ])
                distances, _ = tree.query(receptor_coords)
                within_range = distances <= diffusion_range

                receptor_indices = np.where(receptor_mask)[0]
                color_target = arr_stroke if has_stroke else arr_color
                color_target[receptor_indices[within_range]] = '#4169E1'
                color_target[receptor_indices[~within_range]] = '#40E0D0'
                arr_legend[receptor_indices[within_range]] = 'Receptor (in range)'
                arr_legend[receptor_indices[~within_range]] = 'Receptor (out of range)'

    # Sort for NP mode (render foreground on top)
    if mode == 'np' and legend_arr is not None:
        if 'enabled' not in (diffusion_enabled or []):
            arr_legend = legend_arr[keep_idx]
        group_order_map = {'Neither': 0, 'Receptor (out of range)': 1, 'Receptor (in range)': 2,
                           'Receptor': 3, 'Both': 4, 'Ligand': 5}
        sort_keys = np.array([group_order_map.get(g, 0) for g in arr_legend])
        sort_idx = np.argsort(sort_keys, kind='stable')
        arr_x = arr_x[sort_idx]
        arr_y = arr_y[sort_idx]
        arr_z = arr_z[sort_idx]
        arr_color = arr_color[sort_idx]
        arr_cluster = arr_cluster[sort_idx]
        arr_region = arr_region[sort_idx]
        if has_stroke:
            arr_stroke = arr_stroke[sort_idx]

    # Precompute region full names and divisions
    unique_regions = np.unique(arr_region)
    region_full_name_map = {r: _region_full_name(r, region_descriptions) for r in unique_regions}
    arr_region_full = np.array([region_full_name_map.get(r, r) for r in arr_region])

    def _region_division(region_display):
        if not region_descriptions:
            return ''
        base = region_display.rsplit('-', 1)[0] if region_display.endswith(('-L', '-R')) else region_display
        return region_descriptions.get(base, {}).get('division', '')

    region_division_map = {r: _region_division(r) for r in unique_regions}
    arr_region_div = np.array([region_division_map.get(r, '') for r in arr_region])

    # Build all traces, shapes, and annotations
    all_traces = []
    all_annotations = []
    all_shapes = []

    for slice_idx, z_slice in enumerate(slices):
        row = slice_idx // n_cols + 1
        col = slice_idx % n_cols + 1
        # Get axis references for this subplot
        xaxis = f'x{slice_idx + 1}' if slice_idx > 0 else 'x'
        yaxis = f'y{slice_idx + 1}' if slice_idx > 0 else 'y'

        # Add region boundaries as layout shapes (avoids scatter/scattergl misalignment)
        if show_boundaries and z_slice in region_boundaries:
            for region, hull_points in region_boundaries[z_slice].items():
                is_highlighted = (highlighted_region and
                                  (region == highlighted_region or
                                   region.startswith(f'{highlighted_region}-')))

                # Build SVG path: M x0,y0 L x1,y1 ... Z
                path_parts = [f'M {hull_points[0][0]},{hull_points[0][1]}']
                for px, py in hull_points[1:]:
                    path_parts.append(f'L {px},{py}')
                path_parts.append('Z')
                path_str = ' '.join(path_parts)

                if is_highlighted:
                    all_shapes.append({
                        'type': 'path', 'path': path_str,
                        'fillcolor': 'rgba(34, 197, 94, 0.3)',
                        'line': {'color': '#22c55e', 'width': 2},
                        'layer': 'below', 'xref': xaxis, 'yref': yaxis,
                    })
                else:
                    base_region = region.rsplit('-', 1)[0] if region.endswith(('-L', '-R')) else region
                    hex_color = region_colors.get(base_region, '#F0F0F0')
                    r = int(hex_color[1:3], 16)
                    g = int(hex_color[3:5], 16)
                    b = int(hex_color[5:7], 16)
                    all_shapes.append({
                        'type': 'path', 'path': path_str,
                        'fillcolor': f'rgba({r}, {g}, {b}, 0.3)',
                        'line': {'color': hex_color, 'width': 1},
                        'layer': 'below', 'xref': xaxis, 'yref': yaxis,
                    })

        # Add cells — use numpy boolean mask on pre-extracted arrays (no pandas)
        mask = arr_z == z_slice
        n_pts = mask.sum()
        if n_pts > 0:
            s_x = arr_x[mask]
            s_y = arr_y[mask]
            s_color = arr_color[mask]

            marker_dict = {
                'size': point_size,
                'color': s_color.tolist(),
                'opacity': 0.7,
            }

            if has_stroke:
                marker_dict['line'] = {
                    'color': arr_stroke[mask].tolist(),
                    'width': 2,
                }

            all_traces.append({
                'type': 'scattergl',
                'x': s_x.tolist(),
                'y': s_y.tolist(),
                'mode': 'markers',
                'marker': marker_dict,
                'showlegend': False,
                'hovertemplate': '<b>%{customdata[0]}</b><br>Region: %{customdata[1]} (%{customdata[2]})<br>Division: %{customdata[3]}<extra></extra>',
                'customdata': np.column_stack([arr_cluster[mask], arr_region[mask], arr_region_full[mask], arr_region_div[mask]]).tolist(),
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

    # Build figure as plain dict for speed (bypasses Plotly's per-trace validation
    # which adds ~400ms overhead via add_traces). Dash accepts dicts for dcc.Graph.
    # We use make_subplots only to generate the subplot layout structure.
    fig_dict = fig.to_dict()

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

    # Equalize x/y ranges for 1:1 aspect ratio. Combined with square subplot
    # pixel dimensions (fig_height computed in callers), this enforces 1:1 scaling
    # without scaleanchor (which conflicts with matches + scattergl).
    x_min, x_max = float(arr_x.min()), float(arr_x.max())
    y_min, y_max = float(arr_y.min()), float(arr_y.max())
    x_span = x_max - x_min
    y_span = y_max - y_min
    max_span = max(x_span, y_span) * 1.05  # 5% padding
    x_center = (x_min + x_max) / 2
    y_center = (y_min + y_max) / 2
    half = max_span / 2
    x_range_final = [x_center - half, x_center + half]
    y_range_final = [y_center + half, y_center - half]  # inverted y
    axis_common_x = dict(showticklabels=False, showgrid=False, zeroline=False, matches='x', range=x_range_final)
    axis_common_y = dict(showticklabels=False, showgrid=False, zeroline=False, matches='y', range=y_range_final)

    layout = fig_dict['layout']
    for i in range(1, n_rows * n_cols + 1):
        x_key = f'xaxis{i}' if i > 1 else 'xaxis'
        y_key = f'yaxis{i}' if i > 1 else 'yaxis'
        layout.setdefault(x_key, {}).update(axis_common_x)
        layout.setdefault(y_key, {}).update(axis_common_y)

    layout.update(
        height=fig_height,
        showlegend=False,
        margin=dict(l=10, r=10, t=10, b=10),
        paper_bgcolor='white',
        plot_bgcolor='white',
        annotations=all_annotations,
        shapes=all_shapes,
    )

    fig_dict['data'] = all_traces
    return fig_dict


def create_gene_bar(gene, expr, gene_info, is_ligand=True, use_quantile=False, quantile_value=None):
    """Create a horizontal bar for gene expression with hoverable tooltip.

    Args:
        gene: Gene name
        expr: Expression value (log2CPM)
        gene_info: Dict with gene system info for tooltips
        is_ligand: True for ligand, False for receptor
        use_quantile: If True, display quantile instead of log2CPM
        quantile_value: The quantile percentile (0-100) if use_quantile is True
    """
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

    # Determine display value and bar width
    if use_quantile and quantile_value is not None:
        display_value = f'{quantile_value:.0f}'
        bar_width_pct = quantile_value  # 0-100 directly maps to %
    else:
        display_value = f'{expr:.1f}'
        bar_width_pct = min(100, expr * 10)  # Assuming max ~10 log2CPM

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
            display_value,
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


def create_cell_details(click_data, cells_df, nt_mapping, cluster_expression, region_descriptions, gene_info, use_quantile=False, gene_quantiles=None):
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
    region_division = region_info.get('division', '')

    # Get NT type (direct lookup)
    nt_type = nt_mapping.get(cluster_name, 'Unknown')

    # Get top ligands and receptors for this cluster
    top_ligands = []
    top_receptors = []
    if cluster_name in cluster_expression:
        cluster_expr = cluster_expression[cluster_name]
        for gene, info in cluster_expr.items():
            if info['mean_expr'] >= 1:  # Expression threshold for display
                # Get quantile value if available
                quantile = None
                if gene_quantiles and gene in gene_quantiles:
                    quantile = gene_quantiles[gene].get(cluster_name)

                if info['is_ligand']:
                    top_ligands.append((gene, info['mean_expr'], quantile))
                if info['is_receptor']:
                    top_receptors.append((gene, info['mean_expr'], quantile))

    # Sort by quantile if in quantile mode, otherwise by expression
    if use_quantile:
        top_ligands = sorted(top_ligands, key=lambda x: -(x[2] or 0))[:5]
        top_receptors = sorted(top_receptors, key=lambda x: -(x[2] or 0))[:5]
    else:
        top_ligands = sorted(top_ligands, key=lambda x: -x[1])[:5]
        top_receptors = sorted(top_receptors, key=lambda x: -x[1])[:5]

    # Get region distribution for this cluster
    cluster_regions = cells_df[cells_df['cluster'] == cluster_name]['region'].value_counts().head(5)

    # Header label based on mode
    metric_label = "Quantile" if use_quantile else "log₂(CPM)"

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

        # Region info
        html.Div([
            html.H6("Region", className="text-muted mb-1"),
            html.P([
                html.Span(region, className="fw-bold"),
                html.Span(f" — {region_full_name}", style={'fontSize': '0.85rem'}) if region_full_name and region_full_name != region else None,
            ], className="mb-1", style={'fontSize': '0.9rem'}),
            html.P(f"Division: {region_division}", style={'fontSize': '0.8rem', 'color': '#666'}) if region_division else None,
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
            html.H6(f"Top Ligands [{metric_label}]", className="text-muted mb-1"),
            html.Div([
                create_gene_bar(gene, expr, gene_info, is_ligand=True, use_quantile=use_quantile, quantile_value=quantile)
                for gene, expr, quantile in top_ligands
            ] if top_ligands else [html.P("None expressed", className="text-muted small")]),
        ], className="mb-3"),

        html.Div([
            html.H6(f"Top Receptors [{metric_label}]", className="text-muted mb-1"),
            html.Div([
                create_gene_bar(gene, expr, gene_info, is_ligand=False, use_quantile=use_quantile, quantile_value=quantile)
                for gene, expr, quantile in top_receptors
            ] if top_receptors else [html.P("None expressed", className="text-muted small")]),
        ]),
    ])


def register_callbacks(app, app_data, enable_region_highlight=False, enable_quantile_toggle=False):
    """Register all callbacks for the application."""

    datasets = app_data['datasets']
    default_dataset = app_data['default_dataset']
    nt_mapping = app_data['nt_mapping']
    nt_types = app_data['nt_types']
    np_systems = app_data['np_systems']
    hormone_systems = app_data['hormone_systems']
    gene_info = app_data['gene_info']
    region_descriptions = app_data['region_descriptions']

    def get_dataset(dataset_name):
        """Get dataset bundle by name, falling back to default."""
        if dataset_name and dataset_name in datasets:
            return datasets[dataset_name]
        return datasets[default_dataset]

    @app.callback(
        [
            Output('cluster-controls', 'style'),
            Output('nt-controls', 'style'),
            Output('np-controls', 'style'),
            Output('hormone-controls', 'style'),
            Output('expression-threshold-container', 'style'),
        ],
        Input('viz-mode', 'value'),
    )
    def toggle_controls(mode):
        """Show/hide mode-specific controls."""
        cluster_style = {'display': 'block'} if mode == 'cluster' else {'display': 'none'}
        nt_style = {'display': 'block'} if mode == 'nt' else {'display': 'none'}
        np_style = {'display': 'block'} if mode == 'np' else {'display': 'none'}
        hormone_style = {'display': 'block'} if mode == 'hormone' else {'display': 'none'}
        threshold_style = {'display': 'block'} if mode in ['np', 'hormone'] else {'display': 'none'}
        return cluster_style, nt_style, np_style, hormone_style, threshold_style

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
        Output('np-system-info', 'children'),
        Input('np-system', 'value'),
    )
    def update_np_system_info(system):
        """Display info about selected neuropeptide system."""
        if not system or system not in np_systems:
            return None

        system_data = np_systems[system]
        ligands = sorted(system_data['ligands'])
        # Flatten receptor tuples and deduplicate
        receptors = sorted(set(g for tup in system_data['receptors'] for g in tup))

        # Get functional role from gene_info (look up first ligand)
        func_role = ''
        if ligands and ligands[0] in gene_info:
            func_role = gene_info[ligands[0]][0].get('functional_role', '')

        return [
            html.Div([html.Strong("Ligands: "), ', '.join(ligands)]) if ligands else None,
            html.Div([html.Strong("Receptors: "), ', '.join(receptors)]) if receptors else None,
            html.Div([html.Strong("Function: "), func_role], style={'fontStyle': 'italic'}) if func_role else None,
        ]

    # Region highlight: Click to select region
    if enable_region_highlight:
        @app.callback(
            Output('selected-region', 'data'),
            Input({'type': 'region-btn', 'index': ALL}, 'n_clicks'),
            State('selected-region', 'data'),
            prevent_initial_call=True,
        )
        def select_region(n_clicks_list, current_selection):
            """Handle region button clicks - toggle selection."""
            if not ctx.triggered_id or all(n == 0 for n in n_clicks_list if n):
                raise PreventUpdate

            clicked_region = ctx.triggered_id['index']

            # Toggle: if already selected, deselect; otherwise select
            if current_selection == clicked_region:
                return None
            return clicked_region

        @app.callback(
            Output({'type': 'region-btn', 'index': ALL}, 'style'),
            Input('selected-region', 'data'),
            State({'type': 'region-btn', 'index': ALL}, 'id'),
        )
        def update_button_styles(selected_region, button_ids):
            """Update button styles to show selected state."""
            styles = []
            for btn_id in button_ids:
                region = btn_id['index']
                if region == selected_region:
                    # Selected style - green
                    styles.append({
                        'display': 'inline-block',
                        'padding': '2px 6px',
                        'margin': '2px',
                        'fontSize': '0.75rem',
                        'borderRadius': '4px',
                        'backgroundColor': '#22c55e',
                        'border': '1px solid #16a34a',
                        'cursor': 'pointer',
                        'color': 'white',
                        'fontWeight': '600',
                    })
                else:
                    # Default style
                    styles.append({
                        'display': 'inline-block',
                        'padding': '2px 6px',
                        'margin': '2px',
                        'fontSize': '0.75rem',
                        'borderRadius': '4px',
                        'backgroundColor': '#f0f0f0',
                        'border': '1px solid transparent',
                        'cursor': 'pointer',
                        'color': '#333',
                    })
            return styles

    # Update z-range slider when dataset changes
    @app.callback(
        [
            Output('z-range', 'min'),
            Output('z-range', 'max'),
            Output('z-range', 'marks'),
            Output('z-range', 'value'),
        ],
        Input('dataset-radio', 'value'),
    )
    def update_z_range_slider(dataset_name):
        """Update z-range slider to match the selected dataset's slices."""
        ds = get_dataset(dataset_name)
        slices = ds['slices']
        z_min = slices[0]
        z_max = slices[-1]
        marks = {}
        for i, z in enumerate(slices):
            if i == 0 or i == len(slices) - 1 or i % 4 == 0:
                marks[z] = f'{z:.1f}'
            else:
                marks[z] = ''
        return z_min, z_max, marks, [z_min, z_max]

    # Update subsample slider default when dataset changes
    @app.callback(
        Output('subsample-pct', 'value'),
        Input('dataset-radio', 'value'),
    )
    def update_subsample_default(dataset_name):
        """Scale subsample % to target ~60k points regardless of dataset size."""
        ds = get_dataset(dataset_name)
        n_cells = len(ds['cells_df'])
        default = max(5, min(30, int(60000 / n_cells * 100)))
        # Round to nearest step of 5
        return round(default / 5) * 5

    # Build inputs list for slice-grid callback
    slice_inputs = [
        Input('viz-mode', 'value'),
        Input('cluster-level', 'value'),
        Input('nt-system', 'value'),
        Input('np-system', 'value'),
        Input('hormone-system', 'value'),
        Input('expression-threshold', 'value'),
        Input('display-options', 'value'),
        Input('point-size', 'value'),
        Input('subsample-pct', 'value'),
        Input('diffusion-filter-enabled', 'value'),
        Input('diffusion-range', 'value'),
        Input('rainbow-mode', 'value'),
        Input('dataset-radio', 'value'),
        Input('grid-columns', 'value'),
        Input('z-range', 'value'),
    ]

    if enable_region_highlight:
        slice_inputs.append(Input('selected-region', 'data'))

        @app.callback(
            [Output('slice-grid', 'figure'),
             Output('slice-grid', 'style')],
            slice_inputs,
        )
        def update_slices_with_highlight(mode, cluster_level, nt_system, np_system, hormone_system, threshold, display_options, point_size, subsample_pct, diffusion_enabled, diffusion_range, rainbow_enabled, dataset_name, grid_columns, z_range, selected_region):
            """Update the slice grid visualization with region highlighting."""
            ds = get_dataset(dataset_name)
            n_cols = grid_columns or 2
            # Filter slices by z-range
            z_min, z_max = (z_range or [ds['slices'][0], ds['slices'][-1]])
            filtered_slices = [s for s in ds['slices'] if z_min <= s <= z_max]
            if not filtered_slices:
                filtered_slices = ds['slices']
            import math
            n_rows = math.ceil(len(filtered_slices) / n_cols)
            # Square subplots for 1:1 aspect ratio (ranges are equalized in create_slice_figure)
            cells_df = ds['cells_df']
            subplot_size = 900 / n_cols
            fig_height = min(16000, max(800, int(n_rows * subplot_size)))
            # Default subsample scales with dataset size (target ~60k points)
            default_subsample = max(5, min(30, int(60000 / len(cells_df) * 100)))
            fig = create_slice_figure(
                cells_df=cells_df,
                slices=filtered_slices,
                mode=mode,
                cluster_level=cluster_level,
                nt_system=nt_system,
                np_system=np_system,
                hormone_system=hormone_system,
                threshold=threshold,
                display_options=display_options or [],
                point_size=point_size,
                subsample_pct=subsample_pct or default_subsample,
                diffusion_enabled=diffusion_enabled,
                diffusion_range=diffusion_range or 0.5,
                nt_mapping=nt_mapping,
                np_systems=np_systems,
                hormone_systems=hormone_systems,
                cluster_expression=ds['cluster_expression'],
                cluster_np_expression=ds['cluster_np_expression'],
                region_boundaries=ds['region_boundaries'],
                region_centroids=ds['region_centroids'],
                region_colors=ds['region_colors'],
                highlighted_region=selected_region,
                rainbow_mode='enabled' in (rainbow_enabled or []),
                n_cols=n_cols,
                region_descriptions=region_descriptions,
                fig_height=fig_height,
            )
            return fig, {'height': f'{fig_height}px'}
    else:
        @app.callback(
            [Output('slice-grid', 'figure'),
             Output('slice-grid', 'style')],
            slice_inputs,
        )
        def update_slices(mode, cluster_level, nt_system, np_system, hormone_system, threshold, display_options, point_size, subsample_pct, diffusion_enabled, diffusion_range, rainbow_enabled, dataset_name, grid_columns, z_range):
            """Update the slice grid visualization."""
            ds = get_dataset(dataset_name)
            n_cols = grid_columns or 2
            # Filter slices by z-range
            z_min, z_max = (z_range or [ds['slices'][0], ds['slices'][-1]])
            filtered_slices = [s for s in ds['slices'] if z_min <= s <= z_max]
            if not filtered_slices:
                filtered_slices = ds['slices']
            import math
            n_rows = math.ceil(len(filtered_slices) / n_cols)
            # Square subplots for 1:1 aspect ratio (ranges are equalized in create_slice_figure)
            cells_df = ds['cells_df']
            subplot_size = 900 / n_cols
            fig_height = min(16000, max(800, int(n_rows * subplot_size)))
            # Default subsample scales with dataset size (target ~60k points)
            default_subsample = max(5, min(30, int(60000 / len(cells_df) * 100)))
            fig = create_slice_figure(
                cells_df=cells_df,
                slices=filtered_slices,
                mode=mode,
                cluster_level=cluster_level,
                nt_system=nt_system,
                np_system=np_system,
                hormone_system=hormone_system,
                threshold=threshold,
                display_options=display_options or [],
                point_size=point_size,
                subsample_pct=subsample_pct or default_subsample,
                diffusion_enabled=diffusion_enabled,
                diffusion_range=diffusion_range or 0.5,
                nt_mapping=nt_mapping,
                np_systems=np_systems,
                hormone_systems=hormone_systems,
                cluster_expression=ds['cluster_expression'],
                cluster_np_expression=ds['cluster_np_expression'],
                region_boundaries=ds['region_boundaries'],
                region_centroids=ds['region_centroids'],
                region_colors=ds['region_colors'],
                highlighted_region=None,
                rainbow_mode='enabled' in (rainbow_enabled or []),
                n_cols=n_cols,
                region_descriptions=region_descriptions,
                fig_height=fig_height,
            )
            return fig, {'height': f'{fig_height}px'}

    # Quantile toggle: Expression mode toggle callbacks
    if enable_quantile_toggle:
        @app.callback(
            Output('expr-mode-label', 'children'),
            Input('expr-mode-toggle', 'value'),
        )
        def update_expr_mode_label(use_quantile):
            """Update expression mode label based on toggle."""
            return "Quantile" if use_quantile else "log₂(CPM)"

        @app.callback(
            Output('cell-details', 'children'),
            [
                Input('slice-grid', 'clickData'),
                Input('expr-mode-toggle', 'value'),
                Input('dataset-radio', 'value'),
            ],
        )
        def update_details_experimental(click_data, use_quantile, dataset_name):
            """Update cell details on click (with quantile support)."""
            ds = get_dataset(dataset_name)
            return create_cell_details(
                click_data,
                ds['cells_df'],
                nt_mapping,
                ds['cluster_expression'],
                region_descriptions,
                gene_info,
                use_quantile=use_quantile or False,
                gene_quantiles=ds['gene_quantiles'],
            )
    else:
        @app.callback(
            Output('cell-details', 'children'),
            [
                Input('slice-grid', 'clickData'),
                Input('dataset-radio', 'value'),
            ],
        )
        def update_details(click_data, dataset_name):
            """Update cell details on click."""
            ds = get_dataset(dataset_name)
            return create_cell_details(
                click_data,
                ds['cells_df'],
                nt_mapping,
                ds['cluster_expression'],
                region_descriptions,
                gene_info,
                use_quantile=False,
                gene_quantiles={},
            )
