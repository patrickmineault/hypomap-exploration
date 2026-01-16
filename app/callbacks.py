"""Dash callbacks for HypoMap 3D Viewer interactivity."""

from dash import Input, Output, State, ctx, ALL, MATCH
from dash import html
import dash_bootstrap_components as dbc
import plotly.graph_objects as go
import plotly.express as px
import pandas as pd
import numpy as np
import re
import uuid

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.atlas.region_mapping import (
    HYPOMAP_TO_CCF,
)


def parse_cell_type_with_genes(cell_type_str, gene_descriptions):
    """Parse a cell type string and identify gene names.

    Args:
        cell_type_str: String like "101 ZI Pax6 Gaba" or "C465-344: Notch2.Gjb6.Etnppl.Astrocytes"
        gene_descriptions: Dict of gene_symbol -> {protein_name, description, ...}

    Returns:
        List of (text, is_gene, gene_info) tuples
    """
    if not cell_type_str or not gene_descriptions:
        return [(str(cell_type_str), False, None)]

    cell_type_str = str(cell_type_str)

    # Split by common delimiters (space, dot, underscore, hyphen)
    # But preserve the delimiters for reconstruction
    parts = re.split(r'([\s\.\-_]+)', cell_type_str)

    result = []
    for part in parts:
        # Check if this part is a gene name
        if part in gene_descriptions:
            result.append((part, True, gene_descriptions[part]))
        else:
            result.append((part, False, None))

    return result


# Cache for gene expression data (per dataset)
_gene_expression_cache = {}


def get_h5ad_path(dataset_name):
    """Get h5ad path for a dataset."""
    from src.datasets import get_config
    config = get_config(dataset_name)
    return config.h5ad_path


def get_gene_expression(dataset_name, gene_name, cell_indices):
    """Load expression values for a gene from the h5ad file."""
    cache_key = (dataset_name, gene_name)
    if cache_key in _gene_expression_cache:
        expr_series = _gene_expression_cache[cache_key]
        return expr_series.reindex(cell_indices)

    try:
        import scanpy as sc
        h5ad_path = get_h5ad_path(dataset_name)
        adata = sc.read_h5ad(h5ad_path, backed='r')

        # For human data, genes are in feature_name column
        from src.datasets import get_config
        config = get_config(dataset_name)

        if config.gene_column and config.gene_column in adata.var.columns:
            # Search by gene symbol in feature_name
            gene_mask = adata.var[config.gene_column] == gene_name
            if gene_mask.sum() == 0:
                return None
            gene_idx = np.where(gene_mask)[0][0]
        else:
            # Search by var_names directly
            if gene_name not in adata.var_names:
                return None
            gene_idx = adata.var_names.get_loc(gene_name)

        # Get expression for all cells
        expr = adata.X[:, gene_idx].toarray().flatten()
        expr_series = pd.Series(expr, index=adata.obs_names)

        # Cache it
        _gene_expression_cache[cache_key] = expr_series

        return expr_series.reindex(cell_indices)
    except Exception as e:
        print(f"Error loading gene expression: {e}")
        return None


def get_region_label(region_name, dataset_name):
    """Get short label for a region."""
    if dataset_name == 'human':
        # Human regions are already abbreviated
        return region_name
    # Mouse: use Allen CCF acronyms
    if region_name in HYPOMAP_TO_CCF:
        return HYPOMAP_TO_CCF[region_name][0]
    return region_name.split()[0][:6]


def register_callbacks(app, datasets, gene_descriptions=None, region_descriptions=None, signaling_data=None):
    """Register all callbacks with the Dash app."""

    if gene_descriptions is None:
        gene_descriptions = {}

    if region_descriptions is None:
        region_descriptions = {}

    if signaling_data is None:
        signaling_data = {
            'ligand_to_receptors': {},
            'cluster_expression': {},
            'ligand_genes': set(),
            'receptor_genes': set(),
        }

    ligand_to_receptors = signaling_data['ligand_to_receptors']
    cluster_expression = signaling_data['cluster_expression']
    ligand_genes = signaling_data['ligand_genes']
    receptor_genes = signaling_data['receptor_genes']

    # Camera positions for different views (closer for better viewport fill)
    # Include center to ensure absolute positioning
    CAMERA_VIEWS = {
        'sagittal': dict(eye=dict(x=-1.8, y=0, z=0), up=dict(x=0, y=1, z=0), center=dict(x=0, y=0, z=0)),
        'coronal': dict(eye=dict(x=0, y=0, z=1.8), up=dict(x=0, y=1, z=0), center=dict(x=0, y=0, z=0)),
        'axial': dict(eye=dict(x=0, y=1.8, z=0), up=dict(x=0, y=0, z=-1), center=dict(x=0, y=0, z=0)),
    }

    # Callback to update controls when dataset changes
    @app.callback(
        [Output('region-dropdown', 'options'),
         Output('region-dropdown', 'value'),
         Output('cell-type-slider', 'max'),
         Output('cell-type-slider', 'marks'),
         Output('cell-type-slider', 'value'),
         Output('gene-dropdown', 'options'),
         Output('gene-dropdown', 'value'),
         Output('dataset-description', 'children'),
         Output('current-dataset', 'data')],
        [Input('dataset-dropdown', 'value')]
    )
    def update_controls_for_dataset(dataset_name):
        """Update all controls when dataset changes."""
        if dataset_name not in datasets:
            dataset_name = list(datasets.keys())[0]

        data = datasets[dataset_name]
        regions = data['regions']
        cell_type_levels = data['cell_type_levels']
        marker_genes = data['marker_genes']
        key_receptors = data['key_receptors']
        cells_df = data['cells_df']

        # Region dropdown
        region_options = [{'label': 'All Regions', 'value': 'all'}] + \
                        [{'label': r, 'value': r} for r in sorted(regions)]

        # Cell type slider
        max_idx = len(cell_type_levels) - 1 if cell_type_levels else 0
        marks = {i: get_cell_type_display_name(level)
                for i, level in enumerate(cell_type_levels)} if cell_type_levels else {}
        default_idx = 1 if len(cell_type_levels) > 1 else 0

        # Gene dropdown
        gene_options = [{'label': g, 'value': g} for g in marker_genes + key_receptors]
        default_gene = marker_genes[0] if marker_genes else None

        # Description (short format for new UI)
        description = f"{len(cells_df):,} cells"

        return (
            region_options, 'all',
            max_idx, marks, default_idx,
            gene_options, default_gene,
            description,
            dataset_name
        )

    @app.callback(
        Output('scatter-3d', 'figure'),
        [Input('current-dataset', 'data'),
         Input('region-dropdown', 'value'),
         Input('color-by-radio', 'value'),
         Input('cell-type-slider', 'value'),
         Input('gene-dropdown', 'value'),
         Input('ligand-dropdown', 'value'),
         Input('cell-filter-options', 'value'),
         Input('display-options', 'value'),
         Input('point-size-slider', 'value'),
         Input('z-jitter-slider', 'value'),
         Input('btn-sagittal', 'n_clicks'),
         Input('btn-coronal', 'n_clicks'),
         Input('btn-axial', 'n_clicks')]
    )
    def update_scatter(dataset_name, region_filter, color_by, cell_type_idx, gene, ligand, cell_filters, display_opts, point_size,
                       z_jitter, btn_sag, btn_cor, btn_ax):
        """Update the main 3D scatter plot."""
        if dataset_name not in datasets:
            return create_empty_figure("Dataset not found")

        data = datasets[dataset_name]
        cells_df = data['cells_df']
        region_col = data['region_col']
        cell_type_levels = data['cell_type_levels']
        region_colors = data['region_colors']

        # Determine camera based on which view button was clicked
        # Only set camera when a view button is explicitly clicked, or on initial load/dataset change
        # Otherwise, let Plotly preserve the current view via uirevision
        triggered_id = ctx.triggered_id

        # Triggers that should reset camera to default view
        reset_camera_triggers = ['btn-sagittal', 'btn-coronal', 'btn-axial', 'current-dataset']

        if triggered_id is None:
            # True initial load (first callback execution) - set default sagittal view
            camera = CAMERA_VIEWS['sagittal']
        elif triggered_id == 'btn-coronal':
            camera = CAMERA_VIEWS['coronal']
        elif triggered_id == 'btn-axial':
            camera = CAMERA_VIEWS['axial']
        elif triggered_id == 'btn-sagittal' or triggered_id == 'current-dataset':
            # Sagittal view button or dataset change - reset to sagittal
            camera = CAMERA_VIEWS['sagittal']
        else:
            # Any other trigger (region filter, color mode, etc.) - don't set camera to preserve current view
            camera = None

        # Filter data
        df = cells_df.copy()
        if region_filter and region_filter != 'all' and region_col:
            df = df[df[region_col] == region_filter]

        # Filter to neurons only if checkbox is checked
        if cell_filters and 'neurons_only' in cell_filters:
            # For ABC data: check 'class' column for neurotransmitter types
            if 'class' in df.columns:
                is_neuron = (
                    df['class'].str.contains('Glut', case=False, na=False) |
                    df['class'].str.contains('GABA', case=False, na=False) |
                    df['class'].str.contains('Dopa', case=False, na=False) |
                    df['class'].str.contains('Sero', case=False, na=False) |
                    df['class'].str.contains('Chol', case=False, na=False)
                )
                df = df[is_neuron]
            # For HypoMap data: check Author_Class_Curated or similar
            elif 'Author_Class_Curated' in df.columns:
                is_neuron = df['Author_Class_Curated'].str.contains('Neuron', case=False, na=False)
                df = df[is_neuron]

        # Remove cells without coordinates
        df = df.dropna(subset=['x', 'y', 'z'])

        if len(df) == 0:
            return create_empty_figure("No cells with coordinates")

        # Apply z-jitter if requested
        # Jitter spreads cells within their slice for better visualization
        if z_jitter and z_jitter > 0:
            df = df.copy()  # Don't modify original
            z_values = df['z'].values
            unique_z = np.sort(np.unique(z_values))

            if len(unique_z) > 1:
                # Compute inter-slice spacing (median difference between consecutive slices)
                z_diffs = np.diff(unique_z)
                slice_spacing = np.median(z_diffs)

                # Max jitter is half the slice spacing (so slices don't overlap)
                max_jitter = slice_spacing * 0.45

                # Generate deterministic jitter based on cell index
                # Using hash of index for reproducibility across renders
                rng = np.random.default_rng(seed=42)
                jitter_values = rng.uniform(-1, 1, size=len(df))

                # Scale by jitter amount and apply
                df['z'] = z_values + jitter_values * max_jitter * z_jitter

        # Determine coloring
        if color_by == 'region' and region_col:
            color_col = region_col
            colors = df[region_col].map(lambda r: region_colors.get(r, '#CCCCCC'))
            colorscale = None
            showscale = False
        elif color_by == 'cell_type' and cell_type_levels:
            idx = min(cell_type_idx, len(cell_type_levels) - 1)
            color_col = cell_type_levels[idx]
            unique_types = df[color_col].unique()
            color_map = {t: px.colors.qualitative.Dark24[i % len(px.colors.qualitative.Dark24)]
                        for i, t in enumerate(unique_types)}
            colors = df[color_col].map(color_map)
            colorscale = None
            showscale = False
        elif color_by == 'gene' and gene:
            expr = get_gene_expression(dataset_name, gene, df.index)
            if expr is not None and not expr.isna().all():
                colors = expr.values
                colorscale = 'Viridis'
                showscale = True
            else:
                colors = 'lightgray'
                colorscale = None
                showscale = False
        elif color_by == 'signaling' and ligand and cluster_expression:
            # Signaling mode: color by receptor expression, mark ligand producers
            receptors = ligand_to_receptors.get(ligand, [])

            # Get receptor score for each cell based on cluster
            def get_receptor_score(cluster):
                expr = cluster_expression.get(cluster, {})
                scores = [expr.get(r, 0) for r in receptors if r in expr]
                return max(scores) if scores else 0

            # Get ligand score for each cell
            def get_ligand_score(cluster):
                expr = cluster_expression.get(cluster, {})
                return expr.get(ligand, 0)

            if 'cluster' in df.columns:
                df = df.copy()
                df['receptor_score'] = df['cluster'].map(get_receptor_score)
                df['ligand_score'] = df['cluster'].map(get_ligand_score)
                colors = df['receptor_score'].values
                colorscale = 'Blues'
                showscale = True
            else:
                colors = 'lightgray'
                colorscale = None
                showscale = False
        else:
            if region_col:
                colors = df[region_col].map(lambda r: region_colors.get(r, '#CCCCCC'))
            else:
                colors = 'steelblue'
            colorscale = None
            showscale = False

        # Create 3D scatter
        fig = go.Figure()

        # Add region boundary ellipsoids if enabled
        if 'show_mesh' in (display_opts or []) and region_col:
            for region in df[region_col].unique():
                if region == 'NA' or region == 'Unknown':
                    continue
                region_df = df[df[region_col] == region]
                if len(region_df) < 10:
                    continue

                coords = region_df[['x', 'y', 'z']].values
                center = coords.mean(axis=0)
                cov = np.cov(coords.T)

                eigenvalues, eigenvectors = np.linalg.eigh(cov)
                radii = 2.0 * np.sqrt(np.maximum(eigenvalues, 1))

                u = np.linspace(0, 2 * np.pi, 20)
                v = np.linspace(0, np.pi, 15)
                x_sphere = np.outer(np.cos(u), np.sin(v))
                y_sphere = np.outer(np.sin(u), np.sin(v))
                z_sphere = np.outer(np.ones_like(u), np.cos(v))

                sphere_pts = np.stack([x_sphere.flatten(), y_sphere.flatten(), z_sphere.flatten()])
                scaled_pts = np.diag(radii) @ sphere_pts
                rotated_pts = eigenvectors @ scaled_pts
                ellipsoid_pts = rotated_pts + center.reshape(3, 1)

                region_color = region_colors.get(region, '#CCCCCC')
                fig.add_trace(go.Surface(
                    x=ellipsoid_pts[0].reshape(20, 15),
                    y=ellipsoid_pts[1].reshape(20, 15),
                    z=ellipsoid_pts[2].reshape(20, 15),
                    opacity=0.2,
                    colorscale=[[0, region_color], [1, region_color]],
                    showscale=False,
                    hoverinfo='skip',
                    showlegend=False
                ))

        # Add cells
        marker_dict = dict(
            size=point_size,
            color=colors,
            opacity=0.7,
            line=dict(width=0)
        )
        if colorscale:
            marker_dict['colorscale'] = colorscale
            marker_dict['showscale'] = showscale
            if color_by == 'signaling' and ligand:
                colorbar_title = f'{ligand} receptor %'
            elif gene:
                colorbar_title = gene
            else:
                colorbar_title = ''
            marker_dict['colorbar'] = dict(title=colorbar_title, thickness=15)

        text_data = df[region_col] if region_col else df.index.astype(str)
        if color_by == 'cell_type' and cell_type_levels:
            idx = min(cell_type_idx, len(cell_type_levels) - 1)
            text_data = df[cell_type_levels[idx]]

        fig.add_trace(go.Scatter3d(
            x=df['x'],
            y=df['y'],
            z=df['z'],
            mode='markers',
            marker=marker_dict,
            customdata=df.index.tolist(),
            hovertemplate=(
                f"<b>{color_by.title()}</b>: %{{text}}<br>" +
                "X: %{x:.0f}<br>" +
                "Y: %{y:.0f}<br>" +
                "Z: %{z:.0f}<br>" +
                "<extra></extra>"
            ),
            text=text_data,
            name='Cells'
        ))

        # Add ligand producer overlay when in signaling mode
        if color_by == 'signaling' and ligand and 'ligand_score' in df.columns:
            ligand_cells = df[df['ligand_score'] > 10]  # >10% expressing
            if len(ligand_cells) > 0:
                fig.add_trace(go.Scatter3d(
                    x=ligand_cells['x'],
                    y=ligand_cells['y'],
                    z=ligand_cells['z'],
                    mode='markers',
                    marker=dict(
                        size=point_size + 2,
                        color='red',
                        symbol='diamond',
                        opacity=0.9
                    ),
                    customdata=ligand_cells.index.tolist(),
                    hovertemplate=(
                        f"<b>{ligand} producer</b><br>" +
                        f"Expression: %{{text:.0f}}%<br>" +
                        "<extra></extra>"
                    ),
                    text=ligand_cells['ligand_score'],
                    name=f'{ligand} producers',
                    showlegend=True
                ))

        # Add region labels
        scene_annotations = []
        if 'show_labels' in (display_opts or []) and region_col:
            # For ABC data, split by hemisphere (left/right based on x coordinate = ML axis)
            x_midline = df['x'].median()

            for region in df[region_col].unique():
                if region in ['NA', 'Unknown', 'HY-unassigned']:
                    continue
                region_df = df[df[region_col] == region]
                if len(region_df) < 5:
                    continue

                label = get_region_label(region, dataset_name)

                # Check if region has substantial populations on both sides of midline
                left_df = region_df[region_df['x'] < x_midline]
                right_df = region_df[region_df['x'] >= x_midline]
                left_frac = len(left_df) / len(region_df)
                right_frac = len(right_df) / len(region_df)

                # Bilateral if at least 25% of cells on each side
                is_bilateral = left_frac >= 0.25 and right_frac >= 0.25

                if is_bilateral:
                    # Label each hemisphere separately
                    for hemi_df in [left_df, right_df]:
                        if len(hemi_df) >= 3:
                            centroid = hemi_df[['x', 'y', 'z']].mean()
                            scene_annotations.append(dict(
                                x=centroid['x'],
                                y=centroid['y'],
                                z=centroid['z'],
                                text=label,
                                showarrow=False,
                                bgcolor='rgba(255, 255, 255, 0.8)',
                                borderpad=2,
                                font=dict(size=9, color='black', family='Arial'),
                                opacity=0.85
                            ))
                else:
                    # Midline or unilateral region - single label
                    centroid = region_df[['x', 'y', 'z']].mean()
                    scene_annotations.append(dict(
                        x=centroid['x'],
                        y=centroid['y'],
                        z=centroid['z'],
                        text=label,
                        showarrow=False,
                        bgcolor='rgba(255, 255, 255, 0.8)',
                        borderpad=2,
                        font=dict(size=9, color='black', family='Arial'),
                        opacity=0.85
                    ))

        # Add axis legend
        x_min, x_max = df['x'].min(), df['x'].max()
        y_min, y_max = df['y'].min(), df['y'].max()
        z_min, _ = df['z'].min(), df['z'].max()

        legend_origin = [x_min - (x_max - x_min) * 0.15, y_max + (y_max - y_min) * 0.1, z_min]
        arrow_len = (x_max - x_min) * 0.12

        for axis_name, axis_color, offset in [
            ('LR', 'red', [arrow_len, 0, 0]),      # Left-Right (medial-lateral)
            ('DV', 'green', [0, -arrow_len, 0]),   # Dorsal-Ventral
            ('AP', 'blue', [0, 0, arrow_len])      # Anterior-Posterior
        ]:
            fig.add_trace(go.Scatter3d(
                x=[legend_origin[0], legend_origin[0] + offset[0]],
                y=[legend_origin[1], legend_origin[1] + offset[1]],
                z=[legend_origin[2], legend_origin[2] + offset[2]],
                mode='lines+text',
                line=dict(color=axis_color, width=4),
                text=['', axis_name],
                textposition='middle right',
                textfont=dict(size=10, color=axis_color),
                hoverinfo='skip',
                showlegend=False
            ))

        # Build scene dict - only include camera if a view button was clicked
        scene_dict = dict(
            xaxis=dict(visible=False),
            yaxis=dict(autorange='reversed', visible=False),
            zaxis=dict(visible=False),
            aspectmode='data',
            annotations=scene_annotations
        )
        if camera is not None:
            scene_dict['camera'] = camera

        # Determine uirevision - changes when we want to reset camera, stays constant otherwise
        if triggered_id is None:
            uirevision = 'initial'
        elif triggered_id in reset_camera_triggers:
            uirevision = triggered_id
        else:
            uirevision = 'constant'

        fig.update_layout(
            scene=scene_dict,
            margin=dict(l=0, r=0, t=30, b=0),
            uirevision=uirevision
        )

        return fig

    def create_cell_type_display(cell_type_str):
        """Create display elements for a cell type with highlighted genes.

        Returns:
            Tuple of (elements list, set of gene names found)
        """
        parsed = parse_cell_type_with_genes(cell_type_str, gene_descriptions)

        elements = []
        genes_found = set()

        for text, is_gene, gene_info in parsed:
            if is_gene and gene_info:
                genes_found.add(text)
                elements.append(
                    html.Span(
                        text,
                        style={
                            'color': '#6f42c1',
                            'fontWeight': 'bold',
                        }
                    )
                )
            else:
                elements.append(html.Span(text))

        return elements, genes_found

    @app.callback(
        [Output('cell-prompt', 'style'),
         Output('cell-info', 'style'),
         Output('cell-region', 'children'),
         Output('cell-coords', 'children'),
         Output('cell-types', 'children')],
        [Input('scatter-3d', 'clickData'),
         Input('current-dataset', 'data')]
    )
    def update_selection(click_data, dataset_name):
        """Update the details panel based on single cell click."""
        if dataset_name not in datasets:
            return ({'display': 'block'}, {'display': 'none'}, "", "", [])

        data = datasets[dataset_name]
        cells_df = data['cells_df']
        region_col = data['region_col']
        cell_type_levels = data['cell_type_levels']

        if click_data and click_data.get('points'):
            point = click_data['points'][0]
            point_idx = point.get('customdata')
            if point_idx is not None and point_idx in cells_df.index:
                cell = cells_df.loc[point_idx]

                region_acronym = cell[region_col] if region_col else "N/A"

                # Build region display with full name and description
                if region_acronym in region_descriptions:
                    region_info = region_descriptions[region_acronym]
                    region = html.Div([
                        html.Span(region_acronym, className="fw-bold"),
                        html.Span(f" - {region_info['full_name']}", className="text-muted"),
                        html.P(
                            region_info['description'],
                            className="text-muted small mt-1 mb-0",
                            style={'fontSize': '0.8rem', 'lineHeight': '1.3'}
                        ) if region_info['description'] else None
                    ])
                else:
                    region = region_acronym

                coords = f"X: {cell['x']:.0f}\nY: {cell['y']:.0f}\nZ: {cell['z']:.0f}"

                type_items = []
                all_genes = set()

                for level in cell_type_levels:
                    if level in cell.index:
                        # Format level name for display
                        if '_named' in level:
                            level_name = level.replace('_named', '').replace('C', '')
                        else:
                            level_name = level.capitalize()

                        # Create display with highlighted genes
                        cell_type_elements, genes_found = create_cell_type_display(cell[level])
                        all_genes.update(genes_found)

                        type_items.append(
                            html.Div([
                                html.Span(f"{level_name}: ", className="text-muted"),
                                html.Span(cell_type_elements, className="fw-bold")
                            ], className="small mb-1")
                        )

                # Add gene descriptions section if any genes found
                if all_genes and gene_descriptions:
                    type_items.append(html.Hr(style={'margin': '8px 0'}))
                    type_items.append(
                        html.Div("Marker Genes:", className="text-muted small fw-bold mb-1")
                    )

                    for gene in sorted(all_genes):
                        if gene in gene_descriptions:
                            info = gene_descriptions[gene]
                            protein = info.get('protein_name', '')
                            desc = info.get('description', '')

                            # Create expandable description if long
                            if desc and len(desc) > 150:
                                truncated = desc[:150] + "..."
                                gene_id = f"desc-{gene}-{hash(point_idx) % 10000}"
                                desc_element = html.Div([
                                    html.Span(
                                        truncated,
                                        id={'type': 'gene-desc-short', 'gene': gene_id},
                                        style={'cursor': 'pointer'}
                                    ),
                                    html.Span(
                                        " [more]",
                                        style={'color': '#6f42c1', 'cursor': 'pointer', 'fontSize': '0.7rem'}
                                    ),
                                    dbc.Collapse(
                                        html.Div(desc, style={'marginTop': '4px'}),
                                        id={'type': 'gene-desc-full', 'gene': gene_id},
                                        is_open=False
                                    )
                                ], id={'type': 'gene-desc-container', 'gene': gene_id}, n_clicks=0)
                            elif desc:
                                desc_element = html.Div(desc)
                            else:
                                desc_element = ""

                            type_items.append(
                                html.Div([
                                    html.Span(gene, style={'color': '#6f42c1', 'fontWeight': 'bold'}),
                                    html.Span(f" - {protein}", className="text-muted") if protein else "",
                                    html.Div(
                                        desc_element,
                                        className="text-muted",
                                        style={'fontSize': '0.75rem', 'marginLeft': '8px', 'marginBottom': '4px'}
                                    ) if desc else ""
                                ], className="small mb-1")
                            )

                # Add signaling profile section if we have cluster info
                if 'cluster' in cell.index and cluster_expression:
                    cluster = cell['cluster']
                    cluster_expr = cluster_expression.get(cluster, {})

                    if cluster_expr:
                        # Get top ligands expressed by this cluster
                        ligands_expressed = [
                            (g, pct) for g, pct in cluster_expr.items()
                            if g in ligand_genes and pct > 10
                        ]
                        ligands_expressed.sort(key=lambda x: -x[1])

                        # Get top receptors expressed by this cluster
                        receptors_expressed = [
                            (g, pct) for g, pct in cluster_expr.items()
                            if g in receptor_genes and pct > 10
                        ]
                        receptors_expressed.sort(key=lambda x: -x[1])

                        if ligands_expressed or receptors_expressed:
                            type_items.append(html.Hr(style={'margin': '8px 0'}))
                            type_items.append(
                                html.Div("Signaling Profile:", className="text-muted small fw-bold mb-1")
                            )

                            if ligands_expressed:
                                type_items.append(
                                    html.Div([
                                        html.Span("Ligands: ", className="text-muted small"),
                                        html.Span(
                                            ", ".join([f"{g} ({pct:.0f}%)" for g, pct in ligands_expressed[:5]]),
                                            style={'color': '#dc3545', 'fontSize': '0.8rem'}
                                        )
                                    ], className="mb-1")
                                )

                            if receptors_expressed:
                                type_items.append(
                                    html.Div([
                                        html.Span("Receptors: ", className="text-muted small"),
                                        html.Span(
                                            ", ".join([f"{g} ({pct:.0f}%)" for g, pct in receptors_expressed[:5]]),
                                            style={'color': '#0d6efd', 'fontSize': '0.8rem'}
                                        )
                                    ], className="mb-1")
                                )

                return (
                    {'display': 'none'},
                    {'display': 'block'},
                    region,
                    coords,
                    type_items
                )

        return ({'display': 'block'}, {'display': 'none'}, "", "", [])

    @app.callback(
        Output('cell-type-slider-container', 'style'),
        [Input('color-by-radio', 'value')]
    )
    def toggle_cell_type_slider(color_by):
        if color_by == 'cell_type':
            return {'display': 'block'}
        return {'display': 'none'}

    @app.callback(
        Output('gene-selector-container', 'style'),
        [Input('color-by-radio', 'value')]
    )
    def toggle_gene_selector(color_by):
        if color_by == 'gene':
            return {'display': 'block'}
        return {'display': 'none'}

    @app.callback(
        Output('ligand-selector-container', 'style'),
        [Input('color-by-radio', 'value')]
    )
    def toggle_ligand_selector(color_by):
        if color_by == 'signaling':
            return {'display': 'block'}
        return {'display': 'none'}

    # Callback to expand/collapse gene descriptions
    @app.callback(
        [Output({'type': 'gene-desc-full', 'gene': MATCH}, 'is_open'),
         Output({'type': 'gene-desc-short', 'gene': MATCH}, 'style')],
        [Input({'type': 'gene-desc-container', 'gene': MATCH}, 'n_clicks')],
        [State({'type': 'gene-desc-full', 'gene': MATCH}, 'is_open')],
        prevent_initial_call=True
    )
    def toggle_gene_description(n_clicks, is_open):
        if n_clicks:
            new_is_open = not is_open
            # Hide truncated text when expanded
            style = {'display': 'none'} if new_is_open else {'cursor': 'pointer'}
            return new_is_open, style
        return is_open, {'cursor': 'pointer'}

    # Callback to toggle controls panel collapse
    @app.callback(
        [Output('controls-collapse', 'is_open'),
         Output('collapse-controls-btn', 'children'),
         Output('controls-col', 'width'),
         Output('viz-col', 'width')],
        [Input('collapse-controls-btn', 'n_clicks')],
        [State('controls-collapse', 'is_open')]
    )
    def toggle_controls_collapse(n_clicks, is_open):
        if n_clicks:
            new_is_open = not is_open
        else:
            new_is_open = is_open

        if new_is_open:
            # Expanded: controls gets 2 cols, viz gets 7
            return new_is_open, "◀", 2, 7
        else:
            # Collapsed: controls gets 1 col (just button), viz gets 8
            return new_is_open, "▶", 1, 8


def get_cell_type_display_name(col_name):
    """Convert column name to display name.

    Handles both HypoMap style ('C66_named' -> '66 types')
    and ABC style ('subclass' -> 'Subclass').
    """
    if col_name and '_named' in col_name:
        num = col_name.split('_')[0][1:]
        return f"{num} types"
    # ABC-style: capitalize
    if col_name in ['class', 'subclass', 'supertype', 'cluster']:
        return col_name.capitalize()
    return col_name


def create_empty_figure(message="No data"):
    """Create an empty placeholder figure."""
    fig = go.Figure()
    fig.update_layout(
        xaxis={'visible': False},
        yaxis={'visible': False},
        annotations=[{
            'text': message,
            'xref': 'paper',
            'yref': 'paper',
            'x': 0.5,
            'y': 0.5,
            'showarrow': False,
            'font': {'size': 12, 'color': 'gray'}
        }],
        margin=dict(l=0, r=0, t=0, b=0),
        paper_bgcolor='rgba(0,0,0,0)',
        plot_bgcolor='rgba(0,0,0,0)'
    )
    return fig
