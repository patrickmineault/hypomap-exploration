"""Dash callbacks for HypoMap 3D Viewer interactivity."""

from dash import Input, Output, ctx
from dash import html
import plotly.graph_objects as go
import plotly.express as px
import pandas as pd
import numpy as np

import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

from src.atlas.region_mapping import (
    HYPOMAP_TO_CCF,
)


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


def register_callbacks(app, datasets):
    """Register all callbacks with the Dash app."""

    # Camera positions for different views
    CAMERA_VIEWS = {
        'sagittal': dict(eye=dict(x=-2.5, y=0, z=0), up=dict(x=0, y=1, z=0)),
        'coronal': dict(eye=dict(x=0, y=0, z=2.5), up=dict(x=0, y=1, z=0)),
        'axial': dict(eye=dict(x=0, y=2.5, z=0), up=dict(x=0, y=0, z=-1)),
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

        # Description
        description = f"Interactive 3D visualization of hypothalamic cell types ({len(cells_df)} cells)"

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
         Input('display-options', 'value'),
         Input('point-size-slider', 'value'),
         Input('btn-sagittal', 'n_clicks'),
         Input('btn-coronal', 'n_clicks'),
         Input('btn-axial', 'n_clicks')]
    )
    def update_scatter(dataset_name, region_filter, color_by, cell_type_idx, gene, display_opts, point_size,
                       btn_sag, btn_cor, btn_ax):
        """Update the main 3D scatter plot."""
        if dataset_name not in datasets:
            return create_empty_figure("Dataset not found")

        data = datasets[dataset_name]
        cells_df = data['cells_df']
        region_col = data['region_col']
        cell_type_levels = data['cell_type_levels']
        region_colors = data['region_colors']

        # Determine which view button was clicked
        triggered_id = ctx.triggered_id if ctx.triggered_id else 'btn-sagittal'
        if triggered_id == 'btn-coronal':
            camera = CAMERA_VIEWS['coronal']
        elif triggered_id == 'btn-axial':
            camera = CAMERA_VIEWS['axial']
        else:
            camera = CAMERA_VIEWS['sagittal']

        # Filter data
        df = cells_df.copy()
        if region_filter and region_filter != 'all' and region_col:
            df = df[df[region_col] == region_filter]

        # Remove cells without coordinates
        df = df.dropna(subset=['x', 'y', 'z'])

        if len(df) == 0:
            return create_empty_figure("No cells with coordinates")

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
            marker_dict['colorbar'] = dict(title=gene if gene else '', thickness=15)

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

        # Add region labels
        scene_annotations = []
        if 'show_labels' in (display_opts or []) and region_col:
            for region in df[region_col].unique():
                if region in ['NA', 'Unknown']:
                    continue
                region_df = df[df[region_col] == region]
                if len(region_df) < 5:
                    continue
                centroid = region_df[['x', 'y', 'z']].mean()
                label = get_region_label(region, dataset_name)
                scene_annotations.append(dict(
                    x=centroid['x'],
                    y=centroid['y'],
                    z=centroid['z'],
                    text=label,
                    showarrow=False,
                    bgcolor='rgba(255, 255, 255, 0.8)',
                    borderpad=3,
                    font=dict(size=11, color='black', family='Arial'),
                    opacity=0.9
                ))

        # Add axis legend
        x_min, x_max = df['x'].min(), df['x'].max()
        y_min, y_max = df['y'].min(), df['y'].max()
        z_min, _ = df['z'].min(), df['z'].max()

        legend_origin = [x_min - (x_max - x_min) * 0.15, y_max + (y_max - y_min) * 0.1, z_min]
        arrow_len = (x_max - x_min) * 0.12

        for axis_name, axis_color, offset in [
            ('X', 'red', [arrow_len, 0, 0]),
            ('Y', 'green', [0, -arrow_len, 0]),
            ('Z', 'blue', [0, 0, arrow_len])
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

        fig.update_layout(
            scene=dict(
                xaxis=dict(visible=False),
                yaxis=dict(autorange='reversed', visible=False),
                zaxis=dict(visible=False),
                aspectmode='data',
                camera=camera,
                annotations=scene_annotations
            ),
            margin=dict(l=0, r=0, t=30, b=0),
            uirevision=triggered_id if triggered_id in ['btn-sagittal', 'btn-coronal', 'btn-axial'] else 'constant'
        )

        return fig

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

                region = cell[region_col] if region_col else "N/A"
                coords = f"X: {cell['x']:.0f}\nY: {cell['y']:.0f}\nZ: {cell['z']:.0f}"

                type_items = []
                for level in cell_type_levels:
                    if level in cell.index:
                        level_name = level.replace('_named', '').replace('C', '')
                        type_items.append(
                            html.Div([
                                html.Span(f"{level_name}: ", className="text-muted"),
                                html.Span(str(cell[level]), className="fw-bold")
                            ], className="small mb-1")
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


def get_cell_type_display_name(col_name):
    """Convert column name like 'C66_named' to display name."""
    if col_name and '_named' in col_name:
        num = col_name.split('_')[0][1:]
        return f"{num} types"
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
