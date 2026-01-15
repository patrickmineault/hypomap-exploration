"""Dash UI layout components for HypoMap 3D Viewer."""

from dash import dcc, html
import dash_bootstrap_components as dbc


def get_cell_type_display_name(col_name):
    """Convert column name like 'C66_named' to display name like '66 types'."""
    if col_name and '_named' in col_name:
        num = col_name.split('_')[0][1:]  # Extract number from C66_named -> 66
        return f"{num} types"
    return col_name


def get_dataset_display_name(name):
    """Get display name for a dataset."""
    display_names = {
        'mouse': 'Mouse Hypothalamus',
        'human': 'Human Hypothalamus',
    }
    return display_names.get(name, name.title())


def create_layout(datasets, available_datasets, default_dataset):
    """
    Create the main application layout.

    Args:
        datasets: Dictionary of loaded datasets
        available_datasets: List of available dataset names
        default_dataset: Default dataset to show
    """
    default_data = datasets[default_dataset]
    regions = default_data['regions']
    cell_type_levels = default_data['cell_type_levels']
    marker_genes = default_data['marker_genes']
    key_receptors = default_data['key_receptors']

    return dbc.Container([
        # Store current dataset name
        dcc.Store(id='current-dataset', data=default_dataset),

        # Header with dataset selector
        dbc.Row([
            dbc.Col([
                html.H2("HypoMap 3D Atlas Viewer", className="text-center my-3"),
            ], width=8),
            dbc.Col([
                html.Div([
                    html.Label("Dataset:", className="me-2 fw-bold"),
                    dcc.Dropdown(
                        id='dataset-dropdown',
                        options=[
                            {'label': get_dataset_display_name(name), 'value': name}
                            for name in available_datasets
                        ],
                        value=default_dataset,
                        clearable=False,
                        style={'width': '200px', 'display': 'inline-block'}
                    ),
                ], className="d-flex align-items-center justify-content-end mt-3")
            ], width=4),
        ]),

        dbc.Row([
            dbc.Col([
                html.P(
                    id='dataset-description',
                    className="text-center text-muted",
                    children=f"Interactive 3D visualization of hypothalamic cell types ({len(default_data['cells_df'])} cells)"
                )
            ])
        ]),

        # Main content: 3 columns
        dbc.Row([
            # Left panel: Controls
            dbc.Col([
                create_controls_panel(regions, cell_type_levels, marker_genes, key_receptors)
            ], width=3, className="bg-light p-3"),

            # Center: 3D visualization
            dbc.Col([
                create_visualization_panel()
            ], width=6),

            # Right panel: Details
            dbc.Col([
                create_details_panel()
            ], width=3, className="bg-light p-3"),
        ], className="vh-75"),

    ], fluid=True)


def create_controls_panel(regions, cell_type_levels, marker_genes, key_receptors):
    """Create the left control panel."""
    return html.Div([
        html.H5("Controls", className="mb-3"),

        # Region filter
        html.Label("Filter by Region", className="fw-bold"),
        dcc.Dropdown(
            id='region-dropdown',
            options=[{'label': 'All Regions', 'value': 'all'}] +
                    [{'label': r, 'value': r} for r in sorted(regions)],
            value='all',
            clearable=False,
            className="mb-3"
        ),

        html.Hr(),

        # Color by selection
        html.Label("Color cells by", className="fw-bold"),
        dcc.RadioItems(
            id='color-by-radio',
            options=[
                {'label': '  Region', 'value': 'region'},
                {'label': '  Cell Type', 'value': 'cell_type'},
                {'label': '  Gene Expression', 'value': 'gene'},
            ],
            value='region',
            className="mb-3",
            labelStyle={'display': 'block', 'marginBottom': '8px', 'cursor': 'pointer'}
        ),

        # Cell type level selector (shown when color by cell type)
        html.Div([
            html.Label("Cell Type Level", className="fw-bold"),
            dcc.Slider(
                id='cell-type-slider',
                min=0,
                max=len(cell_type_levels) - 1 if cell_type_levels else 0,
                step=1,
                marks={i: get_cell_type_display_name(level) for i, level in enumerate(cell_type_levels)} if cell_type_levels else {},
                value=1 if len(cell_type_levels) > 1 else 0,
            ),
        ], id='cell-type-slider-container', className="mb-3"),

        # Gene selection (shown when color by gene)
        html.Div([
            html.Label("Select Gene", className="fw-bold"),
            dcc.Dropdown(
                id='gene-dropdown',
                options=[{'label': g, 'value': g} for g in marker_genes + key_receptors],
                value=marker_genes[0] if marker_genes else None,
                placeholder="Search genes...",
                className="mb-2"
            ),
        ], id='gene-selector-container'),

        html.Hr(),

        # Display options
        html.Label("Display Options", className="fw-bold"),
        dcc.Checklist(
            id='display-options',
            options=[
                {'label': ' Show region boundaries', 'value': 'show_mesh'},
                {'label': ' Show region labels', 'value': 'show_labels'},
            ],
            value=['show_labels'],
            className="mb-3"
        ),

        # Point size
        html.Label("Point Size", className="fw-bold"),
        dcc.Slider(
            id='point-size-slider',
            min=1,
            max=10,
            step=1,
            value=3,
            marks={1: '1', 5: '5', 10: '10'}
        ),

    ])


def create_visualization_panel():
    """Create the center 3D visualization panel."""
    return html.Div([
        # View buttons
        html.Div([
            html.Label("View: ", className="me-2", style={'fontWeight': 'bold'}),
            dbc.ButtonGroup([
                dbc.Button("Sagittal", id="btn-sagittal", color="primary", size="sm", outline=True),
                dbc.Button("Coronal", id="btn-coronal", color="primary", size="sm", outline=True),
                dbc.Button("Axial", id="btn-axial", color="primary", size="sm", outline=True),
            ], size="sm"),
        ], className="d-flex align-items-center justify-content-center mb-2"),

        dcc.Graph(
            id='scatter-3d',
            style={'height': '68vh'},
            config={
                'displayModeBar': True,
                'scrollZoom': True,
                'displaylogo': False
            }
        ),

        # Store for camera state
        dcc.Store(id='camera-store', data={'view': 'sagittal'}),

        # Loading indicator
        dcc.Loading(
            id="loading-indicator",
            type="circle",
            children=html.Div(id="loading-output")
        ),

    ])


def create_details_panel():
    """Create the right details panel for single cell info."""
    return html.Div([
        html.H5("Cell Details", className="mb-3"),

        # Prompt when no cell selected
        html.Div(id='cell-prompt', children=[
            html.P("Click on a cell to view details", className="text-muted"),
        ]),

        # Cell info (populated on click)
        html.Div(id='cell-info', style={'display': 'none'}, children=[
            # Region
            html.Div([
                html.Label("Region", className="fw-bold text-muted small"),
                html.P(id='cell-region', className="mb-2"),
            ]),

            # Coordinates
            html.Div([
                html.Label("Coordinates", className="fw-bold text-muted small"),
                html.P(id='cell-coords', className="mb-2 font-monospace small"),
            ]),

            html.Hr(),

            # Cell type hierarchy
            html.Label("Cell Type Hierarchy", className="fw-bold text-muted small"),
            html.Div(id='cell-types', className="mb-2"),
        ]),
    ])


def create_empty_figure():
    """Create an empty placeholder figure."""
    import plotly.graph_objects as go
    fig = go.Figure()
    fig.update_layout(
        xaxis={'visible': False},
        yaxis={'visible': False},
        annotations=[{
            'text': 'No data selected',
            'xref': 'paper',
            'yref': 'paper',
            'showarrow': False,
            'font': {'size': 14, 'color': 'gray'}
        }],
        margin=dict(l=0, r=0, t=0, b=0)
    )
    return fig
