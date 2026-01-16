"""Dash UI layout components for HypoMap 3D Viewer."""

from dash import dcc, html
import dash_bootstrap_components as dbc


def get_cell_type_display_name(col_name):
    """Convert column name to display name.

    Handles both HypoMap style ('C66_named' -> '66 types')
    and ABC style ('subclass' -> 'Subclass').
    """
    if col_name and '_named' in col_name:
        num = col_name.split('_')[0][1:]  # Extract number from C66_named -> 66
        return f"{num} types"
    # ABC-style: capitalize
    if col_name in ['class', 'subclass', 'supertype', 'cluster']:
        return col_name.capitalize()
    return col_name


def get_dataset_display_name(name):
    """Get display name for a dataset."""
    display_names = {
        'mouse_hypomap': 'Mouse HypoMap',
        'human_hypomap': 'Human HypoMap',
        'mouse_abc': 'Mouse ABC',
    }
    return display_names.get(name, name.replace('_', ' ').title())


def create_layout(datasets, available_datasets, default_dataset, ligand_options=None):
    """
    Create the main application layout.

    Args:
        datasets: Dictionary of loaded datasets
        available_datasets: List of available dataset names
        default_dataset: Default dataset to show
        ligand_options: List of available ligands for signaling mode
    """
    default_data = datasets[default_dataset]
    regions = default_data['regions']
    cell_type_levels = default_data['cell_type_levels']
    marker_genes = default_data['marker_genes']
    key_receptors = default_data['key_receptors']
    ligand_options = ligand_options or []

    return dbc.Container([
        # Store current dataset name
        dcc.Store(id='current-dataset', data=default_dataset),

        # Google Font import for cool title
        html.Link(
            href="https://fonts.googleapis.com/css2?family=Orbitron:wght@700&display=swap",
            rel="stylesheet"
        ),

        # Header with centered title and dataset selector
        html.Div([
            html.H1(
                "HypoMap 3D Atlas",
                style={
                    'fontFamily': '"Orbitron", sans-serif',
                    'fontSize': '2.5rem',
                    'fontWeight': '700',
                    'background': 'linear-gradient(135deg, #667eea 0%, #764ba2 100%)',
                    'WebkitBackgroundClip': 'text',
                    'WebkitTextFillColor': 'transparent',
                    'backgroundClip': 'text',
                    'textAlign': 'center',
                    'margin': '0.5rem 0',
                    'letterSpacing': '2px',
                }
            ),
            html.P(
                id='dataset-description',
                className="text-muted mb-2",
                style={'textAlign': 'center', 'fontSize': '0.9rem'},
                children=f"{len(default_data['cells_df']):,} cells"
            ),
            # Dataset selector - small, top right
            html.Div([
                dcc.Dropdown(
                    id='dataset-dropdown',
                    options=[
                        {'label': get_dataset_display_name(name), 'value': name}
                        for name in available_datasets
                    ],
                    value=default_dataset,
                    clearable=False,
                    style={'width': '150px', 'fontSize': '0.85rem'}
                ),
            ], style={'position': 'absolute', 'top': '10px', 'right': '15px'})
        ], style={'position': 'relative', 'paddingTop': '5px'}),

        # Main content: 3 columns
        dbc.Row([
            # Left panel: Controls (collapsible)
            dbc.Col([
                # Header with collapse button on right
                html.Div([
                    html.Span("Controls", className="fw-bold", style={'fontSize': '0.9rem'}),
                    dbc.Button(
                        "◀",
                        id="collapse-controls-btn",
                        color="light",
                        size="sm",
                        style={'padding': '2px 8px', 'fontSize': '0.7rem'}
                    ),
                ], className="d-flex justify-content-between align-items-center mb-2"),
                dbc.Collapse(
                    create_controls_panel(regions, cell_type_levels, marker_genes, key_receptors, ligand_options),
                    id="controls-collapse",
                    is_open=True,
                ),
            ], id="controls-col", className="bg-light p-2", style={'transition': 'all 0.3s'}),

            # Center: 3D visualization
            dbc.Col([
                create_visualization_panel()
            ], id="viz-col", style={'transition': 'all 0.3s'}),

            # Right panel: Details
            dbc.Col([
                create_details_panel()
            ], width=3, className="bg-light p-3"),
        ], className="vh-75"),

    ], fluid=True)


def create_controls_panel(regions, cell_type_levels, marker_genes, key_receptors, ligand_options=None):
    """Create the left control panel."""
    ligand_options = ligand_options or []

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
                {'label': '  Cell Type', 'value': 'cell_type'},
                {'label': '  Region', 'value': 'region'},
                {'label': '  Gene Expression', 'value': 'gene'},
                {'label': '  Signaling', 'value': 'signaling'},
            ],
            value='cell_type',
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

        # Ligand selection (shown when color by signaling)
        html.Div([
            html.Label("Select Ligand", className="fw-bold"),
            dcc.Dropdown(
                id='ligand-dropdown',
                options=[{'label': g, 'value': g} for g in ligand_options],
                value=ligand_options[0] if ligand_options else None,
                placeholder="Search ligands...",
                className="mb-2"
            ),
            html.P(
                "Shows receptor expression (blue) and ligand producers (red diamonds)",
                className="text-muted small"
            ),
        ], id='ligand-selector-container'),

        html.Hr(),

        # Cell filter
        html.Label("Cell Filter", className="fw-bold"),
        dcc.Checklist(
            id='cell-filter-options',
            options=[
                {'label': ' Neurons only', 'value': 'neurons_only'},
            ],
            value=['neurons_only'],
            className="mb-3"
        ),

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

        html.Div(className="mb-3"),

        # Z-jitter (for MERFISH/spatial data with discrete slices)
        html.Label("Z-Jitter", className="fw-bold"),
        html.P("Spread cells within slices", className="text-muted small mb-1"),
        dcc.Slider(
            id='z-jitter-slider',
            min=0,
            max=1,
            step=0.05,
            value=1,
            marks={0: '0', 0.5: '0.5', 1: '1'},
            tooltip={"placement": "bottom", "always_visible": False}
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
            style={'height': '80vh'},
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
