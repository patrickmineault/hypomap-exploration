"""Layout components for HypoMap Coronal Atlas Viewer.

Note: Custom CSS is in app/assets/custom.css (auto-loaded by Dash).
"""

from dash import html, dcc
import dash_bootstrap_components as dbc


DATASET_LABELS = {
    'mouse_abc': 'Hypothalamus',
    'mouse_abc_extended': 'Extended',
}


def create_left_panel(cell_type_levels, nt_types, np_system_names, hormone_names=None, region_list=None, region_descriptions=None, enable_region_highlight=False, dataset_names=None, default_dataset='mouse_abc', default_subsample=30):
    """Create the left control panel."""
    hormone_names = hormone_names or []
    region_list = region_list or []
    region_descriptions = region_descriptions or {}
    dataset_names = dataset_names or ['mouse_abc']

    return dbc.Card([
        dbc.CardHeader(html.H5("Controls", className="mb-0")),
        dbc.CardBody([
            # Dataset selector (only show if multiple datasets available)
            *([
                html.Div([
                    html.Label("Dataset", className="control-label"),
                    dcc.RadioItems(
                        id='dataset-radio',
                        options=[
                            {'label': f' {DATASET_LABELS.get(name, name)}', 'value': name}
                            for name in dataset_names
                        ],
                        value=default_dataset,
                        labelStyle={'display': 'block', 'marginBottom': '6px', 'fontSize': '0.9rem'},
                        inputStyle={'marginRight': '10px'},
                    ),
                ], className="mb-4"),
                html.Hr(style={'borderColor': 'rgba(0,0,0,0.08)'}),
            ] if len(dataset_names) > 1 else []),
            # Always include a hidden radio for callback compatibility when single dataset
            *([
                dcc.RadioItems(
                    id='dataset-radio',
                    options=[{'label': default_dataset, 'value': default_dataset}],
                    value=default_dataset,
                    style={'display': 'none'},
                ),
            ] if len(dataset_names) <= 1 else []),

            # Mode selector
            html.Div([
                html.Label("Visualization Mode", className="control-label"),
                dcc.RadioItems(
                    id='viz-mode',
                    options=[
                        {'label': ' Neuropeptide System', 'value': 'np'},
                        {'label': ' Neurotransmitter System', 'value': 'nt'},
                    ] + ([{'label': ' Hormone Receptor', 'value': 'hormone'}] if hormone_names else []) + [
                        {'label': ' Cluster Granularity', 'value': 'cluster'},
                    ],
                    value='np',
                    labelStyle={'display': 'block', 'marginBottom': '8px', 'fontSize': '0.9rem'},
                    inputStyle={'marginRight': '10px'},
                ),
            ], className="mb-4"),

            # Cluster mode controls
            html.Div(
                id='cluster-controls',
                children=[
                    html.Label("Cell Type Level", className="control-label"),
                    dcc.Dropdown(
                        id='cluster-level',
                        options=[{'label': level.capitalize(), 'value': level}
                                 for level in cell_type_levels],
                        value='subclass',
                        clearable=False,
                        className="modern-dropdown",
                    ),
                ],
                className="mb-4",
            ),

            # NT mode controls
            html.Div(
                id='nt-controls',
                children=[
                    html.Label("Neurotransmitter Type", className="control-label"),
                    dcc.Dropdown(
                        id='nt-system',
                        options=[{'label': nt, 'value': nt} for nt in nt_types],
                        value=nt_types[0] if nt_types else None,
                        clearable=False,
                        className="modern-dropdown",
                    ),
                ],
                className="mb-4",
                style={'display': 'none'},
            ),

            # Hormone mode controls (always include dropdown for callback compatibility)
            html.Div(
                id='hormone-controls',
                children=[
                    html.Label("Hormone System", className="control-label"),
                    dcc.Dropdown(
                        id='hormone-system',
                        options=[{'label': name, 'value': name} for name in hormone_names] if hormone_names else [],
                        value=hormone_names[0] if hormone_names else None,
                        clearable=False,
                        className="modern-dropdown",
                    ),
                ],
                className="mb-4",
                style={'display': 'none'},
            ),

            # NP mode controls
            html.Div(
                id='np-controls',
                children=[
                    html.Label("Neuropeptide System", className="control-label"),
                    dcc.Dropdown(
                        id='np-system',
                        options=[{'label': name, 'value': name} for name in np_system_names],
                        value='Orexin' if 'Orexin' in np_system_names else np_system_names[0] if np_system_names else None,
                        clearable=False,
                        className="modern-dropdown",
                    ),
                    # System info display
                    html.Div(id='np-system-info', className="mt-2", style={'fontSize': '0.8rem', 'color': '#666'}),
                    # Diffusion range filter
                    html.Div([
                        dcc.Checklist(
                            id='diffusion-filter-enabled',
                            options=[{'label': ' Filter receptors by diffusion range', 'value': 'enabled'}],
                            value=[],
                            inputStyle={'marginRight': '8px'},
                            style={'fontSize': '0.85rem'},
                        ),
                        html.Div(
                            id='diffusion-range-container',
                            children=[
                                html.Label("Diffusion Range (mm)", className="control-label mt-2"),
                                dcc.Slider(
                                    id='diffusion-range',
                                    min=0.05,
                                    max=1.0,
                                    step=0.05,
                                    value=0.5,
                                    marks={0.1: '0.1', 0.25: '0.25', 0.5: '0.5', 0.75: '0.75', 1.0: '1.0'},
                                    tooltip={"placement": "bottom", "always_visible": False},
                                ),
                            ],
                            style={'display': 'none'},
                        ),
                    ], className="mt-3"),
                    # Rainbow mode checkbox
                    html.Div([
                        dcc.Checklist(
                            id='rainbow-mode',
                            options=[{'label': ' Rainbow mode (fill by cluster)', 'value': 'enabled'}],
                            value=[],
                            inputStyle={'marginRight': '8px'},
                            style={'fontSize': '0.85rem'},
                        ),
                    ], className="mt-3"),
                    # NP mode legend
                    html.Div([
                        html.Label("Legend", className="control-label mt-3"),
                        html.Div([
                            html.Div([
                                html.Span(style={'display': 'inline-block', 'width': '14px', 'height': '14px', 'backgroundColor': '#e63946', 'marginRight': '6px', 'borderRadius': '50%', 'verticalAlign': 'middle'}),
                                html.Span("Ligand", style={'fontSize': '0.85rem'}),
                            ], style={'marginBottom': '4px'}),
                            html.Div([
                                html.Span(style={'display': 'inline-block', 'width': '14px', 'height': '14px', 'backgroundColor': '#3b5dc9', 'marginRight': '6px', 'borderRadius': '50%', 'verticalAlign': 'middle'}),
                                html.Span("Receptor", style={'fontSize': '0.85rem'}),
                            ], style={'marginBottom': '4px'}),
                            html.Div([
                                html.Span(style={'display': 'inline-block', 'width': '14px', 'height': '14px', 'backgroundColor': '#9b5de5', 'marginRight': '6px', 'borderRadius': '50%', 'verticalAlign': 'middle'}),
                                html.Span("Ligand + Receptor", style={'fontSize': '0.85rem'}),
                            ], style={'marginBottom': '4px'}),
                            html.Div([
                                html.Span(style={'display': 'inline-block', 'width': '14px', 'height': '14px', 'backgroundColor': '#48cae4', 'marginRight': '6px', 'borderRadius': '50%', 'verticalAlign': 'middle'}),
                                html.Span("Out of diffusion range", style={'fontSize': '0.85rem'}),
                            ]),
                        ]),
                    ]),
                ],
                className="mb-4",
                style={'display': 'none'},
            ),

            # Shared expression threshold slider (visible in NP and Hormone modes)
            html.Div(
                id='expression-threshold-container',
                children=[
                    html.Label("Expression Threshold", className="control-label"),
                    dcc.Slider(
                        id='expression-threshold',
                        min=0.5,
                        max=5,
                        step=0.5,
                        value=3,
                        marks={0.5: '0.5', 1: '1', 2: '2', 3: '3', 4: '4', 5: '5'},
                        tooltip={"placement": "bottom", "always_visible": False},
                    ),
                ],
                className="mb-4",
                style={'display': 'none'},
            ),

            # Display options
            html.Hr(style={'borderColor': 'rgba(0,0,0,0.08)'}),
            html.Div([
                html.Label("Display Options", className="control-label"),
                dcc.Checklist(
                    id='display-options',
                    options=[
                        {'label': ' Show region boundaries', 'value': 'show_boundaries'},
                        {'label': ' Show region labels', 'value': 'show_labels'},
                    ],
                    value=['show_boundaries', 'show_labels'],
                    labelStyle={'display': 'block', 'marginBottom': '6px', 'fontSize': '0.85rem'},
                    inputStyle={'marginRight': '10px'},
                ),
            ], className="mb-4"),

            # Region highlight selector (click to select)
            *([
                html.Hr(style={'borderColor': 'rgba(0,0,0,0.08)'}),
                html.Div([
                    html.Label("Highlight Region (click to select)", className="control-label"),
                    html.Div(
                        id='region-list',
                        children=[
                            html.Button(
                                region,
                                id={'type': 'region-btn', 'index': region},
                                n_clicks=0,
                                title=f"[{region_descriptions.get(region, {}).get('division', '')}] {region_descriptions.get(region, {}).get('full_name', region)}: {region_descriptions.get(region, {}).get('description', 'No description available')}",
                                className='region-chip',
                                style={
                                    'display': 'inline-block',
                                    'padding': '2px 6px',
                                    'margin': '2px',
                                    'fontSize': '0.75rem',
                                    'borderRadius': '4px',
                                    'backgroundColor': '#f0f0f0',
                                    'border': '1px solid transparent',
                                    'cursor': 'pointer',
                                    'color': '#333',
                                },
                            )
                            for region in region_list
                        ],
                        style={'maxHeight': '120px', 'overflowY': 'auto', 'lineHeight': '1.8'},
                    ),
                    # Store for currently selected region
                    dcc.Store(id='selected-region', data=None),
                ], className="mb-4"),
            ] if enable_region_highlight else []),

            # Point size slider
            html.Div([
                html.Label("Point Size", className="control-label"),
                dcc.Slider(
                    id='point-size',
                    min=1,
                    max=10,
                    step=1,
                    value=5,
                    marks={1: '1', 5: '5', 10: '10'},
                ),
            ], className="mb-4"),

            # Subsample slider
            html.Div([
                html.Label("Subsample %", className="control-label"),
                dcc.Slider(
                    id='subsample-pct',
                    min=5,
                    max=100,
                    step=5,
                    value=default_subsample,
                    marks={5: '5%', 25: '25%', 50: '50%', 100: '100%'},
                ),
            ], className="mb-3"),
        ]),
    ], className="h-100 modern-card")


def create_center_panel(default_slices=None):
    """Create the center visualization panel with coronal slice grid."""
    default_slices = default_slices or []
    z_min = min(default_slices) if default_slices else 0
    z_max = max(default_slices) if default_slices else 10
    # Build marks at actual slice positions, showing label for first, last, and every ~4th
    z_marks = {}
    for i, z in enumerate(default_slices):
        if i == 0 or i == len(default_slices) - 1 or i % 4 == 0:
            z_marks[z] = f'{z:.1f}'
        else:
            z_marks[z] = ''

    return dbc.Card([
        dbc.CardHeader([
            html.Div([
                # Columns slider
                html.Div([
                    html.Label("Columns", style={'fontSize': '0.8rem', 'fontWeight': '600', 'marginRight': '8px', 'whiteSpace': 'nowrap'}),
                    html.Div(
                        dcc.Slider(
                            id='grid-columns',
                            min=1,
                            max=3,
                            step=1,
                            value=2,
                            marks={1: '1', 2: '2', 3: '3'},
                        ),
                        style={'flex': '1', 'minWidth': '80px'},
                    ),
                ], style={'display': 'flex', 'alignItems': 'center', 'flex': '1', 'marginRight': '24px'}),
                # Z-range slider (snaps to actual slice positions)
                html.Div([
                    html.Span("P", style={'fontSize': '0.75rem', 'fontWeight': '600', 'color': '#888', 'marginRight': '6px'}),
                    html.Div(
                        dcc.RangeSlider(
                            id='z-range',
                            min=z_min,
                            max=z_max,
                            step=None,  # snap to marks only
                            value=[z_min, z_max],
                            marks=z_marks,
                            tooltip={"placement": "bottom", "always_visible": True},
                        ),
                        style={'flex': '1', 'minWidth': '200px'},
                    ),
                    html.Span("A", style={'fontSize': '0.75rem', 'fontWeight': '600', 'color': '#888', 'marginLeft': '6px'}),
                ], style={'display': 'flex', 'alignItems': 'center', 'flex': '3'}),
            ], style={'display': 'flex', 'alignItems': 'center', 'width': '100%'}),
        ]),
        dbc.CardBody([
            dcc.Loading(
                id="loading-slices",
                type="circle",
                color="#667eea",
                children=[
                    html.Div(
                        id='slice-grid-container',
                        children=[
                            dcc.Graph(
                                id='slice-grid',
                                config={
                                    'displayModeBar': True,
                                    'modeBarButtonsToRemove': ['lasso2d', 'select2d'],
                                    'displaylogo': False,
                                    'plotGlPixelRatio': 1,  # Avoid WebGL canvas size limit (16384px) on Retina
                                },
                            ),
                        ],
                    ),
                ],
            ),
        ], className="p-2"),
    ], className="h-100 modern-card")


def create_right_panel(enable_quantile_toggle=False):
    """Create the right details panel."""
    return dbc.Card([
        dbc.CardHeader(html.H5("Cluster Details", className="mb-0", style={'fontWeight': '600'})),
        dbc.CardBody([
            # Expression mode toggle
            *([
                html.Div([
                    html.Label("Expression metric:", style={'fontSize': '0.85rem', 'marginRight': '8px'}),
                    dbc.Switch(
                        id='expr-mode-toggle',
                        label='',
                        value=False,
                        style={'display': 'inline-block'},
                    ),
                    html.Span(
                        id='expr-mode-label',
                        children="log₂(CPM)",
                        style={'fontSize': '0.85rem', 'marginLeft': '4px'},
                    ),
                ], style={'marginBottom': '12px', 'display': 'flex', 'alignItems': 'center'}),
                html.Hr(style={'marginTop': '0'}),
            ] if enable_quantile_toggle else []),
            # Dynamic cell details content
            html.Div(
                id='cell-details',
                children=[
                    html.P("Click on a cell to see details", className="text-muted",
                           style={'fontStyle': 'italic'}),
                ],
            ),
        ]),
    ], className="h-100 modern-card", style={'overflowY': 'auto'})


def create_layout(cell_type_levels, nt_types, np_system_names, hormone_names=None, region_list=None, region_descriptions=None, enable_region_highlight=False, enable_quantile_toggle=False, dataset_names=None, default_dataset='mouse_abc', default_slices=None, default_subsample=30):
    """Create the main application layout."""
    return dbc.Container([
        # Header
        dbc.Row([
            dbc.Col([
                html.H1("Mouse ABC HY", className="app-title mb-1"),
                html.P("Hypothalamus Single-Cell Spatial Atlas", className="app-subtitle mb-0"),
            ], width=12),
        ], className="py-4 mb-4"),

        # Main content: three-column layout
        dbc.Row([
            # Left panel - Controls
            dbc.Col(
                create_left_panel(cell_type_levels, nt_types, np_system_names, hormone_names, region_list, region_descriptions, enable_region_highlight, dataset_names, default_dataset, default_subsample),
                width=2,
                className="pe-3",
            ),

            # Center panel - Slice grid
            dbc.Col(
                create_center_panel(default_slices=default_slices),
                width=7,
                className="px-2",
            ),

            # Right panel - Details
            dbc.Col(
                create_right_panel(enable_quantile_toggle),
                width=3,
                className="ps-3",
            ),
        ], className="g-0", style={'height': '90vh'}),

    ], fluid=True, className="px-4 py-3")
