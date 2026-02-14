import marimo

__generated_with = "0.19.4"
app = marimo.App()


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Region Boundary Algorithm Comparison

    This notebook compares different algorithms for computing region boundaries in coronal slices:

    1. **Convex Hull** (current algorithm): Uses `scipy.spatial.ConvexHull` - fast but can include areas outside the actual region
    2. **Concave Hull** (geopandas): Uses `geopandas.concave_hull` - ratio-based control (0=convex, 1=tight)
    3. **Alpha Shape** (alphashape): Uses `alphashape.alphashape` - alpha parameter controls tightness (0=convex, higher=tighter)

    Select a slice to compare the algorithms.
    """)
    return


@app.cell
def _():
    import pandas as pd
    import numpy as np
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    from scipy.spatial import ConvexHull
    import geopandas as gpd
    from shapely.geometry import MultiPoint
    import alphashape
    from pathlib import Path

    DATA_DIR = Path('../data')
    CELLS_PATH = DATA_DIR / 'processed' / 'mouse_abc' / 'cells_with_coords.parquet'
    MIDLINE_X = 5.5
    return (
        CELLS_PATH,
        ConvexHull,
        MIDLINE_X,
        MultiPoint,
        alphashape,
        go,
        gpd,
        make_subplots,
        pd,
    )


@app.cell
def _(CELLS_PATH, pd):
    # Load cell data
    cells_df = pd.read_parquet(CELLS_PATH)
    cells_df = cells_df[cells_df['region'] != 'HY-unassigned'].copy()
    cells_df['z_slice'] = cells_df['z'].round(1)

    available_slices = sorted(cells_df['z_slice'].unique())
    print(f"Loaded {len(cells_df):,} cells")
    print(f"Available slices: {available_slices}")
    return available_slices, cells_df


@app.cell
def _(MIDLINE_X, cells_df, slice_selector):
    # Get slice data and lateralize regions
    _z = float(slice_selector.value)
    slice_df = cells_df[cells_df['z_slice'] == _z].copy()

    # Lateralize regions
    slice_df['region_display'] = slice_df['region']
    for _region in slice_df['region'].unique():
        _region_mask = slice_df['region'] == _region
        _region_cells = slice_df.loc[_region_mask]

        if len(_region_cells) < 5:
            continue

        _left_count = (_region_cells['x'] < MIDLINE_X).sum()
        _right_count = (_region_cells['x'] >= MIDLINE_X).sum()
        _total = len(_region_cells)

        if _left_count > 0.01 * _total and _right_count > 0.01 * _total:
            _left_mask = _region_mask & (slice_df['x'] < MIDLINE_X)
            _right_mask = _region_mask & (slice_df['x'] >= MIDLINE_X)
            slice_df.loc[_left_mask, 'region_display'] = f"{_region}-L"
            slice_df.loc[_right_mask, 'region_display'] = f"{_region}-R"

    regions_in_slice = sorted(slice_df['region_display'].unique())
    print(f"Slice z={_z}: {len(slice_df):,} cells, {len(regions_in_slice)} regions")
    return regions_in_slice, slice_df


@app.cell
def _(
    ConvexHull,
    MultiPoint,
    alpha_slider,
    alphashape,
    gpd,
    ratio_slider,
    regions_in_slice,
    slice_df,
):
    # Compute boundaries using all three algorithms
    convex_boundaries = {}
    concave_boundaries = {}
    alpha_boundaries = {}

    _ratio = ratio_slider.value
    _alpha = alpha_slider.value

    for _region in regions_in_slice:
        _region_df = slice_df[slice_df['region_display'] == _region]
        if len(_region_df) < 3:
            continue

        _points = _region_df[['x', 'y']].values

        # 1. Convex hull (current algorithm)
        try:
            _hull = ConvexHull(_points)
            _hull_points = _points[_hull.vertices].tolist()
            _hull_points.append(_hull_points[0])  # Close polygon
            convex_boundaries[_region] = _hull_points
        except Exception:
            pass

        # 2. Concave hull (geopandas)
        try:
            _multipoint = MultiPoint(_points)
            _gdf = gpd.GeoDataFrame(geometry=[_multipoint])
            _concave = _gdf.concave_hull(ratio=_ratio).iloc[0]

            if _concave.geom_type == 'Polygon':
                _coords = list(_concave.exterior.coords)
                concave_boundaries[_region] = [[c[0], c[1]] for c in _coords]
            elif _concave.geom_type == 'MultiPolygon':
                _largest = max(_concave.geoms, key=lambda g: g.area)
                _coords = list(_largest.exterior.coords)
                concave_boundaries[_region] = [[c[0], c[1]] for c in _coords]
        except Exception:
            pass

        # 3. Alpha shape (alphashape package)
        try:
            _alpha_shape = alphashape.alphashape(_points, _alpha)

            if _alpha_shape.geom_type == 'Polygon':
                _coords = list(_alpha_shape.exterior.coords)
                alpha_boundaries[_region] = [[c[0], c[1]] for c in _coords]
            elif _alpha_shape.geom_type == 'MultiPolygon':
                _largest = max(_alpha_shape.geoms, key=lambda g: g.area)
                _coords = list(_largest.exterior.coords)
                alpha_boundaries[_region] = [[c[0], c[1]] for c in _coords]
        except Exception:
            pass

    print(f"Computed: {len(convex_boundaries)} convex, {len(concave_boundaries)} concave (geopandas), {len(alpha_boundaries)} alpha shapes")

    # Compute area differences
    def _polygon_area(coords):
        _n = len(coords) - 1
        _area = 0.0
        for _i in range(_n):
            _j = (_i + 1) % _n
            _area += coords[_i][0] * coords[_j][1]
            _area -= coords[_j][0] * coords[_i][1]
        return abs(_area) / 2.0

    area_comparison = []
    for _region in convex_boundaries:
        _convex_area = _polygon_area(convex_boundaries[_region])
        _concave_area = _polygon_area(concave_boundaries[_region]) if _region in concave_boundaries else _convex_area
        _alpha_area = _polygon_area(alpha_boundaries[_region]) if _region in alpha_boundaries else _convex_area

        area_comparison.append({
            'region': _region,
            'convex_area': _convex_area,
            'concave_area': _concave_area,
            'alpha_area': _alpha_area,
            'concave_reduction': 100 * (1 - _concave_area / _convex_area) if _convex_area > 0 else 0,
            'alpha_reduction': 100 * (1 - _alpha_area / _convex_area) if _convex_area > 0 else 0,
        })
    return (
        alpha_boundaries,
        area_comparison,
        concave_boundaries,
        convex_boundaries,
    )


@app.cell
def _(mo):
    # Concave hull ratio parameter (geopandas)
    ratio_slider = mo.ui.slider(
        start=0.0,
        stop=1.0,
        step=0.05,
        value=0.3,
        label='Concave hull ratio (geopandas)'
    )

    # Alpha parameter (alphashape)
    alpha_slider = mo.ui.slider(
        start=0.0,
        stop=10.0,
        step=0.5,
        value=2.0,
        label='Alpha (alphashape)'
    )

    mo.md(f"""
    **Algorithm Parameters:**

    {mo.hstack([ratio_slider, alpha_slider])}

    - **Concave hull ratio** (geopandas): 0.0=convex, 1.0=tightest
    - **Alpha** (alphashape): 0=convex, higher=tighter (try 1-5)
    """)
    return alpha_slider, ratio_slider


@app.cell(hide_code=True)
def _(available_slices, mo):
    # Slice selector
    slice_selector = mo.ui.dropdown(
        options=[str(s) for s in available_slices],
        value='8.0',
        label='Select z slice'
    )
    mo.md(f"**Select slice:** {slice_selector}")
    return (slice_selector,)


@app.cell(hide_code=True)
def _(
    alpha_boundaries,
    concave_boundaries,
    convex_boundaries,
    go,
    make_subplots,
    slice_df,
    slice_selector,
):
    # Create 3-panel comparison plot
    fig_compare = make_subplots(
        rows=1, cols=3,
        subplot_titles=['Convex Hull (current)', 'Concave Hull (geopandas)', 'Alpha Shape (alphashape)'],
        horizontal_spacing=0.05
    )

    # Get unique regions and assign colors
    _regions = sorted(
        set(convex_boundaries.keys()) |
        set(concave_boundaries.keys()) |
        set(alpha_boundaries.keys())
    )
    _colors = {}
    _color_palette = [
        '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
        '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf',
        '#aec7e8', '#ffbb78', '#98df8a', '#ff9896', '#c5b0d5',
        '#c49c94', '#f7b6d2', '#c7c7c7', '#dbdb8d', '#9edae5'
    ]
    for _i, _r in enumerate(_regions):
        _colors[_r] = _color_palette[_i % len(_color_palette)]

    # Plot cells (same for all three)
    for _col in [1, 2, 3]:
        for _region in _regions:
            _region_df = slice_df[slice_df['region_display'] == _region]
            if len(_region_df) == 0:
                continue

            fig_compare.add_trace(
                go.Scattergl(
                    x=_region_df['x'],
                    y=_region_df['y'],
                    mode='markers',
                    marker=dict(size=2, color=_colors[_region], opacity=0.5),
                    name=_region,
                    showlegend=(_col == 1),
                    legendgroup=_region,
                ),
                row=1, col=_col
            )

    # Plot convex boundaries (col 1)
    for _region, _coords in convex_boundaries.items():
        _x = [c[0] for c in _coords]
        _y = [c[1] for c in _coords]
        fig_compare.add_trace(
            go.Scatter(
                x=_x, y=_y,
                mode='lines',
                line=dict(color=_colors[_region], width=2),
                showlegend=False,
                legendgroup=_region,
            ),
            row=1, col=1
        )

    # Plot concave boundaries (col 2)
    for _region, _coords in concave_boundaries.items():
        _x = [c[0] for c in _coords]
        _y = [c[1] for c in _coords]
        fig_compare.add_trace(
            go.Scatter(
                x=_x, y=_y,
                mode='lines',
                line=dict(color=_colors[_region], width=2),
                showlegend=False,
                legendgroup=_region,
            ),
            row=1, col=2
        )

    # Plot alpha shape boundaries (col 3)
    for _region, _coords in alpha_boundaries.items():
        _x = [c[0] for c in _coords]
        _y = [c[1] for c in _coords]
        fig_compare.add_trace(
            go.Scatter(
                x=_x, y=_y,
                mode='lines',
                line=dict(color=_colors[_region], width=2),
                showlegend=False,
                legendgroup=_region,
            ),
            row=1, col=3
        )

    _z = slice_selector.value
    fig_compare.update_layout(
        title=f'Region Boundary Comparison (z={_z})',
        height=600,
        width=1400,
        legend=dict(
            yanchor="top",
            y=0.99,
            xanchor="left",
            x=1.02,
            font=dict(size=8)
        )
    )

    # Make axes equal and matching
    fig_compare.update_xaxes(title_text='x', scaleanchor='y', scaleratio=1)
    fig_compare.update_yaxes(title_text='y')

    fig_compare.show()
    return


@app.cell
def _(area_comparison, mo, pd):
    # Show area comparison table
    if area_comparison:
        _df = pd.DataFrame(area_comparison)
        _df = _df.sort_values('alpha_reduction', ascending=False)

        _avg_concave = _df['concave_reduction'].mean()
        _avg_alpha = _df['alpha_reduction'].mean()

        _table_rows = '\n'.join([
            f"| {r['region']} | {r['convex_area']:.3f} | {r['concave_reduction']:.1f}% | {r['alpha_reduction']:.1f}% |"
            for _, r in _df.head(15).iterrows()
        ])

        _content = f"""
    ## Area Comparison

    **Average area reduction vs convex hull:**
    - Concave (geopandas): {_avg_concave:.1f}%
    - Alpha shape: {_avg_alpha:.1f}%

    Top 15 regions by alpha shape reduction:

    | Region | Convex Area | Concave Reduction | Alpha Reduction |
    |--------|-------------|-------------------|-----------------|
    {_table_rows}
    """
    else:
        _content = "No area comparison data available"

    mo.md(_content)
    return


@app.cell
def _(mo):
    # Region selector for detailed view
    mo.md("""
    ---
    ## Detailed Region View

    Select a specific region to see a zoomed comparison.
    """)
    return


@app.cell
def _(mo, regions_in_slice):
    region_selector = mo.ui.dropdown(
        options=regions_in_slice,
        value=regions_in_slice[0] if regions_in_slice else None,
        label='Select region'
    )
    mo.md(f"**Region:** {region_selector}")
    return (region_selector,)


@app.cell
def _(
    alpha_boundaries,
    concave_boundaries,
    convex_boundaries,
    go,
    make_subplots,
    region_selector,
    slice_df,
):
    # Detailed view of single region
    _region = region_selector.value
    _region_df = slice_df[slice_df['region_display'] == _region]

    fig_detail = make_subplots(
        rows=1, cols=3,
        subplot_titles=['Convex Hull', 'Concave Hull (geopandas)', 'Alpha Shape'],
        horizontal_spacing=0.06
    )

    # Plot cells
    for _col in [1, 2, 3]:
        fig_detail.add_trace(
            go.Scattergl(
                x=_region_df['x'],
                y=_region_df['y'],
                mode='markers',
                marker=dict(size=4, color='steelblue', opacity=0.6),
                name='Cells',
                showlegend=(_col == 1),
            ),
            row=1, col=_col
        )

    # Convex boundary
    if _region in convex_boundaries:
        _coords = convex_boundaries[_region]
        fig_detail.add_trace(
            go.Scatter(
                x=[c[0] for c in _coords],
                y=[c[1] for c in _coords],
                mode='lines',
                line=dict(color='red', width=3),
                fill='toself',
                fillcolor='rgba(255, 0, 0, 0.1)',
                name='Convex Hull',
                showlegend=True,
            ),
            row=1, col=1
        )

    # Concave boundary
    if _region in concave_boundaries:
        _coords = concave_boundaries[_region]
        fig_detail.add_trace(
            go.Scatter(
                x=[c[0] for c in _coords],
                y=[c[1] for c in _coords],
                mode='lines',
                line=dict(color='green', width=3),
                fill='toself',
                fillcolor='rgba(0, 255, 0, 0.1)',
                name='Concave Hull',
                showlegend=True,
            ),
            row=1, col=2
        )

    # Alpha shape boundary
    if _region in alpha_boundaries:
        _coords = alpha_boundaries[_region]
        fig_detail.add_trace(
            go.Scatter(
                x=[c[0] for c in _coords],
                y=[c[1] for c in _coords],
                mode='lines',
                line=dict(color='purple', width=3),
                fill='toself',
                fillcolor='rgba(128, 0, 128, 0.1)',
                name='Alpha Shape',
                showlegend=True,
            ),
            row=1, col=3
        )

    fig_detail.update_layout(
        title=f'Region: {_region} ({len(_region_df)} cells)',
        height=450,
        width=1200,
    )
    fig_detail.update_xaxes(scaleanchor='y', scaleratio=1)
    fig_detail.show()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ---
    ## Notes

    **Convex Hull (current algorithm)**:
    - Fast and deterministic
    - Always produces a valid polygon
    - Can include significant empty space for irregular shapes (e.g., crescent-shaped regions)

    **Concave Hull (geopandas)**:
    - Uses `geopandas.concave_hull()` with ratio parameter
    - Ratio 0.0 = convex hull, 1.0 = tightest fit
    - Good balance around 0.2-0.4
    - Fast, built on GEOS library

    **Alpha Shape (alphashape)**:
    - Uses the `alphashape` package
    - Alpha 0 = convex hull, higher = tighter
    - More control but can create holes/disconnected polygons at high alpha
    - Try values 1-5 for most regions

    **Recommendation**: Start with geopandas concave_hull (ratio=0.3) for a good balance of fit and stability.
    """)
    return


@app.cell
def _():
    import marimo as mo
    return (mo,)


if __name__ == "__main__":
    app.run()
