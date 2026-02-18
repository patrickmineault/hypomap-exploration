"""Diversity computation utilities: SAT-based 3D box finder."""

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree


def compute_diversity(labels, metric):
    """Compute a diversity metric from an array of integer labels.

    Args:
        labels: 1D array of integer category labels
        metric: one of 'Shannon', 'Simpson', 'Type count'

    Returns:
        float diversity value, or NaN if fewer than 2 labels
    """
    if len(labels) < 2:
        return float('nan')

    counts = np.bincount(labels)
    counts = counts[counts > 0]
    n = counts.sum()

    if metric == "Type count":
        return float(len(counts))

    p = counts / n
    if metric == "Shannon":
        return float(-np.sum(p * np.log2(p)))
    elif metric == "Simpson":
        return float(1.0 - np.sum(p ** 2))
    return float('nan')


def build_local_diversity_grid(cells_df, granularity, metric, radius=0.15,
                                grid_res=0.1):
    """Compute local diversity on a regular grid for each z-slice.

    For each grid point, diversity is computed from all cells within `radius`,
    using the same KDTree radius-search as the multi-slice overview.

    Returns:
        div_grids: dict z_slice -> 2D array (nx, ny) of diversity values (NaN where invalid)
        count_grids: dict z_slice -> 2D array (nx, ny) of 1/0 (valid/invalid)
        grid_geo: dict with x_min, y_min, grid_res, nx, ny
    """
    codes, _uniques = pd.factorize(cells_df[granularity])
    int_labels = codes.astype(np.int32)

    x_vals = cells_df['x'].values
    y_vals = cells_df['y'].values
    z_vals = cells_df['z_slice'].values

    x_min, x_max = x_vals.min(), x_vals.max()
    y_min, y_max = y_vals.min(), y_vals.max()

    nx = int(np.ceil((x_max - x_min) / grid_res)) + 1
    ny = int(np.ceil((y_max - y_min) / grid_res)) + 1

    # Build grid point coordinates
    gx = x_min + np.arange(nx) * grid_res
    gy = y_min + np.arange(ny) * grid_res
    gxx, gyy = np.meshgrid(gx, gy, indexing='ij')
    grid_points = np.column_stack([gxx.ravel(), gyy.ravel()])

    slices = np.unique(z_vals)
    div_grids = {}
    count_grids = {}

    for z in slices:
        mask = z_vals == z
        xy = np.column_stack([x_vals[mask], y_vals[mask]])
        labels = int_labels[mask]

        if len(xy) < 2:
            continue

        tree = cKDTree(xy)
        neighbors = tree.query_ball_point(grid_points, radius)

        div_flat = np.full(len(grid_points), float('nan'))
        for i, nbrs in enumerate(neighbors):
            if len(nbrs) >= 2:
                div_flat[i] = compute_diversity(labels[nbrs], metric)

        div_grid = div_flat.reshape(nx, ny)
        valid_grid = (~np.isnan(div_grid)).astype(np.float64)

        div_grids[z] = div_grid
        count_grids[z] = valid_grid

    grid_geo = {
        'x_min': x_min, 'y_min': y_min,
        'grid_res': grid_res, 'nx': nx, 'ny': ny,
    }
    return div_grids, count_grids, grid_geo


def mask_excluded_regions(div_grids, count_grids, grid_geo, cells_df,
                          exclude_regions):
    """Zero out grid points whose nearest cell belongs to an excluded region.

    Modifies div_grids and count_grids in place. Grid points in excluded
    regions get NaN diversity and 0 count, so they contribute nothing to
    downstream SAT queries.

    Args:
        div_grids: dict z_slice -> (nx, ny) diversity array
        count_grids: dict z_slice -> (nx, ny) valid-point array
        grid_geo: dict with x_min, y_min, grid_res, nx, ny
        cells_df: DataFrame with 'x', 'y', 'z_slice', 'region' columns
        exclude_regions: list of region names to exclude
    """
    if not exclude_regions:
        return

    exclude_set = set(exclude_regions)
    nx, ny = grid_geo['nx'], grid_geo['ny']
    gx = grid_geo['x_min'] + np.arange(nx) * grid_geo['grid_res']
    gy = grid_geo['y_min'] + np.arange(ny) * grid_geo['grid_res']
    gxx, gyy = np.meshgrid(gx, gy, indexing='ij')
    grid_points = np.column_stack([gxx.ravel(), gyy.ravel()])

    x_vals = cells_df['x'].values
    y_vals = cells_df['y'].values
    z_vals = cells_df['z_slice'].values
    regions = cells_df['region'].values

    for z in list(div_grids.keys()):
        mask_z = z_vals == z
        xy_z = np.column_stack([x_vals[mask_z], y_vals[mask_z]])
        regions_z = regions[mask_z]

        if len(xy_z) == 0:
            continue

        tree = cKDTree(xy_z)
        _, nearest_idx = tree.query(grid_points)
        nearest_regions = regions_z[nearest_idx]

        bad_mask = np.array([r in exclude_set for r in nearest_regions])
        bad_mask = bad_mask.reshape(nx, ny)

        div_grids[z][bad_mask] = float('nan')
        count_grids[z][bad_mask] = 0.0


def build_diversity_sats(div_grids, count_grids, nx, ny):
    """Build summed area tables from pre-computed local diversity grids.

    For each z-slice, builds two SATs:
    - sum SAT: cumulative sum of diversity values (NaN treated as 0)
    - count SAT: cumulative sum of valid-point indicators

    Returns:
        sum_sats: dict z_slice -> (nx+1, ny+1) SAT of diversity sums
        count_sats: dict z_slice -> (nx+1, ny+1) SAT of valid counts
    """
    sum_sats = {}
    count_sats = {}

    for z in div_grids:
        div_grid = div_grids[z].copy()
        div_grid[np.isnan(div_grid)] = 0.0

        sum_sat = np.zeros((nx + 1, ny + 1), dtype=np.float64)
        sum_sat[1:, 1:] = np.cumsum(np.cumsum(div_grid, axis=0), axis=1)
        sum_sats[z] = sum_sat

        count_sat = np.zeros((nx + 1, ny + 1), dtype=np.float64)
        count_sat[1:, 1:] = np.cumsum(np.cumsum(count_grids[z], axis=0), axis=1)
        count_sats[z] = count_sat

    return sum_sats, count_sats


def find_max_mean_diversity_box(sum_sats, count_sats, grid_geo, z_slices,
                                 n_slices, volume_um3, symmetry,
                                 midline_x=0.0):
    """Find the 3D box with maximum mean local diversity.

    Instead of computing diversity of the box contents, this finds the box
    where the average pre-computed local diversity is highest.

    Args:
        sum_sats: dict z_slice -> SAT of diversity value sums
        count_sats: dict z_slice -> SAT of valid-point counts
        grid_geo: dict with grid geometry
        z_slices: sorted array of z-slice values
        n_slices: number of consecutive z-slices to span
        volume_um3: target volume in millions of cubic micrometers
        symmetry: 'One-sided' or 'Symmetric'
        midline_x: x coordinate of midline for symmetric mode

    Returns:
        dict with 'best', 'best_z_start_idx', 'best_aspect_bx', 'best_aspect_by'
    """
    grid_res = grid_geo['grid_res']
    x_min = grid_geo['x_min']
    y_min = grid_geo['y_min']
    nx = grid_geo['nx']
    ny = grid_geo['ny']
    slice_thickness = 0.2  # mm

    volume_mm3 = volume_um3 * 0.001  # 1 M um3 = 1e6 um3 = 0.001 mm3

    n_slices = min(n_slices, len(z_slices))

    best_div = -1.0
    best_info = None
    best_z_idx = 0
    best_bx = 0
    best_by = 0

    for zi in range(len(z_slices) - n_slices + 1):
        z_start = z_slices[zi]
        z_end = z_slices[zi + n_slices - 1]
        z_extent = z_end - z_start + slice_thickness

        # Sum SATs across selected slices
        window_sum = np.zeros((nx + 1, ny + 1), dtype=np.float64)
        window_count = np.zeros((nx + 1, ny + 1), dtype=np.float64)
        for zj in range(zi, zi + n_slices):
            z_val = z_slices[zj]
            if z_val in sum_sats:
                window_sum += sum_sats[z_val]
                window_count += count_sats[z_val]

        xy_area_mm2 = volume_mm3 / z_extent
        xy_area_vox = xy_area_mm2 / (grid_res ** 2)

        # Log-spaced aspect ratios
        ratios = np.logspace(np.log10(0.22), np.log10(4.5), 18)
        seen_dims = set()

        for ratio in ratios:
            bx = max(1, min(int(np.round(np.sqrt(xy_area_vox * ratio))), nx))
            by = max(1, min(int(np.round(np.sqrt(xy_area_vox / ratio))), ny))
            dim_key = (bx, by)
            if dim_key in seen_dims:
                continue
            seen_dims.add(dim_key)

            # SAT query for sum and count
            div_sum = (window_sum[bx:, by:]
                       - window_sum[:-bx, by:]
                       - window_sum[bx:, :-by]
                       + window_sum[:-bx, :-by])

            n_valid = (window_count[bx:, by:]
                       - window_count[:-bx, by:]
                       - window_count[bx:, :-by]
                       + window_count[:-bx, :-by])

            if symmetry == "Symmetric":
                mid_ix = int(np.round((midline_x - x_min) / grid_res))
                x_start = mid_ix - bx // 2
                x_start = max(0, min(x_start, div_sum.shape[0] - 1))
                div_sum = div_sum[x_start:x_start + 1, :]
                n_valid = n_valid[x_start:x_start + 1, :]

            # Mean diversity where we have enough valid points
            valid = n_valid >= 3
            if not valid.any():
                continue

            mean_div = np.full_like(div_sum, -1.0)
            np.divide(div_sum, n_valid, out=mean_div, where=valid)
            mean_div[~valid] = -1.0

            best_idx = np.unravel_index(np.argmax(mean_div), mean_div.shape)
            local_best = float(mean_div[best_idx])

            if local_best > best_div:
                best_div = local_best
                if symmetry == "Symmetric":
                    box_ix = x_start
                else:
                    box_ix = best_idx[0]
                box_iy = best_idx[1]

                x0 = x_min + box_ix * grid_res
                y0 = y_min + box_iy * grid_res
                x1 = x0 + bx * grid_res
                y1 = y0 + by * grid_res

                best_info = {
                    'z_start': float(z_start),
                    'z_end': float(z_end),
                    'x_range': (float(x0), float(x1)),
                    'y_range': (float(y0), float(y1)),
                    'diversity': float(local_best),
                    'n_valid_points': int(n_valid[best_idx]),
                    'bx_mm': float(bx * grid_res),
                    'by_mm': float(by * grid_res),
                    'z_extent_mm': float(z_extent),
                }
                best_z_idx = zi
                best_bx = bx
                best_by = by

    return {
        'best': best_info,
        'best_z_start_idx': best_z_idx,
        'best_aspect_bx': best_bx,
        'best_aspect_by': best_by,
    }
