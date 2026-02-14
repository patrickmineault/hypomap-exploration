"""Tests for the local-diversity 3D box finder (hypomap/diversity.py)."""

import numpy as np
import pandas as pd
import pytest

from hypomap.diversity import (
    build_diversity_sats,
    build_local_diversity_grid,
    compute_diversity,
    find_max_mean_diversity_box,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_cells(x, y, z_slice, labels):
    """Build a minimal cells DataFrame."""
    return pd.DataFrame(
        {
            "x": x,
            "y": y,
            "z_slice": z_slice,
            "type": labels,
        }
    )


# ---------------------------------------------------------------------------
# Tests: compute_diversity
# ---------------------------------------------------------------------------


class TestComputeDiversity:
    def test_shannon_uniform(self):
        """Uniform distribution of 4 types -> Shannon = log2(4) = 2.0."""
        labels = np.array([0, 1, 2, 3] * 25)
        assert compute_diversity(labels, "Shannon") == pytest.approx(2.0)

    def test_simpson_uniform(self):
        """Uniform distribution of 4 types -> Simpson = 1 - 4*(1/4)^2 = 0.75."""
        labels = np.array([0, 1, 2, 3] * 25)
        assert compute_diversity(labels, "Simpson") == pytest.approx(0.75)

    def test_type_count(self):
        labels = np.array([0, 0, 1, 1, 2])
        assert compute_diversity(labels, "Type count") == 3.0

    def test_single_type_shannon(self):
        """All same type -> Shannon = 0."""
        labels = np.array([0, 0, 0, 0, 0])
        assert compute_diversity(labels, "Shannon") == pytest.approx(0.0)

    def test_too_few_labels(self):
        assert np.isnan(compute_diversity(np.array([0]), "Shannon"))
        assert np.isnan(compute_diversity(np.array([]), "Shannon"))


# ---------------------------------------------------------------------------
# Tests: build_local_diversity_grid
# ---------------------------------------------------------------------------


class TestBuildLocalDiversityGrid:
    def test_single_slice(self):
        """Should produce a grid of diversity values for one slice."""
        rng = np.random.default_rng(42)
        n = 200
        types = [f"T{i % 4}" for i in range(n)]
        df = _make_cells(
            x=rng.uniform(0, 2, n).tolist(),
            y=rng.uniform(0, 2, n).tolist(),
            z_slice=[1.0] * n,
            labels=types,
        )
        div_grids, count_grids, geo = build_local_diversity_grid(
            df, "type", "Shannon", radius=0.3, grid_res=0.1
        )
        assert 1.0 in div_grids
        assert div_grids[1.0].shape == (geo["nx"], geo["ny"])
        # Should have some valid points
        assert count_grids[1.0].sum() > 0

    def test_two_slices(self):
        """Different z-slices produce separate grids."""
        rng = np.random.default_rng(42)
        n = 100
        types = [f"T{i % 3}" for i in range(n)]
        df = _make_cells(
            x=rng.uniform(0, 2, 2 * n).tolist(),
            y=rng.uniform(0, 2, 2 * n).tolist(),
            z_slice=[1.0] * n + [1.2] * n,
            labels=types * 2,
        )
        div_grids, count_grids, geo = build_local_diversity_grid(
            df, "type", "Shannon", radius=0.3, grid_res=0.1
        )
        assert len(div_grids) == 2
        assert 1.0 in div_grids
        assert 1.2 in div_grids

    def test_diverse_region_has_higher_values(self):
        """A region with mixed types should have higher diversity values."""
        rng = np.random.default_rng(42)

        # Left half: uniform (all type A)
        n_uni = 200
        x_uni = rng.uniform(0.0, 1.8, n_uni)
        y_uni = rng.uniform(0.0, 4.0, n_uni)
        labels_uni = ["A"] * n_uni

        # Right half: diverse (8 types)
        n_div = 200
        x_div = rng.uniform(2.2, 4.0, n_div)
        y_div = rng.uniform(0.0, 4.0, n_div)
        types = ["A", "B", "C", "D", "E", "F", "G", "H"]
        labels_div = [types[i % 8] for i in range(n_div)]

        df = _make_cells(
            x=list(x_uni) + list(x_div),
            y=list(y_uni) + list(y_div),
            z_slice=[5.0] * (n_uni + n_div),
            labels=labels_uni + labels_div,
        )
        div_grids, count_grids, geo = build_local_diversity_grid(
            df, "type", "Shannon", radius=0.3, grid_res=0.1
        )
        grid = div_grids[5.0]
        # Split grid at midpoint
        mid_ix = geo["nx"] // 2
        left_vals = grid[:mid_ix, :]
        right_vals = grid[mid_ix:, :]
        left_mean = np.nanmean(left_vals)
        right_mean = np.nanmean(right_vals)
        assert (
            right_mean > left_mean
        ), f"Right (diverse) mean={right_mean:.3f} should exceed left={left_mean:.3f}"


# ---------------------------------------------------------------------------
# Tests: build_diversity_sats
# ---------------------------------------------------------------------------


class TestBuildDiversitySATs:
    def test_sat_total_matches_grid_sum(self):
        """SAT corner value should equal sum of all grid values."""
        rng = np.random.default_rng(42)
        n = 200
        types = [f"T{i % 4}" for i in range(n)]
        df = _make_cells(
            x=rng.uniform(0, 2, n).tolist(),
            y=rng.uniform(0, 2, n).tolist(),
            z_slice=[1.0] * n,
            labels=types,
        )
        div_grids, count_grids, geo = build_local_diversity_grid(
            df, "type", "Shannon", radius=0.3, grid_res=0.1
        )
        sum_sats, count_sats = build_diversity_sats(
            div_grids, count_grids, geo["nx"], geo["ny"]
        )
        # SAT corner = total sum
        grid = div_grids[1.0].copy()
        grid[np.isnan(grid)] = 0.0
        expected_sum = grid.sum()
        actual_sum = sum_sats[1.0][-1, -1]
        assert actual_sum == pytest.approx(expected_sum, rel=1e-6)

        expected_count = count_grids[1.0].sum()
        actual_count = count_sats[1.0][-1, -1]
        assert actual_count == pytest.approx(expected_count, rel=1e-6)


# ---------------------------------------------------------------------------
# Tests: find_max_mean_diversity_box
# ---------------------------------------------------------------------------


class TestFindMaxMeanDiversityBox:
    def _build_diverse_vs_uniform(self):
        """Create a dataset with a known diverse region and a uniform region.

        Left half (x < 2.0): uniform — all type A (200 cells)
        Right half (x >= 2.5): diverse — equal mix of A-H (200 cells total)
        Single z-slice at 5.0, spread over 4mm x 4mm area.
        """
        rng = np.random.default_rng(42)

        n_uni = 200
        x_uni = rng.uniform(0.0, 1.8, n_uni)
        y_uni = rng.uniform(0.0, 4.0, n_uni)
        labels_uni = ["A"] * n_uni

        n_div = 200
        x_div = rng.uniform(2.2, 4.0, n_div)
        y_div = rng.uniform(0.0, 4.0, n_div)
        types = ["A", "B", "C", "D", "E", "F", "G", "H"]
        labels_div = [types[i % 8] for i in range(n_div)]

        df = _make_cells(
            x=list(x_uni) + list(x_div),
            y=list(y_uni) + list(y_div),
            z_slice=[5.0] * (n_uni + n_div),
            labels=labels_uni + labels_div,
        )
        return df

    def _run_search(
        self,
        df,
        metric="Shannon",
        volume_um3=500,
        n_slices=1,
        symmetry="One-sided",
        radius=0.3,
        midline_x=2.0,
    ):
        div_grids, count_grids, geo = build_local_diversity_grid(
            df, "type", metric, radius=radius, grid_res=0.1
        )
        sum_sats, count_sats = build_diversity_sats(
            div_grids, count_grids, geo["nx"], geo["ny"]
        )
        z_slices = np.array(sorted(div_grids.keys()))
        return find_max_mean_diversity_box(
            sum_sats,
            count_sats,
            geo,
            z_slices,
            n_slices,
            volume_um3,
            symmetry,
            midline_x=midline_x,
        )

    def test_finds_diverse_region_shannon(self):
        """Best box should land in the diverse (right) half for Shannon."""
        df = self._build_diverse_vs_uniform()
        result = self._run_search(df, metric="Shannon", volume_um3=500)

        assert result["best"] is not None
        best = result["best"]
        box_cx = (best["x_range"][0] + best["x_range"][1]) / 2
        assert (
            box_cx >= 2.0
        ), f"Box center x={box_cx:.2f} should be >= 2.0 (diverse region)"
        assert best["diversity"] > 0

    def test_finds_diverse_region_simpson(self):
        """Best box should land in the diverse half for Simpson too."""
        df = self._build_diverse_vs_uniform()
        result = self._run_search(df, metric="Simpson", volume_um3=500)

        best = result["best"]
        assert best is not None
        box_cx = (best["x_range"][0] + best["x_range"][1]) / 2
        assert box_cx >= 2.0

    def test_finds_diverse_region_typecount(self):
        """Best box should land in the diverse half for Type count."""
        df = self._build_diverse_vs_uniform()
        result = self._run_search(df, metric="Type count", volume_um3=500)

        best = result["best"]
        assert best is not None
        box_cx = (best["x_range"][0] + best["x_range"][1]) / 2
        assert box_cx >= 2.0

    def test_multi_slice_window(self):
        """Search across 3 slices should find a result spanning them."""
        rng = np.random.default_rng(789)
        dfs = []
        for z in [4.0, 4.2, 4.4]:
            n = 100
            types = [f"T{i % 4}" for i in range(n)]
            dfs.append(
                _make_cells(
                    x=rng.uniform(0.0, 3.0, n).tolist(),
                    y=rng.uniform(0.0, 3.0, n).tolist(),
                    z_slice=[z] * n,
                    labels=types,
                )
            )
        df = pd.concat(dfs, ignore_index=True)

        result = self._run_search(df, volume_um3=600, n_slices=3)

        best = result["best"]
        assert best is not None
        assert best["z_start"] == 4.0
        assert best["z_end"] == 4.4
        assert best["z_extent_mm"] == pytest.approx(0.6, abs=0.01)

    def test_symmetric_mode_centers_at_midline(self):
        """Symmetric mode should place the box centered at the midline x."""
        rng = np.random.default_rng(101)
        midline = 2.0
        n = 400
        types = [f"T{i % 6}" for i in range(n)]
        df = _make_cells(
            x=rng.uniform(0.0, 4.0, n).tolist(),
            y=rng.uniform(0.0, 4.0, n).tolist(),
            z_slice=[5.0] * n,
            labels=types,
        )
        result = self._run_search(
            df, volume_um3=500, symmetry="Symmetric", midline_x=midline
        )

        best = result["best"]
        assert best is not None
        box_cx = (best["x_range"][0] + best["x_range"][1]) / 2
        assert (
            abs(box_cx - midline) < 0.5
        ), f"Symmetric box center x={box_cx:.2f} should be near midline={midline}"

    def test_no_result_when_too_few_points(self):
        """If too few cells for valid diversity, result should be None."""
        df = _make_cells(
            x=[0.0, 5.0],
            y=[0.0, 5.0],
            z_slice=[1.0, 1.0],
            labels=["A", "B"],
        )
        result = self._run_search(df, volume_um3=1, radius=0.05)
        assert result["best"] is None

    def test_box_coordinates_within_data_bounds(self):
        """The returned box should not extend beyond the data extent."""
        rng = np.random.default_rng(404)
        n = 300
        x_lo, x_hi = 3.0, 6.0
        y_lo, y_hi = 1.0, 4.0
        types = [f"T{i % 5}" for i in range(n)]
        df = _make_cells(
            x=rng.uniform(x_lo, x_hi, n).tolist(),
            y=rng.uniform(y_lo, y_hi, n).tolist(),
            z_slice=[5.0] * n,
            labels=types,
        )
        result = self._run_search(df, volume_um3=500)

        best = result["best"]
        assert best is not None
        assert best["x_range"][0] >= x_lo - 0.01
        assert best["x_range"][1] <= x_hi + 0.2
        assert best["y_range"][0] >= y_lo - 0.01
        assert best["y_range"][1] <= y_hi + 0.2

    def test_oversized_volume_clamps_to_grid(self):
        """Volume larger than the grid should clamp box dims, not return None."""
        rng = np.random.default_rng(505)
        n = 200
        types = [f"T{i % 5}" for i in range(n)]
        df = _make_cells(
            x=rng.uniform(0.0, 2.0, n).tolist(),
            y=rng.uniform(0.0, 2.0, n).tolist(),
            z_slice=[5.0] * n,
            labels=types,
        )
        result = self._run_search(df, volume_um3=10000)

        assert result["best"] is not None
        assert result["best"]["diversity"] > 0

    def test_unit_conversion_sanity(self):
        """10 M um^3 on 1 slice should produce a tiny box (~5 voxels)."""
        rng = np.random.default_rng(606)
        n = 2000
        types = [f"T{i % 4}" for i in range(n)]
        df = _make_cells(
            x=rng.uniform(0.0, 2.0, n).tolist(),
            y=rng.uniform(0.0, 2.0, n).tolist(),
            z_slice=[5.0] * n,
            labels=types,
        )
        result = self._run_search(df, volume_um3=10, radius=0.15)

        best = result["best"]
        assert best is not None
        box_area = best["bx_mm"] * best["by_mm"]
        assert (
            box_area < 0.5
        ), f"Box area {box_area:.2f} mm^2 should be small for 10 M um^3"

    def test_mean_diversity_is_average_not_total(self):
        """The reported diversity should be a mean, not a sum."""
        rng = np.random.default_rng(707)
        n = 400
        types = [f"T{i % 4}" for i in range(n)]
        df = _make_cells(
            x=rng.uniform(0.0, 3.0, n).tolist(),
            y=rng.uniform(0.0, 3.0, n).tolist(),
            z_slice=[5.0] * n,
            labels=types,
        )
        result = self._run_search(df, metric="Shannon", volume_um3=500)

        best = result["best"]
        assert best is not None
        # Shannon of 4 uniform types = 2.0; mean should be in reasonable range
        assert (
            0.5 < best["diversity"] < 3.0
        ), f"Mean diversity {best['diversity']:.3f} should be a reasonable average"
