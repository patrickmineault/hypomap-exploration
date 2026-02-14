import marimo

__generated_with = "0.19.4"
app = marimo.App(width="medium")


@app.cell
def _():
    import marimo as mo
    import pandas as pd
    import numpy as np
    return mo, np, pd


@app.cell
def _(mo):
    mo.md(r"""
    # Coordinate System Comparison

    Compare the two coordinate systems available in the ABC MERFISH atlas:

    - **Reconstructed coordinates** (`reconstructed_coordinates.csv`): MERFISH-reconstructed space.
      Axes: x = mediolateral, y = dorsoventral, z = anteroposterior (discrete slices).
    - **CCF coordinates** (`ccf_coordinates.csv`): Allen Common Coordinate Framework.
      Axes: x = anteroposterior, y = dorsoventral, z = mediolateral.

    For a coronal view, we plot **ML vs DV**:
    - Reconstructed: `x_recon` vs `y_recon`
    - CCF: `z_ccf` vs `y_ccf`
    """)
    return


@app.cell
def _(pd):
    # Load both coordinate files (join on cell_label)
    _recon = pd.read_csv(
        '../data/raw/abc_atlas_cache/metadata/MERFISH-C57BL6J-638850-CCF/20231215/reconstructed_coordinates.csv',
        dtype={'cell_label': str},
    )
    _recon = _recon.rename(columns={'x': 'x_recon', 'y': 'y_recon', 'z': 'z_recon'})

    _ccf = pd.read_csv(
        '../data/raw/abc_atlas_cache/metadata/MERFISH-C57BL6J-638850-CCF/20231215/ccf_coordinates.csv',
        dtype={'cell_label': str},
    )
    _ccf = _ccf.rename(columns={'x': 'x_ccf', 'y': 'y_ccf', 'z': 'z_ccf'})

    coords = _recon.merge(
        _ccf[['cell_label', 'x_ccf', 'y_ccf', 'z_ccf']],
        on='cell_label', how='inner',
    )
    print(f"Joined {len(coords):,} cells with both coordinate systems")
    print(f"\nReconstructed ranges:")
    print(f"  x_recon (ML): {coords.x_recon.min():.2f} – {coords.x_recon.max():.2f}")
    print(f"  y_recon (DV): {coords.y_recon.min():.2f} – {coords.y_recon.max():.2f}")
    print(f"  z_recon (AP): {coords.z_recon.min():.2f} – {coords.z_recon.max():.2f}")
    print(f"\nCCF ranges:")
    print(f"  x_ccf (AP): {coords.x_ccf.min():.2f} – {coords.x_ccf.max():.2f}")
    print(f"  y_ccf (DV): {coords.y_ccf.min():.2f} – {coords.y_ccf.max():.2f}")
    print(f"  z_ccf (ML): {coords.z_ccf.min():.2f} – {coords.z_ccf.max():.2f}")
    return (coords,)


@app.cell
def _(coords, mo):
    _slices = sorted(coords.z_recon.unique())
    slice_dropdown = mo.ui.dropdown(
        options={f"z_recon={z:.1f}": z for z in _slices},
        value=f"z_recon={_slices[len(_slices) // 2]:.1f}",
        label="Coronal slice (reconstructed z)",
    )
    slice_dropdown
    return (slice_dropdown,)


@app.cell
def _(coords, np, slice_dropdown):
    import matplotlib.pyplot as plt
    import matplotlib

    matplotlib.rcParams['figure.dpi'] = 150

    _z = slice_dropdown.value
    _slice_df = coords[coords.z_recon == _z]
    _n = len(_slice_df)

    # Subsample if too many cells for scatter
    _max_pts = 50000
    if _n > _max_pts:
        _slice_df = _slice_df.sample(_max_pts, random_state=42)

    _fig, (_ax1, _ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Plot 1: Reconstructed (x_recon=ML vs y_recon=DV)
    _ax1.scatter(_slice_df.x_recon, _slice_df.y_recon, s=0.3, alpha=0.3, c='steelblue', rasterized=True)
    _ax1.axvline(5.7, color='red', lw=0.8, ls='--', alpha=0.7)
    _ax1.axhline(5.3, color='red', lw=0.8, ls='--', alpha=0.7)
    _ax1.plot(5.7, 5.4, 'r+', ms=12, mew=2)
    _ax1.set_xlabel('x_recon (ML, mm)')
    _ax1.set_ylabel('y_recon (DV, mm)')
    _ax1.set_title(f'Reconstructed coords (z_recon={_z:.1f})')
    _ax1.set_aspect('equal')
    _ax1.invert_yaxis()

    # Plot 2: CCF (z_ccf=ML vs y_ccf=DV)
    _ax2.scatter(_slice_df.z_ccf, _slice_df.y_ccf, s=0.3, alpha=0.3, c='darkorange', rasterized=True)
    _ax2.axvline(5.7, color='red', lw=0.8, ls='--', alpha=0.7)
    _ax2.axhline(5.4, color='red', lw=0.8, ls='--', alpha=0.7)
    _ax2.plot(5.7, 5.4, 'r+', ms=12, mew=2)
    _ax2.set_xlabel('z_ccf (ML, mm)')
    _ax2.set_ylabel('y_ccf (DV, mm)')
    _ax2.set_title(f'CCF coords (same cells)')
    _ax2.set_aspect('equal')
    _ax2.invert_yaxis()

    _fig.suptitle(
        f'Coronal slice z_recon={_z:.1f} — {_n:,} cells '
        f'(red cross = 5.7, 4.4)',
        fontsize=11,
    )
    _fig.tight_layout()

    # Stats
    print(f"Cells in slice: {_n:,}")
    print(f"\nReconstructed center: x_recon={np.median(_slice_df.x_recon):.2f}, y_recon={np.median(_slice_df.y_recon):.2f}")
    print(f"CCF center:           z_ccf={np.median(_slice_df.z_ccf):.2f}, y_ccf={np.median(_slice_df.y_ccf):.2f}")

    plt.gca()
    return


@app.cell
def _(coords, np):
    coords[np.round(coords.z_recon*10) == 64]
    return


if __name__ == "__main__":
    app.run()
