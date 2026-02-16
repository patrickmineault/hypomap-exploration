#!/usr/bin/env python
"""Export a marimo notebook with all its local dependencies into a flat directory.

Usage:
    uv run python scripts/export_notebook.py

Copies notebooks/heterogeneity_map.py (as notebook.py), its local Python imports,
and data files into ./export/ as a completely flat folder. Rewrites import and data
paths, then runs marimo export html-wasm to produce a static WASM build.
"""

import re
import shutil
import subprocess
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent

NOTEBOOK = REPO_ROOT / "notebooks" / "heterogeneity_map.py"

DATA_DIR = REPO_ROOT / "data" / "processed" / "mouse_abc_subcortical"
DATA_FILES = [
    DATA_DIR / "cells_with_coords.parquet",
    DATA_DIR / "neuropeptide_expression.parquet",
    DATA_DIR / "coronal_atlas_regions.json",
]

LOCAL_MODULES = []

OUTPUT_DIR = REPO_ROOT / "export"


def rewrite_notebook(source_text: str) -> str:
    """Rewrite paths for flat directory layout."""
    out = source_text

    # Rewrite data path to use mo.notebook_location() / "public" for WASM compat
    out = re.sub(
        r"Path\(['\"]\.\.\/data\/processed\/mouse_abc_subcortical['\"]\)",
        'mo.notebook_location() / "public"',
        out,
    )

    # Add mo to the data-loading cell's parameters (needed for mo.notebook_location)
    out = out.replace(
        "def _(Path, json, np, pd):",
        "def _(Path, json, mo, np, pd):",
    )

    # Remove the unused hypomap.diversity import block
    out = re.sub(
        r"    from hypomap\.diversity import \(\n(?:        .+\n)*    \)\n",
        "",
        out,
    )

    # Insert a micropip cell before the first @app.cell
    micropip_cell = (
        "@app.cell(hide_code=True)\n"
        "def _():\n"
        "    import micropip\n"
        '    await micropip.install("plotly")\n'
        "    return\n"
        "\n"
        "\n"
    )
    out = out.replace("@app.cell\n", micropip_cell + "@app.cell\n", 1)

    return out


def main():
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # Copy and rewrite notebook
    source = NOTEBOOK.read_text()
    rewritten = rewrite_notebook(source)
    dest = OUTPUT_DIR / "notebook.py"
    dest.write_text(rewritten)
    print(f"  notebook.py (rewritten from {NOTEBOOK.name})")

    # Copy local modules
    for mod in LOCAL_MODULES:
        shutil.copy2(mod, OUTPUT_DIR / mod.name)
        print(f"  {mod.name}")

    # Copy data files into public/ (picked up by marimo WASM export)
    public_dir = OUTPUT_DIR / "public"
    public_dir.mkdir(parents=True, exist_ok=True)
    for data_file in DATA_FILES:
        shutil.copy2(data_file, public_dir / data_file.name)
        print(f"  public/{data_file.name}")

    print(
        f"\nExported {1 + len(LOCAL_MODULES) + len(DATA_FILES)} files to {OUTPUT_DIR}"
    )

    # Export as WASM-based HTML for static hosting
    print("\nRunning marimo export html-wasm...")
    subprocess.run(
        [
            "marimo",
            "export",
            "html-wasm",
            "notebook.py",
            "-o",
            "out/",
            "--mode",
            "edit",
            "--sandbox",
        ],
        cwd=OUTPUT_DIR,
        check=True,
    )
    print("Done.")


if __name__ == "__main__":
    main()
