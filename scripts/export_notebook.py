#!/usr/bin/env python
"""Export a marimo notebook with all its local dependencies into a flat directory.

Usage:
    uv run python scripts/export_notebook.py

Copies notebooks/heterogeneity_map.py (as notebook.py), its local Python imports,
and data files into ../marimo-export/heterogeneity/ as a completely flat folder.
Rewrites import paths and data paths, and prepends PEP 723 inline dependency
metadata so the notebook is fully self-describing.
"""

import re
import shutil
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent

NOTEBOOK = REPO_ROOT / "notebooks" / "heterogeneity_map.py"

DATA_DIR = REPO_ROOT / "data" / "processed" / "mouse_abc_subcortical"
DATA_FILES = [
    DATA_DIR / "cells_with_coords.parquet",
    DATA_DIR / "neuropeptide_expression.parquet",
    DATA_DIR / "coronal_atlas_regions.json",
]

LOCAL_MODULES = [
    REPO_ROOT / "hypomap" / "diversity.py",
]

OUTPUT_DIR = REPO_ROOT.parent / "marimo-export" / "heterogeneity"

# PEP 723 inline script metadata (versions from uv.lock)
PEP723_BLOCK = """\
# /// script
# requires-python = ">=3.11"
# dependencies = [
#     "marimo>=0.19.4",
#     "pandas>=2.3.3",
#     "numpy>=1.26.4,<2",
#     "plotly>=6.5.2",
#     "scipy>=1.17.0",
#     "matplotlib>=3.10.8",
# ]
# ///
"""


def rewrite_notebook(source_text: str) -> str:
    """Rewrite paths and prepend PEP 723 metadata for flat directory layout."""
    out = source_text

    # Rewrite data path: Path('../data/processed/mouse_abc_subcortical') -> Path('.')
    out = re.sub(
        r"Path\(['\"]\.\.\/data\/processed\/mouse_abc_subcortical['\"]\)",
        "Path('.')",
        out,
    )

    # Rewrite local imports: from hypomap.diversity import -> from diversity import
    out = re.sub(
        r"from hypomap\.(\w+) import",
        r"from \1 import",
        out,
    )

    # Prepend PEP 723 block after the first line (import marimo)
    # Marimo notebooks start with "import marimo"
    first_newline = out.index("\n") + 1
    out = out[:first_newline] + "\n" + PEP723_BLOCK + "\n" + out[first_newline:]

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

    # Copy data files
    for data_file in DATA_FILES:
        shutil.copy2(data_file, OUTPUT_DIR / data_file.name)
        print(f"  {data_file.name}")

    print(
        f"\nExported {1 + len(LOCAL_MODULES) + len(DATA_FILES)} files to {OUTPUT_DIR}"
    )


if __name__ == "__main__":
    main()
