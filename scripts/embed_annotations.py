"""Precompute sentence-BERT embeddings for cNMF program annotations.

Usage:
    uv run python scripts/embed_annotations.py --run-name nmf_arh_me_vmh

Output:
    data/processed/mouse_abc/cnmf/annotation_embeddings_{run_name}.npz
"""

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
from sentence_transformers import SentenceTransformer

DATA_DIR = Path(__file__).parent.parent / "data"
CNMF_DIR = DATA_DIR / "processed" / "mouse_abc" / "cnmf"


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--run-name", required=True)
    args = parser.parse_args()

    run_name = args.run_name

    # Load all annotation files
    annot_files = sorted(CNMF_DIR.glob(f"program_annotations_{run_name}*.csv"))
    all_data = {}
    for f in annot_files:
        m = re.search(r"_k(\d+)\.csv$", f.name)
        if m:
            k_val = int(m.group(1))
        elif f.name == f"program_annotations_{run_name}.csv":
            k_val = 50
        else:
            continue
        all_data[k_val] = pd.read_csv(f)

    if not all_data:
        print(f"No annotation files found for {run_name}")
        return

    print(f"Found annotations for K={sorted(all_data.keys())}")

    # Build texts
    model = SentenceTransformer("all-MiniLM-L6-v2")
    results = {}
    for k_val in sorted(all_data):
        df = all_data[k_val]
        labels = [f"{row['program']}: {row['name']}" for _, row in df.iterrows()]
        texts = [f"{row['name']}. {row['description']}" for _, row in df.iterrows()]
        embeddings = model.encode(texts, normalize_embeddings=True)
        results[f"labels_k{k_val}"] = np.array(labels)
        results[f"embeddings_k{k_val}"] = np.array(embeddings)
        print(f"  K={k_val}: {len(labels)} programs -> {embeddings.shape}")

    results["ks"] = np.array(sorted(all_data.keys()))

    out_path = CNMF_DIR / f"annotation_embeddings_{run_name}.npz"
    np.savez(out_path, **results)
    print(f"\nSaved to {out_path}")


if __name__ == "__main__":
    main()
