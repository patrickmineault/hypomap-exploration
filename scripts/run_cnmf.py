"""Run cNMF (consensus Non-negative Matrix Factorization) on hypothalamus scRNA-seq.

Wraps the cNMF pipeline (prepare → factorize → combine → consensus) with argparse.
Input is the CPM expression h5ad from extract_scrna_expression.py.

Steps:
  prepare    - Run cNMF prepare on the h5ad
  factorize  - Run NMF factorization (parallelizable across workers)
  combine    - Combine factorization results
  consensus  - Extract consensus spectra for a chosen K

Usage:
    # Quick trial run (500 cells, K=10, 5 iterations) to test the pipeline:
    python scripts/run_cnmf.py --trial-run

    # Run all steps locally (small test):
    python scripts/run_cnmf.py --step all --k-values 50 --n-iter 5

    # Run just prepare (local):
    python scripts/run_cnmf.py --step prepare --k-values 50

    # Factorize with parallel workers (on VM):
    python scripts/run_cnmf.py --step factorize --worker-index 0 --total-workers 28

    # Combine + consensus:
    python scripts/run_cnmf.py --step combine
    python scripts/run_cnmf.py --step consensus --selected-k 50
"""

import argparse
from pathlib import Path

DATA_DIR = Path(__file__).parent.parent / "data"
PROCESSED_DIR = DATA_DIR / "processed" / "mouse_abc"
CNMF_DIR = PROCESSED_DIR / "cnmf"
INPUT_H5AD = PROCESSED_DIR / "scrna_expression_cpm.h5ad"

# cNMF run name (used as prefix for all output files)
CNMF_NAME = "hypo_cnmf"

# Default parameters
DEFAULT_K_VALUES = [50]
DEFAULT_N_ITER = 100
DEFAULT_N_TOP_GENES = 3000
DEFAULT_SEED = 42


def subsample_h5ad(input_path: Path, output_path: Path, max_cells: int):
    """Create a subsampled h5ad for trial runs."""
    import anndata
    import numpy as np

    print(f"Subsampling {input_path} to {max_cells} cells...")
    adata = anndata.read_h5ad(input_path, backed="r")
    n_total = adata.n_obs
    print(f"  Total cells: {n_total}")

    rng = np.random.default_rng(DEFAULT_SEED)
    idx = np.sort(rng.choice(n_total, size=min(max_cells, n_total), replace=False))
    adata_sub = adata[idx].to_memory()
    adata.file.close()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    adata_sub.write_h5ad(output_path)
    print(f"  Saved {adata_sub.n_obs} cells to {output_path}")
    return output_path


def run_prepare(h5ad_path: Path, k_values: list[int], n_iter: int, n_top_genes: int, seed: int):
    """Run cNMF prepare step."""
    from cnmf import cNMF

    print(f"\nRunning cNMF prepare (K={k_values}, n_iter={n_iter})...")
    print(f"  Input: {h5ad_path}")
    cnmf_obj = cNMF(
        output_dir=str(CNMF_DIR),
        name=CNMF_NAME,
    )
    cnmf_obj.prepare(
        counts_fn=str(h5ad_path),
        components=k_values,
        n_iter=n_iter,
        seed=seed,
        num_highvar_genes=n_top_genes,
    )
    print("  Prepare complete.")


def run_factorize(worker_index: int, total_workers: int):
    """Run cNMF factorize step (one worker)."""
    from cnmf import cNMF

    print(f"Running cNMF factorize (worker {worker_index}/{total_workers})...")
    cnmf_obj = cNMF(
        output_dir=str(CNMF_DIR),
        name=CNMF_NAME,
    )
    cnmf_obj.factorize(worker_i=worker_index, total_workers=total_workers)
    print(f"  Worker {worker_index} complete.")


def run_combine():
    """Combine factorization results."""
    from cnmf import cNMF

    print("Running cNMF combine...")
    cnmf_obj = cNMF(
        output_dir=str(CNMF_DIR),
        name=CNMF_NAME,
    )
    cnmf_obj.combine()
    print("  Combine complete.")


def run_consensus(selected_k: int, density_threshold: float = 0.1):
    """Extract consensus spectra for a chosen K."""
    from cnmf import cNMF

    print(f"Running cNMF consensus (K={selected_k}, density_threshold={density_threshold})...")
    cnmf_obj = cNMF(
        output_dir=str(CNMF_DIR),
        name=CNMF_NAME,
    )
    cnmf_obj.consensus(k=selected_k, density_threshold=density_threshold)
    print("  Consensus complete.")

    # Check outputs
    result_dir = CNMF_DIR / CNMF_NAME
    print(f"\nResults in {result_dir}:")
    for f in sorted(result_dir.iterdir()):
        size_mb = f.stat().st_size / 1e6
        print(f"  {f.name} ({size_mb:.1f} MB)")


def run_all(h5ad_path: Path, k_values: list[int], n_iter: int, n_top_genes: int, seed: int, total_workers: int, selected_k: int):
    """Run the full cNMF pipeline."""
    run_prepare(h5ad_path, k_values, n_iter, n_top_genes, seed)

    print(f"\n--- Factorizing with {total_workers} workers ---")
    for i in range(total_workers):
        run_factorize(worker_index=i, total_workers=total_workers)

    run_combine()
    run_consensus(selected_k)


def main():
    parser = argparse.ArgumentParser(
        description="Run cNMF on hypothalamus scRNA-seq expression data"
    )
    parser.add_argument(
        "--step",
        choices=["prepare", "factorize", "combine", "consensus", "all"],
        default="all",
        help="Which cNMF step to run (default: all)",
    )
    parser.add_argument(
        "--k-values",
        type=int,
        nargs="+",
        default=DEFAULT_K_VALUES,
        help=f"K values (number of components) to test (default: {DEFAULT_K_VALUES})",
    )
    parser.add_argument(
        "--selected-k",
        type=int,
        default=50,
        help="K to use for consensus step (default: 50)",
    )
    parser.add_argument(
        "--n-iter",
        type=int,
        default=DEFAULT_N_ITER,
        help=f"Number of NMF iterations per K (default: {DEFAULT_N_ITER})",
    )
    parser.add_argument(
        "--n-top-genes",
        type=int,
        default=DEFAULT_N_TOP_GENES,
        help=f"Number of high-variance genes to use (default: {DEFAULT_N_TOP_GENES})",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=DEFAULT_SEED,
        help=f"Random seed (default: {DEFAULT_SEED})",
    )
    parser.add_argument(
        "--worker-index",
        type=int,
        default=0,
        help="Worker index for factorize step (default: 0)",
    )
    parser.add_argument(
        "--total-workers",
        type=int,
        default=1,
        help="Total number of workers for factorize step (default: 1)",
    )
    parser.add_argument(
        "--density-threshold",
        type=float,
        default=0.1,
        help="Density threshold for consensus (default: 0.1)",
    )
    parser.add_argument(
        "--trial-run",
        action="store_true",
        help="Quick end-to-end test: 500 cells, K=10, 5 iterations, 1 worker",
    )

    args = parser.parse_args()

    # Override params for trial run
    if args.trial_run:
        print("=== TRIAL RUN: 500 cells, K=10, 5 iterations ===\n")
        args.k_values = [10]
        args.selected_k = 10
        args.n_iter = 5
        args.n_top_genes = 500
        args.total_workers = 1
        args.step = "all"

    # Determine input h5ad
    h5ad_path = INPUT_H5AD
    if args.trial_run:
        trial_h5ad = CNMF_DIR / "expression_trial.h5ad"
        h5ad_path = subsample_h5ad(INPUT_H5AD, trial_h5ad, max_cells=500)

    if args.step == "prepare":
        run_prepare(h5ad_path, args.k_values, args.n_iter, args.n_top_genes, args.seed)
    elif args.step == "factorize":
        run_factorize(args.worker_index, args.total_workers)
    elif args.step == "combine":
        run_combine()
    elif args.step == "consensus":
        run_consensus(args.selected_k, args.density_threshold)
    elif args.step == "all":
        run_all(
            h5ad_path,
            args.k_values,
            args.n_iter,
            args.n_top_genes,
            args.seed,
            args.total_workers,
            args.selected_k,
        )


if __name__ == "__main__":
    main()
