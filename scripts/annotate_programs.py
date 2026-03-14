"""Annotate cNMF gene programs using Claude Opus 4.6 with web search.

For each program, extracts top genes, top cell classes/subclasses, and top
spatial regions, then asks Claude to describe what the program represents.

Supports both 10X scRNA-seq and MERFISH cNMF runs (auto-detected by checking
whether cell IDs match cells_with_coords.parquet).

Usage:
    # Annotate the MERFISH ARH/ME/VMH decomposition:
    uv run python scripts/annotate_programs.py --run-name nmf_arh_me_vmh

    # Annotate the neurons-only 10X decomposition:
    uv run python scripts/annotate_programs.py --run-name hypo_cnmf_neurons

    # Dry run (print context for one program, no API call):
    uv run python scripts/annotate_programs.py --run-name nmf_arh_me_vmh --dry-run --programs P1

    # Control parallelism:
    uv run python scripts/annotate_programs.py --run-name nmf_arh_me_vmh --max-workers 10

Output:
    data/processed/mouse_abc/cnmf/program_annotations_{run_name}.csv
"""

import argparse
import json
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

DATA_DIR = Path(__file__).parent.parent / "data"
PROCESSED_DIR = DATA_DIR / "processed" / "mouse_abc"
CNMF_DIR = PROCESSED_DIR / "cnmf"

# Directories to search for cells_with_coords (MERFISH spatial data)
COORDS_SEARCH_DIRS = [
    DATA_DIR / "processed" / "mouse_abc_extended",
    DATA_DIR / "processed" / "mouse_abc",
]

SYSTEM_PROMPT = """\
You are a neuroscience expert annotating gene programs from consensus Non-negative \
Matrix Factorization (cNMF) of mouse hypothalamus single-cell data \
(Allen Brain Cell Atlas).

Each program is a latent factor capturing co-expressed genes across cells. \
You will receive the top marker genes, the cell types that use this program most, \
and the brain regions where it is most active.

Your job: identify what biological process, cell identity, or signaling pathway \
this program likely represents. Use web search if needed to look up gene functions \
or marker gene combinations."""

USER_PROMPT_TEMPLATE = """\
Annotate cNMF program {program_id}.

## Top 30 genes by spectra z-score
{top_genes}

## Top 10 cell classes by mean usage
{top_classes}

## Top 10 cell subclasses by mean usage
{top_subclasses}

## Top 10 brain regions by mean usage (MERFISH spatial)
{top_regions}

Your ENTIRE response must be a single valid JSON object, nothing else. No markdown, no explanation, no preamble — just the JSON. Fields:
- "name": short descriptor (3-6 words, e.g. "VMH SF-1 glutamatergic identity")
- "description": 1-2 sentence evidence summary referencing specific genes and cell types
- "category": one of "neuronal-excitatory", "neuronal-inhibitory", "neuronal-mixed", \
"neuronal-neuropeptide", "glial", "vascular", "immune", "metabolic", "signaling", \
"activity-dependent", "other"
- "certainty": "high" (classic markers, unambiguous), "medium" (likely but multiple \
interpretations), or "low" (unclear or novel)
- "detailed_information": longer free-form discussion — gene function lookups, reasoning, \
alternative interpretations, literature references, anything relevant"""


def load_cnmf_data(run_name: str, k: int = 50) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load gene spectra scores and usages for a cNMF run."""
    cnmf_run_dir = CNMF_DIR / run_name

    gene_spectra = pd.read_csv(
        cnmf_run_dir / f"{run_name}.gene_spectra_score.k_{k}.dt_0_1.txt",
        sep="\t",
        index_col=0,
    )
    gene_spectra.index = [f"P{i}" for i in gene_spectra.index]

    usages = pd.read_csv(
        cnmf_run_dir / f"{run_name}.usages.k_{k}.dt_0_1.consensus.txt",
        sep="\t",
        index_col=0,
    )
    usages.columns = [f"P{c}" for c in usages.columns]

    return gene_spectra, usages


def detect_data_source(usages: pd.DataFrame) -> str:
    """Detect whether cNMF was run on MERFISH or 10X data.

    Returns 'merfish' or '10x'.
    """
    sample_ids = set(usages.index[:100])
    for coords_dir in COORDS_SEARCH_DIRS:
        coords_path = coords_dir / "cells_with_coords.parquet"
        if coords_path.exists():
            coords = pd.read_parquet(coords_path, columns=["cell_id"])
            coord_ids = set(coords["cell_id"])
            overlap = len(sample_ids & coord_ids)
            if overlap > len(sample_ids) * 0.5:
                return "merfish"
    return "10x"


def load_merfish_cell_metadata(usages: pd.DataFrame) -> pd.DataFrame:
    """Load cell metadata for MERFISH-based cNMF runs.

    Joins usages cell IDs directly against cells_with_coords.parquet.
    """
    for coords_dir in COORDS_SEARCH_DIRS:
        coords_path = coords_dir / "cells_with_coords.parquet"
        if not coords_path.exists():
            continue
        coords = pd.read_parquet(
            coords_path,
            columns=["cell_id", "class", "subclass", "supertype", "cluster", "region"],
        )
        coords = coords[~coords["region"].str.contains("unassigned", case=False)]
        coords = coords.set_index("cell_id")
        overlap = len(set(usages.index) & set(coords.index))
        if overlap > len(usages) * 0.5:
            print(f"  Matched {overlap}/{len(usages)} cells in {coords_path}")
            return coords
    raise FileNotFoundError("Could not find matching cells_with_coords.parquet")


def load_10x_cell_metadata() -> pd.DataFrame:
    """Load 10X cell metadata with taxonomy annotations."""
    from abc_atlas_access.abc_atlas_cache.abc_project_cache import AbcProjectCache

    cache_dir = DATA_DIR / "raw" / "abc_atlas_cache"
    cache = AbcProjectCache.from_cache_dir(cache_dir)
    cache.load_latest_manifest()

    cell_meta = cache.get_metadata_dataframe(
        directory="WMB-10X",
        file_name="cell_metadata",
        dtype={"cell_label": str},
    ).set_index("cell_label")

    membership = cache.get_metadata_dataframe(
        directory="WMB-taxonomy",
        file_name="cluster_to_cluster_annotation_membership",
    )
    taxonomy = membership.pivot_table(
        index="cluster_alias",
        columns="cluster_annotation_term_set_name",
        values="cluster_annotation_term_name",
        aggfunc="first",
    )
    col_map = {c: c.lower() for c in taxonomy.columns
               if c.lower() in ("class", "subclass", "supertype", "cluster")}
    taxonomy = taxonomy.rename(columns=col_map)

    cell_meta = cell_meta.join(taxonomy, on="cluster_alias")
    cell_meta = cell_meta[cell_meta["feature_matrix_label"] == "WMB-10Xv3-HY"]
    return cell_meta


def build_program_context(
    program_id: str,
    gene_spectra: pd.DataFrame,
    usages_norm: pd.DataFrame,
    cell_meta: pd.DataFrame,
    data_source: str,
    cluster_usage: pd.DataFrame | None,
) -> str:
    """Build the context string for a single program."""
    # Top genes
    scores = gene_spectra.loc[program_id].sort_values(ascending=False)
    top_genes = scores.head(30)
    genes_str = ", ".join(f"{g} ({v:.4f})" for g, v in top_genes.items())

    # Top classes/subclasses by mean usage
    merged = usages_norm[[program_id]].join(cell_meta[["class", "subclass"]], how="inner")
    class_usage = merged.groupby("class")[program_id].mean().sort_values(ascending=False).head(10)
    classes_str = "\n".join(f"- {c}: {v:.4f}" for c, v in class_usage.items())
    subclass_usage = merged.groupby("subclass")[program_id].mean().sort_values(ascending=False).head(10)
    subclasses_str = "\n".join(f"- {s}: {v:.4f}" for s, v in subclass_usage.items())

    # Top regions
    if data_source == "merfish":
        # Direct: cell_meta already has region column
        region_merged = usages_norm[[program_id]].join(cell_meta[["region"]], how="inner")
        region_usage = region_merged.groupby("region")[program_id].mean().sort_values(ascending=False).head(10)
        regions_str = "\n".join(f"- {r}: {v:.4f}" for r, v in region_usage.items())
    elif cluster_usage is not None:
        # 10X: bridge via MERFISH cluster-level usage
        merfish_path = PROCESSED_DIR / "cells_with_coords.parquet"
        if merfish_path.exists():
            merfish = pd.read_parquet(merfish_path, columns=["cluster", "region"])
            region_mapped = merfish.join(
                cluster_usage[[program_id]], on="cluster", how="inner"
            )
            region_usage = region_mapped.groupby("region")[program_id].mean().sort_values(ascending=False).head(10)
            regions_str = "\n".join(f"- {r}: {v:.4f}" for r, v in region_usage.items())
        else:
            regions_str = "(MERFISH data not available)"
    else:
        regions_str = "(MERFISH data not available)"

    return USER_PROMPT_TEMPLATE.format(
        program_id=program_id,
        top_genes=genes_str,
        top_classes=classes_str,
        top_subclasses=subclasses_str,
        top_regions=regions_str,
    )


def annotate_program(client, user_prompt: str, program_id: str, model: str = "claude-sonnet-4-20250514") -> dict:
    """Call Claude with web search to annotate one program."""
    messages = [{"role": "user", "content": user_prompt}]
    total_input_tokens = 0
    total_output_tokens = 0
    total_web_searches = 0
    n_turns = 0

    response = client.messages.create(
        model=model,
        max_tokens=4096,
        system=SYSTEM_PROMPT,
        tools=[{"type": "web_search_20250305", "name": "web_search", "max_uses": 3}],
        messages=messages,
    )
    n_turns += 1
    total_input_tokens += response.usage.input_tokens
    total_output_tokens += response.usage.output_tokens
    if hasattr(response.usage, "server_tool_use") and response.usage.server_tool_use:
        total_web_searches += getattr(response.usage.server_tool_use, "web_search_requests", 0)

    # Handle tool use loop (web search may require multiple turns)
    while response.stop_reason == "tool_use":
        messages.append({"role": "assistant", "content": response.content})
        tool_results = []
        for block in response.content:
            if block.type == "tool_use":
                tool_results.append({
                    "type": "tool_result",
                    "tool_use_id": block.id,
                    "content": "",  # server-side tool, result comes automatically
                })
        messages.append({"role": "user", "content": tool_results})
        response = client.messages.create(
            model=model,
            max_tokens=4096,
            system=SYSTEM_PROMPT,
            tools=[{"type": "web_search_20250305", "name": "web_search", "max_uses": 3}],
            messages=messages,
        )
        n_turns += 1
        total_input_tokens += response.usage.input_tokens
        total_output_tokens += response.usage.output_tokens
        if hasattr(response.usage, "server_tool_use") and response.usage.server_tool_use:
            total_web_searches += getattr(response.usage.server_tool_use, "web_search_requests", 0)

    # Extract the final text block
    text = ""
    for block in response.content:
        if block.type == "text":
            text = block.text

    # Extract JSON from response
    text = text.strip()
    start = text.find("{")
    end = text.rfind("}")
    if start != -1 and end != -1:
        text = text[start:end + 1]

    try:
        result = json.loads(text)
    except json.JSONDecodeError:
        print(f"  WARNING: Could not parse JSON for {program_id}, raw text:\n{text[:300]}")
        result = {
            "name": "parse_error",
            "description": text[:200],
            "category": "other",
            "certainty": "low",
        }
    result["program"] = program_id
    result["_input_tokens"] = total_input_tokens
    result["_output_tokens"] = total_output_tokens
    result["_web_searches"] = total_web_searches
    result["_turns"] = n_turns
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--run-name", required=True, help="cNMF run name (e.g. nmf_arh_me_vmh)")
    parser.add_argument("--programs", nargs="+", default=None, help="Specific programs to annotate (e.g. P1 P5). Default: all.")
    parser.add_argument("--dry-run", action="store_true", help="Print context for first program without calling API")
    parser.add_argument("--output", type=str, default=None, help="Output CSV path (default: auto)")
    parser.add_argument("--max-workers", type=int, default=5, help="Max parallel API calls (default: 5)")
    parser.add_argument("--k", type=int, default=50, help="Number of components K (default: 50)")
    parser.add_argument("--model", type=str, default="claude-sonnet-4-20250514", help="Claude model to use (default: claude-sonnet-4-20250514)")
    args = parser.parse_args()

    run_name = args.run_name
    k = args.k

    # Load cNMF data
    print(f"Loading cNMF data for {run_name} (K={k})...")
    gene_spectra, usages = load_cnmf_data(run_name, k=k)
    usages_norm = usages.div(usages.sum(axis=1), axis=0)
    programs = sorted(usages.columns, key=lambda x: int(x[1:]))

    if args.programs:
        programs = [p for p in args.programs if p in usages.columns]
        if not programs:
            print(f"ERROR: None of {args.programs} found in usages columns.")
            sys.exit(1)

    # Auto-detect data source
    data_source = detect_data_source(usages)
    print(f"Detected data source: {data_source}")

    # Load cell metadata
    print("Loading cell metadata...")
    cluster_usage = None
    if data_source == "merfish":
        cell_meta = load_merfish_cell_metadata(usages)
    else:
        cell_meta = load_10x_cell_metadata()
        # Precompute cluster-level usage for MERFISH bridge
        merged_cluster = usages_norm.join(cell_meta[["cluster"]], how="inner")
        cluster_usage = merged_cluster.groupby("cluster")[programs].mean()

    # Build all prompts upfront
    print(f"Building context for {len(programs)} programs...")
    prompts = {}
    for prog in programs:
        prompts[prog] = build_program_context(
            prog, gene_spectra, usages_norm, cell_meta, data_source, cluster_usage,
        )

    if args.dry_run:
        prog = programs[0]
        print(f"\n=== Context for {prog} ===")
        print(prompts[prog])
        return

    # Call Claude in parallel
    import anthropic
    from dotenv import load_dotenv
    load_dotenv()
    client = anthropic.Anthropic()

    output_path = Path(args.output) if args.output else CNMF_DIR / f"program_annotations_{run_name}_k{k}.csv"
    results = []
    n_total = len(programs)
    n_done = 0

    print(f"\nAnnotating {n_total} programs with {args.max_workers} parallel workers...")

    with ThreadPoolExecutor(max_workers=args.max_workers) as executor:
        futures = {
            executor.submit(annotate_program, client, prompts[prog], prog, args.model): prog
            for prog in programs
        }

        for future in as_completed(futures):
            prog = futures[future]
            n_done += 1
            try:
                result = future.result()
                results.append(result)
                print(f"  [{n_done}/{n_total}] {prog}: {result.get('name', '?')} [{result.get('certainty', '?')}]")
            except Exception as e:
                print(f"  [{n_done}/{n_total}] {prog}: ERROR - {e}")
                results.append({
                    "program": prog,
                    "name": "api_error",
                    "description": str(e)[:200],
                    "category": "other",
                    "certainty": "low",
                    "detailed_information": "",
                })

    # Print usage summary
    total_input = sum(r.get("_input_tokens", 0) for r in results)
    total_output = sum(r.get("_output_tokens", 0) for r in results)
    total_searches = sum(r.get("_web_searches", 0) for r in results)
    total_turns = sum(r.get("_turns", 0) for r in results)
    input_cost = total_input * 5 / 1e6  # $5/MTok for Opus 4.6, $3/MTok for Sonnet
    output_cost = total_output * 25 / 1e6  # $25/MTok for Opus 4.6, $15/MTok for Sonnet
    search_cost = total_searches * 10 / 1000  # $10/1000 searches

    print(f"\n=== Usage Summary ===")
    print(f"  Turns: {total_turns} ({total_turns/n_total:.1f} avg/program)")
    print(f"  Input tokens: {total_input:,} ({total_input/n_total:,.0f} avg/program)")
    print(f"  Output tokens: {total_output:,} ({total_output/n_total:,.0f} avg/program)")
    print(f"  Web searches: {total_searches} ({total_searches/n_total:.1f} avg/program)")
    print(f"  Est. cost (Opus 4.6 rates): ${input_cost + output_cost + search_cost:.2f}")
    print(f"    Input: ${input_cost:.2f} | Output: ${output_cost:.2f} | Search: ${search_cost:.2f}")

    # Sort by program number and save (strip internal usage fields)
    results.sort(key=lambda r: int(r["program"][1:]))
    cols = ["program", "name", "description", "category", "certainty", "detailed_information"]
    df = pd.DataFrame(results)[cols]
    df.to_csv(output_path, index=False)
    print(f"\nSaved {len(df)} annotations to {output_path}")


if __name__ == "__main__":
    main()
