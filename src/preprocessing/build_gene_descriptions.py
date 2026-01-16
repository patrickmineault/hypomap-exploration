"""Build gene descriptions lookup table from cluster names.

Extracts gene names embedded in cell type cluster names from mouse_hypomap
and mouse_abc datasets, then queries UniProt for descriptions.

Usage:
    python -m src.preprocessing.build_gene_descriptions
"""

import pandas as pd
import re
import time
import urllib.request
import urllib.parse
from pathlib import Path
from typing import Set

# Paths
DATA_DIR = Path(__file__).parent.parent.parent / "data"
PROCESSED_DIR = DATA_DIR / "processed"
OUTPUT_DIR = PROCESSED_DIR / "mouse_common"

# Non-gene terms to filter out (cell types, brain regions, etc.)
NON_GENE_TERMS = {
    # Cell types
    'astrocytes', 'astrocyte', 'astro', 'neurons', 'neuron', 'neuronal',
    'oligodendrocytes', 'oligodendrocyte', 'oligo', 'opc', 'mol',
    'microglia', 'ependymal', 'tanycytes', 'tanycyte', 'pericytes', 'pericyte',
    'endothelial', 'endo', 'vascular', 'vlmc', 'smc', 'abc', 'immune',
    # Neurotransmitter classes
    'gaba', 'glut', 'glutamatergic', 'gabaergic', 'dopa', 'dopaminergic',
    'sero', 'serotonergic', 'chol', 'cholinergic', 'glyc', 'glycinergic',
    # Brain regions
    'hy', 'vmh', 'dmh', 'arh', 'pvh', 'lha', 'ahn', 'pvpo', 'vmpo', 'mpn',
    'tu', 'me', 'zi', 'stn', 'pstn', 'mm', 'pmv', 'pmd', 'sbpv', 'pvhd',
    'rch', 'sch', 'so', 'avpv', 'mpo', 'lpo', 'ph', 'sum', 'mea', 'cnu',
    'mb', 'th', 'my', 'mge', 'lge', 'cge', 'ctx', 'hpf', 'ob', 'str', 'pal',
    # Generic terms
    'nn', 'cells', 'cell', 'type', 'cluster', 'class', 'subclass', 'supertype',
    'named', 'mixed', 'unassigned', 'other', 'unknown',
    # Common non-gene words that might appear
    'like', 'positive', 'negative', 'high', 'low', 'expressing',
}


def extract_genes_from_hypomap(df: pd.DataFrame) -> Set[str]:
    """Extract gene names from mouse_hypomap cluster columns.

    Cluster names follow pattern: "C465-344: Notch2.Gjb6.Etnppl.Astrocytes"
    Genes are dot-separated before the cell type.
    """
    genes = set()

    # Find cell type columns (C#_named)
    c_cols = [c for c in df.columns if '_named' in c]

    for col in c_cols:
        unique_values = df[col].dropna().unique()
        for val in unique_values:
            # Extract part after colon if present
            if ':' in str(val):
                val = str(val).split(':', 1)[1].strip()

            # Split by dots
            parts = str(val).split('.')
            for part in parts:
                # Clean up the part
                part = part.strip()
                # Skip if it's a non-gene term
                if part.lower() in NON_GENE_TERMS:
                    continue
                # Skip if too short or starts with number
                if len(part) < 2 or part[0].isdigit():
                    continue
                # Skip if all uppercase and very short (likely abbreviation like GABA)
                if part.isupper() and len(part) <= 4:
                    continue
                # Gene names typically start with uppercase
                if part[0].isupper():
                    genes.add(part)

    return genes


def extract_genes_from_abc(df: pd.DataFrame) -> Set[str]:
    """Extract gene names from mouse_abc cluster columns.

    Cluster names follow patterns like:
    - "129 VMH Nr5a1 Glut"
    - "318 Astro-NT NN"
    - "101 ZI Pax6 Gaba"
    """
    genes = set()

    # Cell type columns in ABC
    c_cols = ['class', 'subclass', 'supertype', 'cluster']
    c_cols = [c for c in c_cols if c in df.columns]

    for col in c_cols:
        unique_values = df[col].dropna().unique()
        for val in unique_values:
            # Split by spaces and hyphens
            parts = re.split(r'[\s\-_]+', str(val))
            for part in parts:
                # Clean up
                part = part.strip()
                # Only remove trailing numbers after underscore (like _1, _2)
                # But keep numbers that are part of gene names (Gata3, Nr5a1, Slc32a1)
                part = re.sub(r'_\d+$', '', part)

                # Skip if it's a non-gene term
                if part.lower() in NON_GENE_TERMS:
                    continue
                # Skip if too short or is just a number
                if len(part) < 2 or part.isdigit():
                    continue
                # Skip common region/type abbreviations
                if part.isupper() and len(part) <= 4:
                    continue
                # Gene names typically start with uppercase letter
                # and contain lowercase (unlike pure abbreviations)
                if part[0].isupper() and any(c.islower() for c in part):
                    genes.add(part)

    return genes


def query_uniprot_batch(genes: list, organism: str = "mouse") -> pd.DataFrame:
    """Query UniProt for gene descriptions in batches.

    Args:
        genes: List of gene symbols
        organism: "mouse" or "human"

    Returns:
        DataFrame with gene_symbol, protein_name, description, go_biological_process, uniprot_id
    """
    organism_id = "10090" if organism == "mouse" else "9606"

    results = []
    batch_size = 50

    for i in range(0, len(genes), batch_size):
        batch = genes[i:i+batch_size]
        gene_query = " OR ".join([f"gene:{g}" for g in batch])
        query = f"({gene_query}) AND (organism_id:{organism_id})"

        # UniProt REST API
        url = "https://rest.uniprot.org/uniprotkb/search"
        params = {
            "query": query,
            "fields": "gene_names,protein_name,cc_function,go_p,accession",
            "format": "tsv",
            "size": "500"
        }

        try:
            full_url = f"{url}?{urllib.parse.urlencode(params)}"
            with urllib.request.urlopen(full_url, timeout=30) as response:
                content = response.read().decode('utf-8')

            lines = content.strip().split('\n')
            if len(lines) > 1:
                header = lines[0].split('\t')
                for line in lines[1:]:
                    fields = line.split('\t')
                    if len(fields) >= 5:
                        # Parse gene names (first one is primary)
                        gene_names = fields[0].split()
                        primary_gene = gene_names[0] if gene_names else ""

                        results.append({
                            'gene_symbol': primary_gene,
                            'protein_name': fields[1] if len(fields) > 1 else "",
                            'description': fields[2] if len(fields) > 2 else "",
                            'go_biological_process': fields[3].replace('; ', '; ') if len(fields) > 3 else "",
                            'uniprot_id': fields[4] if len(fields) > 4 else "",
                        })
        except Exception as e:
            print(f"  Warning: UniProt query failed for batch {i//batch_size + 1}: {e}")

        # Rate limiting
        if i + batch_size < len(genes):
            time.sleep(0.5)

        if (i + batch_size) % 100 == 0 or i + batch_size >= len(genes):
            print(f"  Processed {min(i + batch_size, len(genes))}/{len(genes)} genes...")

    return pd.DataFrame(results)


def build_gene_descriptions():
    """Main function to build gene descriptions from cluster names."""

    print("=== Building Gene Descriptions from Cluster Names ===\n")

    all_genes = set()

    # 1. Extract from mouse_hypomap
    hypomap_path = PROCESSED_DIR / "mouse_hypomap" / "cells_with_coords.parquet"
    if hypomap_path.exists():
        print("Loading mouse_hypomap...")
        hypomap_df = pd.read_parquet(hypomap_path)
        hypomap_genes = extract_genes_from_hypomap(hypomap_df)
        print(f"  Extracted {len(hypomap_genes)} candidate genes from mouse_hypomap")
        all_genes.update(hypomap_genes)
    else:
        print(f"  Warning: {hypomap_path} not found, skipping")

    # 2. Extract from mouse_abc
    abc_path = PROCESSED_DIR / "mouse_abc" / "cell_metadata.parquet"
    if abc_path.exists():
        print("Loading mouse_abc...")
        abc_df = pd.read_parquet(abc_path)
        abc_genes = extract_genes_from_abc(abc_df)
        print(f"  Extracted {len(abc_genes)} candidate genes from mouse_abc")
        all_genes.update(abc_genes)
    else:
        print(f"  Warning: {abc_path} not found, skipping")

    print(f"\nTotal unique candidate genes: {len(all_genes)}")

    # 3. Query UniProt for descriptions
    print("\nQuerying UniProt for gene descriptions...")
    gene_list = sorted(list(all_genes))
    uniprot_df = query_uniprot_batch(gene_list, organism="mouse")

    # 4. Filter to only genes we found and deduplicate
    # Keep only entries where gene_symbol matches one of our extracted genes (case-insensitive)
    gene_lower_map = {g.lower(): g for g in all_genes}
    uniprot_df['matched_gene'] = uniprot_df['gene_symbol'].str.lower().map(gene_lower_map)
    uniprot_df = uniprot_df.dropna(subset=['matched_gene'])

    # Use original case from our extraction
    uniprot_df['gene_symbol'] = uniprot_df['matched_gene']
    uniprot_df = uniprot_df.drop(columns=['matched_gene'])

    # Deduplicate by gene_symbol, keeping first (usually best annotated)
    uniprot_df = uniprot_df.drop_duplicates(subset=['gene_symbol'], keep='first')

    # Add source column
    uniprot_df['source'] = 'cluster_marker'

    # Reorder columns
    uniprot_df = uniprot_df[['gene_symbol', 'protein_name', 'description',
                             'go_biological_process', 'source', 'uniprot_id']]

    # 5. Save
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    output_path = OUTPUT_DIR / "gene_descriptions.csv"
    uniprot_df.to_csv(output_path, index=False)

    print(f"\n=== Summary ===")
    print(f"Genes extracted from clusters: {len(all_genes)}")
    print(f"Genes with UniProt descriptions: {len(uniprot_df)}")
    print(f"Saved to: {output_path}")

    # Show some examples
    print(f"\nSample entries:")
    print(uniprot_df.head(10).to_string())

    return uniprot_df


if __name__ == "__main__":
    build_gene_descriptions()
