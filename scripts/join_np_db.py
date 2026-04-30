#!/usr/bin/env python3
"""
join_database.py
================
Reconstructs the final merged database from its decomposed components:

    np_map_lewislab.csv    — gene pairs + binding_rank (mechanically derived)
    extra_pairs.csv        — gene pairs NOT in Lewis Lab (manually curated)
    extra_info.csv         — annotation layer (System, Ligand_Name, etc.)
    processing_enzymes.csv — enzyme disambiguation

Output:
    data/processed/mouse_common/np_map.csv     — primary receptors only
    data/processed/mouse_common/np_map_all.csv — primary + secondary receptors

Usage:
    python3 scripts/join_np_db.py
"""

import csv
import sys
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('--lewislab-csv', default='data/processed/mouse_common/np_map_lewislab.csv',
                    help='Path to np_map_lewislab.csv (from build_np_db.py)')
parser.add_argument('--curated-dir', default='data/generated/mouse_common',
                    help='Directory containing extra_pairs.csv, extra_info.csv, processing_enzymes.csv')
parser.add_argument('--output', default='data/processed/mouse_common/np_map.csv',
                    help='Output path for primary-only np_map.csv')
parser.add_argument('--output-all', default='data/processed/mouse_common/np_map_all.csv',
                    help='Output path for all pairs (primary + secondary)')
args = parser.parse_args()

CURATED = args.curated_dir
os.makedirs(os.path.dirname(args.output), exist_ok=True)

# ── Load gene pairs ────────────────────────────────────────────────────

# Lewis Lab pairs (keyed by (Ligand_Gene, Receptor_Gene))
lewislab = {}
with open(args.lewislab_csv) as f:
    for row in csv.DictReader(f):
        key = (row['Ligand_Gene'], row['Receptor_Gene'])
        lewislab[key] = {
            'databases': row['databases'],
            'n_databases': row['n_databases'],
            'binding_rank': row['binding_rank'],
        }

# Extra pairs
extra = {}
with open(f'{CURATED}/extra_pairs.csv') as f:
    for row in csv.DictReader(f):
        key = (row['Ligand_Gene'], row['Receptor_Gene'])
        extra[key] = {
            'References': row['References'],
            'Source_Note': row['Source_Note'],
        }

# ── Load annotations ──────────────────────────────────────────────────

info_rows = []
with open(f'{CURATED}/extra_info.csv') as f:
    info_rows = list(csv.DictReader(f))

# ── Load processing enzymes ───────────────────────────────────────────

pe_lookup = {}
with open(f'{CURATED}/processing_enzymes.csv') as f:
    for row in csv.DictReader(f):
        key = (row['Ligand_Gene'], row['Receptor_Gene'], row['Ligand_Name'])
        pe_lookup[key] = row['Processing_Enzymes']

# ── Reconstruct ───────────────────────────────────────────────────────

output_rows = []
for info in info_rows:
    lg = info['Ligand_Gene']
    rg = info['Receptor_Gene']  # may have ';' for heteromers
    ln = info['Ligand_Name']

    ll = lewislab.get((lg, rg))
    ep = extra.get((lg, rg))

    # Determine Source
    source_parts = []
    if ll:
        if '2023-Zhao' in ll['databases']:
            source_parts.append('NeuronChat')
        source_parts.append('IUPHAR')
    elif ep:
        for s in ep['Source_Note'].split(';'):
            source_parts.append(s.strip())
        source_parts.append('IUPHAR')

    # Determine In_NeuronChat
    in_nc = 'yes' if (ll and '2023-Zhao' in ll['databases']) else 'no'

    # Determine Binding_Rank
    if ll:
        binding_rank = ll['binding_rank']
    elif ep:
        binding_rank = 'primary'  # extra pairs are manually curated primaries
    else:
        binding_rank = ''

    # Get Processing_Enzymes
    pe = pe_lookup.get((lg, rg, ln), '')

    # Get References (only extra_pairs have explicit references)
    refs = ep['References'] if ep else ''

    output_rows.append({
        'System': info['System'],
        'Ligand_Name': info['Ligand_Name'],
        'Receptor_Name': info['Receptor_Name'],
        'Ligand_Gene': lg,
        'Receptor_Gene': rg,
        'Ligand_Class': info['Ligand_Class'],
        'Functional_Role': info['Functional_Role'],
        'Hypothalamic_Nuclei': info['Hypothalamic_Nuclei'],
        'Notes': info['Notes'],
        'Processing_Enzymes': pe,
        'References': refs,
        'Source': ';'.join(source_parts),
        'In_NeuronChat': in_nc,
        'Binding_Rank': binding_rank,
    })

# Write reconstructed file
fieldnames = ['System', 'Ligand_Name', 'Receptor_Name', 'Ligand_Gene', 'Receptor_Gene',
              'Ligand_Class', 'Functional_Role', 'Hypothalamic_Nuclei', 'Notes',
              'Processing_Enzymes', 'References', 'Source', 'In_NeuronChat', 'Binding_Rank']

primary_rows = [r for r in output_rows if r['Binding_Rank'] != 'secondary']

for path, rows in [(args.output, primary_rows), (args.output_all, output_rows)]:
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for row in rows:
            w.writerow(row)
    print(f"Wrote {path}: {len(rows)} rows")

# ── Summary stats ─────────────────────────────────────────────────────

n_primary = sum(1 for r in output_rows if r['Binding_Rank'] == 'primary')
n_secondary = sum(1 for r in output_rows if r['Binding_Rank'] == 'secondary')
n_nc = sum(1 for r in output_rows if r['In_NeuronChat'] == 'yes')
print(f"  Primary: {n_primary}, Secondary: {n_secondary}")
print(f"  In NeuronChat: {n_nc}")

# Sanity checks
orphaned = []
for info in info_rows:
    key = (info['Ligand_Gene'], info['Receptor_Gene'])
    if key not in lewislab and key not in extra:
        orphaned.append(key)

if orphaned:
    print(f"\n  WARNING: {len(orphaned)} annotation rows have no matching gene pair:")
    for lg, rg in orphaned:
        print(f"    {lg} → {rg}")
else:
    print("  OK: All annotation rows matched to a gene pair source")
