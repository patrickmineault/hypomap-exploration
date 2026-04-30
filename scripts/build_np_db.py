#!/usr/bin/env python3
"""
build_decomposition.py
======================
Builds np_map_lewislab.csv from Lewis Lab GitHub repos, then validates
the three curated CSV files against it.

Files produced (mechanically derived — regenerated each run):
  np_map_lewislab.csv    — Gene pairs from ≥1 Lewis Lab database (excl.
                           Dimitrov-only). Includes binding_rank, heteromer
                           notation, and database provenance.

Files validated (curated — never overwritten):
  extra_pairs.csv        — Gene pairs absent from all Lewis Lab databases.
  extra_info.csv         — Annotation layer (System, Ligand_Name, etc.).
  processing_enzymes.csv — Enzyme disambiguation for pro-peptide precursors.

Join logic (implemented in join_database.py):
    (np_map_lewislab UNION extra_pairs)       — all gene pairs
    LEFT JOIN extra_info                       — annotations
    LEFT JOIN processing_enzymes               — enzyme disambiguation
    = np_map_reconstructed.csv
"""

import csv
import openpyxl

import os
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--lewislab-dir', default='data/raw/lewislab',
                    help='Directory containing Lewis Lab raw files')
parser.add_argument('--curated-dir', default='data/generated/mouse_common',
                    help='Directory containing curated CSV files')
parser.add_argument('--output-dir', default='data/processed/mouse_common',
                    help='Directory for output np_map_lewislab.csv')
args = parser.parse_args()

BASE = args.lewislab_dir
OUT  = args.output_dir
CURATED = args.curated_dir
os.makedirs(OUT, exist_ok=True)

# ════════════════════════════════════════════════════════════════════════
# PART A: Rebuild np_map_lewislab.csv (mechanical, auditable)
# ════════════════════════════════════════════════════════════════════════

def normalize_gene(g):
    g = g.strip()
    if not g: return g
    if '&' in g: return ';'.join(normalize_gene(p) for p in g.split('&'))
    if g == g.upper() and len(g) > 1: return g[0].upper() + g[1:].lower()
    return g[0].upper() + g[1:]

# Parse all databases
databases = {}

# 2023-Zhao
pairs = set()
with open(f'{BASE}/Mouse-2023-Zhao-LR-pairs.tsv') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        name = row['interaction_name'].strip()
        if '_' in name:
            lig, rec = name.split('_', 1)
            pairs.add((normalize_gene(lig), normalize_gene(rec)))
databases['2023-Zhao'] = pairs

# 2022-Dimitrov
pairs = set()
with open(f'{BASE}/Mouse-2022-Dimitrov-LR-pairs.csv') as f:
    for row in csv.DictReader(f):
        pairs.add((normalize_gene(row['source_genesymbol']),
                    normalize_gene(row['target_genesymbol'])))
databases['2022-Dimitrov'] = pairs

# 2020-Jin
pairs = set()
with open(f'{BASE}/Mouse-2020-Jin-LR-pairs.csv') as f:
    for row in csv.DictReader(f):
        lig = normalize_gene(row.get('ligand_symbol',''))
        rec_raw = row.get('receptor_symbol','').strip()
        if lig and rec_raw:
            for r in rec_raw.split('&'):
                pairs.add((lig, normalize_gene(r)))
databases['2020-Jin'] = pairs

# 2020-Shao
pairs = set()
with open(f'{BASE}/Mouse-2020-Shao-LR-pairs.txt') as f:
    for row in csv.DictReader(f, delimiter='\t'):
        pairs.add((normalize_gene(row['ligand_gene_symbol']),
                    normalize_gene(row['receptor_gene_symbol'])))
databases['2020-Shao'] = pairs

# XLSX helper
def parse_xlsx(fname, lig_col, rec_col):
    pairs = set()
    wb = openpyxl.load_workbook(f'{BASE}/{fname}', read_only=True)
    ws = wb.active
    hdr = None
    for row in ws.iter_rows(values_only=True):
        if hdr is None: hdr = row; continue
        lig = normalize_gene(str(row[lig_col]) if row[lig_col] else '')
        rec = normalize_gene(str(row[rec_col]) if row[rec_col] else '')
        if lig and rec: pairs.add((lig, rec))
    wb.close()
    return pairs

databases['2020-Baccin'] = parse_xlsx('Mouse-2020-Baccin-LR-pairs.xlsx', 1, 2)
databases['2020-Cain']   = parse_xlsx('Mouse-2020-Cain-LR-pairs.xlsx', 0, 1)
databases['2019-Sheikh'] = parse_xlsx('Mouse-2019-Sheikh-LR-pairs.xlsx', 0, 1)
databases['2018-Skelly'] = parse_xlsx('Mouse-2018-Skelly-LR-pairs.xlsx', 1, 3)
databases['2016-Yuzwa']  = parse_xlsx('Mouse-2016-Yuzwa-LR-pairs.xlsx', 0, 1)
databases['2016-Ding']   = parse_xlsx('Mouse-2016-Ding-LR-pairs.xlsx', 0, 1)

for name, p in databases.items():
    print(f"  {name}: {len(p)} pairs")

# Deduplicate
db_order = ['2023-Zhao','2022-Dimitrov','2020-Jin','2020-Shao','2020-Baccin',
            '2020-Cain','2019-Sheikh','2018-Skelly','2016-Yuzwa','2016-Ding']

all_pairs = {}
for db_name in db_order:
    for pair in databases[db_name]:
        if pair not in all_pairs:
            all_pairs[pair] = [db_name]
        elif db_name not in all_pairs[pair]:
            all_pairs[pair].append(db_name)

# Exclude Dimitrov-only entries
non_dimitrov_only = {k: v for k, v in all_pairs.items()
                     if not (len(v) == 1 and v[0] == '2022-Dimitrov')}

print(f"\nAll unique pairs: {len(all_pairs)}")
print(f"After excluding Dimitrov-only: {len(non_dimitrov_only)}")

# ── Gene set definitions ──────────────────────────────────────────────
# NP precursor genes (GO terms + manual additions)
all_np_ligands = {
    'Adcyap1','Agrp','Avp','Cartpt','Cck','Cort','Gal','Ghrh','Grp',
    'Hcrt','Nmb','Npff','Npy','Nts','Oxt','Pdyn','Ppy','Prlh',
    'Pyy','Qrfp','Spx','Ucn3','Vgf','Vip',
    'Pomc','Penk','Pnoc','Tac1','Tac2','Crh','Kiss1','Ghrl','Gnrh1','Sst',
    'Sct','Pmch','Calca','Calcb','Pthlh','Pth','Pth2','Agt','Apln','Apela',
    'Edn1','Edn2','Edn3','Gcg','Gast','Npvf','Nps','Npw','Npb','Rln1','Rln3',
    'Insl3','Insl5','Uts2','Uts2b','Adm','Adm2','Nmu','Nms','Iapp','Smim20',
    'Pcsk1n','Prok1','Prok2','Gip','Galp',
    'Trh','Ucn','Ucn2',
}

# NP GPCRs (excludes receptor tyrosine kinases, guanylyl cyclases, RAMPs)
all_np_receptors = {
    'Brs3','Cckbr','Gpr139','Gpr171','Gpr37','Gpr83','Grpr','Kiss1r',
    'Nmbr','Nmur1','Nmur2','Npbwr1','Npffr2','Npsr1','Npy1r','Npy6r','Ntsr1',
    'Prlhr','Gpr50','Gpr173','Gpr107',
    'Hcrtr1','Hcrtr2','Mc3r','Mc4r','Mc5r','Oxtr','Avpr1a','Avpr1b','Avpr2',
    'Crhr1','Crhr2','Tacr1','Tacr2','Tacr3','Oprm1','Oprd1','Oprk1','Oprl1',
    'Sstr1','Sstr2','Sstr3','Sstr4','Sstr5','Galr1','Galr2','Galr3','Ghsr',
    'Gnrhr','Vipr1','Vipr2','Adcyap1r1','Npy2r','Npy4r','Npy5r','Cckar',
    'Gcgr','Glp1r','Glp2r','Gipr','Sctr','Ghrhr','Mchr1','Npffr1',
    'Ntsr2','Rxfp1','Rxfp2','Rxfp3','Rxfp4','Calcr','Calcrl','Trhr','Trhr2',
    'Ednra','Ednrb','Agtr1a','Agtr1b','Agtr2','Mas1','Aplnr','Pth1r','Pth2r',
    'Prokr1','Prokr2','Uts2r','Bdkrb2','Qrfpr',
    'Mrgpra1','Gpr160',
}

# ── Exclusion filters (hardcoded for traceability) ────────────────────

# FILTER 1: Neurotrophins — these are NOT neuropeptides. They signal
# through receptor tyrosine kinases (Ntrk1/2/3) and p75NTR (Ngfr),
# not GPCRs. Fundamentally different signaling class.
EXCLUDED_LIGANDS = {
    'Bdnf',   # Brain-derived neurotrophic factor → Ntrk2
    'Ngf',    # Nerve growth factor → Ntrk1
    'Ntf3',   # Neurotrophin-3 → Ntrk3
    'Ntf5',   # Neurotrophin-4/5 → Ntrk2
}

# FILTER 2: Kininogen — a blood plasma protein, not a neuropeptide
# precursor. Kng1 → Rxfp4 is a database artifact.
EXCLUDED_LIGANDS.add('Kng1')

# FILTER 3: RAMPs (Receptor Activity-Modifying Proteins) are NOT
# standalone receptors. They are single-pass accessory subunits that
# modify the pharmacology of Calcr and Calcrl. Entries like
# "Adcyap1 → Ramp2" are meaningless without specifying the GPCR.
EXCLUDED_RECEPTORS = {'Ramp1', 'Ramp2', 'Ramp3'}

# FILTER 4: Receptor tyrosine kinases and p75NTR — not GPCRs.
EXCLUDED_RECEPTORS.update({'Ntrk1', 'Ntrk2', 'Ntrk3', 'Ngfr'})

# FILTER 5: Specific artifact pairs — pharmacologically wrong pairings
# that propagated across multiple databases via copy-paste. Each one has
# been verified against IUPHAR and primary literature (March 2026).
EXCLUDED_PAIRS = {
    # Gene-name confusion: ghrelin does NOT bind GHRH receptor.
    # Ghrelin → GHSR; GHRH → GHRHR. Similar names, different systems.
    # Verified: IUPHAR GHRHR page lists only GHRH as ligand.
    ('Ghrl', 'Ghrhr'),

    # Class B GPCR cross-reactivity artifact: GHRH does NOT bind VIP
    # receptors. Experimentally disproven (Varga et al., PNAS 2000:
    # "VIP receptors do not recognize hGHRH and its analogs").
    ('Ghrh', 'Vipr1'),
    ('Ghrh', 'Vipr2'),

    # Co-expression confusion: NPY does NOT activate MC4R. AgRP (the
    # inverse agonist at MC4R) is co-expressed with NPY in ARH neurons.
    # NPY acts exclusively via NPY receptors (Y1/Y2/Y4/Y5).
    # Verified: IUPHAR MC4R page lists only POMC-derived peptides + AgRP.
    ('Npy', 'Mc4r'),

    # Different peptide families: Substance P (tachykinin) does NOT bind
    # GRPR (bombesin receptor). Co-expressed in spinal itch circuits but
    # signal through separate receptor systems (NK1R vs GRPR).
    ('Tac1', 'Grpr'),

    # Wrong gene: GIP receptor binds GIP (gene: Gip), not glucagon or
    # any other proglucagon product (gene: Gcg). GIPR achieves "exquisite
    # sequence specificity" (eLife 2021). No cross-binding documented.
    ('Gcg', 'Gipr'),

    # No primary literature support: PTHrP signals through PTH1R, not
    # the prolactin-releasing peptide receptor. Listed on IUPHAR PrRPR
    # page but no published binding data could be found in PubMed.
    ('Pthlh', 'Prlhr'),

    # Phylogenetic homology artifact: orexin receptors share ~30% sequence
    # identity with NPFFR2, but orexin peptides have no established binding
    # at NPFFR2. Not listed on IUPHAR NPFFR2 page. No binding data in
    # PubMed.
    ('Hcrt', 'Npffr2'),

    # Synthetic antagonist cross-reactivity, not peptide binding: NPY
    # itself doesn't bind NPFFR2. The confusion comes from NPY receptor
    # antagonists (BIBP3226, GR231118) having off-target NPFFR2 affinity
    # (Mollereau et al., Br J Pharmacol 2002).
    ('Npy', 'Npffr2'),

    # Ki >30 μM at GPR10 (Bonini et al., 2001) — >30,000× weaker than
    # the cognate ligand PrRP (Ki <1 nM). Listed on IUPHAR but reflects
    # structural homology to NPY receptors, not physiological binding.
    ('Npy', 'Prlhr'),

    # GPR83 is NOT an NPY/PYY/PP receptor. Gomes et al. 2016 explicitly
    # tested NPY: <20% displacement at 10 μM, not significant. GPR83
    # remains effectively orphan (PEN proposal contradicted by Giesecke
    # et al. 2023; FAM237A/B proposed by Li et al. 2023 but unconfirmed).
    ('Npy', 'Gpr83'),
    ('Ppy', 'Gpr83'),
    ('Pyy', 'Gpr83'),

    # ── Tier 3 exclusions: verified no-binding or wrong mechanism ──

    # NOP/ORL1 receptor (Oprl1) only binds nociceptin/orphanin FQ
    # (gene: Pnoc). Classical opioid peptides from Pdyn, Penk, Pomc
    # have no significant affinity. IUPHAR ORL1 page: only nociceptin.
    ('Pdyn', 'Oprl1'),
    ('Penk', 'Oprl1'),
    ('Pomc', 'Oprl1'),

    # UCN2 and UCN3 are CRHR2-selective. UCN2 has EC50 >100 nM at
    # CRHR1 (~1000× weaker than at CRHR2). UCN3 shows even greater
    # selectivity. IUPHAR CRHR1 page does not list UCN2/UCN3.
    ('Ucn2', 'Crhr1'),
    ('Ucn3', 'Crhr1'),

    # Bare Calcr without RAMPs does NOT bind adrenomedullin or
    # adrenomedullin 2. ADM requires Calcrl+RAMP2/3 (AM1/AM2 receptors).
    # We already have the correct heteromeric entries.
    ('Adm', 'Calcr'),
    ('Adm2', 'Calcr'),

    # Bare Calcrl (CLR) without RAMPs does NOT bind amylin. Amylin
    # requires Calcr+RAMP1/2/3 (AMY receptors). Calcrl+RAMPs form
    # CGRP and AM receptors, not amylin receptors.
    ('Iapp', 'Calcrl'),

    # TIP39 (Pth2) is an antagonist at PTH1R, not an agonist. It blocks
    # PTH signaling without activating the receptor. Not a signaling pair.
    # Verified: IUPHAR PTH1R page lists TIP39 as antagonist.
    ('Pth2', 'Pth1r'),

    # Relaxin family: no IUPHAR binding data for these cross-pairs.
    # Rln1 is RXFP1-selective; no documented affinity at RXFP3 or RXFP4.
    ('Rln1', 'Rxfp3'),
    ('Rln1', 'Rxfp4'),
    # INSL3 is RXFP2-selective; no documented affinity at RXFP3 or RXFP4.
    ('Insl3', 'Rxfp3'),
    ('Insl3', 'Rxfp4'),

    # ── Additional exclusions from systematic Tier 3 verification ──

    # Enkephalins (met-enk, leu-enk) are inactive at kappa opioid receptor.
    # IUPHAR OPRK1 page lists dynorphins and neoendorphins only.
    # StatPearls: met-enkephalin is "inactive at kappa receptors."
    ('Penk', 'Oprk1'),

    # Beta-endorphin has only unquantifiably low kappa affinity.
    # IUPHAR: OPRK1 "has low affinity for beta-endorphins." No Ki reported.
    # Not a physiologically relevant interaction.
    ('Pomc', 'Oprk1'),

    # INSL5 does not interact with RXFP1 or RXFP2.
    # IUPHAR-confirmed: INSL5 is RXFP4-selective with no RXFP1/RXFP2 binding.
    # (Lok et al., Nature Sci Rep 2016; IUPHAR GtoPdb v.2023.1)
    ('Insl5', 'Rxfp1'),
    ('Insl5', 'Rxfp2'),

    # PACAP (Adcyap1) is NOT listed on IUPHAR secretin receptor page.
    # 2500-fold lower affinity than secretin — no functional binding.
    ('Adcyap1', 'Sctr'),

    # PP (Ppy) does not meaningfully bind Y1 or Y2 receptors.
    # IUPHAR: PP is Y4-selective; Y1 and Y2 binding essentially absent.
    ('Ppy', 'Npy1r'),
    ('Ppy', 'Npy2r'),
}

# ── Obligate heteromeric receptor complexes ───────────────────────────
# Some GPCR complexes require RAMP subunits to form functional receptors.
# Lewis Lab databases list the GPCR and RAMP as separate "receptor" entries.
# We merge them into the pharmacologically correct heteromer notation.
#
# Nomenclature follows IUPHAR (Alexander et al., 2023):
#   Calcrl + Ramp1 = CGRP receptor        (ligands: Calca→CGRP, Calcb→CGRP-β)
#   Calcrl + Ramp2 = AM1 receptor          (ligands: Adm, Adm2)
#   Calcrl + Ramp3 = AM2 receptor          (ligands: Adm, Adm2)
#   Calcr  + Ramp1 = AMY1 receptor         (ligands: Iapp, Calca, Calcb)
#   Calcr  + Ramp2 = AMY2 receptor         (ligands: Iapp, Calca)
#   Calcr  + Ramp3 = AMY3 receptor         (ligands: Iapp, Calca)
#
# The bare Calcr (calcitonin receptor) also signals without a RAMP,
# but for CGRP/adrenomedullin/amylin family peptides, the relevant
# signaling complex is always the heteromer.
#
# Format: (Ligand_Gene, GPCR_component) → list of heteromer Receptor_Genes
# If a ligand has a Lewis Lab entry for both Calcrl and any RAMP, emit
# the heteromer. The GPCR-only and RAMP-only raw rows are consumed.

# Format: (Ligand_Gene, GPCR_component) → (heteromer_list, keep_bare)
# keep_bare=True means also emit the bare GPCR pair (e.g. Calcr alone
# is a functional calcitonin receptor). keep_bare=False means the bare
# pair is pharmacologically meaningless without the RAMP (e.g. Calcrl
# alone doesn't bind CGRP).
OBLIGATE_HETEROMERS = {
    # CGRP receptor: Calcrl + Ramp1.  Bare Calcrl has no known ligand.
    ('Calca', 'Calcrl'): (['Calcrl;Ramp1'], False),
    ('Calcb', 'Calcrl'): (['Calcrl;Ramp1'], False),
    # AM1/AM2 receptors: Calcrl + Ramp2/Ramp3
    ('Adm',  'Calcrl'):  (['Calcrl;Ramp2', 'Calcrl;Ramp3'], False),
    ('Adm2', 'Calcrl'):  (['Calcrl;Ramp2', 'Calcrl;Ramp3'], False),
    # Amylin receptors: Calcr + Ramp1/2/3.  Bare Calcr IS the
    # calcitonin receptor and also binds amylin with lower affinity.
    ('Iapp', 'Calcr'):   (['Calcr;Ramp1', 'Calcr;Ramp2', 'Calcr;Ramp3'], True),
    # Calcitonin: bare Calcr is the primary calcitonin receptor.
    ('Calca', 'Calcr'):  (['Calcr;Ramp1', 'Calcr;Ramp2', 'Calcr;Ramp3'], True),
    # CGRP-β: requires Calcr+Ramp1 (AMY1). Bare Calcr does NOT bind
    # CGRP-β at physiological concentrations.
    ('Calcb', 'Calcr'):  (['Calcr;Ramp1'], False),
}

# ── Apply filters and heteromer mapping ───────────────────────────────

# Step 1: Basic NP ligand × NP receptor filter (excluding junk)
filtered_pairs = {}
for (lig, rec), dbs in non_dimitrov_only.items():
    if lig in EXCLUDED_LIGANDS:
        continue
    if rec in EXCLUDED_RECEPTORS:
        continue
    if (lig, rec) in EXCLUDED_PAIRS:
        continue
    if lig in all_np_ligands and rec in all_np_receptors:
        filtered_pairs[(lig, rec)] = dbs

print(f"NPP/GPCR pairs after filtering: {len(filtered_pairs)}")

# Step 2: Apply heteromer mapping
lewislab_pairs = {}
consumed = set()  # track raw pairs that become heteromers

for (lig, rec), dbs in filtered_pairs.items():
    if (lig, rec) in OBLIGATE_HETEROMERS:
        het_list, keep_bare = OBLIGATE_HETEROMERS[(lig, rec)]
        # Emit heteromeric receptor entries
        for het_rec in het_list:
            key = (lig, het_rec)
            if key not in lewislab_pairs:
                lewislab_pairs[key] = list(dbs)
            else:
                # Merge database lists
                for d in dbs:
                    if d not in lewislab_pairs[key]:
                        lewislab_pairs[key].append(d)
        # Also keep the bare GPCR pair if it's a functional receptor
        if keep_bare:
            lewislab_pairs[(lig, rec)] = dbs
        else:
            consumed.add((lig, rec))
    else:
        lewislab_pairs[(lig, rec)] = dbs

n_het = sum(1 for k in lewislab_pairs if ';' in k[1])
print(f"After heteromer mapping: {len(lewislab_pairs)} pairs ({n_het} heteromeric)")
print(f"Raw pairs consumed by heteromer mapping: {len(consumed)}")

# ── Primary vs. secondary receptor classification ─────────────────────
# For each ligand gene, define which receptors are "primary" targets
# (IUPHAR-endorsed canonical receptors for at least one bioactive
# peptide from that gene). Everything else is "secondary" (real
# pharmacological cross-reactivity but not the principal target).
#
# Ligands with only ONE receptor in the final set are automatically
# primary — no entry needed here.
#
# Convention: if a precursor gene produces multiple peptides with
# different primary receptors (e.g. Tac1 → SP→NK1 + NKA→NK2), ALL
# those receptors are primary for the gene.

PRIMARY_RECEPTORS = {
    # ── Opioid system ──
    # Pdyn → dynorphins: kappa-selective. Cross-react at mu and delta.
    'Pdyn':    {'Oprk1'},
    # Penk → enkephalins: delta-selective. Met-enk has some mu affinity.
    'Penk':    {'Oprd1', 'Oprm1'},
    # Pomc → beta-endorphin (mu-preferring) + melanocortins (MC3R/MC4R).
    # MC5R is peripheral (sebaceous glands), not a primary hypothalamic target.
    'Pomc':    {'Oprm1', 'Mc3r', 'Mc4r'},

    # ── Tachykinin system ──
    # Tac1 → SP (NK1-selective) + NKA (NK2-selective). NK3 is secondary.
    'Tac1':    {'Tacr1', 'Tacr2'},
    # Tac2 → NKB: NK3-selective. NK1/NK2 are secondary.
    'Tac2':    {'Tacr3'},

    # ── Vasopressin / oxytocin ──
    # AVP → all three AVP receptors are primary. OXT receptor is secondary.
    'Avp':     {'Avpr1a', 'Avpr1b', 'Avpr2'},
    # OXT → OXT receptor is primary. AVP receptors are secondary.
    'Oxt':     {'Oxtr'},

    # ── NPY family ──
    # NPY: Y1, Y2, Y5 are principal CNS receptors. Y4 is PP-preferring.
    # Y6 (Npy6r) is functional in mouse but pseudogene in human.
    'Npy':     {'Npy1r', 'Npy2r', 'Npy5r'},
    # PYY: Y1, Y2 (especially PYY3-36 → Y2). Y4/Y5 secondary.
    'Pyy':     {'Npy1r', 'Npy2r'},
    # PP: Y4-selective; Y6 (Npy6r) is also primary in mouse (functional
    # receptor, PP is endogenous ligand — Cell Metabolism 2013). Y1/Y2/Y5 secondary.
    'Ppy':     {'Npy4r', 'Npy6r'},

    # ── Bombesin family ──
    # GRP → GRPR primary. NMBR and BRS3 are secondary.
    'Grp':     {'Grpr'},
    # NMB → NMBR primary. GRPR and BRS3 are secondary.
    'Nmb':     {'Nmbr'},

    # ── CRH / urocortin ──
    # CRH: CRHR1 primary (10× preference). CRHR2 secondary.
    'Crh':     {'Crhr1'},
    # UCN: binds both CRHR1 and CRHR2 equally — both primary.
    'Ucn':     {'Crhr1', 'Crhr2'},

    # ── Melanocortin (AgRP) ──
    # AgRP: MC3R and MC4R primary (inverse agonist). MC5R secondary.
    'Agrp':    {'Mc3r', 'Mc4r'},

    # ── VIP / PACAP / secretin cross-talk ──
    # PACAP: PAC1, VPAC1, VPAC2 all primary (nanomolar affinity).
    # Secretin receptor is secondary (2500× weaker).
    'Adcyap1': {'Adcyap1r1', 'Vipr1', 'Vipr2'},
    # VIP: VPAC1 and VPAC2 primary. PAC1 secondary (100-1000× weaker).
    # Secretin receptor secondary.
    'Vip':     {'Vipr1', 'Vipr2'},
    # Secretin: SCTR primary. VPAC1/VPAC2 secondary (1000-10,000× weaker).
    'Sct':     {'Sctr'},

    # ── GHRH / ghrelin ──
    # GHRH: GHRHR primary. GHSR allosteric co-agonist (secondary).
    'Ghrh':    {'Ghrhr'},

    # ── PTH family ──
    # PTH: PTH1R primary. PTH2R secondary (partial agonist, species-dependent).
    'Pth':     {'Pth1r'},
    # PTHrP: PTH1R primary. PTH2R secondary.
    'Pthlh':   {'Pth1r'},

    # ── Relaxin family ──
    # Relaxin-1: RXFP1 primary. RXFP2 secondary (pKi 8.8 — high but not cognate).
    'Rln1':    {'Rxfp1'},
    # Relaxin-3: RXFP3 primary. RXFP1 (pKi 7.5-8.0), RXFP2 (pKi 7.0),
    # RXFP4 (pKi 8.8-9.0) are all secondary despite high affinities.
    'Rln3':    {'Rxfp3'},
    # INSL3: RXFP2 primary. RXFP1 weak secondary (pKi 5.7).
    'Insl3':   {'Rxfp2'},
    # INSL5: RXFP4 primary. RXFP3 secondary (antagonist, pKi 7.0).
    # RXFP1/RXFP2 weak secondary.
    'Insl5':   {'Rxfp4'},

    # ── Calcitonin / CGRP / amylin / adrenomedullin ──
    # Calca → CGRP (Calcrl;Ramp1 primary) + calcitonin (Calcr primary).
    # Calcitonin also signals via AMY1/2/3 at lower affinity — secondary.
    'Calca':   {'Calcrl;Ramp1', 'Calcr'},
    # Calcb → CGRP-β: Calcrl;Ramp1 primary. Calcr;Ramp1 secondary.
    'Calcb':   {'Calcrl;Ramp1'},
    # Adm: AM1 (Calcrl;Ramp2) and AM2 (Calcrl;Ramp3) both primary.
    'Adm':     {'Calcrl;Ramp2', 'Calcrl;Ramp3'},
    # Adm2/intermedin: same receptors.
    'Adm2':    {'Calcrl;Ramp2', 'Calcrl;Ramp3'},
    # Amylin: AMY1 (Calcr;Ramp1) is highest-affinity primary.
    # AMY2 and AMY3 are also primary (functional amylin receptors).
    # Bare Calcr is secondary (lower affinity without RAMP).
    'Iapp':    {'Calcr;Ramp1', 'Calcr;Ramp2', 'Calcr;Ramp3'},

    # ── Cortistatin ──
    # Binds all 5 SSTRs with high affinity (primary, like somatostatin).
    # GHSR binding is secondary (weaker, ~10× less than ghrelin).
    'Cort':    {'Sstr1', 'Sstr2', 'Sstr3', 'Sstr4', 'Sstr5'},

    # ── UTS2B ──
    # UTS2R is primary. SSTR5 is secondary (real but not cognate).
    'Uts2b':   {'Uts2r'},

    # ── Endothelins ──
    # ET-1 and ET-2 bind both ETA and ETB with similar affinity — both primary.
    # ET-3 is ETB-selective.
    'Edn1':    {'Ednra', 'Ednrb'},
    'Edn2':    {'Ednra', 'Ednrb'},
    'Edn3':    {'Ednrb'},

    # ── Proglucagon ──
    # Gcg produces glucagon (→Gcgr), GLP-1 (→Glp1r), GLP-2 (→Glp2r).
    # All primary for their respective peptides.
    'Gcg':     {'Gcgr', 'Glp1r', 'Glp2r'},
}

def get_binding_rank(lig, rec):
    """Determine primary/secondary status for a ligand-receptor pair."""
    if lig not in PRIMARY_RECEPTORS:
        # Ligand not in mapping → all its receptors are primary
        # (handles single-receptor ligands and those we haven't classified)
        return 'primary'
    if rec in PRIMARY_RECEPTORS[lig]:
        return 'primary'
    return 'secondary'

# Write np_map_lewislab.csv
with open(f'{OUT}/np_map_lewislab.csv', 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['Ligand_Gene','Receptor_Gene','databases','n_databases','binding_rank'])
    for (lig, rec), dbs in sorted(lewislab_pairs.items()):
        rank = get_binding_rank(lig, rec)
        w.writerow([lig, rec, ';'.join(dbs), len(dbs), rank])

n_primary = sum(1 for (l,r) in lewislab_pairs if get_binding_rank(l,r) == 'primary')
n_secondary = len(lewislab_pairs) - n_primary
print(f"Wrote np_map_lewislab.csv: {len(lewislab_pairs)} rows "
      f"({n_primary} primary, {n_secondary} secondary)")

# ════════════════════════════════════════════════════════════════════════
# PART B: Validate extra_pairs.csv (CURATED — never overwritten)
# ════════════════════════════════════════════════════════════════════════
# extra_pairs.csv contains gene pairs absent from all Lewis Lab databases.
# It is manually maintained. This script reads and validates it.

ep_path = f'{CURATED}/extra_pairs.csv'
if not os.path.exists(ep_path):
    raise FileNotFoundError(
        f"{ep_path} not found. This is a curated file that must exist.")

extra_gene_pairs = {}
with open(ep_path) as f:
    for row in csv.DictReader(f):
        key = (row['Ligand_Gene'], row['Receptor_Gene'])
        extra_gene_pairs[key] = row.get('References', '')

# Warn if any extra pair is actually in Lewis Lab (would be redundant)
redundant = [k for k in extra_gene_pairs if k in lewislab_pairs]
print(f"Validated extra_pairs.csv: {len(extra_gene_pairs)} gene pairs")
if redundant:
    print(f"  WARNING: {len(redundant)} pairs also in Lewis Lab (redundant):")
    for lg, rg in redundant:
        print(f"    {lg} -> {rg}")

# ════════════════════════════════════════════════════════════════════════
# PART C: Validate extra_info.csv (CURATED — never overwritten)
# ════════════════════════════════════════════════════════════════════════
# extra_info.csv is manually maintained. This script only VALIDATES it:
# - Every row must reference a valid gene pair (Lewis Lab or extra_pairs)
# - Reports unannotated pairs that need attention

valid_gene_pairs = set(lewislab_pairs.keys()) | set(extra_gene_pairs.keys())

ei_path = f'{CURATED}/extra_info.csv'
if not os.path.exists(ei_path):
    raise FileNotFoundError(
        f"{ei_path} not found. This is a curated file that must exist.")

ei_rows = []
ei_pairs = set()
ei_invalid = []
with open(ei_path) as f:
    for r in csv.DictReader(f):
        key = (r['Ligand_Gene'], r['Receptor_Gene'])
        ei_rows.append(r)
        ei_pairs.add(key)
        if key not in valid_gene_pairs:
            ei_invalid.append(key)

print(f"Validated extra_info.csv: {len(ei_rows)} rows")
if ei_invalid:
    print(f"  WARNING: {len(ei_invalid)} rows reference invalid gene pairs:")
    for lg, rg in ei_invalid:
        print(f"    {lg} -> {rg}")
unannotated = valid_gene_pairs - ei_pairs
print(f"  Unannotated gene pairs: {len(unannotated)}")

# ════════════════════════════════════════════════════════════════════════
# PART D: Validate processing_enzymes.csv (CURATED — never overwritten)
# ════════════════════════════════════════════════════════════════════════

pe_path = f'{CURATED}/processing_enzymes.csv'
if not os.path.exists(pe_path):
    raise FileNotFoundError(
        f"{pe_path} not found. This is a curated file that must exist.")

pe_rows = []
with open(pe_path) as f:
    pe_rows = list(csv.DictReader(f))
print(f"Validated processing_enzymes.csv: {len(pe_rows)} rows")

# ════════════════════════════════════════════════════════════════════════
# PART E: Final summary
# ════════════════════════════════════════════════════════════════════════

# Re-read extra_info for final counts
ei_pairs_final = set()
ei_count_final = 0
with open(f'{CURATED}/extra_info.csv') as f:
    for r in csv.DictReader(f):
        ei_pairs_final.add((r['Ligand_Gene'], r['Receptor_Gene']))
        ei_count_final += 1

unannotated_final = valid_gene_pairs - ei_pairs_final

print(f"\n=== FINAL SUMMARY ===")
print(f"np_map_lewislab.csv:    {len(lewislab_pairs):>4d} gene pairs (mechanically derived)")
print(f"extra_pairs.csv:        {len(extra_gene_pairs):>4d} gene pairs (manually curated)")
print(f"extra_info.csv:         {ei_count_final:>4d} rows (annotation layer — CURATED)")
print(f"processing_enzymes.csv: {len(pe_rows):>4d} rows (enzyme disambiguation — CURATED)")
print(f"Total gene pairs:       {len(valid_gene_pairs):>4d}")
print(f"Annotated pairs:        {len(ei_pairs_final):>4d}")
print(f"Unannotated pairs:      {len(unannotated_final):>4d}")
if unannotated_final:
    for lg, rg in sorted(unannotated_final):
        print(f"  {lg} -> {rg}")
