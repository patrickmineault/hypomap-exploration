"""Build ligand-receptor list from NeuronChat database + manual additions.

Downloads the NeuronChat mouse interaction database and adds manually curated
neuropeptide pairs that are missing from the database (e.g., Oxt/Oxtr, Avp/Avpr).

Usage:
    python -m src.preprocessing.build_ligand_receptor_list

Output:
    data/processed/mouse_common/ligand_receptor_mouse.csv
"""

import pandas as pd
from pathlib import Path
import urllib.request

# Paths
DATA_DIR = Path(__file__).parent.parent.parent / "data"
OUTPUT_DIR = DATA_DIR / "processed" / "mouse_common"
OUTPUT_FILE = OUTPUT_DIR / "ligand_receptor_mouse.csv"

# NeuronChat database URL
NEURONCHAT_URL = "https://raw.githubusercontent.com/Wei-BioMath/NeuronChatAnalysis2022/main/NeuronChatDB_table/interactionDB_mouse.txt"

# Manual additions: neuropeptide pairs missing from NeuronChat
# Format: (ligand_name, receptor_name, ligand_type, interaction_type, gene_ligand, gene_receptor)
MANUAL_NEUROPEPTIDES = [
    # Oxytocin signaling - critical for social behavior, PVH/SO neurons
    ("Oxt", "Oxtr", "Neuropeptide", "ligand-receptor", "Oxt", "Oxtr"),

    # Vasopressin signaling - critical for stress response, PVH/SO neurons
    ("Avp", "Avpr1a", "Neuropeptide", "ligand-receptor", "Avp", "Avpr1a"),
    ("Avp", "Avpr1b", "Neuropeptide", "ligand-receptor", "Avp", "Avpr1b"),
    ("Avp", "Avpr2", "Neuropeptide", "ligand-receptor", "Avp", "Avpr2"),

    # AgRP - critical for feeding, ARH neurons (if missing from NeuronChat)
    # Note: AgRP signals through melanocortin receptors by antagonizing alpha-MSH
    ("Agrp", "Mc3r", "Neuropeptide", "ligand-receptor", "Agrp", "Mc3r"),
    ("Agrp", "Mc4r", "Neuropeptide", "ligand-receptor", "Agrp", "Mc4r"),

    # POMC-derived peptides (alpha-MSH) - satiety signaling
    ("Pomc", "Mc3r", "Neuropeptide", "ligand-receptor", "Pomc", "Mc3r"),
    ("Pomc", "Mc4r", "Neuropeptide", "ligand-receptor", "Pomc", "Mc4r"),

    # Orexin/Hypocretin - arousal, LHA neurons
    ("Hcrt", "Hcrtr1", "Neuropeptide", "ligand-receptor", "Hcrt", "Hcrtr1"),
    ("Hcrt", "Hcrtr2", "Neuropeptide", "ligand-receptor", "Hcrt", "Hcrtr2"),

    # MCH - sleep/feeding, LHA neurons
    ("Pmch", "Mchr1", "Neuropeptide", "ligand-receptor", "Pmch", "Mchr1"),
    ("Pmch", "Mchr2", "Neuropeptide", "ligand-receptor", "Pmch", "Mchr2"),

    # Kisspeptin - reproduction, AVPV/ARH neurons
    ("Kiss1", "Kiss1r", "Neuropeptide", "ligand-receptor", "Kiss1", "Kiss1r"),

    # GHRH - growth hormone releasing hormone
    ("Ghrh", "Ghrhr", "Neuropeptide", "ligand-receptor", "Ghrh", "Ghrhr"),

    # Galanin - feeding, arousal
    ("Gal", "Galr1", "Neuropeptide", "ligand-receptor", "Gal", "Galr1"),
    ("Gal", "Galr2", "Neuropeptide", "ligand-receptor", "Gal", "Galr2"),
    ("Gal", "Galr3", "Neuropeptide", "ligand-receptor", "Gal", "Galr3"),

    # CART - feeding regulation (receptor unclear, included as ligand)
    ("Cartpt", "Unknown", "Neuropeptide", "ligand-receptor", "Cartpt", "Cartpt"),

    # Neuromedin S - circadian, feeding
    ("Nms", "Nmur1", "Neuropeptide", "ligand-receptor", "Nms", "Nmur1"),
    ("Nms", "Nmur2", "Neuropeptide", "ligand-receptor", "Nms", "Nmur2"),

    # Prokineticin 2 - circadian output from SCH
    ("Prok2", "Prokr1", "Neuropeptide", "ligand-receptor", "Prok2", "Prokr1"),
    ("Prok2", "Prokr2", "Neuropeptide", "ligand-receptor", "Prok2", "Prokr2"),

    # BDNF - neurotrophin
    ("Bdnf", "Ntrk2", "Neuropeptide", "ligand-receptor", "Bdnf", "Ntrk2"),
]

# Hormone receptors (ligands are circulating, receptors in hypothalamus)
MANUAL_HORMONE_RECEPTORS = [
    # Leptin receptor - metabolic sensing in ARH
    ("Lep", "Lepr", "Hormone", "ligand-receptor", "Lep", "Lepr"),

    # Ghrelin receptor - hunger signal
    ("Ghrl", "Ghsr", "Hormone", "ligand-receptor", "Ghrl", "Ghsr"),

    # Insulin receptor - metabolic sensing
    ("Ins", "Insr", "Hormone", "ligand-receptor", "Ins", "Insr"),

    # Glucocorticoid receptor - stress response (ligand is corticosterone)
    ("Cort", "Nr3c1", "Hormone", "ligand-receptor", "Cort", "Nr3c1"),

    # Estrogen receptor - Cyp19a1 (aromatase) synthesizes estradiol
    ("Estrogen", "Esr1", "Hormone", "ligand-receptor", "Cyp19a1", "Esr1"),

    # Melatonin receptor - Aanat is rate-limiting enzyme in melatonin synthesis
    ("Melatonin", "Mtnr1a", "Hormone", "ligand-receptor", "Aanat", "Mtnr1a"),
    ("Melatonin", "Mtnr1b", "Hormone", "ligand-receptor", "Aanat", "Mtnr1b"),

    # Endocannabinoid receptors - Dagla makes 2-AG, Napepld makes anandamide
    ("2-AG", "Cnr1", "Hormone", "ligand-receptor", "Dagla", "Cnr1"),
    ("Anandamide", "Cnr1", "Hormone", "ligand-receptor", "Napepld", "Cnr1"),
]

# Histamine - not in NeuronChat, Hdc marks histaminergic neurons (TMN)
MANUAL_HISTAMINE = [
    ("Histamine", "Hrh1", "Amine", "ligand-receptor", "Hdc", "Hrh1"),
    ("Histamine", "Hrh2", "Amine", "ligand-receptor", "Hdc", "Hrh2"),
    ("Histamine", "Hrh3", "Amine", "ligand-receptor", "Hdc", "Hrh3"),
]


def download_neuronchat_db() -> pd.DataFrame:
    """Download and parse NeuronChat mouse interaction database."""
    print(f"Downloading NeuronChat database from:\n  {NEURONCHAT_URL}")

    # Download the file
    with urllib.request.urlopen(NEURONCHAT_URL) as response:
        content = response.read().decode('utf-8')

    # Parse tab-separated file
    lines = content.strip().split('\n')
    header = lines[0].split()

    data = []
    for line in lines[1:]:
        parts = line.split()
        if len(parts) >= 4:
            # Extract interaction_name (e.g., "Vip_Vipr1") and parse ligand/receptor
            interaction_name = parts[1] if parts[0].isdigit() else parts[0]
            ligand_type = parts[2] if parts[0].isdigit() else parts[1]
            interaction_type = parts[3] if parts[0].isdigit() else parts[2]
            lig_contributor = parts[4] if parts[0].isdigit() else parts[3]

            # Parse ligand and receptor from interaction name
            if '_' in interaction_name:
                ligand_name, receptor_name = interaction_name.split('_', 1)
            else:
                continue

            # Get gene names (handle complex multi-gene ligands)
            gene_ligand = lig_contributor
            gene_receptor = receptor_name

            data.append({
                'ligand_name': ligand_name,
                'receptor_name': receptor_name,
                'ligand_type': ligand_type,
                'interaction_type': interaction_type,
                'gene_ligand': gene_ligand,
                'gene_receptor': gene_receptor,
            })

    df = pd.DataFrame(data)
    print(f"  Loaded {len(df)} interactions from NeuronChat")
    return df


def add_manual_entries(df: pd.DataFrame) -> pd.DataFrame:
    """Add manually curated ligand-receptor pairs."""
    columns = ['ligand_name', 'receptor_name', 'ligand_type',
               'interaction_type', 'gene_ligand', 'gene_receptor']

    # Combine all manual entries
    all_manual = (
        MANUAL_NEUROPEPTIDES +
        MANUAL_HORMONE_RECEPTORS +
        MANUAL_HISTAMINE
    )
    manual_df = pd.DataFrame(all_manual, columns=columns)

    # Check which manual entries are already in the database
    existing = set(zip(df['gene_ligand'], df['gene_receptor']))
    new_entries = []
    skipped = []

    for _, row in manual_df.iterrows():
        key = (row['gene_ligand'], row['gene_receptor'])
        if key not in existing:
            new_entries.append(row)
        else:
            skipped.append(f"{row['gene_ligand']}-{row['gene_receptor']}")

    if skipped:
        print(f"  Skipped {len(skipped)} entries (already in NeuronChat): {skipped[:5]}...")

    if new_entries:
        new_df = pd.DataFrame(new_entries)
        df = pd.concat([df, new_df], ignore_index=True)
        print(f"  Added {len(new_entries)} manual entries")

    return df


def main():
    """Build the ligand-receptor list."""
    print("=== Building Ligand-Receptor List ===\n")

    # 1. Download NeuronChat database
    df = download_neuronchat_db()

    # 2. Add manual entries
    print("\nAdding manual neuropeptide entries...")
    df = add_manual_entries(df)

    # 3. Summary by type
    print("\nSummary by ligand type:")
    for ltype, count in df['ligand_type'].value_counts().items():
        print(f"  {ltype}: {count}")

    # 4. Save
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUTPUT_FILE, index=False)
    print(f"\nSaved {len(df)} interactions to {OUTPUT_FILE}")

    # 5. List neuropeptides
    neuropeptides = df[df['ligand_type'] == 'Neuropeptide']
    ligands = sorted(neuropeptides['gene_ligand'].unique())
    receptors = sorted(neuropeptides['gene_receptor'].unique())
    print(f"\nNeuropeptide ligands ({len(ligands)}): {ligands}")
    print(f"Neuropeptide receptors ({len(receptors)}): {receptors}")


if __name__ == "__main__":
    main()
