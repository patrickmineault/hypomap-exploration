"""Extract hypothalamus connectivity matrix from Allen Mouse Brain Connectivity Atlas.

This script downloads connectivity data from the Allen Brain Institute API and extracts
a connectivity matrix for hypothalamic regions.

Requires: pandas, requests

Reference: Oh SW, et al. (2014). A mesoscale connectome of the mouse brain. Nature 508, 207-214.
"""

from pathlib import Path
from typing import Dict, List

import pandas as pd
import requests

# Output paths
DATA_DIR = Path(__file__).parent.parent.parent / "data"
PROCESSED_DIR = DATA_DIR / "processed" / "mouse_common"
OUTPUT_FILE = PROCESSED_DIR / "hypothalamus_connectivity.csv"

# Allen Brain API endpoints
BASE_URL = "http://api.brain-map.org/api/v2"
STRUCTURE_GRAPH_ID = 1  # Adult mouse structure graph

# Hypothalamus region acronyms (Allen CCF)
HYPOTHALAMUS_ACRONYMS = [
    # Major hypothalamic nuclei
    "ARH",  # Arcuate hypothalamic nucleus
    "VMH",  # Ventromedial hypothalamic nucleus
    "DMH",  # Dorsomedial nucleus of hypothalamus
    "PVH",  # Paraventricular hypothalamic nucleus
    "LHA",  # Lateral hypothalamic area
    "SCH",  # Suprachiasmatic nucleus
    "SO",  # Supraoptic nucleus
    "AHN",  # Anterior hypothalamic nucleus
    "PH",  # Posterior hypothalamic nucleus
    "TU",  # Tuberal nucleus
    "ME",  # Median eminence
    "ZI",  # Zona incerta (technically subthalamic but connected)
    # Preoptic regions
    "MPO",  # Medial preoptic area
    "MPN",  # Medial preoptic nucleus
    "LPO",  # Lateral preoptic area
    "AVPV",  # Anteroventral periventricular nucleus
    "MEPO",  # Median preoptic nucleus
    # Mammillary regions
    "MM",  # Medial mammillary nucleus
    "LM",  # Lateral mammillary nucleus
    "SUM",  # Supramammillary nucleus
    "TM",  # Tuberomammillary nucleus
    "PMd",  # Dorsal premammillary nucleus
    "PMv",  # Ventral premammillary nucleus
    # Periventricular regions
    "PVa",  # Periventricular hypothalamic nucleus, anterior
    "PVi",  # Periventricular hypothalamic nucleus, intermediate
    "PVpo",  # Periventricular hypothalamic nucleus, preoptic
]


def get_structures() -> pd.DataFrame:
    """Fetch structure information from Allen API."""
    url = f"{BASE_URL}/data/query.json"
    params = {
        "criteria": f"model::Structure,rma::criteria,[graph_id$eq{STRUCTURE_GRAPH_ID}],rma::options[num_rows$eqall]"
    }

    print("Fetching structure ontology...")
    response = requests.get(url, params=params)
    response.raise_for_status()

    data = response.json()
    structures = pd.DataFrame(data["msg"])
    print(f"  Found {len(structures)} structures")

    return structures


def get_hypothalamus_structure_ids(structures: pd.DataFrame) -> Dict[str, int]:
    """Get structure IDs for hypothalamic regions."""
    hypo_structures = structures[structures["acronym"].isin(HYPOTHALAMUS_ACRONYMS)]

    # Also try to get child structures of HY (hypothalamus)
    hy_id = structures[structures["acronym"] == "HY"]["id"].values
    if len(hy_id) > 0:
        hy_id = hy_id[0]
        # Get all structures with HY as parent
        hy_children = structures[
            structures["structure_id_path"].str.contains(f"/{hy_id}/")
        ]
        hypo_structures = pd.concat([hypo_structures, hy_children]).drop_duplicates(
            subset="id"
        )

    print(f"  Found {len(hypo_structures)} hypothalamic structures")

    return dict(zip(hypo_structures["acronym"], hypo_structures["id"]))


def get_experiments() -> pd.DataFrame:
    """Fetch connectivity experiments from Mouse Connectivity Atlas."""
    url = f"{BASE_URL}/data/query.json"

    # Get all connectivity experiments (Product ID 5 = Mouse Connectivity)
    params = {
        "criteria": "model::SectionDataSet,rma::criteria,products[id$eq5],rma::options[num_rows$eq500]"
    }

    print("Fetching connectivity experiments...")
    response = requests.get(url, params=params)
    response.raise_for_status()

    data = response.json()
    if not data["msg"]:
        return pd.DataFrame()

    experiments = pd.DataFrame(data["msg"])
    print(f"  Found {len(experiments)} experiments")

    return experiments


def get_projection_unionizes(
    experiment_ids: List[int], structure_ids: List[int]
) -> pd.DataFrame:
    """Fetch projection density data for experiments into structures."""
    all_unionizes = []

    print(f"Fetching projection data for {len(experiment_ids)} experiments...")

    for i, exp_id in enumerate(experiment_ids[:50]):  # Limit to first 50 for speed
        if (i + 1) % 10 == 0:
            print(f"  Progress: {i + 1}/{min(50, len(experiment_ids))}")

        url = f"{BASE_URL}/data/query.json"
        structure_id_str = ",".join(map(str, structure_ids))
        params = {
            "criteria": f"model::ProjectionStructureUnionize,rma::criteria,[section_data_set_id$eq{exp_id}],[structure_id$in{structure_id_str}],rma::options[num_rows$eqall]"
        }

        try:
            response = requests.get(url, params=params)
            response.raise_for_status()
            data = response.json()
            if data["msg"]:
                unionizes = pd.DataFrame(data["msg"])
                unionizes["experiment_id"] = exp_id
                all_unionizes.append(unionizes)
        except Exception as e:
            print(f"    Warning: Failed to fetch data for experiment {exp_id}: {e}")
            continue

    if all_unionizes:
        return pd.concat(all_unionizes, ignore_index=True)
    return pd.DataFrame()


def build_connectivity_matrix(
    unionizes: pd.DataFrame, id_to_acronym: Dict[int, str]
) -> pd.DataFrame:
    """Build a structure-to-structure connectivity matrix from projection data."""

    # Pivot: rows = experiments/injection sites, columns = target structures
    pivot = unionizes.pivot_table(
        index="experiment_id",
        columns="structure_id",
        values="projection_density",
        aggfunc="mean",
    )

    # Rename columns to acronyms
    pivot.columns = [id_to_acronym.get(c, str(c)) for c in pivot.columns]

    # Compute mean projection density from each injection site
    # This gives a coarse connectivity profile
    summary = pivot.describe().T
    summary["mean_projection"] = pivot.mean()

    return pivot, summary


def main():
    """Main function to extract hypothalamus connectivity."""
    print("=" * 60)
    print("Extracting Hypothalamus Connectivity from Allen Brain Atlas")
    print("=" * 60)
    print()

    # Step 1: Get structure ontology
    structures = get_structures()

    # Step 2: Get hypothalamus structure IDs
    hypo_id_map = get_hypothalamus_structure_ids(structures)

    if not hypo_id_map:
        print("Error: No hypothalamic structures found!")
        return

    print(f"\nHypothalamic regions found: {list(hypo_id_map.keys())}")

    hypo_ids = list(hypo_id_map.values())
    acronym_map = {v: k for k, v in hypo_id_map.items()}

    # Step 3: Get connectivity experiments
    experiments = get_experiments()

    if len(experiments) == 0:
        print("Warning: No experiments found. Creating structure list only.")

        # At minimum, save the structure information
        PROCESSED_DIR.mkdir(parents=True, exist_ok=True)
        structures_df = pd.DataFrame(
            [
                {
                    "acronym": k,
                    "structure_id": v,
                    "name": (
                        structures[structures["id"] == v]["name"].values[0]
                        if len(structures[structures["id"] == v]) > 0
                        else k
                    ),
                }
                for k, v in hypo_id_map.items()
            ]
        )
        structures_file = PROCESSED_DIR / "hypothalamus_structures.csv"
        structures_df.to_csv(structures_file, index=False)
        print(f"\nSaved structure list to: {structures_file}")
        return

    # Step 4: Get projection data
    unionizes = get_projection_unionizes(experiments["id"].tolist(), hypo_ids)

    if len(unionizes) == 0:
        print("Warning: No projection data retrieved.")
        return

    print(f"\nRetrieved {len(unionizes)} projection records")

    # Step 5: Build connectivity matrix
    connectivity_df, summary_df = build_connectivity_matrix(unionizes, acronym_map)

    # Ensure output directory exists
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

    # Save full data
    connectivity_df.to_csv(OUTPUT_FILE)
    print(f"\nSaved connectivity matrix to: {OUTPUT_FILE}")
    print(f"  Shape: {connectivity_df.shape}")

    # Save summary
    summary_file = PROCESSED_DIR / "hypothalamus_connectivity_summary.csv"
    summary_df.to_csv(summary_file)
    print(f"Saved summary statistics to: {summary_file}")

    # Save structure metadata
    structures_df = pd.DataFrame(
        [{"acronym": k, "structure_id": v} for k, v in hypo_id_map.items()]
    )
    structures_file = PROCESSED_DIR / "hypothalamus_structures.csv"
    structures_df.to_csv(structures_file, index=False)
    print(f"Saved structure list to: {structures_file}")

    print("\n" + "=" * 60)
    print("Done!")
    print("=" * 60)


if __name__ == "__main__":
    main()
