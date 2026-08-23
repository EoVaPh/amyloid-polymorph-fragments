import requests
from collections import defaultdict


url = "https://ff54g8ykd7.execute-api.eu-central-1.amazonaws.com/prod/structures"
output_file = "amyloid_explorer_families.txt"


def get_structures():
    '''Fetch all structures from the Amyloid Explorer API.'''

    response = requests.get(url)
    response.raise_for_status()

    data = response.json()

    return data["structures"]


def group_by_protein(structures):
    '''Group PDB structures by protein name.'''

    families = defaultdict(set)

    for structure in structures:
        protein = structure.get("Protein")
        pdb_id = structure.get("Name")

        if not protein or not pdb_id:
            continue

        families[protein].add(pdb_id.lower())

    return families


def save_families(families, output_file):
    '''Save protein families and their PDB IDs to a text file.'''

    with open(output_file, "w", encoding="utf-8") as f:

        for protein in sorted(families):
            f.write(f">{protein}\n")

            for pdb_id in sorted(families[protein]):
                f.write(f"{pdb_id}\n")

            f.write("\n")


def main():
    print("Fetching data from Amyloid Explorer...")

    structures = get_structures()
    print(f"Number of structures retrieved: {len(structures)}")

    families = group_by_protein(structures)
    print(f"Number of protein families found: {len(families)}")

    save_families(families, output_file)
    print(f"Results saved to: {output_file}")


if __name__ == "__main__":
    main()