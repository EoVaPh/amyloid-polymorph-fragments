# Download structures using entries from Amyloid Atlas.

import urllib.request
import re
import os

# Sawaya, Michael R., et al.
# "The expanding amyloid family: Structure, stability, function, and pathogenesis."
# Cell 184.19 (2021): 4857-4873.
amyloid_atlas_file = open('DBs/Amyloid Atlas 2026.html', 'r',
                          errors='ignore')
lines = amyloid_atlas_file.readlines()
amyloid_atlas_file.close()

fibril_pdbs = []

for line in lines:
    match = None

    try:
        match = re.search(
            r'<a\s+href="https://www\.rcsb\.org/structure/([A-Za-z0-9]+)"',
            line
        ).group(1)
    except:
        continue

    fibril_pdbs.append(match.strip())

cnt = 0

os.makedirs('CIFs', exist_ok=True)

print(len(fibril_pdbs), 'entries in Amyloid Atlas found.')

for pdb in fibril_pdbs:
    cif_path = os.path.join('CIFs', pdb + '.cif')
    pdb_path = os.path.join('CIFs', pdb + '.pdb')

    if os.path.exists(cif_path):
        print(f'{pdb}: CIF already exists')
        cnt += 1
        print(cnt)
        continue

    try:
        urllib.request.urlretrieve(
            'http://files.rcsb.org/download/' + pdb + '.cif',
            'CIFs/' + pdb + '.cif'
        )

        if os.path.exists(pdb_path):
            os.remove(pdb_path)
            print(f'{pdb}: old PDB removed')

    except Exception:
        print('Failed to download ' + pdb + ' in mmCIF format.')

        try:

            urllib.request.urlretrieve(
                f'https://files.rcsb.org/download/{pdb}.pdb',
                pdb_path
            )

            print(f'{pdb}: PDB downloaded')

        except Exception:
            print(f'{pdb}: neither CIF nor PDB is available')


    cnt += 1
    print(cnt)

# Check for missed entries from Amyloid Explorer.

# Kyriazis, Vasilis, et al.
# "Amyloid Explorer: a global atlas of amyloid fibril structures and thermodynamic principles."
# bioRxiv (2025): 2025-10.
AmyloidExplorer_file = open('DBs/Amyloid Explorer.html', 'r')
text = AmyloidExplorer_file.read()
AmyloidExplorer_file.close()

# Check which structures are still missing in the pool directory.
substring = 'strct='
num_chars = 4
pattern = rf"{re.escape(substring)}(.{{{num_chars}}})"

AmyloidExplorer_pdbs = re.findall(pattern, text)

path = 'CIFs'
pdbs = os.listdir(path)

miss_pdbs = []

for pdb in AmyloidExplorer_pdbs:
    # Check all possible names of a file with the structure.
    cif_exists = (pdb.lower() + '.cif' in pdbs or pdb.upper() + '.cif' in pdbs)

    if not cif_exists:
        miss_pdbs.append(pdb)

print('Missed structure files from Amyloid Explorer:', miss_pdbs)


import os
import urllib.request
import json


url = (
    'https://ff54g8ykd7.execute-api.eu-central-1.amazonaws.com'
    '/prod/structures'
)

with urllib.request.urlopen(url) as response:
    data = json.load(response)

structures = data["structures"]

print(f"Found {len(structures)} structures in Amyloid Explorer.")

pdb_folder = "CIFs"
os.makedirs(pdb_folder, exist_ok=True)

already_exists = 0
downloaded = 0
not_found = 0

for structure in structures:

    pdb = structure["Name"].lower()

    cif_path = os.path.join(pdb_folder, pdb + ".cif")
    pdb_path = os.path.join(pdb_folder, pdb + ".pdb")

    if os.path.exists(cif_path) or os.path.exists(pdb_path):
        print(f"{pdb}: already exists")
        already_exists += 1
        continue

    cif_url = (f"https://files.rcsb.org/download/{pdb}.cif")

    try:
        urllib.request.urlretrieve(cif_url, cif_path)
        print(f"{pdb}: CIF downloaded")
        downloaded += 1

    except Exception:
        if os.path.exists(cif_path):
            os.remove(cif_path)

        print(f"{pdb}: CIF not found")
        not_found += 1

print()
print("Finished.")
print(f"Already existed: {already_exists}")
print(f"Downloaded CIF:  {downloaded}")
print(f"CIF not found:   {not_found}")
