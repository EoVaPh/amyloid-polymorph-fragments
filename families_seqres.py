from pathlib import Path
from Bio import SeqIO


fasta_file = Path("all_seqres_from_cifs.txt")
families_file = Path("amyloid_explorer_families.txt")

output_folder = Path("families_seqres_from_cifs")
output_folder.mkdir(exist_ok=True)


families = {}

current_family = None

with open(families_file, "r", encoding="utf-8") as file:

    for line in file:
        line = line.strip()

        if not line:
            continue

        if line.startswith(">"):
            current_family = line[1:].strip()
            families[current_family] = []

        else:
            families[current_family].append(line)

records = list(SeqIO.parse(fasta_file, "fasta"))

for number, (family_name, structures) in enumerate(families.items(), start=1):

    structures = set(structures)
    selected_records = []

    for record in records:
        structure_name = record.id[:4]

        if structure_name in structures:
            selected_records.append(record)

    safe_filename = family_name.replace("/", "_")
    output_file = output_folder / f"{safe_filename}.fasta"
    SeqIO.write(selected_records, output_file, "fasta")

    print(f"[{number}/{len(families)}]"
          f"{family_name}: {len(selected_records)} sequences")