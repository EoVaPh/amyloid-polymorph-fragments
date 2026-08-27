from pathlib import Path
import subprocess
from io import StringIO
from Bio import SeqIO
 
mafft = Path(r"C:\Users\User\Documents\bioinf\smtb\MAFFT\mafft-win\mafft.bat")
input_folder = Path("families_seqres_from_cifs")
output_folder = Path("MAFFT_seqres_from_cifs_multiple_alignment")
output_folder.mkdir(parents=True, exist_ok=True)

for fasta_file in input_folder.glob("*.fasta"):

    output_file = output_folder / f"{fasta_file.stem}.txt"

    print(f"Aligning: {fasta_file.name}")

    command = [
        "cmd.exe",
        "/C",
        str(mafft),
        "--auto",
        str(fasta_file)
    ]

    result = subprocess.run(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        encoding="utf-8"
    )

    if result.returncode != 0:
        print(f"ERROR: {fasta_file.name}")
        print(result.stderr)
        continue

    # Read MAFFT alignment from stdout
    alignment = StringIO(result.stdout)

    # Save in the required format
    with open(output_file, "w", encoding="utf-8") as out:

        for record in SeqIO.parse(alignment, "fasta"):

            sequence_id = record.id
            sequence = str(record.seq).replace("\n", "")

            out.write(f"{sequence_id} {sequence}\n")

    print(f"Saved: {output_file}")

print("\nAll alignments completed.")