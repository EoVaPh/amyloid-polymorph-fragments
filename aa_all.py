from pathlib import Path


def count_seqres_aa(pdb_id: str, folder: str, target_aa: str) -> list:
    '''Find a SEQRES file by PDB ID, count number of amino acid residues, and
    return them as a list where 1 is a coincidence and 0 is not.'''

    folder = Path(folder)
    file_path = folder / f'{pdb_id}_seq.txt'

    if not file_path.exists():
        raise FileNotFoundError(f'File not found: {file_path}')

    with open(file_path, 'r') as file:
        sequence = file.read().strip()

    number_aa = list()

    for aa in sequence:
        if aa == target_aa:
            number_aa.append(1)
        else:
            number_aa.append(0)

    return number_aa


def calculate_all_seqres_aa(folder: str, output_file: str, target_aa: str):
    '''Calculate IDR disorder propensities for all SEQRES files in a folder
    and save the results to one output file.'''

    folder = Path(folder)
    output_file = Path(output_file)

    seq_files = sorted(folder.glob("*_seq.txt"))
    total = len(seq_files)

    print(f"Found {total} SEQRES files.")

    with open(output_file, "w") as out:
        for i, file_path in enumerate(seq_files, start=1):

            pdb_id = file_path.stem.removesuffix("_seq")
            print(f"[{i}/{total}] Processing {pdb_id}...", end=" ")

            try:
                count = count_seqres_aa(pdb_id, folder, target_aa)
            except Exception as e:
                print(f"Error processing {pdb_id}: {e}")
                continue

            out.write(f">{pdb_id}\n")
            out.write(" ".join(f"{value}" for value in count) + "\n")

            print("OK")

    print(f"Results saved to: {output_file}")


calculate_all_seqres_aa(folder='extracted_chains', output_file='Ps.txt',
                        target_aa='P')
