from pathlib import Path
from aiupred import AIUPred

predictor = AIUPred()

def calculate_seqres_idr(pdb_id: str, folder: str) -> list:
    '''Find a SEQRES file by PDB ID, calculate AIUPred disorder
    propensities, and return them as a list.'''

    folder = Path(folder)
    file_path = folder / f'{pdb_id}_seq.txt'

    if not file_path.exists():
        raise FileNotFoundError(f'File not found: {file_path}')

    with open(file_path, 'r') as file:
        sequence = file.read().strip()

    disorder_propensities = predictor.predict_disorder(sequence)

    return disorder_propensities.tolist()


def calculate_all_seqres_idr(folder: str, output_file: str):
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
                idr = calculate_seqres_idr(pdb_id, folder)
            except Exception as e:
                print(f"Error processing {pdb_id}: {e}")
                continue

            out.write(f">{pdb_id}\n")
            out.write(" ".join(f"{value:.3f}" for value in idr) + "\n")

            print("OK")

    print(f"Results saved to: {output_file}")


calculate_all_seqres_idr(folder='extracted_chains', output_file='idrs.txt')