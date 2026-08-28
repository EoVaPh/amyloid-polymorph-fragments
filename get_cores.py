from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm


def read_alignment(filename: Path) -> list:
    '''Read alignment in format:

    PDB ID  ----AAAAAAA...
    PDB ID  -----AAAAAA...'''

    alignment = []

    with open(filename, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            parts = line.split(maxsplit=1)

            if len(parts) != 2:
                continue

            structure_id = parts[0]
            sequence = parts[1]

            alignment.append((structure_id, sequence))

    return alignment


def get_family_structure_names(filename: Path) -> list:
    '''Get structure names for family from alignment file
    in format:

    PDB_ID sequence
    PDB_ID sequence'''

    structure_names = []

    with open(filename, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()

            if not line:
                continue

            parts = line.split(maxsplit=1)

            if len(parts) != 2:
                continue

            structure_names.append(parts[0])

    return structure_names


def make_residue_mapping(sequence, structure_angles):
    '''Map alignment positions to residue numbers.
    A gap gets None and does not consume a residue number.'''

    residue_numbers = sorted(structure_angles.keys())

    mapping = {}
    residue_index = 0

    for alignment_position, aa in enumerate(sequence):

        # Gap: no residue exists at this position
        if aa == "-":
            mapping[alignment_position] = None
            continue

        # No more residues in the structure
        if residue_index >= len(residue_numbers):
            mapping[alignment_position] = None
            continue

        mapping[alignment_position] = (residue_numbers[residue_index])
        residue_index += 1

    return mapping


def find_core(alignment, threshold=0.875):
    '''Find the longest continuous alignment region where
    at least `threshold` of sequences contain a residue.

    Gaps are allowed if their fraction does not exceed
    (1 - threshold).'''

    if not alignment:
        return []

    alignment_length = len(alignment[0][1])
    number_of_sequences = len(alignment)

    best_start = None
    best_end = None
    best_length = 0

    current_start = None
    current_length = 0

    for position in range(alignment_length):

        number_of_residues = sum(
            sequence[position] != "-"
            for _, sequence in alignment
        )

        fraction = (number_of_residues / number_of_sequences)

        if fraction >= threshold:

            if current_start is None:
                current_start = position
                current_length = 1

            else:
                current_length += 1

            if current_length > best_length:

                best_length = current_length
                best_start = current_start
                best_end = position

        else:

            current_start = None
            current_length = 0

    if best_start is None:
        return []

    return list(range(best_start, best_end + 1))


alignment_folder = Path('MAFFT_seqres_from_cifs_multiple_alignment')
alignment_files = sorted(alignment_folder.rglob("*.txt"))
output_file = Path("cores.txt")
count = 0

with open(output_file, "w", encoding="utf-8") as out:
    for alignment_file in alignment_files:
        alignment = read_alignment(alignment_file)
        core_positions = find_core(alignment, threshold=0.93)

        if not core_positions:
            print(f"{alignment_file.name}: core not found")
            count += 1
            continue

        # Take the sequence of the first structure
        # to represent the core.
        reference_sequence = alignment[0][1]
        core = "".join(reference_sequence[position] for position in core_positions)

        family_name = alignment_file.stem

        out.write(f">{family_name}\n"
                  f"{core}\n")

print(f'For {count} families core was not found')