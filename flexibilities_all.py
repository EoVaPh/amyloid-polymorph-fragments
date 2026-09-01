from pathlib import Path
import math
from Bio import SeqIO
from Bio.SeqUtils import ProtParamData


input_file = Path("all_seqres_from_cifs.txt")
output_file = Path("flexibilities_new.txt")

flex = ProtParamData.Flex
window_size = 9
weights = [0.25, 0.4375, 0.625, 0.8125, 1.0]


def calculate_window_flexibility(window: str) -> float:
    '''Calculate flexibility for one 9-residue window.
    The formula is identical to Biopython 1.87.
    If the window contains X, return NaN.'''

    # Unknown amino acid
    if "X" in window:
        return math.nan

    score = 0.0

    for j in range(window_size // 2):
        front = window[j]
        back = window[window_size - j - 1]
        score += (flex[front] + flex[back]) * weights[j]

    middle = window[window_size // 2 + 1]
    score += flex[middle]

    return score / 5.25


def calculate_flexibilities(sequence: str) -> list[float]:
    '''Calculate flexibility for every possible 9-residue window.
    Unlike Biopython 1.87, the last possible window is included.'''

    scores = []

    for start in range(len(sequence) - window_size + 1):
        window = sequence[start:start + window_size]
        score = calculate_window_flexibility(window)
        scores.append(score)

    return scores


with open(output_file, "w", encoding="utf-8") as out:
    for record in SeqIO.parse(input_file, "fasta"):

        sequence_id = record.id
        sequence = str(record.seq)

        if len(sequence) < window_size:
            print(f"{sequence_id}: sequence is too short "
                  f"({len(sequence)} aa), skipped")
            continue

        flexibilities = calculate_flexibilities(sequence)
        out.write(f">{sequence_id}\n")
        values = []

        for value in flexibilities:
            if math.isnan(value):
                values.append("NaN")
            else:
                values.append(f"{value:.3f}")

        out.write(" ".join(values))
        out.write("\n")

        print(f"{sequence_id}: "
              f"{len(sequence)} aa -> "
              f"{len(flexibilities)} windows")

print(f"\nResults saved to: {output_file}")