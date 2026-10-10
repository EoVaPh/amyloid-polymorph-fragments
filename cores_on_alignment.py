from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from Bio import AlignIO

from core_phi_psi_functions import read_alignment


alignment_folder = Path("MAFFT_clusters_chains_final")
core_sequences_file = Path("clusters_chain_cores.txt")
output_folder = Path("cores_on_alignment")
output_folder.mkdir(exist_ok=True)


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


def read_core_sequences(filename):
    '''Read the core sequence from a file in fasta-like format.'''
    
    core_sequences = {}
    family_name = None
    sequence = []

    with open(filename, "r") as file:
        for line in file:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if family_name is not None:
                    core_sequences[family_name] = "".join(sequence)
                family_name = line[1:].split()[0]
                sequence = []
            else:
                sequence.append(line)

    if family_name is not None:
        core_sequences[family_name] = "".join(sequence)

    return core_sequences


def plot_alignment(alignment, core_sequence, output_file):
    if not alignment:
        raise ValueError("The alignment file is empty.")

    alignment_length = len(alignment[0][1])

    if any(len(sequence) != alignment_length
           for _, sequence in alignment):
        raise ValueError("Sequences have different alignment lengths.")

    core_positions = set(find_core(alignment, 0.875))

    fig, ax = plt.subplots(
        figsize=(16, max(3, len(alignment) * 0.28))
    )

    for row, (pdbid, sequence) in enumerate(alignment):
        y = len(alignment) - row - 1

        segments = []
        colors = []
        widths = []

        for position, aa in enumerate(sequence):
            in_core = position in core_positions
            is_gap = aa == "-"

            segments.append([
                (position, y),
                (position + 1, y)
            ])

            colors.append("#D95F59" if in_core else "#555555")
            widths.append(0.7 if is_gap else 3.0)

        collection = LineCollection(
            segments,
            colors=colors,
            linewidths=widths,
            capstyle="butt"
        )

        ax.add_collection(collection)

        ax.text(
            -0.5, y, pdbid,
            ha="right",
            va="center",
            fontsize=8
        )

    for position, aa in zip(core_positions, core_sequence):
        ax.text(
            position + 0.5,
            len(alignment) - 0.35,
            aa,
            ha="center",
            va="bottom",
            fontsize=8,
            color="#D95F59",
            fontweight="bold"
        )

    ax.set_xlim(0, alignment_length)
    ax.set_ylim(-1, len(alignment) + 1)
    ax.set_xlabel("Alignment position")
    ax.set_yticks([])

    ax.spines["left"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    fig.tight_layout()
    fig.savefig(output_file, dpi=600, bbox_inches="tight")
    plt.close(fig)


core_sequences = read_core_sequences(core_sequences_file)
alignment_files = sorted(alignment_folder.rglob("*.txt"))

for alignment_file in alignment_files:
    family_name = alignment_file.stem
    print(f"Processing {family_name}...")

    alignment = read_alignment(alignment_file)
    if not alignment:
        print(f"  Alignment is empty: "
                f"{alignment_file}")
        continue

    if family_name not in core_sequences:
        print(f"  Core sequence not found for {family_name}")
        continue

    plot_alignment(alignment, core_sequences[family_name], output_folder / f"{family_name}.png")