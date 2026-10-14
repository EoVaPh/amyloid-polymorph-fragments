from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

from core_phi_psi_functions import (read_alignment, read_angles, collect_core_angles)


alignment_folder = Path("MAFFT_clusters_chains_final")
angles_file = Path("all_phi_psi_for_cifs.txt")
core_sequences_file = Path("clusters_chain_cores.txt")
output_folder = Path("core_angles_boxplots")
output_folder.mkdir(parents=True, exist_ok=True)


def read_core_sequences(filename):
    '''Read core sequences from a FASTA-like file.'''

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


def plot_core_boxplots(family_name, core_sequence, core_phis, core_psis, output_file):
    '''Plot phi and psi distributions for every core position.'''

    if not core_phis or not core_psis:
        print(f"  No angle data for {family_name}")
        return

    if len(core_phis) != len(core_psis):
        raise ValueError(f"{family_name}: phi and psi data have different lengths.")

    if len(core_sequence) != len(core_phis):
        raise ValueError(f"{family_name}: core sequence length ({len(core_sequence)}) "
                         f"does not match the number of core positions ({len(core_phis)}).")

    phi_data = []
    psi_data = []
    phi_positions = []
    psi_positions = []

    aa_colors = {"A": "#FFD92F", "V": "#FFD92F", "I": "#FFD92F", "L": "#FFD92F", "M": "#FFD92F",
                "K": "#D7191C", "R": "#D7191C",
                "D": "#2C7BB6", "E": "#2C7BB6",
                "F": "#F28E2B", "W": "#F28E2B", "Y": "#F28E2B", "H": "#F28E2B",
                "S": "#2CA25F", "T": "#2CA25F", "N": "#2CA25F", "Q": "#2CA25F", "C": "#2CA25F",
                "G": "#9467BD",
                "P": "#8C564B"}

    for i, (phis, psis) in enumerate(zip(core_phis, core_psis)):
        phis = [np.degrees(value) for value in phis if value is not None and np.isfinite(value)]
        psis = [np.degrees(value) for value in psis if value is not None and np.isfinite(value)]
        if phis:
            phi_data.append(phis)
            phi_positions.append(i + 1 - 0.18)
        if psis:
            psi_data.append(psis)
            psi_positions.append(i + 1 + 0.18)

    fig, ax = plt.subplots(figsize=(max(12, len(core_sequence) * 0.55), 6))

    if phi_data:
        phi_boxes = ax.boxplot(
            phi_data,
            positions=phi_positions,
            widths=0.3,
            patch_artist=True,
            manage_ticks=False,
            showfliers=False
        )

        for box in phi_boxes["boxes"]:
            box.set_facecolor("#BB5EB6")
            box.set_alpha(0.8)

    if psi_data:
        psi_boxes = ax.boxplot(
            psi_data,
            positions=psi_positions,
            widths=0.3,
            patch_artist=True,
            manage_ticks=False,
            showfliers=False
        )

        for box in psi_boxes["boxes"]:
            box.set_facecolor("#46B89E")
            box.set_alpha(0.8)

    ax.set_xticks(range(1, len(core_sequence) + 1))
    ax.set_xticklabels(list(core_sequence), fontsize=8)
    for label, aa in zip(ax.get_xticklabels(), core_sequence):
        label.set_color(aa_colors.get(aa.upper(), "#555555"))

    ax.set_xlabel("Core position (amino acid; position number)")
    ax.set_ylabel("Ramachandran angle (°)")
    ax.set_ylim(-180, 180)
    ax.set_yticks([-180, -90, 0, 90, 180])

    ax.axhline(0, color="gray", linewidth=0.7, alpha=0.5)

    ax.set_title(f"{family_name}: phi and psi distributions")

    ax.legend(
        handles=[
            Patch(facecolor="#BB5EB6", label="φ"),
            Patch(facecolor="#46B89E", label="ψ")
        ],
        loc="upper left",
        bbox_to_anchor=(1.01, 1),
        frameon=False
    )

    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    fig.tight_layout()
    fig.savefig(output_file, dpi=600, bbox_inches="tight")
    plt.close(fig)


def main():
    core_sequences = read_core_sequences(core_sequences_file)
    alignment_files = sorted(alignment_folder.rglob("*.txt"))

    for alignment_file in alignment_files:
        family_name = alignment_file.stem
        print(f"Processing {family_name}...")

        alignment = read_alignment(alignment_file)
        if not alignment:
            print("  Empty alignment, skipped.")
            continue

        angles = read_angles(alignment, angles_file)

        _, core_phis, core_psis = collect_core_angles(alignment, angles, threshold=0.865)
        core_sequence = core_sequences.get(family_name)
        if core_sequence is None:
            print("  Core sequence not found, skipped.")
            continue

        plot_core_boxplots(family_name, core_sequence, core_phis, core_psis, output_folder / f"{family_name}.png")


if __name__ == "__main__":
    main()