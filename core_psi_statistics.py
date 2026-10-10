from pathlib import Path
from collections import Counter

import numpy as np
import matplotlib.pyplot as plt

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

from matplotlib.colors import LinearSegmentedColormap

from adjustText import adjust_text

from core_phi_psi_functions import read_alignment, read_angles, collect_core_angles

import matplotlib

matplotlib.rcParams['figure.dpi'] = 300
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rc('font', family='STIXGeneral')
matplotlib.rc('font', weight='ultralight')


alignment_folder = Path("mafft_clusters_chains_final")
core_sequences_file = Path("clusters_chain_cores.txt")
angles_file = Path("all_phi_psi_for_cifs.txt")
alignment_files = sorted(alignment_folder.rglob("*.txt"))

core_threshold = 0.865
positive_threshold = 0.8


def read_core_sequences(filename: Path) -> dict:
    '''Read core sequences in FASTA-like format:
    >family
    core sequence'''

    sequences = {}
    current_name = None

    with open(filename, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                current_name = line[1:].strip()
                sequences[current_name] = ""

            elif current_name is not None:
                sequences[current_name] += line

    return sequences


def collect_amino_acid_statistics(core_sequence, core_psis, threshold):
    '''Count amino acids at core positions where
    MORE than `threshold` of available psi cosines
    are positive and total occurrences of each 
    amino acid in core.'''

    n_positions = len(core_sequence)

    if len(core_psis) != n_positions:
        raise ValueError(f"Core psi position count mismatch: "
                         f"{len(core_psis)} vs {n_positions}")

    positive_counts = np.zeros(n_positions, dtype=int)
    total_counts = np.zeros(n_positions, dtype=int)

    for position, psis in enumerate(core_psis):
        for psi in psis:
            if psi is None:
                continue

            cosine = np.cos(psi)
            if not np.isfinite(cosine):
                continue

            total_counts[position] += 1
            if cosine > 0:
                positive_counts[position] += 1

    fractions = np.divide(positive_counts,
                          total_counts,
                          out=np.full(n_positions, np.nan, dtype=float),
                          where=total_counts > 0)

    total_statistics = Counter()
    turn_statistics = Counter()

    for position, aa in enumerate(core_sequence):
        total_statistics[aa] += 1
        if fractions[position] > threshold:
            turn_statistics[aa] += 1

    return total_statistics, turn_statistics


def plot_amino_acid_statistics(statistics, ax, panel_label="A"):
    """Plot amino acid frequency in turns."""

    amino_acids = list("ACDEFGHIKLMNPQRSTVWY")
    counts = [statistics[aa] for aa in amino_acids]

    bars = ax.bar(amino_acids, counts, width=0.7)

    ax.set_xlabel("Amino acid", fontsize=14)
    ax.set_ylabel("Frequency in turns", fontsize=14)

    ax.set_title(f"{panel_label}. Amino acid frequency in turns", fontsize=18, pad=12)

    max_count = np.nanmax(counts)

    ax.set_ylim(0, max_count * 1.15 if max_count > 0 else 1)

    ax.tick_params(axis="both", labelsize=12)

    for bar, count in zip(bars, counts):
        ax.text(bar.get_x() + bar.get_width() / 2,
                bar.get_height(),
                f"{count:.3f}",
                ha="center",
                va="bottom",
                fontsize=10)


def plot_amino_acid_mass_polarity(statistics, ax, panel_label="B"):

    amino_acid_masses = {
        "A": 89.09, "C": 121.16, "D": 133.10, "E": 147.13,
        "F": 165.19, "G": 75.07, "H": 155.16, "I": 131.17,
        "K": 146.19, "L": 131.17, "M": 149.21, "N": 132.12,
        "P": 115.13, "Q": 146.15, "R": 174.20, "S": 105.09,
        "T": 119.12, "V": 117.15, "W": 204.23, "Y": 181.19
    }

    amino_acid_polarity = {
        "A": 8.1, "C": 5.5, "D": 13.0, "E": 12.3,
        "F": 5.2, "G": 9.0, "H": 10.4, "I": 5.2,
        "K": 11.3, "L": 4.9, "M": 5.7, "N": 11.6,
        "P": 8.0, "Q": 10.5, "R": 10.5, "S": 9.2,
        "T": 8.6, "V": 5.9, "W": 5.4, "Y": 6.2
    }

    amino_acids = list(amino_acid_masses.keys())

    masses = np.array([amino_acid_masses[aa] for aa in amino_acids])
    polarities = np.array([amino_acid_polarity[aa] for aa in amino_acids])
    frequencies = np.array([statistics[aa] for aa in amino_acids])

    ax.scatter(masses,
               polarities,
               frequencies,
               s=20,
               color="#a53860")

    z_min = frequencies.min()
    ax.scatter(masses,
               polarities,
               np.full_like(frequencies, z_min),
               s=10,
               alpha=0.75)

    y_min = polarities.min()
    ax.scatter(masses,
               np.full_like(polarities, y_min),
               frequencies,
               s=10,
               alpha=0.75)

    x_min = masses.min()
    ax.scatter(np.full_like(masses, x_min),
               polarities,
               frequencies,
               s=10,
               alpha=0.75)

    dx = (max(masses) - min(masses)) * 0.015
    dy = (max(polarities) - min(polarities)) * 0.015
    dz = (max(frequencies) - min(frequencies)) * 0.015

    for aa, mass, polarity, frequency in zip(amino_acids, masses, polarities, frequencies):
        ax.text(mass + dx,
                polarity + dy,
                frequency + dz,
                aa,
                fontsize=11)

    ax.set_xlim(min(masses) - (max(masses) - min(masses)) * 0.01,
                max(masses) + (max(masses) - min(masses)) * 0.01)
    ax.set_ylim(min(polarities) - (max(polarities) - min(polarities)) * 0.01,
                max(polarities) + (max(polarities) - min(polarities)) * 0.01)
    ax.set_zlim(min(frequencies) - (max(frequencies) - min(frequencies)) * 0.01,
                max(frequencies) + (max(frequencies) - min(frequencies)) * 0.01)

    ax.set_xlabel("Molecular mass (Da)", labelpad=10, fontsize=14)
    ax.set_ylabel("Polarity (Grantham)", labelpad=10, fontsize=14)
    ax.set_zlabel("Frequency in turns", labelpad=10, fontsize=14)

    ax.set_title(f"{panel_label}. Amino acid properties and frequency", fontsize=18, pad=12)

    ax.tick_params(axis="both", labelsize=10)

    ax.view_init(elev=25,azim=45)


def prepare_pca_data(turn_frequency):
    amino_acids = list("ACDEFGHIKLMNPQRSTVWY")
    masses = {
        "A": 89.09, "C": 121.16, "D": 133.10, "E": 147.13,
        "F": 165.19, "G": 75.07, "H": 155.16, "I": 131.17,
        "K": 146.19, "L": 131.17, "M": 149.21, "N": 132.12,
        "P": 115.13, "Q": 146.15, "R": 174.20, "S": 105.09,
        "T": 119.12, "V": 117.15, "W": 204.23, "Y": 181.19
    }

    hydrophobicity = {
        "A": 1.8, "C": 2.5, "D": -3.5, "E": -3.5,
        "F": 2.8, "G": -0.4, "H": -3.2, "I": 4.5,
        "K": -3.9, "L": 3.8, "M": 1.9, "N": -3.5,
        "P": -1.6, "Q": -3.5, "R": -4.5, "S": -0.8,
        "T": -0.7, "V": 4.2, "W": -0.9, "Y": -1.3
    }

    polarity = {
        "A": 8.1, "C": 5.5, "D": 13.0, "E": 12.3,
        "F": 5.2, "G": 9.0, "H": 10.4, "I": 5.2,
        "K": 11.3, "L": 4.9, "M": 5.7, "N": 11.6,
        "P": 8.0, "Q": 10.5, "R": 10.5, "S": 9.2,
        "T": 8.6, "V": 5.9, "W": 5.4, "Y": 6.2
    }

    aromaticity = {aa: int(aa in "FWY") for aa in amino_acids}
    charge = {aa: 1 if aa in "KR" else -1 if aa in "DE" else 0 for aa in amino_acids}

    X = [[turn_frequency[aa], masses[aa], hydrophobicity[aa], aromaticity[aa], charge[aa], polarity[aa]] for aa in amino_acids]
    
    return amino_acids, np.array(X, dtype=float)


def pca_plot(turn_frequency, ax_pca, ax_table, panel_label_pca="C",
             panel_label_table="D"):

    amino_acids, X = prepare_pca_data(turn_frequency)
    X_scaled = StandardScaler().fit_transform(X)
    pca = PCA()
    X_pca = pca.fit_transform(X_scaled)
    turn_values = np.array([turn_frequency[aa] for aa in amino_acids])
    norm = plt.Normalize(vmin=turn_values.min(), vmax=turn_values.max())
    cmap = LinearSegmentedColormap.from_list("custom_cmap", ["#366FA8", "#4E2D7C", "#AC3232"])

    for aa, x, y in zip(amino_acids, X_pca[:, 0], X_pca[:, 1]):
        color = cmap(norm(turn_frequency[aa]))
        ax_pca.scatter(x,
                       y,
                       s=60,
                       color=color)

        if aa in "PH":
            offset = (5, -5)
        else:
            offset = (5, 5)

        ax_pca.annotate(aa,
                       (x, y),
                       xytext=offset,
                       textcoords="offset points",
                       fontsize=11)

    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])

    cbar = ax_pca.figure.colorbar(sm, ax=ax_pca, pad=0.02)
    cbar.set_label("Frequency of occurrence in turns", fontsize=13)
    cbar.ax.tick_params(labelsize=10)

    ax_pca.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0] * 100:.1f}%)", fontsize=14)
    ax_pca.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1] * 100:.1f}%)", fontsize=14)
    ax_pca.set_title(f"{panel_label_pca}. Principal component analysis", fontsize=18, pad=12)

    ax_pca.tick_params(axis="both", labelsize=11)

    variables = ["Frequency",
                 "Molecular mass",
                 "Hydrophobicity",
                 "Aromaticity",
                 "Charge",
                 "Polarity"]

    loadings = pca.components_.T

    im = ax_table.imshow(loadings,
                         aspect="auto",
                         cmap="coolwarm",
                         vmin=-1,
                         vmax=1)

    ax_table.set_xticks(range(6))
    ax_table.set_xticklabels([f"PC{i + 1}\n"
                              f"({pca.explained_variance_ratio_[i] * 100:.1f}%)" for i in range(6)],
                              fontsize=11)

    ax_table.set_yticks(range(6))
    ax_table.set_yticklabels(variables, fontsize=12)

    for i in range(6):
        for j in range(6):
            ax_table.text(j,
                          i,
                          f"{loadings[i, j]:.2f}",
                          ha="center",
                          va="center",
                          fontsize=11)

    ax_table.figure.colorbar(im,
                             ax=ax_table,
                             pad=0.02,
                             label="Loading")

    ax_table.set_xlabel("Principal component", fontsize=14)
    ax_table.set_title(f"{panel_label_table}. PCA loadings", fontsize=18, pad=12)


core_sequences = read_core_sequences(core_sequences_file)

total_amino_acid_statistics = Counter()
turn_amino_acid_statistics = Counter()

processed_families = 0
skipped_families = 0

for alignment_file in alignment_files:
    family_name = alignment_file.stem
    print(f"Processing {family_name}...")

    if family_name not in core_sequences:
        print(f"  Core sequence not found "
              f"for {family_name}")
        skipped_families += 1
        continue

    core_sequence = core_sequences[family_name]

    alignment = read_alignment(alignment_file)

    if not alignment:
        print(f"  Alignment is empty: "
              f"{alignment_file}")
        skipped_families += 1
        continue

    angles = read_angles(alignment, angles_file)
    core, core_phis, core_psis = collect_core_angles(alignment, angles, threshold=core_threshold)

    if not core:
        print(f"  Core was not found "
              f"for {family_name}")
        skipped_families += 1
        continue

    if len(core_sequence) != len(core):
        print(f"  WARNING: core length mismatch "
              f"for {family_name}: "
              f"sequence = {len(core_sequence)}, "
              f"core positions = {len(core)}")
        skipped_families += 1
        continue

    total_statistics, turn_statistics = collect_amino_acid_statistics(core_sequence=core_sequence, core_psis=core_psis, threshold=positive_threshold)
    total_amino_acid_statistics.update(total_statistics)
    turn_amino_acid_statistics.update(turn_statistics)
    processed_families += 1

print(f"Processed families: "
      f"{processed_families}")
print(f"Skipped families: "
      f"{skipped_families}")

amino_acids = list("ACDEFGHIKLMNPQRSTVWY")

turn_frequency = {aa: (turn_amino_acid_statistics[aa] / total_amino_acid_statistics[aa] 
                       if total_amino_acid_statistics[aa] > 0
                       else np.nan) for aa in amino_acids}

fig = plt.figure(figsize=(16, 13))

gs = fig.add_gridspec(2,
                      2,
                      wspace=0.28,
                      hspace=0.32)

ax_A = fig.add_subplot(gs[0, 0])
plot_amino_acid_statistics(statistics=turn_frequency, ax=ax_A, panel_label="A")

ax_B = fig.add_subplot(gs[0, 1], projection="3d")
plot_amino_acid_mass_polarity(statistics=turn_frequency, ax=ax_B, panel_label="B")
pos = ax_B.get_position()
ax_B.set_position([pos.x0 - 0.02, 
                   pos.y0 - 0.02,
                   pos.width + 0.04,
                   pos.height + 0.04])

ax_C = fig.add_subplot(gs[1, 0])
ax_D = fig.add_subplot(gs[1, 1])
pca_plot(turn_frequency=turn_frequency,
         ax_pca=ax_C,
         ax_table=ax_D,
         panel_label_pca="C",
         panel_label_table="D")

fig.savefig(
    "amino_acid_combined_figure.png",
    dpi=600,
    bbox_inches="tight")

plt.close(fig)