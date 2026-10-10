from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

from core_phi_psi_functions import (read_alignment, read_angles, collect_core_angles, plot_scatter)


alignment_folder = Path("MAFFT_clusters_chains_final")
angles_file = Path("all_phi_psi_for_cifs.txt")
output_folder = Path("ramachandran_plots")
output_folder.mkdir(exist_ok=True)


def calculate_trig_std(angle_lists, trig_function):
    '''Calculate the standard deviation of sin(angle) or cos(angle)
    for each core position.'''

    std_values = []

    for angles in angle_lists:
        angles = np.array([np.nan if x is None else x for x in angles], dtype=float)
        angles = angles[~np.isnan(angles)]

        if len(angles) == 0:
            std_values.append(np.nan)
            continue

        trig_values = trig_function(angles)
        std_values.append(np.std(trig_values))

    return std_values


def calculate_trig_variability(core_phis, core_psis):
    '''Calculate combined variability of sin(phi), cos(phi),
    sin(psi), and cos(psi) for every core position.'''

    sin_phi_std = calculate_trig_std(core_phis, np.sin)
    cos_phi_std = calculate_trig_std(core_phis, np.cos)
    sin_psi_std = calculate_trig_std(core_psis, np.sin)
    cos_psi_std = calculate_trig_std(core_psis, np.cos)

    variability = []

    for i in range(len(core_phis)):
        values = [sin_phi_std[i], cos_phi_std[i], sin_psi_std[i], cos_psi_std[i]]

        if not all(np.isfinite(x) for x in values):
            variability.append(np.nan)
            continue

        value = np.sqrt(np.mean(np.square(values)))
        variability.append(value)

    return variability


def get_top_variable_positions(core, core_phis, core_psis, top_n=4):
    '''Return the top positions with the highest combined
    trigonometric variability.'''

    variability = calculate_trig_variability(core_phis, core_psis)
    valid_positions = [i for i, value in enumerate(variability) if np.isfinite(value)]
    top_positions = sorted(valid_positions, key=lambda i: variability[i], reverse=True)[:top_n]

    return [{"position": i,
             "aa": core[i],
             "variability": variability[i],
             "phis": core_phis[i],
             "psis": core_psis[i]} for i in top_positions]


def plot_ramachandran_family(family_name, top_positions, output_dir):
    '''Plot Ramachandran angles for the n most variable
    core positions of one family.'''

    output_dir = Path(output_dir)

    fig, ax = plt.subplots(figsize=(9, 9))
    colors = ["#366FA8", "#AC3232", "#4C956C", "#D49A36"]

    for i, item in enumerate(top_positions):
        phis = item["phis"]
        psis = item["psis"]
        valid_angles = [(phi, psi) for phi, psi in zip(phis, psis) if phi is not None
                                                                    and psi is not None
                                                                    and np.isfinite(phi)
                                                                    and np.isfinite(psi)]

        if not valid_angles:
            continue

        angles = np.degrees(np.array(valid_angles, dtype=float))
        position = item["position"] + 1
        aa = item["aa"]
        ax.scatter(
            angles[:, 0],
            angles[:, 1],
            s=35,
            alpha=0.75,
            color=colors[i % len(colors)],
            label=f"Position {position} ({aa})"
        )

    ax.set_xlim(-180, 180)
    ax.set_ylim(-180, 180)

    ax.set_xticks(np.arange(-180, 181, 60))
    ax.set_yticks(np.arange(-180, 181, 60))

    ax.set_xlabel(r"$\phi$ (°)", fontsize=14)
    ax.set_ylabel(r"$\psi$ (°)", fontsize=14)

    ax.set_title(f"Ramachandran plot: {family_name}", fontsize=14)

    ax.axhline(0, color="black", linewidth=0.6, alpha=0.5)
    ax.axvline(0, color="black", linewidth=0.6, alpha=0.5)

    ax.grid(alpha=0.2)
    ax.set_aspect("equal", adjustable="box")

    ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1), fontsize=10)
    fig.tight_layout()

    fig.savefig(
        output_dir / f"{family_name}_ramachandran.png",
        dpi=300,
        bbox_inches="tight")

    plt.close(fig)


alignment_files = sorted(alignment_folder.rglob("*.txt"))

for alignment_file in alignment_files:
    family_name = alignment_file.stem
    print(f"Processing {family_name}...")

    alignment = read_alignment(alignment_file)
    if not alignment:
        print(f"  Alignment is empty: "
              f"{alignment_file}")
        continue

    angles = read_angles(alignment, angles_file)
    core, core_phis, core_psis = collect_core_angles(alignment, angles, threshold=0.865)
    if not core:
        print(f"  Core was not found for "
              f"{family_name}")
        continue

    top_positions = get_top_variable_positions(core, core_phis, core_psis, top_n=4)
    plot_ramachandran_family(family_name, top_positions, output_dir=output_folder)