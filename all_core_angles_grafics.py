from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

from core_phi_psi_functions import (read_alignment, read_angles, collect_core_angles, plot_scatter)


alignment_folder = Path("MAFFT_clusters_chains_final")
angles_file = Path("all_phi_psi_for_cifs.txt")

output_folder = Path("angles_chain_scatter")
output_folder.mkdir(exist_ok=True)


alignment_files = sorted(alignment_folder.rglob("*.txt"))

for alignment_file in alignment_files:

    family_name = alignment_file.stem
    print(f"Processing {family_name}...")

    family_output_folder = (output_folder / family_name)
    family_output_folder.mkdir(parents=True, exist_ok=True)

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

    core_cos_phis = [[np.cos(phi) for phi in phis if phi is not None] for phis in core_phis]
    core_cos_psis = [[np.cos(psi) for psi in psis if psi is not None] for psis in core_psis]
    core_sin_phis = [[np.sin(phi) for phi in phis if phi is not None] for phis in core_phis]
    core_sin_psis = [[np.sin(psi) for psi in psis if psi is not None] for psis in core_psis]

    positions = np.arange(len(core))
    core_residues = [core[i] for i in range(len(core))]

    plot_scatter(core_cos_phis,
                positions,
                r"cos($\phi$)",
                f"{family_name}: cos($\\phi$) angles in the core",
                family_output_folder
                / f"{family_name}_phi_cosine.png")

    plot_scatter(core_cos_psis,
                positions,
                r"cos($\psi$)",
                f"{family_name}: cos($\\psi$) angles in the core",
                family_output_folder
                / f"{family_name}_psi_cosine.png")

    plot_scatter(core_sin_phis,
                positions,
                r"sin($\phi$)",
                f"{family_name}: sin($\\phi$) angles in the core",
                family_output_folder
                / f"{family_name}_phi_sine.png")

    plot_scatter(core_sin_psis,
                positions,
                r"sin($\psi$)",
                f"{family_name}: sin($\\psi$) angles in the core",
                family_output_folder
                / f"{family_name}_psi_sine.png")

print("Done.")