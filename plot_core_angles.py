"""
Plot phi and psi angles of core residues for each aligned family.

This script reads alignment files from ``MAFFT_clusters_chains_final``,
retrieves phi and psi angles from ``phi_psi.txt``, identifies core positions,
and writes scatter plots of sine and cosine transformed angles to
``angles_chain_scatter``.
"""

from pathlib import Path

import numpy as np

from core_phi_psi_functions import (
    collect_core_angles,
    plot_scatter,
    read_alignment,
    read_angles,
)

#: Similarity threshold used to identify core positions.
CORE_THRESHOLD = 0.865

#: Directory containing input alignment files.
ALIGNMENT_FOLDER = Path("MAFFT_clusters_chains_final")

#: File containing phi and psi angle data.
ANGLES_FILE = Path("phi_psi.txt")

#: Root directory where output scatter plots are written.
OUTPUT_FOLDER = Path("angles_chain_scatter")


def process_alignment_file(
    alignment_file: Path,
    angles_file: Path,
    output_folder: Path,
) -> None:
    """Process one alignment file and generate core-angle scatter plots.

    :param alignment_file: Path to the alignment file.
    :type alignment_file: pathlib.Path
    :param angles_file: Path to the phi/psi angle file.
    :type angles_file: pathlib.Path
    :param output_folder: Root directory for output plots.
    :type output_folder: pathlib.Path
    :return: None
    :rtype: None
    """
    family_name = alignment_file.stem
    print(f"Processing {family_name}...")

    family_output_folder = output_folder / family_name
    family_output_folder.mkdir(parents=True, exist_ok=True)

    alignment = read_alignment(alignment_file)
    if not alignment:
        print(f"  Alignment is empty: {alignment_file}")
        return

    angles = read_angles(alignment, angles_file)
    core, core_phis, core_psis = collect_core_angles(
        alignment,
        angles,
        threshold=CORE_THRESHOLD,
    )

    if not core:
        print(f"  Core was not found for {family_name}")
        return

    core_cos_phis = [
        [np.cos(phi) for phi in phis if phi is not None]
        for phis in core_phis
    ]
    core_cos_psis = [
        [np.cos(psi) for psi in psis if psi is not None]
        for psis in core_psis
    ]
    core_sin_phis = [
        [np.sin(phi) for phi in phis if phi is not None]
        for phis in core_phis
    ]
    core_sin_psis = [
        [np.sin(psi) for psi in psis if psi is not None]
        for psis in core_psis
    ]

    positions = np.arange(len(core))

    plot_scatter(
        core_cos_phis,
        positions,
        r"cos($\phi$)",
        rf"{family_name}: cos($\phi$) angles in the core",
        family_output_folder / f"{family_name}_phi_cosine.png",
    )

    plot_scatter(
        core_cos_psis,
        positions,
        r"cos($\psi$)",
        rf"{family_name}: cos($\psi$) angles in the core",
        family_output_folder / f"{family_name}_psi_cosine.png",
    )

    plot_scatter(
        core_sin_phis,
        positions,
        r"sin($\phi$)",
        rf"{family_name}: sin($\phi$) angles in the core",
        family_output_folder / f"{family_name}_phi_sine.png",
    )

    plot_scatter(
        core_sin_psis,
        positions,
        r"sin($\psi$)",
        rf"{family_name}: sin($\psi$) angles in the core",
        family_output_folder / f"{family_name}_psi_sine.png",
    )


def main() -> None:
    """Run the core-angle plotting workflow.

    :return: None
    :rtype: None
    """
    OUTPUT_FOLDER.mkdir(parents=True, exist_ok=True)

    alignment_files = sorted(ALIGNMENT_FOLDER.rglob("*.txt"))

    for alignment_file in alignment_files:
        process_alignment_file(
            alignment_file,
            ANGLES_FILE,
            OUTPUT_FOLDER,
        )

    print("Done.")


if __name__ == "__main__":
    main()
