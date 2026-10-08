"""Utilities for reading alignments and analyzing core phi/psi angles."""

from pathlib import Path
from typing import (
    Dict,
    List,
    Optional,
    Sequence,
    Tuple,
    TypedDict,
)

import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np

Alignment = List[Tuple[str, str]]


class ResidueData(TypedDict):
    """Angle data for one residue."""

    chain: str
    aa: str
    phi: Optional[float]
    psi: Optional[float]


StructureAngles = Dict[int, ResidueData]
AnglesByStructure = Dict[str, StructureAngles]
ResidueMapping = Dict[int, Optional[int]]


def read_alignment(filename: Path) -> Alignment:
    """Read a multiple sequence alignment file.

    The expected format is::

        PDB_ID  ----AAAAAAA...
        PDB_ID  -----AAAAAA...

    :param filename: Path to the alignment file.
    :type filename: pathlib.Path
    :return: A list of ``(structure_id, sequence)`` pairs.
    :rtype: list[tuple[str, str]]
    """
    alignment: Alignment = []

    with filename.open("r", encoding="utf-8") as file_handle:
        for line in file_handle:
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


def read_angles(
    alignment: Alignment,
    filename: Path,
) -> AnglesByStructure:
    """Read phi and psi angles for structures in an alignment.

    The expected format is::

        >PDB ID
        chain  number  amino_acid  phi  psi

    Missing angles are represented by the string ``None``.

    :param alignment: Alignment records returned by
        :func:`read_alignment`.
    :type alignment: list[tuple[str, str]]
    :param filename: Path to the angle file.
    :type filename: pathlib.Path
    :return: Mapping from structure ID to residue number and angle data.
    :rtype: AnglesByStructure
    """
    needed_structures = {structure_id for structure_id, _ in alignment}

    angles: AnglesByStructure = {}
    current_structure: Optional[str] = None

    with filename.open("r", encoding="utf-8") as file_handle:
        for line in file_handle:
            line = line.strip()

            if not line:
                continue

            if line.endswith((".pdb", ".cif")):
                current_structure = line.lstrip("> ").strip()[:4]

                if current_structure in needed_structures:
                    angles[current_structure] = {}

                continue

            if current_structure is None:
                continue

            if current_structure not in needed_structures:
                continue

            parts = line.split()

            if len(parts) != 5:
                continue

            chain = parts[0]
            residue_number = int(parts[1])
            aa = parts[2]
            phi = (
                None
                if parts[3] == "None"
                else float(parts[3])
            )
            psi = (
                None
                if parts[4] == "None"
                else float(parts[4])
            )

            angles[current_structure][residue_number] = {
                "chain": chain,
                "aa": aa,
                "phi": phi,
                "psi": psi,
            }

    return angles


def make_residue_mapping(
    sequence: str,
    structure_angles: StructureAngles,
) -> ResidueMapping:
    """Map alignment positions to residue numbers.

    Gaps receive ``None`` and do not consume a residue number.

    :param sequence: Aligned sequence for one structure.
    :type sequence: str
    :param structure_angles: Angle data keyed by residue number.
    :type structure_angles: StructureAngles
    :return: Mapping from alignment position to residue number.
    :rtype: ResidueMapping
    """
    residue_numbers = sorted(structure_angles.keys())
    mapping: ResidueMapping = {}
    residue_index = 0

    for alignment_position, aa in enumerate(sequence):
        if aa == "-":
            mapping[alignment_position] = None
            continue

        if residue_index >= len(residue_numbers):
            mapping[alignment_position] = None
            continue

        mapping[alignment_position] = residue_numbers[residue_index]
        residue_index += 1

    return mapping


def find_core(
    alignment: Alignment,
    threshold: float = 0.875,
) -> List[int]:
    """Find the longest continuous core region in an alignment.

    A position is included when at least ``threshold`` of the sequences
    contain a residue at that position. Gaps are allowed when their
    fraction does not exceed ``1 - threshold``.

    :param alignment: Alignment records returned by
        :func:`read_alignment`.
    :type alignment: list[tuple[str, str]]
    :param threshold: Minimum fraction of sequences that must contain a
        residue at a position for it to be considered part of the core.
    :type threshold: float
    :return: Alignment positions in the longest core region.
    :rtype: list[int]
    """
    if not alignment:
        return []

    alignment_length = len(alignment[0][1])
    number_of_sequences = len(alignment)

    best_start: Optional[int] = None
    best_end: Optional[int] = None
    best_length = 0

    current_start: Optional[int] = None
    current_length = 0

    for position in range(alignment_length):
        number_of_residues = sum(
            sequence[position] != "-"
            for _, sequence in alignment
        )
        fraction = number_of_residues / number_of_sequences

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

    if best_start is None or best_end is None:
        return []

    return list(range(best_start, best_end + 1))


def collect_core_angles(
    alignment: Alignment,
    angles: AnglesByStructure,
    threshold: float = 0.875,
) -> Tuple[
    str,
    List[List[Optional[float]]],
    List[List[Optional[float]]],
]:
    """Collect phi and psi angles for the core region.

    The core is determined once for the entire family. A missing residue,
    whether due to a gap or missing angle data, contributes ``None`` for
    both phi and psi.

    :param alignment: Alignment records returned by
        :func:`read_alignment`.
    :type alignment: list[tuple[str, str]]
    :param angles: Angle data returned by :func:`read_angles`.
    :type angles: AnglesByStructure
    :param threshold: Minimum fraction of sequences that must contain a
        residue at a position for it to be considered part of the core.
    :type threshold: float
    :return: A tuple containing the core sequence, phi values, and psi
        values. Each angle list contains one inner list per core
        position.
    :rtype: tuple[str, list[list[float or None]], list[list[float or None]]]
    """
    if not alignment:
        return "", [], []

    core_positions = find_core(alignment, threshold)

    if not core_positions:
        return "", [], []

    residue_mappings: Dict[
        str,
        Optional[ResidueMapping],
    ] = {}

    for structure_id, sequence in alignment:
        if structure_id not in angles:
            residue_mappings[structure_id] = None
            continue

        residue_mappings[structure_id] = make_residue_mapping(
            sequence,
            angles[structure_id],
        )

    reference_sequence = alignment[0][1]
    core = "".join(
        reference_sequence[position]
        for position in core_positions
    )

    core_phis: List[List[Optional[float]]] = []
    core_psis: List[List[Optional[float]]] = []

    for position in core_positions:
        phis: List[Optional[float]] = []
        psis: List[Optional[float]] = []

        for structure_id, _ in alignment:
            mapping = residue_mappings.get(structure_id)

            if mapping is None:
                phis.append(None)
                psis.append(None)
                continue

            if position not in mapping:
                phis.append(None)
                psis.append(None)
                continue

            residue_number = mapping[position]

            if residue_number is None:
                phis.append(None)
                psis.append(None)
                continue

            data = angles[structure_id].get(residue_number)

            if data is None:
                phis.append(None)
                psis.append(None)
                continue

            phis.append(data.get("phi"))
            psis.append(data.get("psi"))

        core_phis.append(phis)
        core_psis.append(psis)

    return core, core_phis, core_psis


def calculate_cos_std(
    phi_angles: List[List[Optional[float]]],
) -> List[float]:
    """Calculate the standard deviation of cosine-transformed phi angles.

    For each core position, this function removes missing values, computes
    the cosine of each phi angle, and returns the standard deviation of
    those cosine values.

    :param phi_angles: Phi angles for each core position.
    :type phi_angles: list[list[float or None]]
    :return: Standard deviation of cosine values for each core position.
        Positions with no valid angles receive ``numpy.nan``.
    :rtype: list[float]
    """
    std_values: List[float] = []

    for angles in phi_angles:
        angles_array = np.array(
            [np.nan if angle is None else angle for angle in angles],
            dtype=float,
        )
        angles_array = angles_array[~np.isnan(angles_array)]

        if len(angles_array) == 0:
            std_values.append(np.nan)
            continue

        cos_values = np.cos(angles_array)
        std = np.std(cos_values)
        std_values.append(std)

    return std_values


def plot_scatter(
    data: Sequence[Sequence[float]],
    positions: Sequence[int],
    ylabel: str,
    title: str,
    filename: Path,
) -> None:
    """Plot a jittered scatter plot for core positions.

    :param data: Values to plot. Each inner sequence contains the values
        for one core position.
    :type data: sequence[sequence[float]]
    :param positions: Tick labels for the core positions.
    :type positions: sequence[int]
    :param ylabel: Label for the y-axis.
    :type ylabel: str
    :param title: Plot title.
    :type title: str
    :param filename: Output image filename.
    :type filename: pathlib.Path
    :return: None
    :rtype: None
    """
    plt.figure(figsize=(16, 6))

    font_path = (
        r"C:\\Users\\User\\AppData\\Local\\Microsoft\\Windows\\Fonts"
        r"\\MavenPro-VariableFont_wght.ttf"
    )
    fm.fontManager.addfont(font_path)
    plt.rcParams["font.family"] = "Maven Pro"

    rng = np.random.default_rng(42)
    jitter_width = 0.15

    for position, values in enumerate(data):
        plt.axvline(
            position,
            color="lightgray",
            alpha=0.6,
            linewidth=0.8,
        )

        negative = [value for value in values if value < 0]
        positive = [value for value in values if value > 0]
        zero = [value for value in values if value == 0]

        negative_x = position + rng.uniform(
            -jitter_width,
            jitter_width,
            len(negative),
        )
        plt.scatter(
            negative_x,
            negative,
            color="#254441",
            alpha=0.8,
            edgecolors="none",
            s=20,
        )

        positive_x = position + rng.uniform(
            -jitter_width,
            jitter_width,
            len(positive),
        )
        plt.scatter(
            positive_x,
            positive,
            color="#c9184a",
            alpha=0.8,
            edgecolors="none",
            s=20,
        )

        zero_x = position + rng.uniform(
            -jitter_width,
            jitter_width,
            len(zero),
        )
        plt.scatter(
            zero_x,
            zero,
            color="black",
            alpha=0.8,
            edgecolors="none",
            s=20,
        )

    plt.xlabel("Residue number in core", fontsize=16)
    plt.ylabel(ylabel, fontsize=16)
    plt.title(title, fontsize=22)

    plt.xticks(positions, rotation=90)
    plt.ylim(-1.05, 1.05)

    plt.tick_params(axis="x", labelsize=14)
    plt.tick_params(axis="y", labelsize=16)

    plt.tight_layout()
    plt.savefig(
        filename,
        dpi=300,
        bbox_inches="tight",
    )
    plt.close()
