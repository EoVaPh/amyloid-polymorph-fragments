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


def read_angles(alignment: list, filename: Path) -> dict:
    '''Read angle file:

    >PDB ID
    chain   number  aminoacid   phi psi'''


    needed_structures = {structure_id for structure_id, _ in alignment}

    angles = {}

    current_structure = None

    with open(filename, "r", encoding="utf-8") as f:

        for line in f:
            line = line.strip()

            if not line:
                continue

            if line.endswith((".pdb", ".cif")):

                current_structure = (
                    line.lstrip("> ").strip()[:4]
                )

                if current_structure in needed_structures:
                    angles[current_structure] = {}

                continue

            if current_structure not in needed_structures:
                continue

            parts = line.split()

            if len(parts) != 5:
                continue

            chain = parts[0]
            residue_number = int(parts[1])
            aa = parts[2]

            phi = (None if parts[3] == "None" else float(parts[3]))

            psi = (None if parts[4] == "None" else float(parts[4]))

            angles[current_structure][residue_number] = {
                "chain": chain,
                "aa": aa,
                "phi": phi,
                "psi": psi
            }

    return angles


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

        mapping[alignment_position] = (
            residue_numbers[residue_index]
        )

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


def collect_core_angles(alignment, angles, threshold=0.875):
    '''
    Find the core and collect phi/psi values.

    A residue missing because of a gap contributes
    None for both phi and psi.

    Returns
    -------
    core : str
        Core sequence from the first sequence.

    core_phis : list
        Phi values for every core position.

    core_psis : list
        Psi values for every core position.
    '''
    
    if not alignment:
        return "", [], []

    # Core is determined once for the entire family
    core_positions = find_core(alignment, threshold)

    if not core_positions:
        return "", [], []

    # Create residue mappings for all structures
    residue_mappings = {}

    for structure_id, sequence in alignment:

        if structure_id not in angles:
            # Keep the structure in the alignment.
            # Its angles will simply be None.
            residue_mappings[structure_id] = None
            continue

        residue_mappings[structure_id] = (make_residue_mapping(sequence, angles[structure_id]))

    # Core sequence from the first alignment sequence
    reference_sequence = alignment[0][1]

    core = "".join(reference_sequence[position] for position in core_positions)

    core_phis = []
    core_psis = []

    # Process every core alignment position
    for position in core_positions:

        phis = []
        psis = []

        # Go through EVERY structure in the family
        for structure_id, sequence in alignment:

            mapping = residue_mappings.get(structure_id)

            # No angle data for this structure
            if mapping is None:
                phis.append(None)
                psis.append(None)
                continue

            # Gap in this particular sequence
            if position not in mapping:
                phis.append(None)
                psis.append(None)
                continue

            residue_number = mapping[position]

            # Mapping explicitly contains None
            if residue_number is None:
                phis.append(None)
                psis.append(None)
                continue

            # Get phi/psi for this actual residue
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


def calculate_cos_std(phi_angles):
    '''For each core position:

        1. remove None / NaN;
        2. calculate cos(angle);
        3. calculate std(cos(angle)).

    Returns a list of std values for all core positions.'''
    
    std_values = []

    for angles in phi_angles:
        angles = np.array([np.nan if x is None else x for x in angles], dtype=float)
        angles = angles[~np.isnan(angles)]

        if len(angles) == 0:
            std_values.append(np.nan)
            continue

        cos_values = np.cos(angles)
        std = np.std(cos_values)
        std_values.append(std)

    return std_values


def plot_scatter(data, positions, ylabel, title, filename):
    
    plt.figure(figsize=(16, 6))

    font_path = r"C:\\Users\\User\\AppData\\Local\\Microsoft\\Windows\\Fonts\\MavenPro-VariableFont_wght.ttf"
    fm.fontManager.addfont(font_path)
    plt.rcParams["font.family"] = "Maven Pro"

    rng = np.random.default_rng(42)
    jitter_width = 0.15

    for position, values in enumerate(data):

        plt.axvline(
            position,
            color="lightgray",
            alpha=0.6,
            linewidth=0.8
        )

        positive = [value for value in values if value > 0]
        negative = [value for value in values if value < 0]
        zero = [value for value in values if value == 0]

        negative_x = (
            position
            + rng.uniform(
                -jitter_width,
                jitter_width,
                len(negative)
            )
        )

        plt.scatter(
            negative_x,
            negative,
            color="#254441",
            alpha=0.8,
            edgecolors='none',
            s=20
        )

        positive_x = (
            position
            + rng.uniform(
                -jitter_width,
                jitter_width,
                len(positive)
            )
        )

        plt.scatter(
            positive_x,
            positive,
            color="#c9184a",
            alpha=0.8,
            edgecolors='none',
            s=20
        )

        zero_x = (
            position
            + rng.uniform(
                -jitter_width,
                jitter_width,
                len(zero)
            )
        )

        plt.scatter(
            zero_x,
            zero,
            color="black",
            alpha=0.8,
            edgecolors='none',
            s=20
        )

    plt.xlabel("Residue number in core", fontsize=16)
    plt.ylabel(ylabel, fontsize=16)
    plt.title(title, fontsize=22)

    plt.xticks(positions, rotation=90)
    plt.ylim(-1.05, 1.05)

    plt.tick_params(axis='x', labelsize=14)
    plt.tick_params(axis='y', labelsize=16)

    plt.tight_layout()
    plt.savefig(
        filename,
        dpi=300,
        bbox_inches="tight"
    )

    plt.close()