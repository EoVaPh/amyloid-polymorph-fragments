from pathlib import Path
from itertools import combinations

import numpy as np
from Bio.Align import PairwiseAligner


# ============================================================
# SETTINGS
# ============================================================

families_file = Path('clusters_renamed.txt')

chains_folder = Path("extracted_chains")

output_file = Path("pairwise_tm_scores.txt")


# ============================================================
# SEQUENCE ALIGNMENT
# ============================================================

def align_seqs(seq_1: str, seq_2: str) -> tuple:
    '''Do global pairwise alignment of two amino acid sequences.'''

    aligner = PairwiseAligner(
        open_gap_score=-3,
        extend_gap_score=-2,
        mismatch_score=-1,
        left_gap_score=-1,
        right_gap_score=-1
    )

    alignment = aligner.align(seq_1, seq_2)[0]

    return alignment[0], alignment[1]


# ============================================================
# READ SEQRES
# ============================================================

def read_seq(seq_file_path: str) -> str:
    '''Read a SEQRES sequence from a text file.'''

    seq_file = open(seq_file_path, 'r')

    seq = seq_file.read().strip()

    seq_file.close()

    return seq


# ============================================================
# READ STRUCTURAL CHAIN
# ============================================================

def read_chain(chain_file_path: str) -> tuple:
    '''
    Read structural chain.

    File format:

    x y z chain residue_number amino_acid
    '''

    sequence = []
    coordinates = []
    residue_numbers = []

    with open(chain_file_path, "r", encoding="utf-8") as file:

        for line in file:

            line = line.strip()

            if not line:
                continue

            parts = line.split()

            x = float(parts[0])
            y = float(parts[1])
            z = float(parts[2])

            residue_number = int(parts[4])
            amino_acid = parts[5]

            coordinates.append([x, y, z])
            sequence.append(amino_acid)
            residue_numbers.append(residue_number)

    return (
        "".join(sequence),
        np.array(coordinates, dtype=float),
        residue_numbers
    )


# ============================================================
# READ FAMILIES
# ============================================================

def read_families(file_path: str) -> dict:
    '''
    Read families from:

    >Family name
    pdbid1
    pdbid2
    ...
    '''

    families = {}

    current_family = None

    with open(file_path, 'r', encoding='utf-8') as file:

        for line in file:

            line = line.strip()

            if not line:
                continue

            if line.startswith('>'):

                current_family = line[1:].strip()

                families[current_family] = []

            elif current_family is not None:

                pdbid = line.lower()

                if pdbid not in families[current_family]:

                    families[current_family].append(pdbid)

    return families


# ============================================================
# CHAIN → SEQRES
# ============================================================

def mapping_chain_to_seqres(seqres: str, chain: str) -> tuple:
    '''
    Align a structural chain to its corresponding SEQRES sequence
    and map chain positions to SEQRES positions.
    '''

    aligned_seqres, aligned_chain = align_seqs(seqres, chain)

    chain_to_seqres = {}

    seqres_position = 0
    chain_position = 0

    for aa_seqres, aa_chain in zip(
        aligned_seqres,
        aligned_chain
    ):

        current_seqres_position = None
        current_chain_position = None

        if aa_seqres != '-':

            current_seqres_position = seqres_position

            seqres_position += 1

        if aa_chain != '-':

            current_chain_position = chain_position

            chain_position += 1

        if (
            current_seqres_position is not None
            and current_chain_position is not None
        ):

            chain_to_seqres[current_chain_position] = (
                current_seqres_position
            )

    return (
        aligned_seqres,
        aligned_chain,
        chain_to_seqres
    )


# ============================================================
# SEQRES → COMMON ALIGNMENT
# ============================================================

def mapping_seqres_to_common_alignment(
        aligned_seq: str
) -> dict:
    '''
    Map SEQRES positions to positions
    in the common alignment.
    '''

    seqres_to_common_alignment = {}

    seqres_position = 0

    for alignment_position, aa in enumerate(aligned_seq):

        if aa != '-':

            seqres_to_common_alignment[
                seqres_position
            ] = alignment_position

            seqres_position += 1

    return seqres_to_common_alignment


# ============================================================
# PROJECT CHAIN TO COMMON ALIGNMENT
# ============================================================

def project_chain_to_common_alignment(
        chain: str,
        chain_to_seqres: dict,
        seqres_to_common_alignment: dict,
        alignment_length: int
) -> list:
    '''
    Project a structural chain onto the common SEQRES alignment.
    '''

    result = ['-'] * alignment_length

    for chain_position, seqres_position in chain_to_seqres.items():

        if seqres_position not in seqres_to_common_alignment:
            continue

        alignment_position = (
            seqres_to_common_alignment[seqres_position]
        )

        result[alignment_position] = chain[chain_position]

    return result


# ============================================================
# CHAIN POSITION → COMMON ALIGNMENT POSITION
# ============================================================

def project_chain_positions_to_common_alignment(
        chain_to_seqres: dict,
        seqres_to_common_alignment: dict
) -> dict:
    '''
    Map chain positions to common alignment positions.

    Result:

    common_alignment_position → chain_position
    '''

    chain_to_common = {}

    for chain_position, seqres_position in chain_to_seqres.items():

        if seqres_position not in seqres_to_common_alignment:
            continue

        alignment_position = (
            seqres_to_common_alignment[seqres_position]
        )

        chain_to_common[alignment_position] = chain_position

    return chain_to_common


# ============================================================
# KABSCH ALIGNMENT
# ============================================================

def kabsch_distances(
        coordinates_1: np.ndarray,
        coordinates_2: np.ndarray
) -> np.ndarray:
    '''
    Superimpose two sets of Cα coordinates using the Kabsch
    algorithm and return pairwise distances after alignment.
    '''

    center_1 = coordinates_1.mean(axis=0)
    center_2 = coordinates_2.mean(axis=0)

    x = coordinates_1 - center_1
    y = coordinates_2 - center_2

    covariance = x.T @ y

    U, S, Vt = np.linalg.svd(covariance)

    correction = np.eye(3)

    if np.linalg.det(Vt.T @ U.T) < 0:

        correction[-1, -1] = -1

    rotation = (
        Vt.T
        @ correction
        @ U.T
    )

    x_rotated = x @ rotation

    distances = np.linalg.norm(
        x_rotated - y,
        axis=1
    )

    return distances


# ============================================================
# TM-SCORE
# ============================================================

def calculate_d0(length: int) -> float:
    '''
    Calculate the TM-score normalization parameter d0.
    '''

    if length <= 15:

        return 0.5

    return 1.24 * (length - 15) ** (1 / 3) - 1.8


def calculate_tm_score(
        distances: np.ndarray,
        normalization_length: int
) -> float:
    '''
    Calculate TM-score using the supplied normalization length.
    '''

    if len(distances) == 0:
        return np.nan

    d0 = calculate_d0(normalization_length)

    tm_score = np.mean(
        1.0 / (
            1.0 + (distances / d0) ** 2
        )
    )

    # Convert from mean over aligned residues
    # to normalization by the complete sequence length.
    tm_score *= len(distances) / normalization_length

    return tm_score


# ============================================================
# SYMMETRIC TM-SCORE
# ============================================================

def calculate_symmetric_tm_score(
        distances: np.ndarray,
        length_1: int,
        length_2: int
) -> tuple:
    '''
    Calculate two directional TM-scores and their mean.

    TM_1 uses length_1 as normalization length.
    TM_2 uses length_2 as normalization length.

    The final score is:

        TM_sym = (TM_1 + TM_2) / 2
    '''

    tm_1 = calculate_tm_score(
        distances,
        length_1
    )

    tm_2 = calculate_tm_score(
        distances,
        length_2
    )

    tm_symmetric = (
        tm_1 + tm_2
    ) / 2

    return tm_1, tm_2, tm_symmetric


# ============================================================
# TM-SCORE FOR ONE PAIR
# ============================================================

def calculate_pair_tm_score(
        pdbid_1: str,
        pdbid_2: str
) -> float:
    '''
    Calculate symmetric TM-score for one pair of structures.

    Sequence mapping:

        chain
          ↓
        SEQRES
          ↓
        common SEQRES alignment
          ↓
        common experimental positions
          ↓
        Cα coordinates
          ↓
        Kabsch
          ↓
        symmetric TM-score
    '''

    # --------------------------------------------------------
    # File paths
    # --------------------------------------------------------

    seqres_file_1 = (
        chains_folder / f"{pdbid_1}_seq.txt"
    )

    seqres_file_2 = (
        chains_folder / f"{pdbid_2}_seq.txt"
    )

    chain_file_1 = (
        chains_folder / f"{pdbid_1}.txt"
    )

    chain_file_2 = (
        chains_folder / f"{pdbid_2}.txt"
    )

    # --------------------------------------------------------
    # Read sequences
    # --------------------------------------------------------

    seqres_1 = read_seq(seqres_file_1)
    seqres_2 = read_seq(seqres_file_2)

    chain_1, coords_1, residue_numbers_1 = (
        read_chain(chain_file_1)
    )

    chain_2, coords_2, residue_numbers_2 = (
        read_chain(chain_file_2)
    )

    # --------------------------------------------------------
    # Chain → SEQRES mapping
    # --------------------------------------------------------

    (
        aligned_seqres_1,
        aligned_chain_1,
        chain_to_seqres_1
    ) = mapping_chain_to_seqres(
        seqres_1,
        chain_1
    )

    (
        aligned_seqres_2,
        aligned_chain_2,
        chain_to_seqres_2
    ) = mapping_chain_to_seqres(
        seqres_2,
        chain_2
    )

    # --------------------------------------------------------
    # Pairwise alignment of the two SEQRES sequences
    # --------------------------------------------------------

    common_seqres_1, common_seqres_2 = align_seqs(
        seqres_1,
        seqres_2
    )

    alignment_length = len(common_seqres_1)

    # --------------------------------------------------------
    # SEQRES → common alignment
    # --------------------------------------------------------

    seqres_to_common_1 = (
        mapping_seqres_to_common_alignment(
            common_seqres_1
        )
    )

    seqres_to_common_2 = (
        mapping_seqres_to_common_alignment(
            common_seqres_2
        )
    )

    # --------------------------------------------------------
    # Chain → common alignment
    # --------------------------------------------------------

    chain_to_common_1 = (
        project_chain_positions_to_common_alignment(
            chain_to_seqres_1,
            seqres_to_common_1
        )
    )

    chain_to_common_2 = (
        project_chain_positions_to_common_alignment(
            chain_to_seqres_2,
            seqres_to_common_2
        )
    )

    # --------------------------------------------------------
    # Find positions where both structures have
    # experimental residues
    # --------------------------------------------------------

    common_positions = sorted(
        set(chain_to_common_1)
        &
        set(chain_to_common_2)
    )

    # --------------------------------------------------------
    # Check number of common residues
    # --------------------------------------------------------

    if len(common_positions) < 3:

        return np.nan

    # --------------------------------------------------------
    # Extract corresponding Cα coordinates
    # --------------------------------------------------------

    matched_coords_1 = np.array([
        coords_1[
            chain_to_common_1[position]
        ]
        for position in common_positions
    ])

    matched_coords_2 = np.array([
        coords_2[
            chain_to_common_2[position]
        ]
        for position in common_positions
    ])

    # --------------------------------------------------------
    # Kabsch superposition
    # --------------------------------------------------------

    distances = kabsch_distances(
        matched_coords_1,
        matched_coords_2
    )

    # --------------------------------------------------------
    # Symmetric TM-score
    # --------------------------------------------------------

    tm_1, tm_2, tm_symmetric = (
        calculate_symmetric_tm_score(
            distances,
            len(seqres_1),
            len(seqres_2)
        )
    )

    return tm_symmetric


# ============================================================
# MAIN
# ============================================================

def main():

    families = read_families(
        families_file
    )

    with open(
        output_file,
        "w",
        encoding="utf-8"
    ) as output:

        for family_name, pdbids in families.items():

            print(
                f"Processing family: {family_name}"
            )

            output.write(
                f">{family_name}\n"
            )

            # ------------------------------------------------
            # Every pair within the family
            # ------------------------------------------------

            for pdbid_1, pdbid_2 in combinations(
                pdbids,
                2
            ):

                try:

                    tm_score = calculate_pair_tm_score(
                        pdbid_1,
                        pdbid_2
                    )

                    if np.isnan(tm_score):

                        print(
                            f"  {pdbid_1} vs {pdbid_2}: "
                            f"not enough common residues"
                        )

                        continue

                    output.write(
                        f"{pdbid_1}_{pdbid_2} "
                        f"{tm_score:.3f}\n"
                    )

                except Exception as error:

                    print(
                        f"  ERROR: "
                        f"{pdbid_1} vs {pdbid_2}: "
                        f"{error}"
                    )

            output.write("\n")

    print()
    print(
        f"Results saved to: {output_file}"
    )


# ============================================================
# RUN
# ============================================================

if __name__ == "__main__":
    main()