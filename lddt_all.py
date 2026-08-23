from pathlib import Path
from itertools import combinations
from typing import List, Tuple, Sequence
import math

from Bio.Align import PairwiseAligner


input_folder = Path('extracted_chains')
families_file = Path('amyloid_explorer_families.txt')
output_file = Path('lddts.txt')


def read_seq(seq_file_path: str) -> str:
    """Read a SEQRES sequence from a text file."""

    seq_file = open(seq_file_path, 'r')
    seq = seq_file.read().strip()
    seq_file.close()

    return seq


def read_chain(chain_file_path: str) -> str:
    """Read a structural chain from a text file."""

    chain_file = open(chain_file_path, 'r')
    residues = chain_file.readlines()
    chain_file.close()

    seq = ''

    for residue in residues:
        resnum = int(residue.split()[-2])
        seq += residue.split()[-1]

    return seq


def get_positions(pdbid: str) -> list:
    """Read atom coordinates from a chain file."""

    chain_file = open(input_folder / f'{pdbid}.txt', 'r')

    residues = chain_file.readlines()
    chain_file.close()

    positions = []

    for residue in residues:

        residue_tokens = residue.strip().split()

        positions.append(
            (
                float(residue_tokens[0]),
                float(residue_tokens[1]),
                float(residue_tokens[2])
            )
        )

    return positions


def read_families(file_path: Path) -> dict:
    '''Read families from a file of the form:

    >Family name
    pdbid1
    pdbid2
    pdbid3

    >Another family
    pdbid4
    ...'''

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


def align_seqs(seq_1: str, seq_2: str) -> tuple:
    '''Do global pairwise alignment of two amino acid sequences.'''

    aligner = PairwiseAligner(open_gap_score=-3,
                            extend_gap_score=-2,
                            mismatch_score=-1,
                            left_gap_score=-1,
                            right_gap_score=-1)

    alignment = aligner.align(seq_1, seq_2)[0]

    return alignment[0], alignment[1]


def mapping_seqres_to_common_alignment(aligned_seq: str) -> dict:
    '''Map SEQRES positions to positions in the common alignment.'''

    seqres_to_common_alignment = {}
    seqres_position = 0

    for alignment_position, aa in enumerate(aligned_seq):
        if aa != '-':
            seqres_to_common_alignment[seqres_position] = alignment_position
            seqres_position += 1

    return seqres_to_common_alignment


def mapping_chain_to_seqres(seqres: str, chain: str) -> tuple:
    '''Align a structural chain to its corresponding SEQRES
    and map chain positions to SEQRES positions.'''

    aligned_seqres, aligned_chain = align_seqs(seqres, chain)
    chain_to_seqres = {}

    seqres_position = 0
    chain_position = 0

    for aa_seqres, aa_chain in zip(aligned_seqres, aligned_chain):

        current_seqres_position = None
        current_chain_position = None

        if aa_seqres != '-':
            current_seqres_position = seqres_position
            seqres_position += 1

        if aa_chain != '-':
            current_chain_position = chain_position
            chain_position += 1

        if (current_seqres_position is not None and current_chain_position is not None):
            chain_to_seqres[current_chain_position] = current_seqres_position

    return aligned_seqres, aligned_chain, chain_to_seqres


def project_chain_to_common_alignment(chain: str, chain_to_seqres: dict, seqres_to_common_alignment: dict, alignment_length: int) -> list:
    '''Project a structural chain onto the common SEQRES alignment.'''

    result = ['-'] * alignment_length

    for (chain_position, seqres_position) in chain_to_seqres.items():
        if seqres_position not in seqres_to_common_alignment:
            continue

        alignment_position = (seqres_to_common_alignment[seqres_position])

        result[alignment_position] = (chain[chain_position])

    return result


def process_pair(pdbid_1: str, pdbid_2: str, sequences: dict, chains: dict) -> dict:
    '''Process one pair of structures.
    1. Align two SEQRES sequences.
    2. Map each SEQRES to the common alignment.
    3. Map each chain to its SEQRES.
    4. Project both chains onto the common SEQRES alignment.'''

    seq_1 = sequences[pdbid_1]
    seq_2 = sequences[pdbid_2]

    chain_1 = chains[pdbid_1]
    chain_2 = chains[pdbid_2]

    aligned_seqres_1, aligned_seqres_2 = align_seqs(seq_1, seq_2)

    seqres_1_to_common = (mapping_seqres_to_common_alignment(aligned_seqres_1))
    seqres_2_to_common = (mapping_seqres_to_common_alignment(aligned_seqres_2))

    common_alignment_length = len(aligned_seqres_1)

    chain_1_seqres_alignment, chain_1_chain_alignment, chain_1_to_seqres = mapping_chain_to_seqres(seq_1, chain_1)
    chain_2_seqres_alignment, chain_2_chain_alignment, chain_2_to_seqres = mapping_chain_to_seqres(seq_2,chain_2)

    projected_chain_1 = (project_chain_to_common_alignment(chain_1,
                                                            chain_1_to_seqres,
                                                            seqres_1_to_common,
                                                            common_alignment_length))

    projected_chain_2 = (project_chain_to_common_alignment(chain_2,
                                                            chain_2_to_seqres,
                                                            seqres_2_to_common,
                                                            common_alignment_length))

    return {'chain_alignment': (projected_chain_1,
                                projected_chain_2)}


def strip_alignment(alignment_1, alignment_2) -> tuple:
    '''Remove leading and trailing positions containing only gaps.'''

    assert len(alignment_1) == len(alignment_2)

    N = len(alignment_1)

    begin_n = 0
    end_n = N - 1

    for n in range(N):
        if (alignment_1[n] != '-' or alignment_2[n] != '-'):
            begin_n = n
            break

    for n in reversed(range(N)):
        if (alignment_1[n] != '-' or alignment_2[n] != '-'):
            end_n = n
            break

    return (''.join(alignment_1[begin_n:end_n + 1]),
            ''.join(alignment_2[begin_n:end_n + 1]))


def get_len_longest_shared_region(aligned_chain_1: str, aligned_chain_2: str) -> int:
    '''Find the length of the longest continuous region
    where the two chains contain the same amino acids.'''

    assert len(aligned_chain_1) == len(aligned_chain_2)

    len_shared_region = 0
    len_longest_shared_region = 0

    N = len(aligned_chain_1)

    for n in range(N):
        if aligned_chain_1[n] != '-':
            if aligned_chain_1[n] == aligned_chain_2[n]:
                len_shared_region += 1
                len_longest_shared_region = max(
                    len_longest_shared_region,
                    len_shared_region)
            else:
                len_shared_region = 0

    return len_longest_shared_region


def num_substitutions(seq_1: str, seq_2: str) -> int:
    '''Count amino acid substitutions between equal-length sequences.'''

    assert len(seq_1) == len(seq_2)

    N = len(seq_1)

    count = 0
    for n in range(N):
        if seq_1[n] != seq_2[n]:
            count += 1

    return count


def get_shared_regions(w: int, aligned_chain_1: str, aligned_chain_2: str, positions_1: list, positions_2: list) -> tuple:
    '''Find all identical continuous regions of length w.

    Windows containing gaps are excluded.
    Windows containing amino acid substitutions are excluded.'''

    assert len(aligned_chain_1) == len(aligned_chain_2)

    N = len(aligned_chain_1)

    i_1 = -1
    i_2 = -1

    aa_regions_1 = []
    pos_regions_1 = []

    aa_regions_2 = []
    pos_regions_2 = []

    for i in range(
        6,
        N - w + 1 - 6
    ):

        if aligned_chain_1[i] != '-':
            i_1 += 1

        if aligned_chain_2[i] != '-':
            i_2 += 1

        region_sequence_1 = []
        region_positions_1 = []

        region_sequence_2 = []
        region_positions_2 = []

        for j in range(w):

            if (
                aligned_chain_1[i + j] != '-'
                and aligned_chain_2[i + j] != '-'
            ):

                region_sequence_1 += (
                    aligned_chain_1[i + j]
                )

                region_positions_1.append(
                    positions_1[i_1 + j]
                )

                region_sequence_2 += (
                    aligned_chain_2[i + j]
                )

                region_positions_2.append(
                    positions_2[i_2 + j]
                )

            else:

                break

        if len(region_sequence_1) == w:

            if num_substitutions(
                region_sequence_1,
                region_sequence_2
            ) <= 0:

                aa_regions_1.append(
                    region_sequence_1
                )

                pos_regions_1.append(
                    region_positions_1
                )

                aa_regions_2.append(
                    region_sequence_2
                )

                pos_regions_2.append(
                    region_positions_2
                )

    return (
        aa_regions_1,
        pos_regions_1,
        aa_regions_2,
        pos_regions_2
    )


# ============================================================
# LDDT
# ============================================================

def calc_lddt(
    reference: List[Tuple[float, float, float]],
    model: List[Tuple[float, float, float]],
    inclusion_radius: float = 15.0,
    thresholds: Sequence[float] = (
        0.5,
        1.0,
        2.0,
        4.0
    )
) -> float:
    """
    Calculate the Local Distance Difference Test (LDDT) score
    between two structures represented as lists of corresponding
    3D coordinates.

    Parameters
    ----------
    reference : list of tuples
        Reference coordinates.

    model : list of tuples
        Model coordinates, same order and length as reference.

    inclusion_radius : float
        Include reference pairs whose distance is <= this value.

    thresholds : sequence of floats
        Distance-difference thresholds used by LDDT.

    Returns
    -------
    float
        The LDDT score in the range [0, 1].
    """

    if len(reference) != len(model):

        raise ValueError(
            "reference and model must have the same length"
        )

    if len(reference) < 2:

        raise ValueError(
            "at least two points are required"
        )

    if not thresholds:

        raise ValueError(
            "thresholds must not be empty"
        )

    def distance(
        a: Tuple[float, float, float],
        b: Tuple[float, float, float]
    ) -> float:

        return math.sqrt(
            (a[0] - b[0]) ** 2
            + (a[1] - b[1]) ** 2
            + (a[2] - b[2]) ** 2
        )

    n = len(reference)

    total_score = 0.0
    considered_points = 0

    for i in range(n):

        ref_i = reference[i]
        model_i = model[i]

        neighbor_count = 0
        preserved_count = 0

        for j in range(n):

            if i == j:
                continue

            d_ref = distance(
                ref_i,
                reference[j]
            )

            if d_ref > inclusion_radius:
                continue

            d_model = distance(
                model_i,
                model[j]
            )

            diff = abs(
                d_model - d_ref
            )

            for threshold in thresholds:

                if diff < threshold:
                    preserved_count += 1

            neighbor_count += 1

        # Ignore points with no neighbours
        # inside the inclusion radius.
        if neighbor_count == 0:
            continue

        point_score = (
            preserved_count
            / (neighbor_count * len(thresholds))
        )

        total_score += point_score
        considered_points += 1

    if considered_points == 0:
        return 0.0

    return total_score / considered_points


# ============================================================
# Main
# ============================================================

def main():

    print('Reading families...')

    families = read_families(
        families_file
    )

    print(
        f'Found {len(families)} families.'
    )

    # --------------------------------------------------------
    # Check available structures
    # --------------------------------------------------------

    valid_families = {}

    for family_name, pdb_ids in families.items():

        valid_pdb_ids = []

        for pdbid in pdb_ids:

            seq_file = (
                input_folder
                / f'{pdbid}_seq.txt'
            )

            chain_file = (
                input_folder
                / f'{pdbid}.txt'
            )

            if (
                seq_file.exists()
                and chain_file.exists()
            ):

                valid_pdb_ids.append(
                    pdbid
                )

            else:

                print(
                    f'WARNING: files not found '
                    f'for {pdbid}'
                )

        # Only families with at least
        # two structures are useful.
        if len(valid_pdb_ids) >= 2:

            valid_families[
                family_name
            ] = valid_pdb_ids

    # --------------------------------------------------------
    # Count total pairs
    # --------------------------------------------------------

    total_pairs = sum(
        len(list(combinations(pdb_ids, 2)))
        for pdb_ids in valid_families.values()
    )

    print(
        f'Families with at least two structures: '
        f'{len(valid_families)}'
    )

    print(
        f'Total unique pairs: {total_pairs}'
    )

    # --------------------------------------------------------
    # Cache sequences, chains and coordinates
    # --------------------------------------------------------

    print(
        'Reading structure files...'
    )

    all_pdb_ids = sorted(
        {
            pdbid
            for pdb_ids in valid_families.values()
            for pdbid in pdb_ids
        }
    )

    sequences = {}
    chains = {}
    positions = {}

    for number, pdbid in enumerate(
        all_pdb_ids,
        start=1
    ):

        sequences[pdbid] = read_seq(
            input_folder
            / f'{pdbid}_seq.txt'
        )

        chains[pdbid] = read_chain(
            input_folder
            / f'{pdbid}.txt'
        )

        positions[pdbid] = get_positions(
            pdbid
        )

        print(
            f'Loaded '
            f'[{number}/{len(all_pdb_ids)}] '
            f'{pdbid}'
        )

    # --------------------------------------------------------
    # Calculate LDDT
    # --------------------------------------------------------

    processed_pairs = 0
    successful_pairs = 0

    with open(
        output_file,
        'w',
        encoding='utf-8'
    ) as output:

        for family_name, pdb_ids in (
            valid_families.items()
        ):

            print()

            print(
                f'Processing family: '
                f'{family_name} '
                f'({len(pdb_ids)} structures)'
            )

            output.write(
                f'>{family_name}\n'
            )

            for pdbid_1, pdbid_2 in combinations(
                pdb_ids,
                2
            ):

                processed_pairs += 1

                print(
                    f'[{processed_pairs}/{total_pairs}] '
                    f'{family_name}: '
                    f'{pdbid_1} / {pdbid_2}'
                )

                try:

                    # ------------------------------------------------
                    # Align the two structures
                    # ------------------------------------------------

                    results = process_pair(
                        pdbid_1,
                        pdbid_2,
                        sequences,
                        chains
                    )

                    aligned_chain_1, aligned_chain_2 = (
                        strip_alignment(
                            results[
                                'chain_alignment'
                            ][0],
                            results[
                                'chain_alignment'
                            ][1]
                        )
                    )

                    # ------------------------------------------------
                    # Check longest shared region
                    # ------------------------------------------------

                    longest_shared_region = (
                        get_len_longest_shared_region(
                            aligned_chain_1,
                            aligned_chain_2
                        )
                    )

                    if longest_shared_region < 7:

                        continue

                    # ------------------------------------------------
                    # Find common windows of length 6
                    # ------------------------------------------------

                    (
                        aa_regions_1,
                        pos_regions_1,
                        aa_regions_2,
                        pos_regions_2
                    ) = get_shared_regions(
                        6,
                        aligned_chain_1,
                        aligned_chain_2,
                        positions[pdbid_1],
                        positions[pdbid_2]
                    )

                    # ------------------------------------------------
                    # Calculate LDDT for every window
                    # ------------------------------------------------

                    lddt_values = []

                    for region_number in range(
                        len(pos_regions_1)
                    ):

                        lddt = calc_lddt(
                            pos_regions_1[
                                region_number
                            ],
                            pos_regions_2[
                                region_number
                            ]
                        )

                        lddt_values.append(
                            lddt
                        )

                    # ------------------------------------------------
                    # No results
                    # ------------------------------------------------

                    if not lddt_values:

                        continue

                    successful_pairs += 1

                    # ------------------------------------------------
                    # Write results
                    # ------------------------------------------------

                    output.write(
                        f'{pdbid_1}_{pdbid_2}\n'
                    )

                    output.write(
                        ' '.join(
                            f'{value:.3f}'
                            for value in lddt_values
                        )
                        + '\n'
                    )

                except Exception as e:

                    print(
                        f'ERROR: '
                        f'{family_name}: '
                        f'{pdbid_1} / {pdbid_2}: '
                        f'{e}'
                    )

            output.write('\n')

    # --------------------------------------------------------
    # Finish
    # --------------------------------------------------------

    print()
    print('Finished.')

    print(
        f'Total pairs processed: '
        f'{processed_pairs}'
    )

    print(
        f'Pairs with LDDT results: '
        f'{successful_pairs}'
    )

    print(
        f'Results saved to: '
        f'{output_file}'
    )


if __name__ == '__main__':
    main()