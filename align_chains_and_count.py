from pathlib import Path
import numpy as np
from scipy.spatial.transform import Rotation

from Bio.Align import PairwiseAligner

import math
from typing import List, Tuple, Sequence

from scipy.stats import linregress

from matplotlib import pyplot as plot
import matplotlib

matplotlib.rcParams['figure.dpi'] = 300
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rc('font', family='STIXGeneral')
matplotlib.rc('font', weight='ultralight')


def read_seq(seq_file_path: str) -> str:
    '''Read a SEQRES sequence from a text file.'''

    seq_file = open(seq_file_path, 'r')
    seq = seq_file.read().strip()
    seq_file.close()

    return seq


def read_chain(chain_file_path: str) -> str:
    '''Read a structural chain from a text file.'''

    chain_file = open(chain_file_path, 'r')
    residues = chain_file.readlines()
    chain_file.close()

    seq = ''

    for residue in residues:
        resnum = int(residue.split()[-2])
        seq += residue.split()[-1]

    return seq


def read_families(file_path: str) -> dict:
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

    aligner = PairwiseAligner(open_gap_score = -3,
                              extend_gap_score = -2,
                              mismatch_score = -1,
                              left_gap_score = -1,
                              right_gap_score = -1)
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
    '''Align a structural chain to its corresponding SEQRES sequence
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
            chain_to_seqres[current_chain_position] = (current_seqres_position)

    return aligned_seqres, aligned_chain, chain_to_seqres


def project_chain_to_common_alignment(chain: str, chain_to_seqres: dict, seqres_to_common_alignment: dict, alignment_length: int) -> list:
    '''Project a structural chain onto the common SEQRES alignment.'''

    result = ['-'] * alignment_length

    for chain_position, seqres_position in chain_to_seqres.items():
        if seqres_position not in seqres_to_common_alignment:
            continue

        alignment_position = seqres_to_common_alignment[seqres_position]

        result[alignment_position] = chain[chain_position]

    return result


def process_pair(pdbid_1: str, pdbid_2: str) -> dict:
    '''Process a pair of structures:
    1. Align two SEQRES sequences
    2. Map each SEQRES to the common alignment
    3. Map each chain to it's SEQRES
    4. Project both chains onto the common SEQRES alignment'''

    seq_1 = read_seq('extracted_chains/' + pdbid_1 + '_seq.txt')
    seq_2 = read_seq('extracted_chains/' + pdbid_2 + '_seq.txt')

    chain_1 = read_chain('extracted_chains/' + pdbid_1 + '.txt')
    chain_2 = read_chain('extracted_chains/' + pdbid_2 + '.txt')

    aligned_seqres_1, aligned_seqres_2 = align_seqs(seq_1, seq_2)

    seqres_1_to_common = mapping_seqres_to_common_alignment(aligned_seqres_1)
    seqres_2_to_common = mapping_seqres_to_common_alignment(aligned_seqres_2)
    common_alignment_length = len(aligned_seqres_1)

    chain_1_seqres_alignment, chain_1_chain_alignment, chain_1_to_seqres = mapping_chain_to_seqres(seq_1, chain_1)
    chain_2_seqres_alignment, chain_2_chain_alignment, chain_2_to_seqres = mapping_chain_to_seqres(seq_2, chain_2)

    projected_chain_1 = project_chain_to_common_alignment(chain_1, chain_1_to_seqres, seqres_1_to_common, common_alignment_length)
    projected_chain_2 = project_chain_to_common_alignment(chain_2, chain_2_to_seqres, seqres_2_to_common, common_alignment_length)

    return {
        'seqres_alignment': (
            aligned_seqres_1,
            aligned_seqres_2
        ),

        'chain_1': (
            chain_1_seqres_alignment,
            chain_1_chain_alignment
        ),

        'chain_2': (
            chain_2_seqres_alignment,
            chain_2_chain_alignment
        ),

        'chain_alignment': (
            projected_chain_1,
            projected_chain_2
        ),

        'chain_1_to_seqres': chain_1_to_seqres,
        'chain_2_to_seqres': chain_2_to_seqres
    }


def strip_alignment(alignment_1, alignment_2) -> tuple:
    assert len(alignment_1) == len(alignment_2)

    N = len(alignment_1)

    for n in range(0, N):
        if alignment_1[n] != '-' or alignment_2[n] != '-':
            begin_n = n
            break

    for n in reversed(range(0, N)):
        if alignment_1[n] != '-' or alignment_2[n] != '-':
            end_n = n
            break

    return ''.join(alignment_1[begin_n:end_n+1]),\
           ''.join(alignment_2[begin_n:end_n+1])


def get_len_longest_shared_region(aligned_chain_1: str,
                                  aligned_chain_2: str) -> int:
    assert len(aligned_chain_1) == len(aligned_chain_2)

    len_shared_region, len_longest_shared_region = 0, 0

    N = len(aligned_chain_1)

    for n in range(N):
        if aligned_chain_1[n] != '-':
            if aligned_chain_1[n] == aligned_chain_2[n]:
                len_shared_region += 1
                len_longest_shared_region = max(len_longest_shared_region,
                                                len_shared_region)
            else:
                len_shared_region = 0

    return len_longest_shared_region


def get_positions(pdbid: str) -> list:
    chain_file = open('extracted_chains/' + pdbid + '.txt', 'r')
    residues = chain_file.readlines()
    chain_file.close()

    positions = []

    for residue in residues:
        residue_tokens = residue.strip().split()
        positions.append((
            float(residue_tokens[0]),
            float(residue_tokens[1]),
            float(residue_tokens[2])
        ))

    return positions


def num_substitutions(seq_1: str, seq_2: str) -> int:
    assert len(seq_1) == len(seq_2)

    N = len(seq_1)

    count = 0

    for n in range(N):
        if seq_1[n] != seq_2[n]:
            count += 1

    return count


def get_shared_regions(w: int, aligned_chain_1: str, aligned_chain_2: str,
                               positions_1: list, positions_2: list) -> tuple:
    '''Find all shared continuous regions of the same length and return their
       sequences and coordinates.'''

    assert len(aligned_chain_1) == len(aligned_chain_2)

    N = len(aligned_chain_1)

    i_1, i_2 = -1, -1

    aa_regions_1, pos_regions_1, aa_regions_2, pos_regions_2 = [], [], [], []

    for i in range(6, N - w + 1 - 6):
        if aligned_chain_1[i] != '-':
            i_1 += 1

        if aligned_chain_2[i] != '-':
            i_2 += 1

        region_sequence_1, region_positions_1 = [], []
        region_sequence_2, region_positions_2 = [], []

        for j in range(w):
            if aligned_chain_1[i + j] != '-' and aligned_chain_2[i + j] != '-':
                region_sequence_1 += aligned_chain_1[i + j]
                region_positions_1.append(positions_1[i_1 + j])
                region_sequence_2 += aligned_chain_2[i + j]
                region_positions_2.append(positions_2[i_2 + j])
            else:
                break

        if len(region_sequence_1) == w:
            if num_substitutions(region_sequence_1, region_sequence_2) <= 0:
                aa_regions_1.append(region_sequence_1)
                pos_regions_1.append(region_positions_1)
                aa_regions_2.append(region_sequence_2)
                pos_regions_2.append(region_positions_2)

    return aa_regions_1, pos_regions_1, aa_regions_2, pos_regions_2



def normalize_pair(
    pdbid_1: str,
    pdbid_2: str
) -> tuple:
    """
    Make pair independent of structure order.

    12gb_2lmn and 2lmn_12gb become the same pair.
    """

    return tuple(sorted((
        pdbid_1.lower(),
        pdbid_2.lower()
    )))


def read_pair_values_file(file_path: str) -> dict:
    """
    Read file of the form:

    >Family
    pdbid1_pdbid2
    value1 value2 value3 ...

    Returns:

    {
        ('pdbid1', 'pdbid2'): [value1, value2, ...],
        ...
    }
    """

    values_by_pair = {}

    with open(file_path, 'r', encoding='utf-8') as file:
        lines = [
            line.strip()
            for line in file
            if line.strip()
        ]

    i = 0

    while i < len(lines):

        if lines[i].startswith('>'):
            i += 1
            continue

        pair_name = lines[i]

        if i + 1 >= len(lines):
            raise ValueError(
                f'No values found for pair {pair_name}'
            )

        values = [
            float(x)
            for x in lines[i + 1].split()
        ]

        pdbid_1, pdbid_2 = pair_name.split('_', 1)

        pair_key = normalize_pair(
            pdbid_1,
            pdbid_2
        )

        values_by_pair[pair_key] = values

        i += 2

    return values_by_pair


def read_seqres_idr(pdb_id: str, idr_file: str) -> list:
    '''Read pre-calculated IDR propensities for a PDB ID from a file.'''

    with open(idr_file, 'r') as file:
        lines = file.readlines()

    for i, line in enumerate(lines):
        if line.strip() == f'>{pdb_id}':
            if i + 1 >= len(lines):
                raise ValueError(f'No IDR values found for {pdb_id}')

            return [float(x) for x in lines[i + 1].split()]

    raise ValueError(f'IDR values not found for {pdb_id}')


def select_atomseq_idr(seqres_aligned: str, chain_aligned: str, aiupred_values: list) -> list:
    '''Having pre-calculated IDR propensities for a SEQRES, transfer them to the
    chain in accordance with SEQRES-chain alignment.'''

    assert len(seqres_aligned) == len(chain_aligned)

    IDR = []

    seqres_index = 0

    for i in range(len(seqres_aligned)):
        if seqres_aligned[i] != "-":
            p = aiupred_values[seqres_index]
            seqres_index += 1
        if chain_aligned[i] != "-":
            IDR.append(p)

    return IDR


def find_mismatches(aligned_chain_1: str, aligned_chain_2: str) -> str:
    assert len(aligned_chain_1) == len(aligned_chain_2)

    N = len(aligned_chain_1)

    mismatches = ''

    for n in range(N):
        if aligned_chain_1[n] != '-' and aligned_chain_2[n] != '-' and \
           aligned_chain_1[n] != aligned_chain_2[n]:
                mismatches += '*'
        else:
            mismatches += ' '

    return mismatches



families = read_families('amyloid_explorer_families.txt')

lddt_data = read_pair_values_file('lddts.txt')
rmsd_data = read_pair_values_file('rmsds.txt')

verbose = False

OUTPUT_DIR = Path('idr_lddt_rmsd_graphics')
OUTPUT_DIR.mkdir(exist_ok=True)

for family_number, (family_name, pdb_ids) in enumerate(
    families.items(),
    start=1
):

    print()
    print('=' * 80)
    print(
        f'[{family_number}/{len(families)}] {family_name}'
    )
    print(
        f'Structures: {len(pdb_ids)}'
    )
    print('=' * 80)

    # --------------------------------------------------------
    # Arrays for current family
    # --------------------------------------------------------

    rmsd_length_6 = []
    lddt_length_6 = []
    idr_length_6 = []

    # --------------------------------------------------------
    # All unique pairs in this family
    # --------------------------------------------------------

    total_pairs = len(pdb_ids) * (len(pdb_ids) - 1) // 2
    pair_counter = 0

    for pdbid_1 in pdb_ids:

        for pdbid_2 in pdb_ids:

            if pdbid_1 >= pdbid_2:
                continue

            pair_counter += 1

            print(
                f'[{pair_counter}/{total_pairs}] '
                f'{pdbid_1} / {pdbid_2}'
            )

            # ------------------------------------------------
            # Get precalculated LDDT and RMSD
            # ------------------------------------------------

            pair_key = normalize_pair(
                pdbid_1,
                pdbid_2
            )

            if pair_key not in lddt_data:
                print(
                    f'  WARNING: LDDT not found for '
                    f'{pdbid_1} / {pdbid_2}'
                )
                continue

            if pair_key not in rmsd_data:
                print(
                    f'  WARNING: RMSD not found for '
                    f'{pdbid_1} / {pdbid_2}'
                )
                continue

            pair_lddt = lddt_data[pair_key]
            pair_rmsd = rmsd_data[pair_key]

            # ------------------------------------------------
            # Alignment
            # ------------------------------------------------

            results = process_pair(
                pdbid_1,
                pdbid_2
            )

            aligned_chain_1, aligned_chain_2 = strip_alignment(
                results['chain_alignment'][0],
                results['chain_alignment'][1]
            )

            len_longest_shared_region = \
                get_len_longest_shared_region(
                    aligned_chain_1,
                    aligned_chain_2
                )

            if len_longest_shared_region < 7:
                continue

            # ------------------------------------------------
            # Verbose output
            # ------------------------------------------------

            if verbose:

                print("SEQRES alignment")
                print()
                print(results['seqres_alignment'][0])
                print(results['seqres_alignment'][1])
                print()

                print("CHAIN 1 -> SEQRES 1")
                print()
                print(results['chain_1'][0])
                print(results['chain_1'][1])
                print()

                print("CHAIN 2 -> SEQRES 2")
                print()
                print(results['chain_2'][0])
                print(results['chain_2'][1])
                print()

                print("FINAL CHAIN ALIGNMENT")

            print()
            print(aligned_chain_1)
            print(aligned_chain_2)
            print(
                find_mismatches(
                    aligned_chain_1,
                    aligned_chain_2
                )
            )
            print()

            # ------------------------------------------------
            # Coordinates
            # ------------------------------------------------

            positions_1 = get_positions(pdbid_1)
            positions_2 = get_positions(pdbid_2)

            # ------------------------------------------------
            # Shared windows
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
                positions_1,
                positions_2
            )

            # ------------------------------------------------
            # Check number of windows
            # ------------------------------------------------

            if len(pair_lddt) != len(pos_regions_1):

                print(
                    f'WARNING: LDDT/window mismatch: '
                    f'{len(pair_lddt)} vs '
                    f'{len(pos_regions_1)}'
                )

                continue

            if len(pair_rmsd) != len(pos_regions_1):

                print(
                    f'WARNING: RMSD/window mismatch: '
                    f'{len(pair_rmsd)} vs '
                    f'{len(pos_regions_1)}'
                )

                continue

            # ------------------------------------------------
            # IDR
            # ------------------------------------------------

            seqres_idr_1 = read_seqres_idr(
                pdbid_1,
                'idrs.txt'
            )

            seqres_idr_2 = read_seqres_idr(
                pdbid_2,
                'idrs.txt'
            )

            chain_idr_1 = select_atomseq_idr(
                results['chain_1'][0],
                results['chain_1'][1],
                seqres_idr_1
            )

            chain_idr_2 = select_atomseq_idr(
                results['chain_2'][0],
                results['chain_2'][1],
                seqres_idr_2
            )

            idr_by_position_1 = dict(
                zip(
                    positions_1,
                    chain_idr_1
                )
            )

            idr_by_position_2 = dict(
                zip(
                    positions_2,
                    chain_idr_2
                )
            )

            # ------------------------------------------------
            # Process windows
            # ------------------------------------------------

            for r in range(
                len(pos_regions_1)
            ):

                # --------------------------------------------
                # RMSD from file
                # --------------------------------------------

                rmsd_length_6.append(
                    pair_rmsd[r]
                )

                # --------------------------------------------
                # LDDT from file
                # --------------------------------------------

                lddt_length_6.append(
                    1 - pair_lddt[r]
                )

                # --------------------------------------------
                # IDR
                # --------------------------------------------

                window_pair_idr = []

                for pos_1, pos_2 in zip(
                    pos_regions_1[r],
                    pos_regions_2[r]
                ):

                    idr_1 = idr_by_position_1[pos_1]
                    idr_2 = idr_by_position_2[pos_2]

                    mean_idr_residue = (
                        idr_1 + idr_2
                    ) / 2

                    window_pair_idr.append(
                        mean_idr_residue
                    )

                mean_idr_window = np.mean(
                    window_pair_idr
                )

                idr_length_6.append(
                    mean_idr_window
                )

    # ========================================================
    # PLOTS FOR CURRENT FAMILY
    # ========================================================

    if len(idr_length_6) < 2:

        print(
            'Not enough data to build graphs.'
        )

        continue

    # --------------------------------------------------------
    # Family filename
    # --------------------------------------------------------

    safe_family_name = ''.join(
        c if c.isalnum() or c in '-_.'
        else '_'
        for c in family_name
    )

    # ========================================================
    # IDR vs LDDT
    # ========================================================

    slope, intercept, r, p, se = linregress(
        idr_length_6,
        lddt_length_6
    )

    idr_range = np.arange(
        min(idr_length_6),
        max(idr_length_6),
        0.01
    )

    figure = plot.figure(
        figsize=(8, 6)
    )

    plot.plot(
        idr_range,
        intercept + slope * idr_range,
        color='#450920'
    )

    plot.scatter(
        idr_length_6,
        lddt_length_6,
        color='#a53860',
        alpha=0.1
    )

    plot.xticks(fontsize=12)
    plot.yticks(fontsize=12)

    plot.xlabel(
        'IDR',
        fontsize=16
    )

    plot.ylabel(
        '1 – LDDT',
        fontsize=16
    )

    plot.title(
        f'{family_name} family screened using window of length 6',
        fontsize=16
    )

    plot.tight_layout()

    plot.savefig(
        OUTPUT_DIR / f'{safe_family_name}_idr_lddt.png',
        dpi=300
    )

    plot.close(figure)

    # ========================================================
    # IDR vs RMSD
    # ========================================================

    slope, intercept, r, p, se = linregress(
        idr_length_6,
        rmsd_length_6
    )

    idr_range = np.arange(
        min(idr_length_6),
        max(idr_length_6),
        0.01
    )

    figure = plot.figure(
        figsize=(8, 6)
    )

    plot.plot(
        idr_range,
        intercept + slope * idr_range,
        color='#450920'
    )

    plot.scatter(
        idr_length_6,
        rmsd_length_6,
        color='#a53860',
        alpha=0.1
    )

    plot.xticks(fontsize=12)
    plot.yticks(fontsize=12)

    plot.xlabel(
        'IDR',
        fontsize=16
    )

    plot.ylabel(
        'RMSD',
        fontsize=16
    )

    plot.title(
        f'{family_name} family screened using window of length 6',
        fontsize=16
    )

    plot.tight_layout()

    plot.savefig(
        OUTPUT_DIR / f'{safe_family_name}_idr_rmsd.png',
        dpi=300
    )

    plot.close(figure)

    print(
        f'Saved graphs for {family_name}'
    )