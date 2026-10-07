from typing import List, Tuple
from tmtools import tm_align


def read_families(families_file_path: str) -> List[List]:
    '''Read all families from a families file.'''

    families_file = open(families_file_path, 'r')
    families_lines = families_file.readlines()
    families_file.close()

    families = list()

    for line in families_lines:
        line_content = line.strip()

        if line_content.startswith('>'):
            families.append(list())
        else:
            families[-1].append(line_content)

    return families


def read_chain_coords_seq(ID: str) -> Tuple[List, str]:
    '''Read coordinates of alpha-carbons of a chain and the aa sequence.'''

    chain_file = open('extracted_chains/' + ID + '.txt', 'r')
    lines = chain_file.readlines()
    chain_file.close()

    coords = list()
    seq = ''

    for line in lines:
        tokens = line.strip().split()
        coords.append((float(tokens[0]), float(tokens[1]), float(tokens[2])))
        seq += tokens[5]

    return coords, seq


def calculate_TM_score_from_IDs(ID_1: str, ID_2: str) -> float:
    '''Calculate TM score for amyloid chains with given PDB IDs.'''

    coords_1, seq_1 = read_chain_coords_seq(ID_1)
    coords_2, seq_2 = read_chain_coords_seq(ID_2)

    # Perform TM-align and get results.
    tm_alignment = tm_align(coords_1, coords_2, seq_1, seq_2)

    # The TM score (normalized by the length of the first structure by default).
    tm_score = min(tm_alignment.tm_norm_chain1, tm_alignment.tm_norm_chain2)

    return tm_score


families = read_families('clusters_renamed.txt')

print(len(families), 'families read.')

TM_scores_file = open('TM_scores.txt', 'w')

cnt = 0
for family in families:
    cnt += 1
    print('#### Family ' + str(cnt) + '. ####')
    print('#### ' + str(len(family)) + ' structures. ####')

    family_size = len(family)

    for i in range(family_size):
        print(i + 1, '/', len(family))
        for j in range(i+1, family_size):
            ID_i, ID_j = family[i], family[j]
            TM_scores_file.write(
                str(ID_i) + ' ' + str(ID_j) + ' ' + \
                str(calculate_TM_score_from_IDs(ID_i, ID_j)) + '\n'
            )

TM_scores_file.close()
