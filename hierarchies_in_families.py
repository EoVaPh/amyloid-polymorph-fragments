from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
from scipy.spatial.distance import squareform
import numpy as np
from matplotlib import pyplot as plot
import matplotlib
from pathlib import Path

matplotlib.rcParams['figure.dpi'] = 300
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rc('font', family='STIXGeneral')
matplotlib.rc('font', weight='ultralight')

rmsd_file = open('all_rmsd.txt', 'r')
rmsd_lines = rmsd_file.readlines()
rmsd_file.close()

def cluster(family_name: str, distances: list):
    # Extract all unique object names
    objects = sorted(set([row[0] for row in distances] + \
                         [row[1] for row in distances]))
    n = len(objects)
    obj_to_idx = {obj: i for i, obj in enumerate(objects)}

    # Initialise an n x n matrix (diagonal stays 0)
    D = np.zeros((n, n))

    for name1, name2, dist in distances:
        i = obj_to_idx[name1]
        j = obj_to_idx[name2]
        D[i, j] = dist
        D[j, i] = dist   # assuming distances are symmetric; remove if not

    condensed = squareform(D)

    Z = linkage(condensed, method='average')

    # Plot dendrogram to decide a cut height
    plot.figure(figsize=(10, 5))
    dendrogram(Z, labels=objects, leaf_rotation=90)
    plot.ylabel("RMSD (Å)")
    plot.tight_layout()
    plot.savefig('figures/dendrogram_' + family_name.replace('/', '_') + '.png')
    plot.close()

    # Cut at your desired distance threshold
    max_dist = 0.2 * max(Z[:, 2])
    cluster_labels = fcluster(Z, t=max_dist, criterion='distance')

    # Group names by cluster label
    clusters = {}
    for name, lbl in zip(objects, cluster_labels):
        clusters.setdefault(lbl, []).append(name)

    representatives = []

    # Print the sets of names
    for cid, members in clusters.items():
        print(f"Cluster {cid}: {set(members)}")
        cluster_ids = set(members)
        longest_id, max_length = None, 0
        for id in cluster_ids:
            if (sequences_path / (id + '.pdb_sequences.fa')).is_file():
                ext = 'pdb'
            else:
                ext = 'cif'
            length = len(open(sequences_path / (id + '.' + ext + '_sequences.fa'), 'r').readlines()[3].strip())
            if length > max_length:
                max_length = length
                longest_id = id
        representatives.append(longest_id)

    return representatives


representatives_file = open('families_representatives.txt', 'w')

sequences_path = Path('sequences')

distances = []
for line in rmsd_lines:
    if '>' in line:
        if len(distances) > 0:
            print(family_name)

            family_representatives = cluster(family_name, distances)
            representatives_file.write('>' + family_name)

            for representative in family_representatives:
                representatives_file.write(representative + '\n')

        distances = []
        family_name = line.strip()[1:]
    else:
        tokens = line.split()
        tokens[2] = float(tokens[2])
        distances.append(tokens)

# Process the last family.
if len(distances) > 0:
    print(family_name)

    family_representatives = cluster(family_name, distances)
    representatives_file.write('>' + family_name + '\n')

    for representative in family_representatives:
        representatives_file.write(representative + '\n')

representatives_file.close()
