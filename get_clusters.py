""" import os

os.system(r'C:\\Users\\User\\Documents\\bioinf\\smtb\\mmseqs\\mmseqs.bat '
    r'easy-cluster all_seqres_from_cifs.txt cluster_results cluster_tmp '
    r'--min-seq-id 0.3 -c 0.4 --cov-mode 0 --threads 8') """


from pathlib import Path
from collections import defaultdict
import os
from Bio.Align import PairwiseAligner


MMSEQS = Path(r"C:\Users\User\Documents\bioinf\smtb\mmseqs\mmseqs.bat")

seqres_file = Path("all_seqres_from_cifs.txt")

FIRST_OUTPUT = Path("cluster_results")
FIRST_TMP = Path("cluster_tmp")

short_input = Path("short_sequences.txt")
SHORT_OUTPUT = Path("short_clusters")
SHORT_TMP = Path("short_cluster_tmp")
SHORT_CLUSTER_FILE = Path("short_clusters_cluster.tsv")


def sequences_are_similar(name1, name2, threshold):

    seq1 = sequences[name1]
    seq2 = sequences[name2]

    alignments = aligner.align(seq1, seq2)

    if len(alignments) == 0:
        return False

    alignment = alignments[0]
    score = alignment.score

    return score > threshold


def cluster_short_sequences(names, threshold):

    if not names:
        return []

    graph = {name: set() for name in names}

    for i in range(len(names)):
        for j in range(i + 1, len(names)):

            name1 = names[i]
            name2 = names[j]

            seq1 = sequences[name1]
            seq2 = sequences[name2]

            alignment = aligner.align(seq1, seq2)[0]
            score = alignment.score

            print(f"{name1} vs {name2}: "
                  f"score = {score}")

            if score > threshold:
                graph[name1].add(name2)
                graph[name2].add(name1)

    visited = set()
    result_clusters = []

    for name in names:
        if name in visited:
            continue

        cluster = []
        stack = [name]
        visited.add(name)

        while stack:
            current = stack.pop()
            cluster.append(current)

            for neighbor in graph[current]:
                if neighbor not in visited:
                    visited.add(neighbor)
                    stack.append(neighbor)

        result_clusters.append(cluster)

    return result_clusters


sequences = {}

with open(seqres_file, "r", encoding="utf8") as file:

    current_name = None
    for line in file:
        line = line.strip()

        if not line:
            continue

        if line.startswith(">"):
            current_name = line[1:].strip()

        else:
            sequences[current_name] = line

os.system(r'C:\\Users\\User\\Documents\\bioinf\\smtb\\mmseqs\\mmseqs.bat easy-cluster all_seqres_from_cifs.txt cluster_results cluster_tmp --min-seq-id 0.3 -c 0.4 --cov-mode 0 --threads 8')

clusters = defaultdict(list)

with open("cluster_results_cluster.tsv", "r", encoding="utf8") as file:

    for line in file:
        line = line.strip()

        if not line:
            continue

        representative, member = line.split("\t")
        clusters[representative].append(member)

print(f"Number of clusters after first clustering: "
      f"{len(clusters)}")


normal_clusters = []
short_sequences = []

for representative, members in clusters.items():
    short_members = [name for name in members if len(sequences[name]) < 11]
    short_members_for_clustering = []

    for name in short_members:
        sequence = sequences[name]
        if set(sequence) == {"X"}:
            continue

        short_members_for_clustering.append(name)

    if short_members_for_clustering:
        print(f"Short sequences from cluster {representative}: "
              f"{len(short_members_for_clustering)}")
        short_sequences.extend(short_members_for_clustering)

    remaining_members = [name for name in members if name not in short_members_for_clustering]
    if remaining_members:
        normal_clusters.append(remaining_members)


print()
print(f"Short sequences for pairwise clustering: "
      f"{len(short_sequences)}")


aligner = PairwiseAligner(open_gap_score = -3,
                          extend_gap_score = -2,
                          mismatch_score = -1,
                          left_gap_score = -1,
                          right_gap_score = -1)


if short_sequences:
    short_clusters = cluster_short_sequences(short_sequences, 4)
else:
    short_clusters = []


final_clusters = []
final_clusters.extend(normal_clusters)
final_clusters.extend(short_clusters)

with open("clusters_final.txt", "w", encoding="utf8") as file:
    for number, members in enumerate(final_clusters, start=1):
        file.write(f">Cluster {number}\n")

        for name in members:
            file.write(f"{name} {sequences[name]}\n")

        file.write("\n")

print()
print(f"Normal clusters: {len(normal_clusters)}")
print(f"New short-sequence clusters: {len(short_clusters)}")
print(f"Total final clusters: {len(final_clusters)}")