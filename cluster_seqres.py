"""Cluster sequences with MMseqs2 and recluster short sequences.

This script reads sequences from ``seqres.txt``, runs MMseqs2 easy-cluster,
separates clusters containing only short sequences, reclusters those short
sequences with stricter identity and coverage thresholds, and writes the final
cluster assignments to ``clusters_final.txt``.
"""

from __future__ import annotations

import subprocess
from collections import defaultdict
from pathlib import Path


def read_sequences(seqres_path: Path) -> dict[str, str]:
    """Read sequences from ``seqres_path``.

    :param seqres_path: Path to the input sequence file.
    :type seqres_path: pathlib.Path
    :return: Mapping from sequence name to sequence.
    :rtype: dict[str, str]
    """
    sequences = {}
    with seqres_path.open("r", encoding="utf8") as file:
        current_name = None
        for line in file:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                current_name = line[1:].strip()
            else:
                sequences[current_name] = line

    return sequences


def run_mmseqs_cluster(
    input_fasta: Path,
    output_prefix: Path,
    tmp_prefix: Path,
    min_seq_id: float,
    coverage: float,
    cov_mode: int = 0,
    threads: int = 8,
) -> None:
    """Run MMseqs2 ``easy-cluster``.

    :param input_fasta: Path to the input FASTA file.
    :type input_fasta: pathlib.Path
    :param output_prefix: Prefix for MMseqs2 output files.
    :type output_prefix: pathlib.Path
    :param tmp_prefix: Prefix for MMseqs2 temporary files.
    :type tmp_prefix: pathlib.Path
    :param min_seq_id: Minimum sequence identity threshold.
    :type min_seq_id: float
    :param coverage: Coverage threshold.
    :type coverage: float
    :param cov_mode: Coverage mode.
    :type cov_mode: int
    :param threads: Number of CPU threads.
    :type threads: int
    :return: None
    :rtype: None
    """
    subprocess.run(
        [
            "mmseqs",
            "easy-cluster",
            str(input_fasta),
            str(output_prefix),
            str(tmp_prefix),
            "--min-seq-id",
            str(min_seq_id),
            "-c",
            str(coverage),
            "--cov-mode",
            str(cov_mode),
            "--threads",
            str(threads),
        ],
        check=True,
    )


def read_clusters(cluster_tsv: Path) -> dict[str, list[str]]:
    """Read MMseqs2 cluster TSV output.

    :param cluster_tsv: Path to the MMseqs2 cluster TSV file.
    :type cluster_tsv: pathlib.Path
    :return: Mapping from representative sequence name to cluster members.
    :rtype: dict[str, list[str]]
    """
    clusters = defaultdict(list)
    with cluster_tsv.open("r", encoding="utf8") as file:
        for line in file:
            line = line.strip()

            if not line:
                continue

            representative, member = line.split("\t")
            clusters[representative].append(member)

    return clusters


def write_short_sequences(
    names: list[str],
    sequences: dict[str, str],
    output_path: Path,
) -> None:
    """Write selected sequences to a FASTA file.

    :param names: Sequence names to write.
    :type names: list[str]
    :param sequences: Mapping from sequence name to sequence.
    :type sequences: dict[str, str]
    :param output_path: Path to the output FASTA file.
    :type output_path: pathlib.Path
    :return: None
    :rtype: None
    """
    with output_path.open("w", encoding="utf8") as file:
        for name in names:
            file.write(f">{name}\n")
            file.write(f"{sequences[name]}\n")


def write_final_clusters(
    clusters: list[list[str]],
    output_path: Path,
) -> None:
    """Write final cluster assignments.

    :param clusters: Final cluster assignments.
    :type clusters: list[list[str]]
    :param output_path: Path to the output cluster file.
    :type output_path: pathlib.Path
    :return: None
    :rtype: None
    """
    with output_path.open("w", encoding="utf8") as file:
        for number, members in enumerate(clusters):
            file.write(f">Cluster {number}\n")
            for name in members:
                file.write(f"{name}\n")
            file.write("\n")


def main() -> None:
    """Run the clustering workflow.

    .. note::
       Clusters in which every sequence is shorter than 11 residues are
       reclustered with higher identity and coverage thresholds.

    :return: None
    :rtype: None
    """
    sequences = read_sequences(Path("seqres.txt"))

    run_mmseqs_cluster(
        input_fasta=Path("seqres.txt"),
        output_prefix=Path("cluster_results"),
        tmp_prefix=Path("cluster_tmp"),
        min_seq_id=0.3,
        coverage=0.4,
    )

    clusters = read_clusters(Path("cluster_results_cluster.tsv"))

    print(f"Number of clusters after first clustering: {len(clusters)}")

    normal_clusters = []
    short_sequences = []

    for representative, members in clusters.items():
        lengths = [len(sequences[name]) for name in members]

        if all(length < 11 for length in lengths):
            print(
                f"Short cluster: {representative} "
                f"({len(members)} sequences)"
            )
            short_sequences.extend(members)
        else:
            normal_clusters.append(members)

    print()
    print(f"Short sequences for re-clustering: {len(short_sequences)}")

    if short_sequences:
        write_short_sequences(
            short_sequences,
            sequences,
            Path("short_sequences.txt"),
        )

        print()
        print("Created file: short_sequences.txt")

        run_mmseqs_cluster(
            input_fasta=Path("short_sequences.txt"),
            output_prefix=Path("short_clusters"),
            tmp_prefix=Path("short_cluster_tmp"),
            min_seq_id=0.5,
            coverage=0.5,
        )

        short_clusters = read_clusters(Path("short_clusters_cluster.tsv"))
    else:
        short_clusters = {}

    final_clusters = []
    final_clusters.extend(normal_clusters)
    final_clusters.extend(short_clusters.values())

    write_final_clusters(final_clusters, Path("clusters_final.txt"))

    print()
    print(f"Normal clusters: {len(normal_clusters)}")
    print(f"New short-sequence clusters: {len(short_clusters)}")
    print(f"Total final clusters: {len(final_clusters)}")


if __name__ == "__main__":
    main()
