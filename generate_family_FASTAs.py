"""Generate one FASTA file per family from its structure chains.

This module reads family definitions and a FASTA file of chain sequences. For
each family, it selects records whose structure identifiers match the family
and writes the selected records to a separate FASTA file.
"""

from __future__ import annotations

from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def parse_families(families_file: Path) -> dict[str, list[str]]:
    """Parse family names and their structure identifiers.

    Family names are given on lines that start with ``>``. Each subsequent
    non-empty line is treated as a structure identifier, and only its first
    four characters are kept.

    :param families_file: Path to the family definition file.
    :type families_file: pathlib.Path
    :return: Mapping from family names to structure identifiers.
    :rtype: dict[str, list[str]]
    :raises ValueError: If a structure identifier appears before a family name.
    """
    families: dict[str, list[str]] = {}
    current_family: str | None = None

    with families_file.open("r", encoding="utf-8") as file:
        for line in file:
            line = line.strip()

            if not line:
                continue

            if line.startswith(">"):
                current_family = line[1:].strip()
                families[current_family] = []
            else:
                if current_family is None:
                    raise ValueError(
                        "Structure identifier found before family name."
                    )
                families[current_family].append(line[:4])

    return families


def load_fasta_records(fasta_file: Path) -> list[SeqRecord]:
    """Load all FASTA records from a file.

    :param fasta_file: Path to the FASTA file containing chain sequences.
    :type fasta_file: pathlib.Path
    :return: FASTA records in the order in which they appear in the file.
    :rtype: list[Bio.SeqRecord.SeqRecord]
    """
    return list(SeqIO.parse(fasta_file, "fasta"))


def select_records_for_family(
    records: list[SeqRecord],
    structures: set[str],
) -> list[SeqRecord]:
    """Select records whose identifiers belong to a family.

    A record is selected when the first four characters of its identifier are
    present in ``structures``.

    :param records: All FASTA records to filter.
    :type records: list[Bio.SeqRecord.SeqRecord]
    :param structures: Structure identifiers associated with the family.
    :type structures: set[str]
    :return: Records that belong to the family.
    :rtype: list[Bio.SeqRecord.SeqRecord]
    """
    selected_records: list[SeqRecord] = []

    for record in records:
        structure_name = record.id[:4]

        if structure_name in structures:
            selected_records.append(record)

    return selected_records


def write_family_fasta(
    family_name: str,
    records: list[SeqRecord],
    output_folder: Path,
) -> Path:
    """Write selected records to a family-specific FASTA file.

    :param family_name: Name of the family.
    :type family_name: str
    :param records: Records to write.
    :type records: list[Bio.SeqRecord.SeqRecord]
    :param output_folder: Directory where the FASTA file will be created.
    :type output_folder: pathlib.Path
    :return: Path to the written FASTA file.
    :rtype: pathlib.Path
    """
    safe_filename = family_name.replace("/", "_")
    output_file = output_folder / f"{safe_filename}.fasta"
    SeqIO.write(records, output_file, "fasta")
    return output_file


def main() -> None:
    """Generate one FASTA file per family.

    :return: ``None``.
    :rtype: None
    """
    fasta_file = Path("chains.txt")
    families_file = Path("clusters_renamed.txt")
    output_folder = Path("family_FASTAs")
    output_folder.mkdir(exist_ok=True)

    families = parse_families(families_file)
    records = load_fasta_records(fasta_file)

    for number, (family_name, structures) in enumerate(
        families.items(), start=1
    ):
        structure_set = set(structures)
        selected_records = select_records_for_family(
            records, structure_set
        )
        write_family_fasta(family_name, selected_records, output_folder)

        print(
            f"[{number}/{len(families)}] "
            f"{family_name}: {len(selected_records)} sequences"
        )


if __name__ == "__main__":
    main()
