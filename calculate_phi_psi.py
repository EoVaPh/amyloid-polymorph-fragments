"""Calculate phi and psi dihedral angles from mmCIF files."""

import traceback
from pathlib import Path

from Bio.PDB import MMCIFParser
from Bio.PDB.vectors import calc_dihedral


AMINO_ACIDS = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def get_amino_acid(code):
    """Return the one-letter code for a three-letter amino acid code.

    :param code: Three-letter amino acid code.
    :type code: str
    :return: One-letter amino acid code, or ``"X"`` if unknown.
    :rtype: str
    """
    return AMINO_ACIDS.get(code.upper(), "X")


def get_phi_psi(filename, output_file, chain_id):
    """Calculate phi and psi dihedral angles for a chain.

    The function reads a structure from an mmCIF file and appends the
    calculated angles to a common output file. Alternative residue
    locations are handled by selecting the first available residue
    variant. Alternative atom locations are also handled by selecting
    the first available atom variant.

    Non-canonical amino acids are represented as ``"X"``. Residues
    without a CA atom are skipped.

    :param filename: Path to the input mmCIF file.
    :type filename: pathlib.Path
    :param output_file: Path to the common output file.
    :type output_file: pathlib.Path
    :param chain_id: ID of the chain for which angles are calculated.
    :type chain_id: str
    :raises ValueError: If the file format is unsupported or the chain
        is not found.
    :return: ``None``.
    :rtype: None
    """
    filename = Path(filename)
    output_file = Path(output_file)

    if filename.suffix.lower() in [".cif", ".mmcif"]:
        parser = MMCIFParser(QUIET=True, auth_chains=False)
    else:
        raise ValueError(f"Unsupported file format: {filename}")

    structure = parser.get_structure("structure", filename)
    model = next(structure.get_models())

    if chain_id not in model:
        raise ValueError(
            f"Chain '{chain_id}' not found in {filename.name}"
        )
    chain = model[chain_id]

    print(f"Processing: {filename.name}")
    print(f"Chain: {chain.id}")

    def get_residue_variant(residue):
        """Return the first available variant of a disordered residue.

        :param residue: Residue to inspect.
        :type residue: Bio.PDB.Residue.Residue
        :return: Selected residue variant, or ``None`` if none exists.
        :rtype: Bio.PDB.Residue.Residue or None
        """
        if residue.is_disordered() == 2:
            variants = list(residue.child_dict.values())
            if variants:
                return variants[0]
            return None

        return residue

    def get_atom(residue, atom_name):
        """Return an atom from a residue, selecting the first variant.

        :param residue: Residue to inspect.
        :type residue: Bio.PDB.Residue.Residue or None
        :param atom_name: Name of the atom to retrieve.
        :type atom_name: str
        :return: Selected atom, or ``None`` if it is unavailable.
        :rtype: Bio.PDB.Atom.Atom or None
        """
        if residue is None:
            return None
        if atom_name not in residue:
            return None

        atom = residue[atom_name]

        if atom.is_disordered():
            variants = list(atom.child_dict.values())
            if variants:
                return variants[0]
            return None

        return atom

    residues = []

    for residue in chain.get_residues():
        selected_residue = get_residue_variant(residue)

        if selected_residue is None:
            continue
        if get_atom(selected_residue, "CA") is None:
            continue

        residues.append(selected_residue)

    with open(output_file, "a", encoding="utf-8") as file:
        file.write(f">{filename.name}\n")

        for i, residue in enumerate(residues):
            amino_acid = get_amino_acid(residue.resname)
            residue_number = residue.id[1]

            phi = None
            psi = None

            if i > 0:
                previous_residue = residues[i - 1]
                previous_c = get_atom(previous_residue, "C")
                current_n = get_atom(residue, "N")
                current_ca = get_atom(residue, "CA")
                current_c = get_atom(residue, "C")

                if all(
                    [
                        previous_c is not None,
                        current_n is not None,
                        current_ca is not None,
                        current_c is not None,
                    ]
                ):
                    phi = calc_dihedral(
                        previous_c.get_vector(),
                        current_n.get_vector(),
                        current_ca.get_vector(),
                        current_c.get_vector(),
                    )

            if i < len(residues) - 1:
                next_residue = residues[i + 1]
                current_n = get_atom(residue, "N")
                current_ca = get_atom(residue, "CA")
                current_c = get_atom(residue, "C")
                next_n = get_atom(next_residue, "N")

                if all(
                    [
                        current_n is not None,
                        current_ca is not None,
                        current_c is not None,
                        next_n is not None,
                    ]
                ):
                    psi = calc_dihedral(
                        current_n.get_vector(),
                        current_ca.get_vector(),
                        current_c.get_vector(),
                        next_n.get_vector(),
                    )

            file.write(
                f"{chain.id}\t"
                f"{residue_number}\t"
                f"{amino_acid}\t"
                f"{phi}\t"
                f"{psi}\n"
            )

        file.write("\n")


def check_sequence(angles_file, structure_name, chain_file):
    """Check whether an angles file sequence matches a chain file.

    :param angles_file: Path to the file containing phi/psi angles.
    :type angles_file: pathlib.Path
    :param structure_name: Structure name to locate in the angles file.
    :type structure_name: str
    :param chain_file: Path to the file containing the expected sequence.
    :type chain_file: pathlib.Path
    :return: ``True`` if the sequences match, otherwise ``False``.
    :rtype: bool
    """
    angles_seq = ""
    inside_structure = False

    with open(angles_file, "r", encoding="utf-8") as file:
        for line in file:
            line = line.strip()

            if not line:
                continue
            if line.startswith(">"):
                current_name = line[1:].strip()
                inside_structure = current_name == structure_name
                continue

            if not inside_structure:
                continue

            parts = line.split()

            if len(parts) < 3:
                continue

            angles_seq += parts[2]

    atom_seq = ""

    with open(chain_file, "r", encoding="utf-8") as file:
        for line in file:
            parts = line.split()

            if not parts:
                continue

            atom_seq += parts[-1]

    angles_seq = angles_seq.strip().upper()
    atom_seq = atom_seq.strip().upper()

    if angles_seq != atom_seq:
        print(f"\nSequence mismatch: {structure_name}")
        print(f"angles_seq: {angles_seq}")
        print(f"atom_seq:   {atom_seq}")
        print(f"angles length: {len(angles_seq)}")
        print(f"atom length:   {len(atom_seq)}")

        min_length = min(len(angles_seq), len(atom_seq))

        for i in range(min_length):
            if angles_seq[i] != atom_seq[i]:
                print(
                    f"First difference at position "
                    f"{i + 1}: "
                    f"{angles_seq[i]} != {atom_seq[i]}"
                )
                break

        return False

    return True


def main():
    """Run the phi/psi calculation workflow.

    :return: ``None``.
    :rtype: None
    """
    structure_folder = Path("CIFs")
    sequences_folder = Path("extracted_chains")
    output_file = Path("phi_psi.txt")

    output_file.write_text("", encoding="utf-8")

    for structure_file in sorted(structure_folder.iterdir()):
        if (
            not structure_file.is_file()
            or structure_file.suffix.lower() not in [".cif", ".mmcif"]
        ):
            continue

        try:
            sequence_file = (
                sequences_folder / f"{structure_file.stem}.txt"
            )

            with open(sequence_file, "r", encoding="utf-8") as file:
                first_line = file.readline().strip().split()
            chain_id = first_line[3]

            get_phi_psi(structure_file, output_file, chain_id)

            if sequence_file.exists():
                result = check_sequence(
                    output_file, structure_file.name, sequence_file
                )

                if result:
                    print(
                        f"Sequence check True for "
                        f"{structure_file}"
                    )
                else:
                    print(
                        f"Sequence check False for "
                        f"{structure_file}"
                    )
            else:
                print(f"Sequence file not found: {sequence_file}")

        except Exception as exc:
            print(
                f"ERROR processing "
                f"{structure_file.name}: {exc}"
            )
            traceback.print_exc()

    print("Completed")


if __name__ == "__main__":
    main()
