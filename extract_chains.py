"""
Extract the longest chain from PDB/mmCIF structures.

For each structure in an input directory, this module extracts:

1. The CA atom coordinates of the longest chain. The longest chain is
   the chain with the most amino acid CA atoms in the first model.
2. The SEQRES sequence of that same chain.

Outputs are written to an output directory as:

``<stem>.txt``
    ``x y z chain_id res_num res_name`` per line.
``<stem>_seq.txt``
    One-letter SEQRES sequence, without a trailing newline.

Output file names contain only lowercase letters.

Notes
-----
* Only ``ATOM``/``HETATM`` residues that map to a known amino acid
  three-letter code, including modified residues such as ``MSE`` and
  ``SEP``, are considered. Ions and ligands named ``CA`` (calcium) are
  therefore correctly ignored.
* The chain-ID lookup against SEQRES uses the ``label_asym_id`` from
  ``_pdbx_poly_seq_scheme``, which is consistent with using
  ``MMCIFParser(auth_chains=False)``.
"""

import os
import re
import warnings

from Bio import BiopythonWarning, SeqIO
from Bio.Data import PDBData
from Bio.PDB import MMCIFParser, PDBParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.PDB.Residue import DisorderedResidue


# ---------------------------------------------------------------------------
# Constants.
# ---------------------------------------------------------------------------

CIF_EXTS = {".cif", ".mmcif", ".mcif"}
PDB_EXTS = {".pdb", ".ent"}

# Extended three-letter to one-letter mapping, including modified
# residues (MSE, SEP, PTR, ...) and the "UNK" placeholder.
_3TO1_EXTENDED = {
    k.upper(): v for k, v in PDBData.protein_letters_3to1_extended.items()
}
_3TO1_EXTENDED.setdefault("UNK", "X")


# ---------------------------------------------------------------------------
# Small helper functions.
# ---------------------------------------------------------------------------

def _as_list(value):
    """Normalize an MMCIF2Dict value into a list.

    :param value: Value to normalize.
    :type value: object or None
    :returns: Normalized list.
    :rtype: list
    """
    if value is None:
        return []
    if isinstance(value, list):
        return value
    return [value]


def _safe_chain_id(chain_id):
    """Return display and filename-safe chain IDs.

    :param chain_id: Chain identifier.
    :type chain_id: str or None
    :returns: Tuple of display chain ID and filename-safe chain ID.
    :rtype: tuple[str, str]
    """
    display = (chain_id or "").strip() or "_"
    safe = re.sub(r"[^A-Za-z0-9_-]", "_", display)
    return display, safe


def _lowercase_stem(stem):
    """Return a filename stem containing only lowercase letters.

    :param stem: Original filename stem.
    :type stem: str
    :returns: Lowercase filename stem.
    :rtype: str
    """
    return stem.lower()


def _residue_number_string(residue):
    """Return the residue number string.

    The string is the sequence ID followed by the insertion code if an
    insertion code is present.

    :param residue: Residue object.
    :type residue: Bio.PDB.Residue.Residue
    :returns: Residue number string.
    :rtype: str
    """
    _, seq_id, ins_code = residue.id
    base = str(seq_id)
    if ins_code and ins_code.strip():
        return base + ins_code.strip()
    return base


def _residue_to_one_letter(resname):
    """Map a three-letter residue name to a one-letter code.

    Unknown names fall back to ``'X'``.

    :param resname: Three-letter residue name.
    :type resname: str
    :returns: One-letter residue code.
    :rtype: str
    """
    return _3TO1_EXTENDED.get(resname.upper(), "X")


def _is_amino_acid(residue):
    """Return ``True`` if the residue maps to a known amino acid.

    :param residue: Residue object.
    :type residue: Bio.PDB.Residue.Residue
    :returns: ``True`` if the residue is an amino acid, otherwise
        ``False``.
    :rtype: bool
    """
    return residue.get_resname().upper() in _3TO1_EXTENDED


# ---------------------------------------------------------------------------
# SEQRES extraction.
# ---------------------------------------------------------------------------

def _get_seqres_mmcif(path, chain_id):
    """Extract SEQRES for a chain from an mmCIF file.

    :param path: Path to the mmCIF file.
    :type path: str
    :param chain_id: Chain identifier.
    :type chain_id: str
    :returns: SEQRES sequence.
    :rtype: str
    """
    cif_dict = MMCIF2Dict(path)

    asym_ids = _as_list(cif_dict.get("_pdbx_poly_seq_scheme.asym_id"))
    seq_ids = _as_list(cif_dict.get("_pdbx_poly_seq_scheme.seq_id"))
    mon_ids = _as_list(cif_dict.get("_pdbx_poly_seq_scheme.mon_id"))

    rows = []
    for asym_id, seq_id, mon_id in zip(asym_ids, seq_ids, mon_ids):
        if asym_id != chain_id:
            continue
        try:
            key = int(seq_id)
        except (TypeError, ValueError):
            key = float("inf")
        rows.append((key, seq_id, mon_id))

    rows.sort(key=lambda r: r[0])

    sequence = []
    seen_positions = set()
    for _, seq_id, mon_id in rows:
        if seq_id in seen_positions:
            continue
        seen_positions.add(seq_id)
        sequence.append(_residue_to_one_letter(mon_id))

    return "".join(sequence)


def _get_seqres_pdb(path, chain_id):
    """Extract SEQRES for a chain from a PDB file.

    :param path: Path to the PDB file.
    :type path: str
    :param chain_id: Chain identifier.
    :type chain_id: str
    :returns: SEQRES sequence.
    :rtype: str
    """
    for record in SeqIO.parse(path, "pdb-seqres"):
        # Biopython sets record.id to e.g. "1ABC:A".
        if ":" in record.id:
            rec_chain = record.id.split(":", 1)[1]
        else:
            rec_chain = record.annotations.get("chain", "")
        if rec_chain == chain_id:
            return str(record.seq)
    return ""


def get_seqres(path, chain_id):
    """Return the SEQRES sequence for a chain.

    Returns ``''`` if the sequence is unavailable.

    :param path: Path to the structure file.
    :type path: str
    :param chain_id: Chain identifier.
    :type chain_id: str
    :returns: SEQRES sequence.
    :rtype: str
    :raises ValueError: If the file extension is unsupported.
    """
    ext = os.path.splitext(path)[1].lower()
    if ext in CIF_EXTS:
        return _get_seqres_mmcif(path, chain_id)
    if ext in PDB_EXTS:
        return _get_seqres_pdb(path, chain_id)
    raise ValueError(f"Unsupported file extension: {ext}")


# ---------------------------------------------------------------------------
# CA extraction.
# ---------------------------------------------------------------------------

def _iter_chain_ca(chain):
    """Yield CA information for amino acid residues in ``chain``.

    For every amino acid residue in ``chain`` that has a usable CA
    atom, yield ``(ca_atom, residue_number_string, one_letter_code)``.

    Disordered residues (altloc residues) and disordered CA atoms
    (altloc atoms) are handled. Entries that have no CA in any altloc
    are skipped.

    :param chain: Chain object.
    :type chain: Bio.PDB.Chain.Chain
    :yields: Tuple of CA atom, residue number string, and one-letter
        code.
    :rtype: iterator
    """
    for residue in chain:
        if not _is_amino_acid(residue):
            continue

        working_residue = residue
        if isinstance(residue, DisorderedResidue):
            candidates = [
                r for r in residue.disordered_get_list() if "CA" in r
            ]
            if not candidates:
                continue
            working_residue = candidates[0]

        if "CA" not in working_residue:
            continue

        ca = working_residue["CA"]
        try:
            if ca.is_disordered():
                ca = ca.disordered_get_list()[0]
        except AttributeError:
            pass

        # Sanity check: ensure the atom is actually an alpha carbon.
        element = getattr(ca, "element", None)
        if element is not None and element.strip().upper() not in ("C", ""):
            continue

        yield (
            ca,
            _residue_number_string(residue),
            _residue_to_one_letter(working_residue.get_resname()),
        )


# ---------------------------------------------------------------------------
# Main entry point.
# ---------------------------------------------------------------------------

def extract_longest_chain_to_files(input_dir, output_dir):
    """Process all PDB/mmCIF structures in ``input_dir``.

    For each structure:

    * Parse the first model.
    * Select the chain with the most amino acid CA atoms.
    * Write CA coordinates and SEQRES to ``output_dir``.

    Output file names contain only lowercase letters.

    :param input_dir: Directory containing PDB/mmCIF files.
    :type input_dir: str
    :param output_dir: Directory for output files.
    :type output_dir: str
    """
    warnings.simplefilter("ignore", BiopythonWarning)
    os.makedirs(output_dir, exist_ok=True)

    cif_parser = MMCIFParser(QUIET=True, auth_chains=False)
    pdb_parser = PDBParser(QUIET=True)

    for filename in sorted(os.listdir(input_dir)):
        path = os.path.join(input_dir, filename)
        if not os.path.isfile(path):
            continue

        stem, ext = os.path.splitext(filename)
        stem_lower = _lowercase_stem(stem)
        ext_lower = ext.lower()

        if ext_lower in CIF_EXTS:
            parser = cif_parser
        elif ext_lower in PDB_EXTS:
            parser = pdb_parser
        else:
            continue

        try:
            structure = parser.get_structure(stem_lower, path)
        except Exception as exc:
            print(f"Skipping {filename}: {exc}")
            continue

        if len(structure) == 0:
            print(f"Skipping {filename}: no models")
            continue

        model = structure[0]

        best_chain_id = None
        best_ca_atoms = []
        best_res_nums = []
        best_res_names = []

        for chain in model:
            ca_atoms, res_nums, res_names = [], [], []
            for ca, res_num, res_name in _iter_chain_ca(chain):
                ca_atoms.append(ca)
                res_nums.append(res_num)
                res_names.append(res_name)

            if len(ca_atoms) > len(best_ca_atoms):
                best_chain_id = chain.id
                best_ca_atoms = ca_atoms
                best_res_nums = res_nums
                best_res_names = res_names

        if best_chain_id is None or not best_ca_atoms:
            print(f"Skipping {filename}: no CA atoms found")
            continue

        chain_id_display, _ = _safe_chain_id(best_chain_id)

        # --- CA coordinate output -----------------------------------------
        out_filename = f"{stem_lower}.txt"
        out_path = os.path.join(output_dir, out_filename)

        with open(out_path, "w") as out_f:
            for ca, res_num, res_name in zip(
                best_ca_atoms, best_res_nums, best_res_names
            ):
                x, y, z = ca.coord
                out_f.write(
                    f"{x:.3f} {y:.3f} {z:.3f} "
                    f"{chain_id_display} {res_num} {res_name}\n"
                )

        print(f"Wrote {out_filename} ({len(best_ca_atoms)} CA atoms)")

        # --- SEQRES output ------------------------------------------------
        seq_filename = f"{stem_lower}_seq.txt"
        seq_path = os.path.join(output_dir, seq_filename)

        try:
            seq = get_seqres(path, chain_id_display)
        except Exception as exc:
            print(f"  Could not extract SEQRES for {filename}: {exc}")
            continue

        if not seq:
            print(
                f"  No SEQRES found for chain {chain_id_display} "
                f"in {filename}"
            )

        with open(seq_path, "w") as seq_file:
            seq_file.write(seq or "")

        print(f"Wrote {seq_filename} ({len(seq or '')} residues)")


if __name__ == "__main__":
    extract_longest_chain_to_files("CIFs", "extracted_chains")
