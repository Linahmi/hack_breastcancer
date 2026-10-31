#!/usr/bin/env python3
"""RED/GREEN position constraints for the trastuzumab scFv from PDB 1N8Z.

RED   = any scFv position with a heavy atom within CUTOFF of any HER2 heavy
        atom. Direct paratope. Do not mutate.
GREEN = Kabat CDR position that is NOT in contact. Safe to scan.

Positions are 1-based on the 242 aa parent scFv (VH 1-120, linker 121-135,
VL 136-242), never PDB numbering. The mapping from structure to parent index
is established by exact sequence alignment and asserted residue by residue,
so a numbering slip fails loudly instead of producing a plausible wrong answer.

Waters are excluded from the antigen by default. Crystallographic HOH in a
2.52 A structure is solvent, not epitope, and counting it inflates the
paratope: with waters included, position 201 (framework) scores 3.52 A on a
single HOH while its nearest HER2 protein atom is 5.45 A away. The first
committed red_green.csv was produced with waters included; --include-waters
reproduces it exactly, so the difference is auditable rather than asserted.

Usage:
    python3 structure/contacts.py 1N8Z.cif > structure/red_green.csv
    python3 structure/contacts.py 1N8Z.cif --sweep           # cutoff robustness
    python3 structure/contacts.py 1N8Z.cif --include-waters   # legacy behaviour

Requires numpy only. No Biopython, no network.
"""

import sys
from difflib import SequenceMatcher
from pathlib import Path

import numpy as np

CUTOFF = 4.5  # Angstrom, any heavy atom to any heavy atom

# 1N8Z chain roles, established by sequence (not assumed from chain letters).
CHAIN_VH, CHAIN_VL, CHAIN_AG = "B", "A", "C"

# Kabat CDRs, in parent 1-based numbering.
CDRS = {
    "CDR-H1": (31, 35),
    "CDR-H2": (50, 66),
    "CDR-H3": (99, 109),
    "CDR-L1": (159, 169),
    "CDR-L2": (185, 191),
    "CDR-L3": (224, 232),
}

THREE2ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V", "MSE": "M",
}


def parent_sequence():
    """The 242 aa parent, from the repo's single source of truth."""
    fasta = Path(__file__).resolve().parent.parent / "data" / "parent.fasta"
    return "".join(
        ln.strip() for ln in fasta.read_text().splitlines() if not ln.startswith(">")
    )


def parse_mmcif_atoms(path):
    """Minimal mmCIF _atom_site reader. Returns a list of row dicts."""
    cols, rows, in_loop, reading = [], [], False, False
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if s.startswith("_atom_site."):
                cols.append(s.split(".", 1)[1].split()[0])
                in_loop = True
                continue
            if in_loop and not reading:
                if s.startswith(("ATOM", "HETATM")):
                    reading = True
                else:
                    continue
            if reading:
                if not s.startswith(("ATOM", "HETATM")):
                    break
                f = s.split()
                if len(f) == len(cols):
                    rows.append(dict(zip(cols, f)))
    if not rows:
        raise SystemExit(f"no _atom_site rows parsed from {path}")
    return rows


def chain_residues(rows, chain, keep=None):
    """Ordered residues of one auth chain, model 1, altloc . or A.

    By default keeps standard amino acids only. `keep` may name extra
    comp_ids to retain (e.g. {"HOH"}) for the legacy comparison.

    Returns [(auth_seq_id, one_letter, coords_array, atoms), ...] in file order.
    """
    seen, order = {}, []
    extra = keep or set()
    for r in rows:
        if r["auth_asym_id"] != chain or r["pdbx_PDB_model_num"] != "1":
            continue
        if r["label_alt_id"] not in (".", "A"):
            continue
        comp = r["label_comp_id"]
        if comp not in THREE2ONE and comp not in extra:
            continue  # waters, sugars, ions
        if r["type_symbol"] == "H":
            continue
        key = (int(r["auth_seq_id"]), r["pdbx_PDB_ins_code"], comp)
        if key not in seen:
            seen[key] = {"aa": THREE2ONE.get(comp, "X"), "xyz": [], "atoms": {}}
            order.append(key)
        xyz = (float(r["Cartn_x"]), float(r["Cartn_y"]), float(r["Cartn_z"]))
        seen[key]["xyz"].append(xyz)
        seen[key]["atoms"][r["label_atom_id"]] = xyz
    return [
        (k[0], seen[k]["aa"], np.asarray(seen[k]["xyz"], float), seen[k]["atoms"])
        for k in order
    ]


def map_to_parent(residues, parent, lo, hi, label):
    """Align a chain's residues onto parent[lo:hi]; assert exact identity.

    Returns {parent_position_1based: residue_tuple}.
    """
    seq = "".join(aa for _, aa, _, _ in residues)
    target = parent[lo:hi]
    blocks = [b for b in SequenceMatcher(None, seq, target, autojunk=False)
              .get_matching_blocks() if b.size]
    if not blocks:
        raise SystemExit(f"{label}: no alignment to parent[{lo}:{hi}]")
    best = max(blocks, key=lambda b: b.size)
    expected = hi - lo
    if best.size != expected:
        raise SystemExit(
            f"{label}: aligned only {best.size}/{expected} residues to "
            f"parent[{lo}:{hi}] - numbering or chain assignment is wrong"
        )
    mapping = {}
    for off in range(best.size):
        pos = lo + best.b + off + 1  # 1-based parent position
        res = residues[best.a + off]
        assert res[1] == parent[pos - 1], (
            f"{label}: structure has {res[1]} at parent position {pos}, "
            f"parent has {parent[pos - 1]}"
        )
        mapping[pos] = res
    return mapping


def min_distances(scfv, antigen_atoms):
    """Minimum heavy-atom distance from each scFv position to the antigen."""
    out = {}
    for pos, (_, _, xyz, _) in scfv.items():
        d = np.sqrt(((xyz[:, None, :] - antigen_atoms[None, :, :]) ** 2).sum(-1))
        out[pos] = float(d.min())
    return out


def cdr_of(pos):
    for name, (a, b) in CDRS.items():
        if a <= pos <= b:
            return name
    return "framework"


def main():
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    cif = sys.argv[1]
    sweep = "--sweep" in sys.argv
    with_waters = "--include-waters" in sys.argv

    parent = parent_sequence()
    if len(parent) != 242:
        raise SystemExit(f"parent is {len(parent)} aa, expected 242")

    rows = parse_mmcif_atoms(cif)
    vh = chain_residues(rows, CHAIN_VH)
    vl = chain_residues(rows, CHAIN_VL)
    ag = chain_residues(rows, CHAIN_AG, keep={"HOH"} if with_waters else None)
    if not (vh and vl and ag):
        raise SystemExit("missing one of the expected chains B (VH), A (VL), C (HER2)")

    # VH -> parent 1..120, VL -> parent 136..242. Linker 121-135 is not in
    # the crystal (the construct is a Fab) and has no coordinates.
    scfv = {}
    scfv.update(map_to_parent(vh, parent, 0, 120, "chain B / VH"))
    scfv.update(map_to_parent(vl, parent, 135, 242, "chain A / VL"))

    ag_atoms = np.concatenate([xyz for _, _, xyz, _ in ag], axis=0)
    dmin = min_distances(scfv, ag_atoms)

    if sweep:
        print("cutoff_A\tn_red\tred_positions", file=sys.stderr)
        for c in (4.0, 4.25, 4.5, 4.75, 5.0):
            red = sorted(p for p, d in dmin.items() if d <= c)
            print(f"{c}\t{len(red)}\t{red}", file=sys.stderr)
        core = sorted(p for p, d in dmin.items() if d <= 4.0)
        print(f"\ncore (<=4.0 A, stable across the range): {len(core)} residues",
              file=sys.stderr)
        return

    red = sorted(p for p, d in dmin.items() if d <= CUTOFF)
    cdr_positions = {p for a, b in CDRS.values() for p in range(a, b + 1)}
    green = sorted(cdr_positions - set(red))

    # GREEN must never touch the linker or a conserved cysteine.
    for p in green:
        assert not (121 <= p <= 135), f"GREEN position {p} is in the linker"
        assert parent[p - 1] != "C", f"GREEN position {p} is a cysteine"

    print("position,aa,class,min_dist_A")
    for p in red:
        print(f"{p},{parent[p - 1]},RED,{dmin[p]:.2f}")
    for p in green:
        print(f"{p},{parent[p - 1]},GREEN,")

    print(
        f"# {len(red)} RED, {len(green)} GREEN at {CUTOFF} A; "
        f"{len(scfv)}/242 parent positions had coordinates",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
