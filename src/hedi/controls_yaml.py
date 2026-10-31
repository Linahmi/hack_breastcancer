"""Emit Boltz-2 input YAMLs for the three controls in results/controls.csv.

Only parent_wt has been cofolded so far (outputs/parent_quick, and that was
a deliberate smoke test). These YAMLs are what the other two need.

The reason this is a script and not three hand-written files: the scFv
disulfide constraints MOVE between controls. The parent has its cysteines at
22/96 (VH) and 158/223 (VL); reverse_vl_vh puts VL first, so the same four
cysteines land at 23/88 and 144/218. Copy-pasting the parent's constraint
block into that YAML would restrain four positions that are not cysteines,
and Boltz would either error or silently distort the control. Every bond
emitted here is checked to land on a C in its own sequence.

    python src/hedi/controls_yaml.py            # writes src/hedi/inputs/
    python src/hedi/controls_yaml.py --print

NOTE: no `msa: empty` is emitted. outputs/parent_quick used single-sequence
mode for speed and its PROVENANCE.md says not to use its numbers. Controls
that are meant to be compared against designs need the same MSA treatment as
the designs.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from scripts.parent import LINKER_SEQ, PARENT, VH, VL  # noqa: E402

CONTROLS_CSV = REPO / "results" / "controls.csv"
OUT_DIR = Path(__file__).resolve().parent / "inputs"

# HER2 domain IV, P04626 511-652 (src/aida/her2_domainIV.fasta on the aida
# branch). Chain B of every cofold.
HER2_DIV = (
    "CHQLCARGHCWGPGPTQCVNCSQFLRGQECVEECRVLQGLPREYVNARHCLPCHPECQPQ"
    "NGSVTCFGPEADQCVACAHYKDPPFCVARCPSGVKPDLSYMPIWKFPDEEGACQPCPINC"
    "THSCVDLDDKGCPAEQRASPLT"
)

# The 10 domain IV disulfides, as used in src/aida/inputs/parent_smoketest.yaml.
HER2_BONDS = [
    (1, 10), (5, 18), (21, 30), (34, 50), (53, 66),
    (57, 74), (77, 86), (90, 113), (116, 124), (120, 132),
]

# The parent scFv's two intradomain disulfides, 1-based on the 242 aa parent.
PARENT_BONDS = [(22, 96), (158, 223)]


def remap_parent_position(name: str, pos: int) -> int:
    """Map a 1-based parent position onto the same residue in a control.

    parent_wt and cdrh3_scrambled keep parent numbering -- the scramble
    permutes CDR-H3 (97-109), which contains no cysteine. reverse_vl_vh
    reorders the domains, so every position shifts.
    """
    if name in ("parent_wt", "cdrh3_scrambled"):
        return pos
    if name == "reverse_vl_vh":
        vl_len = VL[1] - VL[0] + 1
        if VL[0] <= pos <= VL[1]:          # VL moves to the front
            return pos - VL[0] + 1
        if VH[0] <= pos <= VH[1]:          # VH moves after VL + linker
            return pos + vl_len + len(LINKER_SEQ)
        return pos - VH[1] + vl_len        # linker
    raise ValueError(f"no position mapping defined for control {name!r}")


def scfv_bonds(name: str, seq: str) -> list[tuple[int, int]]:
    """The control's own disulfide pairs, verified to land on cysteines."""
    bonds = [
        (remap_parent_position(name, a), remap_parent_position(name, b))
        for a, b in PARENT_BONDS
    ]
    for a, b in bonds:
        for p in (a, b):
            if not 1 <= p <= len(seq):
                raise ValueError(f"{name}: bond position {p} outside 1..{len(seq)}")
            if seq[p - 1] != "C":
                raise ValueError(
                    f"{name}: position {p} is {seq[p-1]}, not C -- the disulfide "
                    "remapping is wrong, do not run this YAML"
                )
    found = sorted(i for i, c in enumerate(seq, 1) if c == "C")
    if sorted(p for pair in bonds for p in pair) != found:
        raise ValueError(
            f"{name}: bonds {bonds} do not cover exactly the cysteines {found}"
        )
    return bonds


def verify_her2() -> None:
    for a, b in HER2_BONDS:
        for p in (a, b):
            if HER2_DIV[p - 1] != "C":
                raise ValueError(f"HER2 bond position {p} is {HER2_DIV[p-1]}, not C")
    found = sorted(i for i, c in enumerate(HER2_DIV, 1) if c == "C")
    if sorted(p for pair in HER2_BONDS for p in pair) != found:
        raise ValueError("HER2 bonds do not cover exactly its 20 cysteines")


def build_yaml(name: str, seq: str) -> str:
    bonds = scfv_bonds(name, seq)
    lines = [
        "version: 1",
        "sequences:",
        "  - protein:",
        "      id: A",
        f"      sequence: {seq}",
        "  - protein:",
        "      id: B",
        f"      sequence: {HER2_DIV}",
        "constraints:",
    ]
    for chain, pairs in (("A", bonds), ("B", HER2_BONDS)):
        for a, b in pairs:
            lines += [
                "  - bond:",
                f"      atom1: [{chain}, {a}, SG]",
                f"      atom2: [{chain}, {b}, SG]",
            ]
    return "\n".join(lines) + "\n"


def main(argv: list[str]) -> int:
    verify_her2()
    if not CONTROLS_CSV.exists():
        print(f"missing {CONTROLS_CSV}; run src/hedi/controls.py first")
        return 1

    rows = list(csv.DictReader(open(CONTROLS_CSV)))
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for r in rows:
        name, seq = r["name"], r["sequence"]
        text = build_yaml(name, seq)
        bonds = scfv_bonds(name, seq)
        print(f"{name:18s} {len(seq)} aa  scFv disulfides {bonds}")
        if "--print" in argv:
            print(text)
        else:
            (OUT_DIR / f"{name}.yaml").write_text(text)

    if "--print" not in argv:
        print(f"\nwrote {len(rows)} YAML(s) to {OUT_DIR.relative_to(REPO)}")
        print(
            "\nparent_wt is included for a like-for-like rerun: outputs/parent_quick\n"
            "used msa:empty / recycling 1 / 1 sample and its own PROVENANCE.md says\n"
            "its numbers are not a usable baseline. Run all three the same way."
        )
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
