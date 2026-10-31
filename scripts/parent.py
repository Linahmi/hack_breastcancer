"""Single source of truth for the parent sequence and its domain boundaries.

Everyone imports from here. Nobody retypes the sequence, nobody hardcodes
an index. Boundaries are computed from the linker, not assumed.

    from scripts.parent import PARENT, VH, LINKER, VL, REGIONS
"""

from pathlib import Path

LINKER_SEQ = "GGGGS" * 3  # fixed by the rules, 15 aa

_fasta = Path(__file__).resolve().parent.parent / "data" / "parent.fasta"
PARENT = "".join(
    line.strip() for line in _fasta.read_text().splitlines() if not line.startswith(">")
)

_i = PARENT.find(LINKER_SEQ)
if _i == -1:
    raise ValueError("linker (GGGGS)x3 not found in parent sequence")

# 1-based, inclusive, on the 242 aa scFv. This is THE numbering convention.
VH = (1, _i)                                  # 1 .. end of VH
LINKER = (_i + 1, _i + len(LINKER_SEQ))       # the 15 linker residues
VL = (_i + len(LINKER_SEQ) + 1, len(PARENT))  # start of VL .. 242


def region(pos: int) -> str:
    """Region of a 1-based position on the parent."""
    if VH[0] <= pos <= VH[1]:
        return "VH"
    if LINKER[0] <= pos <= LINKER[1]:
        return "linker"
    if VL[0] <= pos <= VL[1]:
        return "VL"
    raise ValueError(f"position {pos} outside 1..{len(PARENT)}")


REGIONS = {"VH": VH, "linker": LINKER, "VL": VL}

if __name__ == "__main__":
    print(f"parent length   {len(PARENT)}")
    print(f"VH              {VH[0]}-{VH[1]}  ({VH[1] - VH[0] + 1} aa)")
    print(f"linker          {LINKER[0]}-{LINKER[1]}  ({LINKER[1] - LINKER[0] + 1} aa)")
    print(f"VL              {VL[0]}-{VL[1]}  ({VL[1] - VL[0] + 1} aa)")
    print(f"headroom to 260 {260 - len(PARENT)} aa")
