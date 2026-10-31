#!/usr/bin/env python3
"""Build results/positions.csv: every parent position, binned by distance to HER2.

Reuses contacts.py so there is exactly one implementation of the chain mapping
and the distance calculation. Waters are excluded, like contacts.py: the
antigen is HER2 protein atoms only.

    contact       <= 4.5 A to the nearest HER2 heavy atom
    second_shell  4.5 A < d <= 8.0 A
    free          d > 8.0 A

min_distance_A is written for contact and second_shell rows only, so a blank
means "not in the shell" rather than "distance unknown".

Linker positions 121-135 have no coordinates: the crystal is a Fab, the linker
exists only in the scFv construct. They are reported free with no distance.

imgt and delta_sasa are in the repo header but not computed here; they stay
blank rather than being guessed.

Usage:
    python3 structure/positions_table.py 1N8Z.cif > results/positions.csv
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from contacts import (  # noqa: E402
    CDRS,
    CHAIN_AG,
    CHAIN_VH,
    CHAIN_VL,
    CUTOFF,
    chain_residues,
    map_to_parent,
    min_distances,
    parent_sequence,
    parse_mmcif_atoms,
)

SECOND_SHELL = 8.0
LINKER = range(121, 136)


def region_of(pos):
    for name, (a, b) in CDRS.items():
        if a <= pos <= b:
            return name
    return "linker" if pos in LINKER else "FR"


def main():
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)

    parent = parent_sequence()
    if len(parent) != 242:
        raise SystemExit(f"parent is {len(parent)} aa, expected 242")

    rows = parse_mmcif_atoms(sys.argv[1])
    vh = chain_residues(rows, CHAIN_VH)
    vl = chain_residues(rows, CHAIN_VL)
    ag = chain_residues(rows, CHAIN_AG)  # protein only: no HOH
    if not (vh and vl and ag):
        raise SystemExit("missing one of the expected chains B (VH), A (VL), C (HER2)")

    scfv = {}
    scfv.update(map_to_parent(vh, parent, 0, 120, "chain B / VH"))
    scfv.update(map_to_parent(vl, parent, 135, 242, "chain A / VL"))

    dmin = min_distances(scfv, np.concatenate([xyz for _, _, xyz, _ in ag], axis=0))

    for p in set(dmin) & set(LINKER):
        raise SystemExit(f"linker position {p} unexpectedly has coordinates")

    print("position,wt_aa,region,imgt,category,min_distance_A,delta_sasa")
    tally = {"contact": 0, "second_shell": 0, "free": 0}
    for pos in range(1, 243):
        d = dmin.get(pos)
        if d is not None and d <= CUTOFF:
            cat, dist = "contact", f"{d:.2f}"
        elif d is not None and d <= SECOND_SHELL:
            cat, dist = "second_shell", f"{d:.2f}"
        else:
            cat, dist = "free", ""
        tally[cat] += 1
        print(f"{pos},{parent[pos - 1]},{region_of(pos)},,{cat},{dist},")

    print(
        f"# {tally['contact']} contact, {tally['second_shell']} second_shell, "
        f"{tally['free']} free; {len(scfv)}/242 positions had coordinates",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
