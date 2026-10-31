#!/usr/bin/env python3
"""242x242 CA-CA distance matrix for the parent scFv, from PDB 1N8Z.

M[i, j] is the CA-CA distance in Angstrom between parent positions i+1 and
j+1, so the array indexes directly by parent position minus one. Positions
with no coordinates in the crystal are NaN rows and columns: that is exactly
the (GGGGS)x3 linker, 121-135, which is absent because 1N8Z is a Fab, not an
scFv. The array stays 242x242 so no caller has to re-derive an offset.

Usage:
    python3 structure/ca_distances.py 1N8Z.cif structure/ca_distances.npy

Load with:
    import numpy as np
    D = np.load("structure/ca_distances.npy")
    D[21, 95]            # positions 22 and 96, the VH disulfide pair
    np.nanmin(D[32])     # closest CA to position 33

Requires numpy only.
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from contacts import (  # noqa: E402
    chain_residues,
    map_to_parent,
    parent_sequence,
    parse_mmcif_atoms,
)

CHAIN_VH, CHAIN_VL = "B", "A"
LINKER = (121, 135)

# Expected disulfide pairs on the parent, as a sanity check on the mapping.
DISULFIDES = [(22, 96), (158, 223)]


def build(cif):
    parent = parent_sequence()
    n = len(parent)
    if n != 242:
        raise SystemExit(f"parent is {n} aa, expected 242")

    rows = parse_mmcif_atoms(cif)
    scfv = {}
    scfv.update(map_to_parent(chain_residues(rows, CHAIN_VH), parent, 0, 120,
                              "chain B / VH"))
    scfv.update(map_to_parent(chain_residues(rows, CHAIN_VL), parent, 135, 242,
                              "chain A / VL"))

    ca = np.full((n, 3), np.nan)
    missing_ca = []
    for pos, (_, _, _, atoms) in scfv.items():
        if "CA" in atoms:
            ca[pos - 1] = atoms["CA"]
        else:
            missing_ca.append(pos)
    if missing_ca:
        print(f"warning: no CA atom for positions {missing_ca}", file=sys.stderr)

    diff = ca[:, None, :] - ca[None, :, :]
    D = np.sqrt((diff ** 2).sum(-1))
    np.fill_diagonal(D, 0.0)
    # A NaN coordinate already propagates to NaN here; make the diagonal of a
    # missing position NaN too, so "has coordinates" is testable as ~isnan.
    no_coord = np.isnan(ca).any(axis=1)
    D[no_coord, :] = np.nan
    D[:, no_coord] = np.nan
    return parent, D, no_coord


def main():
    if len(sys.argv) < 3:
        raise SystemExit(__doc__)
    cif, out = sys.argv[1], sys.argv[2]
    parent, D, no_coord = build(cif)

    missing = sorted(int(i) + 1 for i in np.where(no_coord)[0])
    expected = list(range(LINKER[0], LINKER[1] + 1))
    print(f"shape           {D.shape}")
    print(f"with coords     {(~no_coord).sum()}/242")
    print(f"NaN positions   {missing}")
    if missing != expected:
        print(f"WARNING: expected NaN only for the linker {expected}", file=sys.stderr)
    else:
        print("                (exactly the linker 121-135, as expected)")

    print("\ndisulfide sanity check (CA-CA, expect roughly 5-7 A):")
    ok = True
    for a, b in DISULFIDES:
        d = D[a - 1, b - 1]
        wt = f"{parent[a - 1]}{a}-{parent[b - 1]}{b}"
        flag = "OK" if 4.0 <= d <= 8.0 else "OUT OF RANGE"
        if flag != "OK":
            ok = False
        print(f"  {wt:<12} {d:5.2f} A   {flag}")
    if not ok:
        raise SystemExit("disulfide distance out of range: mapping is suspect")

    finite = D[~np.isnan(D)]
    print(f"\nfinite entries  {finite.size} "
          f"(min {finite.min():.2f}, max {finite.max():.2f} A)")

    np.save(out, D)
    print(f"saved           {out}")


if __name__ == "__main__":
    main()
