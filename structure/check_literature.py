#!/usr/bin/env python3
"""Validate structure/literature_mutations.csv against the parent sequence.

The point of this check is that a numbering translation error is silent: a
wrong Kabat->parent offset still produces a plausible-looking row. If the
published wild-type residue does not match PARENT[position-1], the
translation is wrong and the row must not be used.

Also cross-references each position against red_green.csv, so it is visible
whether a published hotspot landed on a RED (contact) or GREEN (mutable)
position.

Usage:
    python3 structure/check_literature.py
Exits non-zero on any mismatch.
"""

import csv
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from contacts import parent_sequence  # noqa: E402


def load_classes():
    """position -> RED/GREEN from red_green.csv, if present."""
    path = HERE / "red_green.csv"
    if not path.exists():
        return {}
    out = {}
    with path.open(newline="") as fh:
        for row in csv.DictReader(fh):
            out[int(row["position"])] = row["class"]
    return out


def main():
    parent = parent_sequence()
    classes = load_classes()
    path = HERE / "literature_mutations.csv"

    rows = list(csv.DictReader(path.open(newline="")))
    if not rows:
        raise SystemExit(f"{path} has no rows")

    bad = 0
    print(f"{'pos':>5} {'wt':>3} {'mut':>4} {'parent':>7} {'class':>6}  verdict")
    for row in rows:
        pos = int(row["position"])
        wt, mut = row["wt_aa"], row["mut_aa"]
        if not 1 <= pos <= len(parent):
            print(f"{pos:>5} {wt:>3} {mut:>4} {'-':>7} {'-':>6}  POSITION OUT OF RANGE")
            bad += 1
            continue
        actual = parent[pos - 1]
        cls = classes.get(pos, "-")
        if actual != wt:
            print(f"{pos:>5} {wt:>3} {mut:>4} {actual:>7} {cls:>6}  "
                  f"MISMATCH: numbering translation is wrong")
            bad += 1
        else:
            print(f"{pos:>5} {wt:>3} {mut:>4} {actual:>7} {cls:>6}  ok")

    # Sanity rules that apply to any mutation we would actually build.
    for row in rows:
        pos = int(row["position"])
        if 121 <= pos <= 135:
            print(f"  WARNING: position {pos} is in the linker")
            bad += 1
        if parent[pos - 1] == "C":
            print(f"  WARNING: position {pos} is a conserved cysteine")
            bad += 1

    print(f"\n{len(rows)} rows, {bad} problem(s)")
    if classes:
        red = [int(r["position"]) for r in rows
               if classes.get(int(r["position"])) == "RED"]
        green = [int(r["position"]) for r in rows
                 if classes.get(int(r["position"])) == "GREEN"]
        print(f"on RED positions:   {sorted(red)}")
        print(f"on GREEN positions: {sorted(green)}")
    if bad:
        raise SystemExit(1)
    print("all wild-type residues match the parent sequence")


if __name__ == "__main__":
    main()
