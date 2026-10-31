"""Format and developability validator for candidate scFv designs.

    python scripts/validate.py results/candidates.csv

Exit code 0 if every row passes, 1 otherwise. Run this before anything
goes near a GPU, and once more on the final 15 before submission.
"""

import csv
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from scripts.parent import PARENT, LINKER_SEQ  # noqa: E402

MIN_MUTATIONS = 5
MAX_LENGTH = 260
VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")


def mutations(seq: str) -> list[str]:
    """Mutations vs parent, as H31N style strings. Only valid for equal length."""
    if len(seq) != len(PARENT):
        return []
    return [
        f"{wt}{i}{mt}" for i, (wt, mt) in enumerate(zip(PARENT, seq), 1) if wt != mt
    ]


def check(name: str, seq: str) -> list[str]:
    errs = []
    seq = seq.strip().upper()

    bad = set(seq) - VALID_AA
    if bad:
        errs.append(f"invalid residues: {''.join(sorted(bad))}")
        return errs

    if len(seq) > MAX_LENGTH:
        errs.append(f"length {len(seq)} > {MAX_LENGTH}")
    if len(seq) < len(PARENT) - 20:
        errs.append(f"length {len(seq)} suspiciously short")

    # linker present exactly once, and in the middle third
    hits = [m.start() for m in re.finditer(LINKER_SEQ, seq)]
    if len(hits) != 1:
        errs.append(f"linker (GGGGS)x3 found {len(hits)}x, expected once")
    elif not (len(seq) * 0.3 < hits[0] < len(seq) * 0.7):
        errs.append("linker not in the middle - VH/VL order may be wrong")

    # mutation count (only checkable at parent length; flag otherwise)
    if len(seq) == len(PARENT):
        n = len(mutations(seq))
        if n < MIN_MUTATIONS:
            errs.append(f"only {n} mutations, need >= {MIN_MUTATIONS}")
    else:
        errs.append(f"length {len(seq)} != parent 242 - count mutations by hand")

    # developability
    if seq.count("C") != PARENT.count("C"):
        errs.append(f"cysteine count {seq.count('C')} != parent {PARENT.count('C')}")
    for motif, label in [(r"N[^P][ST]", "N-glycosylation sequon"), (r"NG", "deamidation NG"), (r"DG", "isomerisation DG")]:
        new = len(re.findall(motif, seq)) - len(re.findall(motif, PARENT))
        if new > 0:
            errs.append(f"{new} new {label}")

    return errs


def main(path: str) -> int:
    rows = list(csv.DictReader(open(path)))
    if not rows:
        print("no rows")
        return 1

    failed = 0
    for row in rows:
        errs = check(row.get("name", "?"), row["sequence"])
        if errs:
            failed += 1
            print(f"FAIL  {row.get('name', '?')}")
            for e in errs:
                print(f"        {e}")

    print(f"\n{len(rows) - failed}/{len(rows)} passed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else "results/candidates.csv"))
