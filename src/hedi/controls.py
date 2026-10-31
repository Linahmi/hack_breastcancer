"""Generate results/controls.csv -- the three cofold controls.

These get cofolded against HER2 domain IV BEFORE any design, so that the
ipSAE numbers coming out of the pipeline can be calibrated against known
answers instead of read in a vacuum:

  parent_wt         positive control. Must score well. If it does not, the
                    pipeline is wrong (bad target trim, wrong chain order,
                    broken PAE) and no design score means anything.
  reverse_vl_vh     VL-(GGGGS)3-VH instead of VH-(GGGGS)3-VL. Same residues,
                    same composition, same linker -- only the domain order
                    differs. A scFv does tolerate both orders in reality, so
                    this is a sensitivity probe, not a guaranteed negative:
                    it says whether the score responds to domain topology at
                    all, and how much of the signal is just "this is an
                    antibody-shaped thing near HER2".
  cdrh3_scrambled   negative control. CDR-H3 permuted, composition identical.
                    CDR-H3 dominates the trastuzumab paratope, so this should
                    score clearly worse. If it does not, ipSAE is not reading
                    the interface and the whole ranking is suspect.

All boundaries come from scripts/parent.py; nothing here hardcodes an index.

    python src/hedi/controls.py            # writes results/controls.csv
    python src/hedi/controls.py --print    # show them without writing
"""

from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from scripts.parent import LINKER, LINKER_SEQ, PARENT, VH, VL  # noqa: E402

OUT_PATH = REPO / "results" / "controls.csv"

# Mirrors scripts/validate.py so we do not hand the team a control that trips
# their gate for an avoidable reason.
LIABILITY_MOTIFS = [(r"N[^P][ST]", "N-glycosylation sequon"), (r"NG", "deamidation NG"), (r"DG", "isomerisation DG")]
MIN_MUTATIONS = 5


def segments() -> dict[str, str]:
    """VH / linker / VL as sequences, sliced with parent.py's 1-based spans."""
    vh = PARENT[VH[0] - 1:VH[1]]
    linker = PARENT[LINKER[0] - 1:LINKER[1]]
    vl = PARENT[VL[0] - 1:VL[1]]
    assert linker == LINKER_SEQ, f"linker slice {linker!r} != {LINKER_SEQ!r}"
    assert vh + linker + vl == PARENT, "segments do not reassemble the parent"
    return {"VH": vh, "linker": linker, "VL": vl}


def find_cdrh3() -> tuple[int, int, str]:
    """Locate CDR-H3 by its flanking motifs. Returns (start, end, seq), 1-based inclusive.

    CDR-H3 runs from just after the conserved VH cysteine (the second of the
    VH pair) to just before the J-segment W-G-x-G motif. Both anchors are
    found by search, so this stays correct if the numbering convention moves.
    """
    vh = segments()["VH"]
    cys = vh.rfind("C")
    if cys == -1:
        raise ValueError("no cysteine found in VH; cannot anchor CDR-H3")
    m = re.search(r"WG.G", vh[cys:])
    if m is None:
        raise ValueError("no W-G-x-G motif found after the VH cysteine")
    start0 = cys + 1          # 0-based, inclusive
    end0 = cys + m.start()    # 0-based, exclusive
    seq = PARENT[start0:end0]
    if not 5 <= len(seq) <= 25:
        raise ValueError(f"CDR-H3 came out {len(seq)} aa ({seq!r}); anchors look wrong")
    return start0 + 1, end0, seq


def _hamming(a: str, b: str) -> int:
    return sum(1 for x, y in zip(a, b) if x != y) + abs(len(a) - len(b))


def _new_liabilities(seq: str) -> list[str]:
    """Motifs present more often than in the parent, as validate.py counts them."""
    out = []
    for motif, label in LIABILITY_MOTIFS:
        delta = len(re.findall(motif, seq)) - len(re.findall(motif, PARENT))
        if delta > 0:
            out.append(f"{delta} new {label}")
    return out


def scramble_cdrh3() -> tuple[str, int, str]:
    """Permute CDR-H3 deterministically. Returns (sequence, seed, permuted_cdr).

    Composition is preserved exactly, so cysteine count, length and pI are
    untouched and the only variable is CDR-H3 order. We search seeds for a
    permutation that also differs from the parent in at least MIN_MUTATIONS
    positions and introduces no new liability motif -- otherwise the control
    would fail scripts/validate.py for a reason unrelated to its purpose.
    """
    import random

    start, end, cdr = find_cdrh3()
    for seed in range(10_000):
        rng = random.Random(seed)
        chars = list(cdr)
        rng.shuffle(chars)
        permuted = "".join(chars)
        if _hamming(permuted, cdr) < MIN_MUTATIONS:
            continue
        seq = PARENT[:start - 1] + permuted + PARENT[end:]
        assert len(seq) == len(PARENT), "scramble changed the length"
        assert sorted(seq) == sorted(PARENT), "scramble changed the composition"
        if _new_liabilities(seq):
            continue
        return seq, seed, permuted
    raise RuntimeError("no acceptable CDR-H3 permutation found in 10000 seeds")


def build() -> list[dict]:
    seg = segments()
    start, end, cdr = find_cdrh3()
    scrambled, seed, permuted = scramble_cdrh3()

    rows = [
        {
            "name": "parent_wt",
            "sequence": PARENT,
            "role": "positive",
            "note": f"trastuzumab scFv parent, unmodified; VH 1-{VH[1]}, linker {LINKER[0]}-{LINKER[1]}, VL {VL[0]}-{VL[1]}",
        },
        {
            "name": "reverse_vl_vh",
            "sequence": seg["VL"] + seg["linker"] + seg["VH"],
            "role": "topology probe",
            "note": f"VL-(GGGGS)3-VH; linker moves to {len(seg['VL']) + 1}-{len(seg['VL']) + len(seg['linker'])}",
        },
        {
            "name": "cdrh3_scrambled",
            "sequence": scrambled,
            "role": "negative",
            "note": f"CDR-H3 {start}-{end} {cdr}->{permuted} (seed {seed}); composition identical",
        },
    ]
    for r in rows:
        r["length"] = len(r["sequence"])
        r["n_mut_vs_parent"] = _hamming(r["sequence"], PARENT)
    return rows


FIELDS = ["name", "sequence", "role", "length", "n_mut_vs_parent", "note"]


def main(argv: list[str]) -> int:
    rows = build()
    start, end, cdr = find_cdrh3()

    print(f"parent      {len(PARENT)} aa   VH 1-{VH[1]}  linker {LINKER[0]}-{LINKER[1]}  VL {VL[0]}-{VL[1]}")
    print(f"CDR-H3      {start}-{end}  {cdr}  ({len(cdr)} aa)\n")
    for r in rows:
        print(f"{r['name']:18s} {r['length']:4d} aa  {r['n_mut_vs_parent']:3d} diffs vs parent  [{r['role']}]")
        print(f"                   {r['note']}")

    if "--print" in argv:
        for r in rows:
            print(f"\n>{r['name']}\n{r['sequence']}")
        return 0

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_PATH, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS)
        w.writeheader()
        w.writerows(rows)
    print(f"\nwrote {OUT_PATH.relative_to(REPO)}  ({len(rows)} rows)")
    print(
        "\nNote: parent_wt has 0 mutations, so scripts/validate.py will FAIL it\n"
        "on its >=5 mutation rule. That is correct -- these are controls, not\n"
        "submissions. Do not run controls.csv through the submission gate."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
