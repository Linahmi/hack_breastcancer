"""Append the five-mutation jointly disfavoured anti-design control."""

import csv
import itertools
import sys
from pathlib import Path

import numpy as np


REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO))

from scripts.parent import PARENT  # noqa: E402
from scripts.validate import check, mutations  # noqa: E402


CANON = "ACDEFGHIKLMNPQRSTVWY"
AA_INDEX = {aa: i for i, aa in enumerate(CANON)}
HARD_EXCLUSIONS = set(range(121, 136)) | {22, 96, 158, 223}
CDR_RANGES = {
    "H1": (31, 35),
    "H2": (50, 66),
    "H3": (99, 109),
    "L1": (159, 169),
    "L2": (185, 191),
    "L3": (224, 232),
}


def percentile_disfavor(scores: np.ndarray) -> np.ndarray:
    """Most negative score gets 1; most positive gets 0."""
    order = np.argsort(scores, kind="mergesort")
    ranks = np.empty(len(scores), dtype=float)
    ranks[order] = np.arange(len(scores))
    return 1.0 - ranks / (len(scores) - 1)


def main() -> int:
    path = REPO / "results/controls.csv"
    rows = list(csv.DictReader(path.open(newline="")))
    if len(rows) != 3:
        raise SystemExit(f"expected 3 existing controls, found {len(rows)}")
    if any(row["name"] == "ctrl_antidesign" for row in rows):
        raise SystemExit("ctrl_antidesign already exists")

    esm = np.load(REPO / "src/lina/esm2_logprobs.npy")
    ablang = np.load(REPO / "src/lina/ablang2_logprobs.npy")
    records = []
    for loop, (start, end) in CDR_RANGES.items():
        for position in range(start, end + 1):
            if position in HARD_EXCLUSIONS:
                continue
            wt_index = AA_INDEX[PARENT[position - 1]]
            for aa_index, mutant in enumerate(CANON):
                if aa_index == wt_index or mutant == "C":
                    continue
                records.append(
                    {
                        "loop": loop,
                        "position": position,
                        "mutant": mutant,
                        "esm_score": float(esm[position - 1, aa_index] - esm[position - 1, wt_index]),
                        "ablang_score": float(
                            ablang[position - 1, aa_index] - ablang[position - 1, wt_index]
                        ),
                    }
                )

    esm_disfavor = percentile_disfavor(
        np.asarray([record["esm_score"] for record in records])
    )
    ablang_disfavor = percentile_disfavor(
        np.asarray([record["ablang_score"] for record in records])
    )
    by_loop = {loop: [] for loop in CDR_RANGES}
    for record, esm_rank, ablang_rank in zip(records, esm_disfavor, ablang_disfavor):
        if record["esm_score"] >= 0 or record["ablang_score"] >= 0:
            continue
        record["joint_disfavor"] = float((esm_rank + ablang_rank) / 2)
        by_loop[record["loop"]].append(record)
    for loop in by_loop:
        by_loop[loop].sort(
            key=lambda record: (
                -record["joint_disfavor"],
                record["position"],
                record["mutant"],
            )
        )
        by_loop[loop] = by_loop[loop][:10]

    best = None
    for loops in itertools.combinations(CDR_RANGES, 5):
        for selected in itertools.product(*(by_loop[loop] for loop in loops)):
            positions = [record["position"] for record in selected]
            if len(set(positions)) != 5:
                continue
            sequence = list(PARENT)
            for record in selected:
                sequence[record["position"] - 1] = record["mutant"]
            sequence = "".join(sequence)
            errors = check("ctrl_antidesign", sequence)
            if errors:
                continue
            objective = sum(record["joint_disfavor"] for record in selected)
            if best is None or objective > best[0]:
                best = (objective, selected, sequence)

    if best is None:
        raise RuntimeError("no valid five-loop anti-design found")
    _, selected, sequence = best
    mutation_list = ";".join(mutations(sequence))
    note = (
        f"{mutation_list} - This is a deliberate negative control for the submission - "
        "if it binds as well as our designed variants, our scores carried no signal."
    )
    row = {
        "name": "ctrl_antidesign",
        "sequence": sequence,
        "role": "negative",
        "length": len(sequence),
        "n_mut_vs_parent": len(mutations(sequence)),
        "note": note,
    }
    errors = check(row["name"], sequence)
    if errors:
        raise ValueError(errors)

    with path.open("a", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["name", "sequence", "role", "length", "n_mut_vs_parent", "note"],
        )
        writer.writerow(row)

    print("ctrl_antidesign validation: PASS")
    print(note)
    print("position wt mut loop esm2_score ablang2_score joint_disfavor")
    for record in sorted(selected, key=lambda item: item["position"]):
        position = record["position"]
        print(
            f"{position:8d} {PARENT[position-1]:2s} {record['mutant']:3s} "
            f"{record['loop']:4s} {record['esm_score']:.6f} "
            f"{record['ablang_score']:.6f} {record['joint_disfavor']:.6f}"
        )
    print(
        "cysteines: "
        + ", ".join(f"{position}={sequence[position-1]}" for position in (22, 96, 158, 223))
    )
    print(f"controls rows: {len(rows) + 1}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
