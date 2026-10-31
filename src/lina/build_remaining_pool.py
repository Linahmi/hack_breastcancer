"""Select the remaining 28 pool designs and prepare their wave assignment."""

import csv
import sys
from pathlib import Path

import numpy as np


REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO))

from src.lina.build_wave01 import (  # noqa: E402
    build_choice_table,
    enumerate_arm,
    germline_map,
    hamming,
    rationale,
    sample_aggressive,
)


REMAINING_QUOTAS = {
    "consensus_green": 8,
    "contact_test": 8,
    "germline_revert": 6,
    "aggressive": 6,
}


def main() -> int:
    candidates_path = REPO / "results/candidates.csv"
    existing = list(csv.DictReader(candidates_path.open(newline="")))
    if len(existing) != 8:
        raise SystemExit(f"expected exactly 8 wave-01 designs, found {len(existing)}")

    esm = np.load(REPO / "src/lina/esm2_logprobs.npy")
    ablang = np.load(REPO / "src/lina/ablang2_logprobs.npy")
    choices, consensus_positions = build_choice_table(esm, ablang)
    if len(consensus_positions) != 35:
        raise ValueError(len(consensus_positions))

    with (REPO / "structure/red_green.csv").open(newline="") as handle:
        structure_rows = list(csv.DictReader(handle))
    class_by_position = {int(row["position"]): row["class"] for row in structure_rows}
    green_positions = {
        position for position in choices if class_by_position.get(position) == "GREEN"
    }
    red_positions = {
        position for position in choices if class_by_position.get(position) == "RED"
    }
    germline = germline_map()
    germline_choices = {}
    for position, position_choices in choices.items():
        target = germline.get(position)
        if target is None:
            continue
        matching = [choice for choice in position_choices if choice.mutant == target]
        if target != choice_wt(position) and matching:
            germline_choices[position] = matching

    arm_candidates = {
        "consensus_green": enumerate_arm(green_positions, choices),
        "contact_test": enumerate_arm(red_positions, choices),
        "germline_revert": enumerate_arm(germline_choices, germline_choices),
        "aggressive": sample_aggressive(set(choices), choices),
    }

    chosen_sequences = [row["sequence"] for row in existing]
    selected = []
    arm_indices = {arm: 2 for arm in REMAINING_QUOTAS}
    remaining = dict(REMAINING_QUOTAS)

    slot_order = []
    while any(count > 0 for count in remaining.values()):
        for arm in ("consensus_green", "contact_test", "germline_revert", "aggressive"):
            if remaining[arm] > 0:
                slot_order.append(arm)
                remaining[arm] -= 1

    for slot, arm in enumerate(slot_order, start=1):
        available = [
            design
            for design in arm_candidates[arm]
            if design.sequence not in chosen_sequences
        ]
        if not available:
            raise RuntimeError(f"no candidates left for {arm}")

        ranked = sorted(
            available,
            key=lambda design: (
                min(hamming(design.sequence, sequence) for sequence in chosen_sequences),
                design.spread if arm == "aggressive" else design.quality,
                design.quality,
            ),
            reverse=True,
        )
        shortlist = ranked[:4]
        chosen = max(
            shortlist,
            key=lambda design: (
                min(hamming(design.sequence, sequence) for sequence in chosen_sequences),
                design.spread if arm == "aggressive" else design.quality,
                design.quality,
            ),
        )
        min_distance = min(
            hamming(chosen.sequence, sequence) for sequence in chosen_sequences
        )
        arm_indices[arm] += 1
        name = f"{arm}_{arm_indices[arm]:02d}"
        selected.append(
            {
                "name": name,
                "sequence": chosen.sequence,
                "rationale": rationale(arm, chosen),
                "tool": arm,
                "slot": slot,
                "shortlist_size": len(shortlist),
                "min_hamming_at_selection": min_distance,
            }
        )
        chosen_sequences.append(chosen.sequence)
        print(
            f"slot {slot:02d} {arm:17s} shortlist={len(shortlist)} "
            f"chosen_min_hamming={min_distance:2d} {name}"
        )

    if len(selected) != 28:
        raise ValueError(f"expected 28 selections, got {len(selected)}")

    # The first three dispatch waves are balanced 2/2/2/2 across the arms.
    # Four designs remain after three full eight-row waves because 36 total
    # designs is not divisible by eight; keep those four explicit in wave 05.
    for index, row in enumerate(selected):
        row["wave"] = 2 + index // 8

    plan_path = REPO / "src/lina/remaining_plan.csv"
    with plan_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "wave",
                "slot",
                "name",
                "sequence",
                "rationale",
                "tool",
                "shortlist_size",
                "min_hamming_at_selection",
            ],
        )
        writer.writeheader()
        writer.writerows(selected)

    with candidates_path.open("a", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["name", "sequence", "rationale", "tool"]
        )
        writer.writerows(
            {
                "name": row["name"],
                "sequence": row["sequence"],
                "rationale": row["rationale"],
                "tool": row["tool"],
            }
            for row in selected
        )

    print(f"prepared {plan_path.relative_to(REPO)} with {len(selected)} rows")
    print("candidates total: 36")
    return 0


def choice_wt(position: int) -> str:
    from scripts.parent import PARENT

    return PARENT[position - 1]


if __name__ == "__main__":
    raise SystemExit(main())
