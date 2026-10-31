"""Build eight diverse ESM2-only designs and wave 07."""

import csv
import itertools
import sys
from collections import Counter
from pathlib import Path

import numpy as np


REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO))

from scripts.parent import PARENT  # noqa: E402
from scripts.validate import check, mutations  # noqa: E402


CANON = "ACDEFGHIKLMNPQRSTVWY"
HARD_EXCLUSIONS = set(range(121, 136)) | {22, 96, 158, 223}
CDR_RANGES = {
    "H1": (31, 35),
    "H2": (50, 66),
    "H3": (99, 109),
    "L1": (159, 169),
    "L2": (185, 191),
    "L3": (224, 232),
}

# Forty distinct positions: five per design, no position reused within the arm.
# The four lowest-margin eligible positions (35, 84, 88, 224) are omitted.
DESIGN_SPECS = (
    ((52, "N"), (149, "R"), (156, "L"), (177, "E"), (218, "E")),
    ((76, "S"), (119, "V"), (136, "S"), (157, "S"), (241, "V")),
    ((138, "V"), (155, "Q"), (195, "A"), (227, "S"), (239, "L")),
    ((1, "Q"), (48, "M"), (153, "T"), (225, "G"), (240, "T")),
    ((9, "P"), (78, "M"), (212, "N"), (220, "D"), (232, "Y")),
    ((31, "S"), (68, "A"), (120, "G"), (139, "L"), (215, "A")),
    ((100, "S"), (162, "S"), (186, "S"), (213, "V"), (235, "S")),
    ((27, "Y"), (91, "S"), (144, "A"), (150, "P"), (230, "A")),
)


def cdr_for(position: int) -> str | None:
    for name, (start, end) in CDR_RANGES.items():
        if start <= position <= end:
            return name
    return None


def flagged_positions(scores: np.ndarray) -> set[int]:
    flagged = set()
    scored_positions = (
        set(range(1, 121)) | set(range(136, len(PARENT) + 1))
    ) - HARD_EXCLUSIONS
    for position in scored_positions:
        wild_type = CANON.index(PARENT[position - 1])
        if any(
            scores[position - 1, index] > scores[position - 1, wild_type]
            for index in range(len(CANON))
            if index != wild_type
        ):
            flagged.add(position)
    return flagged


def main() -> int:
    esm2 = np.load(REPO / "src/lina/esm2_logprobs.npy")
    ablang2 = np.load(REPO / "src/lina/ablang2_logprobs.npy")
    distances = np.load(REPO / "structure/ca_distances.npy")

    esm2_flagged = flagged_positions(esm2)
    ablang2_flagged = flagged_positions(ablang2)
    both_flagged = esm2_flagged & ablang2_flagged
    eligible = esm2_flagged - ablang2_flagged
    if (len(esm2_flagged), len(both_flagged), len(eligible)) != (79, 35, 44):
        raise ValueError(
            "unexpected model position counts: "
            f"ESM2={len(esm2_flagged)}, both={len(both_flagged)}, "
            f"ESM2-only={len(eligible)}"
        )

    candidates_path = REPO / "results/candidates.csv"
    with candidates_path.open(newline="") as handle:
        existing = list(csv.DictReader(handle))
    existing_names = {row["name"] for row in existing}
    planned_names = {f"esm2_only_{index:02d}" for index in range(1, 9)}
    collisions = existing_names & planned_names
    if collisions:
        raise SystemExit(f"refusing to append existing designs: {sorted(collisions)}")

    previously_touched = {
        position
        for row in existing
        for position, (wild_type, mutant) in enumerate(
            zip(PARENT, row["sequence"]), start=1
        )
        if wild_type != mutant
    }
    selected_positions = {
        position for design in DESIGN_SPECS for position, _ in design
    }
    if len(selected_positions) != 40:
        raise ValueError("ESM2-only positions are not unique across the arm")
    if selected_positions & previously_touched:
        raise ValueError(
            f"ESM2-only arm reuses pool positions: "
            f"{sorted(selected_positions & previously_touched)}"
        )

    rows = []
    for index, design in enumerate(DESIGN_SPECS, 1):
        positions = tuple(position for position, _ in design)
        if not set(positions) <= eligible:
            raise ValueError(f"non-ESM2-only position in design {index}: {positions}")

        cdr_counts = Counter(filter(None, (cdr_for(position) for position in positions)))
        if any(count > 2 for count in cdr_counts.values()):
            raise ValueError(f"CDR cap failure in design {index}: {cdr_counts}")

        for left, right in itertools.combinations(positions, 2):
            distance = float(distances[left - 1, right - 1])
            if np.isfinite(distance) and distance < 10.0:
                raise ValueError(
                    f"CA distance failure in design {index}: "
                    f"{left}-{right}={distance:.3f} A"
                )

        sequence = list(PARENT)
        for position, mutant in design:
            wild_type_index = CANON.index(PARENT[position - 1])
            mutant_index = CANON.index(mutant)
            non_wild_type = [
                aa_index
                for aa_index in range(len(CANON))
                if aa_index != wild_type_index
            ]
            best_esm2 = max(
                non_wild_type, key=lambda aa_index: esm2[position - 1, aa_index]
            )
            if mutant_index != best_esm2:
                raise ValueError(
                    f"{PARENT[position - 1]}{position}{mutant} is not the "
                    "top ESM2 non-WT residue"
                )
            if esm2[position - 1, mutant_index] <= esm2[position - 1, wild_type_index]:
                raise ValueError(f"ESM2 does not prefer mutation at {position}")
            if any(
                ablang2[position - 1, aa_index]
                > ablang2[position - 1, wild_type_index]
                for aa_index in non_wild_type
            ):
                raise ValueError(f"AbLang2 also flags position {position}")
            sequence[position - 1] = mutant

        sequence_text = "".join(sequence)
        name = f"esm2_only_{index:02d}"
        errors = check(name, sequence_text)
        if errors:
            raise ValueError(f"{name}: {errors}")
        mutation_list = ";".join(mutations(sequence_text))
        rows.append(
            {
                "name": name,
                "sequence": sequence_text,
                "rationale": (
                    f"{mutation_list} - ESM2-only, AbLang2 disagrees; tests whether "
                    "requiring model consensus was justified or merely conservative"
                ),
                "tool": "esm2_only",
            }
        )

    with candidates_path.open("a", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["name", "sequence", "rationale", "tool"]
        )
        writer.writerows(rows)

    wave_path = REPO / "results/wave_07.csv"
    with wave_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["name", "sequence"])
        writer.writeheader()
        writer.writerows(
            {"name": row["name"], "sequence": row["sequence"]} for row in rows
        )

    print(
        f"ESM2 flagged {len(esm2_flagged)} allowed positions; "
        f"both models flagged {len(both_flagged)}; ESM2-only={len(eligible)}"
    )
    print(f"appended {len(rows)} designs; candidates total={len(existing) + len(rows)}")
    print(f"wrote {wave_path.relative_to(REPO)} with {len(rows)} rows")
    for row in rows:
        print(f"{row['name']}: {row['rationale']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
