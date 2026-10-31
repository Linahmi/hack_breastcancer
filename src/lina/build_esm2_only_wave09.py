"""Build higher-mutation ESM2-only designs to break the count confound."""

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
DESIGN_SPECS = (
    ((52, "N"), (139, "L"), (149, "R"), (156, "L"), (177, "E"), (218, "E"), (230, "A")),
    ((9, "P"), (68, "A"), (76, "S"), (100, "S"), (119, "V"), (136, "S"), (157, "S"), (241, "V")),
    ((27, "Y"), (120, "G"), (138, "V"), (150, "P"), (155, "Q"), (186, "S"), (195, "A"), (227, "S"), (239, "L")),
    ((1, "Q"), (31, "S"), (48, "M"), (78, "M"), (91, "S"), (153, "T"), (225, "G"), (235, "S"), (240, "T")),
)


def cdr_for(position: int) -> str | None:
    for name, (start, end) in CDR_RANGES.items():
        if start <= position <= end:
            return name
    return None


def main() -> int:
    candidates_path = REPO / "results/candidates.csv"
    with candidates_path.open(newline="") as handle:
        existing = list(csv.DictReader(handle))
    by_name = {row["name"]: row for row in existing}

    baseline_rows = [by_name[f"esm2_only_{index:02d}"] for index in range(1, 9)]
    allowed_mutations = {
        (int(mutation[1:-1]), mutation[-1])
        for row in baseline_rows
        for mutation in mutations(row["sequence"])
    }
    if len({position for position, _ in allowed_mutations}) != 40:
        raise ValueError("baseline ESM2-only arm does not contain 40 positions")

    esm2 = np.load(REPO / "src/lina/esm2_logprobs.npy")
    ablang2 = np.load(REPO / "src/lina/ablang2_logprobs.npy")
    distances = np.load(REPO / "structure/ca_distances.npy")
    rows = []
    selected_positions = []

    for index, spec in enumerate(DESIGN_SPECS, 9):
        positions = tuple(position for position, _ in spec)
        if len(positions) not in (7, 8, 9):
            raise ValueError(f"bad mutation count for esm2_only_{index:02d}")
        if not set(spec) <= allowed_mutations:
            raise ValueError(f"mutation outside baseline ESM2-only set: {spec}")
        if set(positions) & HARD_EXCLUSIONS:
            raise ValueError(f"hard exclusion mutated: {positions}")

        cdr_counts = Counter(filter(None, (cdr_for(position) for position in positions)))
        if any(count > 2 for count in cdr_counts.values()):
            raise ValueError(f"CDR cap failure: {cdr_counts}")
        for left, right in itertools.combinations(positions, 2):
            distance = float(distances[left - 1, right - 1])
            if np.isfinite(distance) and distance < 10.0:
                raise ValueError(f"CA distance {left}-{right}={distance:.3f} A")

        sequence = list(PARENT)
        for position, mutant in spec:
            wild_type_index = CANON.index(PARENT[position - 1])
            mutant_index = CANON.index(mutant)
            if esm2[position - 1, mutant_index] <= esm2[position - 1, wild_type_index]:
                raise ValueError(f"ESM2 does not prefer {position}{mutant}")
            if any(
                ablang2[position - 1, aa_index]
                > ablang2[position - 1, wild_type_index]
                for aa_index in range(len(CANON))
                if aa_index != wild_type_index
            ):
                raise ValueError(f"AbLang2 also flags position {position}")
            sequence[position - 1] = mutant

        name = f"esm2_only_{index:02d}"
        sequence_text = "".join(sequence)
        errors = check(name, sequence_text)
        if errors:
            raise ValueError(f"{name}: {errors}")
        mutation_list = ";".join(mutations(sequence_text))
        rows.append(
            {
                "name": name,
                "sequence": sequence_text,
                "rationale": (
                    f"{mutation_list} - ESM2-only, AbLang2 disagrees; mutation "
                    "count deliberately raised to break the confound Aida identified "
                    "between ESM2-only position choice and the five-mutation minimum"
                ),
                "tool": "esm2_only",
            }
        )
        selected_positions.append(set(positions))

    if len(set().union(*selected_positions)) != sum(map(len, selected_positions)):
        raise ValueError("wave-09 designs reuse positions across the four rows")
    if set(by_name) & {row["name"] for row in rows}:
        raise ValueError("wave-09 name collision")
    if {row["sequence"] for row in existing} & {row["sequence"] for row in rows}:
        raise ValueError("wave-09 sequence collision")

    with candidates_path.open("a", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["name", "sequence", "rationale", "tool"]
        )
        writer.writerows(rows)

    wave_path = REPO / "results/wave_09.csv"
    with wave_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["name", "sequence"])
        writer.writeheader()
        writer.writerows(
            {"name": row["name"], "sequence": row["sequence"]} for row in rows
        )

    print(f"appended {len(rows)} designs; candidates total={len(existing) + len(rows)}")
    print(f"wrote {wave_path.relative_to(REPO)} with {len(rows)} rows")
    for row in rows:
        print(f"{row['name']}: {row['rationale']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
