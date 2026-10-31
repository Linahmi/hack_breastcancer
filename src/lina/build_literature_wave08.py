"""Build four matched D102W pairs and two pass-supported literature designs."""

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


PAIR_SOURCES = (
    "consensus_green_03",
    "consensus_green_07",
    "consensus_green_10",
    "consensus_green_09",
)
PASSING_CONSENSUS = {
    "w01_consensus_green_02",
    "consensus_green_03",
    "consensus_green_04",
    "consensus_green_05",
    "consensus_green_06",
    "consensus_green_07",
    "consensus_green_09",
    "consensus_green_10",
}
HOTSPOT_EXCLUDED_SPECS = (
    ((54, "G"), (102, "S"), (107, "L"), (164, "I"), (190, "Q")),
    ((32, "Y"), (55, "G"), (107, "L"), (164, "I"), (190, "Q")),
)
HARD_EXCLUSIONS = set(range(121, 136)) | {22, 96, 158, 223}
ALANINE_HOTSPOTS = {50, 99, 105, 226}
CDR_RANGES = {
    "H1": (31, 35),
    "H2": (50, 66),
    "H3": (99, 109),
    "L1": (159, 169),
    "L2": (185, 191),
    "L3": (224, 232),
}


def cdr_for(position: int) -> str | None:
    for name, (start, end) in CDR_RANGES.items():
        if start <= position <= end:
            return name
    return None


def validate_design(
    name: str,
    sequence: str,
    distances: np.ndarray,
    *,
    enforce_ca_floor: bool,
) -> None:
    errors = check(name, sequence)
    if errors:
        raise ValueError(f"{name}: {errors}")
    positions = tuple(
        index
        for index, (wild_type, mutant) in enumerate(
            zip(PARENT, sequence), start=1
        )
        if wild_type != mutant
    )
    if set(positions) & HARD_EXCLUSIONS:
        raise ValueError(f"{name}: hard exclusion mutated")
    cdr_counts = Counter(filter(None, (cdr_for(position) for position in positions)))
    if any(count > 2 for count in cdr_counts.values()):
        raise ValueError(f"{name}: CDR cap failure {cdr_counts}")
    for left, right in itertools.combinations(positions, 2):
        distance = float(distances[left - 1, right - 1])
        if enforce_ca_floor and np.isfinite(distance) and distance < 10.0:
            raise ValueError(
                f"{name}: CA distance {left}-{right}={distance:.3f} A"
            )


def main() -> int:
    candidates_path = REPO / "results/candidates.csv"
    with candidates_path.open(newline="") as handle:
        existing = list(csv.DictReader(handle))
    by_name = {row["name"]: row for row in existing}
    distances = np.load(REPO / "structure/ca_distances.npy")

    rows = []
    for index, source_name in enumerate(PAIR_SOURCES, 1):
        source = by_name[source_name]
        if source["tool"] != "consensus_green":
            raise ValueError(f"{source_name} is not consensus_green")
        source_sequence = source["sequence"]
        if source_sequence[101] != "S" or PARENT[101] != "D":
            raise ValueError(f"{source_name} does not carry D102S")
        sequence = source_sequence[:101] + "W" + source_sequence[102:]
        name = f"lit_d102w_pair_{index:02d}"
        # Exact matching takes priority: all passing D102S source backgrounds
        # predate the 10 A heuristic and carry at least one sub-10 A pair.
        validate_design(name, sequence, distances, enforce_ca_floor=False)
        mutation_list = ";".join(mutations(sequence))
        rows.append(
            {
                "name": name,
                "sequence": sequence,
                "rationale": (
                    f"{mutation_list} - matched D102W pair for {source_name}; "
                    "identical background with only D102S changed to D102W; "
                    "Gerstner et al., J Mol Biol 2002, "
                    "doi:10.1016/S0022-2836(02)00677-0"
                ),
                "tool": "literature",
            }
        )

    passing_rows = [by_name[name] for name in PASSING_CONSENSUS]
    passing_mutations = {
        mutation
        for row in passing_rows
        for mutation in mutations(row["sequence"])
    }
    for index, spec in enumerate(HOTSPOT_EXCLUDED_SPECS, 1):
        positions = {position for position, _ in spec}
        if positions & ALANINE_HOTSPOTS:
            raise ValueError(f"hot spot included: {positions & ALANINE_HOTSPOTS}")
        sequence = list(PARENT)
        for position, mutant in spec:
            mutation = f"{PARENT[position - 1]}{position}{mutant}"
            if mutation not in passing_mutations:
                raise ValueError(f"not supported by passing consensus: {mutation}")
            sequence[position - 1] = mutant
        sequence_text = "".join(sequence)
        name = f"lit_hotspot_excluded_{index:02d}"
        validate_design(name, sequence_text, distances, enforce_ca_floor=True)
        mutation_list = ";".join(mutations(sequence_text))
        rows.append(
            {
                "name": name,
                "sequence": sequence_text,
                "rationale": (
                    f"{mutation_list} - literature hot-spot-excluded design; "
                    "every filler mutation occurred in a passing consensus_green "
                    "design; excludes R50, W99, Y105 and H226 alanine-scan hot spots"
                ),
                "tool": "literature",
            }
        )

    existing_names = set(by_name)
    existing_sequences = {row["sequence"] for row in existing}
    for row in rows:
        if row["name"] in existing_names:
            raise ValueError(f"duplicate name: {row['name']}")
        if row["sequence"] in existing_sequences:
            raise ValueError(f"duplicate sequence: {row['name']}")

    with candidates_path.open("a", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["name", "sequence", "rationale", "tool"]
        )
        writer.writerows(rows)

    wave_path = REPO / "results/wave_08.csv"
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
