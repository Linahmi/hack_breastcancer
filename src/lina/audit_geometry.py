"""Write a machine-readable geometry audit for every candidate."""

import csv
import itertools
import sys
from collections import Counter
from pathlib import Path

import numpy as np


REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO))

from scripts.parent import PARENT  # noqa: E402


CDR_RANGES = {
    "H1": (31, 35),
    "H2": (50, 66),
    "H3": (99, 109),
    "L1": (159, 169),
    "L2": (185, 191),
    "L3": (224, 232),
}
CA_COVARIATE_THRESHOLD_A = 10.0


def cdr_for(position: int) -> str | None:
    for name, (start, end) in CDR_RANGES.items():
        if start <= position <= end:
            return name
    return None


def main() -> int:
    distances = np.load(REPO / "structure/ca_distances.npy")
    if distances.shape != (len(PARENT), len(PARENT)):
        raise ValueError(f"unexpected CA matrix shape: {distances.shape}")

    wave_by_name: dict[str, str] = {}
    for wave_path in sorted((REPO / "results").glob("wave_[0-9][0-9].csv")):
        with wave_path.open(newline="") as handle:
            for row in csv.DictReader(handle):
                if row["name"] in wave_by_name:
                    raise ValueError(f"candidate appears in multiple waves: {row['name']}")
                wave_by_name[row["name"]] = wave_path.stem

    with (REPO / "results/candidates.csv").open(newline="") as handle:
        candidates = list(csv.DictReader(handle))

    output_rows = []
    for row in candidates:
        sequence = row["sequence"]
        if len(sequence) != len(PARENT):
            raise ValueError(f"{row['name']}: length {len(sequence)}")
        positions = [
            index
            for index, (wild_type, mutant) in enumerate(
                zip(PARENT, sequence), start=1
            )
            if wild_type != mutant
        ]
        mutation_list = ";".join(
            f"{PARENT[position - 1]}{position}{sequence[position - 1]}"
            for position in positions
        )

        cdr_counts = Counter(filter(None, (cdr_for(position) for position in positions)))
        max_cdr_count = max(cdr_counts.values(), default=0)
        cdr_cap_pass = max_cdr_count <= 2

        pair_distances = []
        below_threshold = []
        unresolved = []
        for left, right in itertools.combinations(positions, 2):
            distance = float(distances[left - 1, right - 1])
            pair_distances.append(distance)
            if not np.isfinite(distance):
                unresolved.append(f"{left}-{right}")
            elif distance < CA_COVARIATE_THRESHOLD_A:
                below_threshold.append(f"{left}-{right}:{distance:.3f}")

        finite_distances = [distance for distance in pair_distances if np.isfinite(distance)]
        min_distance = min(finite_distances, default=float("nan"))
        output_rows.append(
            {
                "name": row["name"],
                "wave": wave_by_name.get(row["name"], ""),
                "tool": row["tool"],
                "mutations": mutation_list,
                "mutation_count": len(positions),
                "max_mutations_in_one_cdr": max_cdr_count,
                "cdr_cap_at_most_2": str(cdr_cap_pass).lower(),
                "min_ca_distance_A": (
                    "nan" if not np.isfinite(min_distance) else f"{min_distance:.3f}"
                ),
                "pairs_below_10A": len(below_threshold),
                "pairs_without_ca": len(unresolved),
                "below_10A_pairs": ";".join(below_threshold),
                "pairs_without_ca_detail": ";".join(unresolved),
            }
        )

    output_path = REPO / "results/geometry_audit.csv"
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output_rows[0]))
        writer.writeheader()
        writer.writerows(output_rows)

    print(f"wrote {output_path.relative_to(REPO)}: {len(output_rows)} covariate rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
