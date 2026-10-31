"""Build four measured-mutation designs without changing prior designs.

The literature mutations are used exactly as reported. Consensus fillers are
restricted to positions and residues preferred over wild type by both ESM2 and
AbLang2, and every design obeys the current CDR and CA-distance spreading rule.
"""

import csv
import sys
from pathlib import Path

import numpy as np


REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO))

from scripts.parent import PARENT  # noqa: E402
from scripts.validate import check, mutations  # noqa: E402
from src.lina.build_wave01 import (  # noqa: E402
    MutationChoice,
    build_choice_table,
    respects_spreading_rule,
    sequence_for,
)


LITERATURE_MUTATIONS = {
    (50, "A"),
    (99, "A"),
    (102, "W"),
    (105, "A"),
    (226, "A"),
}

# Each tuple is (two exact literature mutations, three both-model consensus
# fillers). The four designs collectively cover all five literature rows.
DESIGN_SPECS = (
    (((50, "A"), (102, "W")), ((55, "G"), (168, "L"), (190, "Q"))),
    (((50, "A"), (99, "A")), ((54, "G"), (164, "I"), (169, "S"))),
    (((50, "A"), (105, "A")), ((55, "G"), (164, "I"), (190, "Q"))),
    (((50, "A"), (226, "A")), ((54, "G"), (102, "S"), (190, "Q"))),
)


def load_literature() -> dict[tuple[int, str], dict[str, str]]:
    path = REPO / "structure/literature_mutations.csv"
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))

    by_mutation = {}
    for row in rows:
        position = int(row["position"])
        if row["wt_aa"] != PARENT[position - 1]:
            raise ValueError(
                f"literature numbering mismatch at {position}: "
                f"{row['wt_aa']} != {PARENT[position - 1]}"
            )
        by_mutation[(position, row["mut_aa"])] = row
    if set(by_mutation) != LITERATURE_MUTATIONS:
        raise ValueError(f"unexpected literature mutation set: {set(by_mutation)}")
    return by_mutation


def short_citation(row: dict[str, str]) -> str:
    if "Kelley RF" in row["source"]:
        return (
            "Kelley & O'Connell, Biochemistry 1993, "
            "doi:10.1021/bi00078a005 (>=100-fold affinity loss for alanine)"
        )
    if "Gerstner RB" in row["source"]:
        return (
            "Gerstner et al., J Mol Biol 2002, "
            "doi:10.1016/S0022-2836(02)00677-0 (~3-fold improved affinity)"
        )
    raise ValueError(f"no short citation for {row['source']}")


def main() -> int:
    literature = load_literature()
    choices, _ = build_choice_table(
        np.load(REPO / "src/lina/esm2_logprobs.npy"),
        np.load(REPO / "src/lina/ablang2_logprobs.npy"),
    )

    rows = []
    for index, (literature_specs, consensus_specs) in enumerate(DESIGN_SPECS, 1):
        selected = []
        for position, mutant in literature_specs:
            if (position, mutant) not in literature:
                raise ValueError(f"missing literature mutation {position}{mutant}")
            selected.append(MutationChoice(position, mutant, 0.0))

        for position, mutant in consensus_specs:
            matching = [choice for choice in choices.get(position, ()) if choice.mutant == mutant]
            if not matching:
                raise ValueError(
                    f"{PARENT[position - 1]}{position}{mutant} is not preferred "
                    "over wild type by both models"
                )
            selected.append(matching[0])

        selected_tuple = tuple(selected)
        positions = tuple(choice.position for choice in selected_tuple)
        if len(positions) not in (5, 6) or len(set(positions)) != len(positions):
            raise ValueError(f"bad mutation count/positions: {positions}")
        if not respects_spreading_rule(positions):
            raise ValueError(f"spreading-rule failure: {positions}")

        sequence = sequence_for(selected_tuple)
        name = f"literature_{index:02d}"
        errors = check(name, sequence)
        if errors:
            raise ValueError(f"{name}: {errors}")

        mutation_list = ";".join(mutations(sequence))
        citations = "; ".join(
            f"{PARENT[position - 1]}{position}{mutant}: "
            f"{short_citation(literature[(position, mutant)])}"
            for position, mutant in literature_specs
        )
        filler_list = ", ".join(
            f"{PARENT[position - 1]}{position}{mutant}"
            for position, mutant in consensus_specs
        )
        rows.append(
            {
                "name": name,
                "sequence": sequence,
                "rationale": (
                    f"{mutation_list} - measured literature mutations [{citations}]; "
                    f"both-model consensus fillers: {filler_list}; max 2 per CDR and "
                    "all mutated CAs >=10 A apart"
                ),
                "tool": "literature",
            }
        )

    candidates_path = REPO / "results/candidates.csv"
    with candidates_path.open(newline="") as handle:
        existing = list(csv.DictReader(handle))
    existing_names = {row["name"] for row in existing}
    collisions = existing_names & {row["name"] for row in rows}
    if collisions:
        raise SystemExit(f"refusing to append existing designs: {sorted(collisions)}")

    with candidates_path.open("a", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["name", "sequence", "rationale", "tool"]
        )
        writer.writerows(rows)

    for row in rows:
        print(f"{row['name']}: {row['rationale']}")
    print(f"appended {len(rows)} literature designs; candidates total {len(existing) + len(rows)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
