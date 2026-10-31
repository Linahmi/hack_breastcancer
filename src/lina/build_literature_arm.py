"""Rebuild the literature arm and write the priority wave-05 handoff.

Two designs use the affinity-improving D102W result with four both-model
consensus fillers. Two matched designs use consensus mutations while excluding
every published alanine-scan hot spot. Every design obeys the current CDR and
CA-distance spreading rule.
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


EXPECTED_LITERATURE_MUTATIONS = {
    (50, "A"),
    (99, "A"),
    (102, "W"),
    (105, "A"),
    (226, "A"),
}
AFFINITY_IMPROVING_MUTATION = (102, "W")
ALANINE_HOTSPOT_MUTATIONS = {
    (50, "A"),
    (99, "A"),
    (105, "A"),
    (226, "A"),
}
ALANINE_HOTSPOT_POSITIONS = {position for position, _ in ALANINE_HOTSPOT_MUTATIONS}

# Each tuple is (design mode, exact literature mutations, both-model consensus
# mutations). The consensus-only designs deliberately mirror the filler sets of
# the D102W designs so the literature arm contains two interpretable pairs.
DESIGN_SPECS = (
    (
        "d102w",
        (AFFINITY_IMPROVING_MUTATION,),
        ((55, "G"), (106, "Y"), (164, "I"), (190, "Q")),
    ),
    (
        "d102w",
        (AFFINITY_IMPROVING_MUTATION,),
        ((55, "G"), (107, "L"), (168, "L"), (190, "Q")),
    ),
    (
        "hotspot_avoid",
        (),
        ((32, "Y"), (55, "G"), (106, "Y"), (164, "I"), (190, "Q")),
    ),
    (
        "hotspot_avoid",
        (),
        ((32, "Y"), (55, "G"), (107, "L"), (168, "L"), (190, "Q")),
    ),
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
    if set(by_mutation) != EXPECTED_LITERATURE_MUTATIONS:
        raise ValueError(f"unexpected literature mutation set: {set(by_mutation)}")
    for mutation in ALANINE_HOTSPOT_MUTATIONS:
        if "loss" not in by_mutation[mutation]["reported_effect"]:
            raise ValueError(f"alanine hot spot is not recorded as loss: {mutation}")
    if (
        "improved affinity"
        not in by_mutation[AFFINITY_IMPROVING_MUTATION]["reported_effect"]
    ):
        raise ValueError("D102W is not recorded as affinity-improving")
    return by_mutation


def short_citation(row: dict[str, str]) -> str:
    if "Kelley RF" in row["source"]:
        return (
            "Kelley & O'Connell, Biochemistry 1993, "
            "doi:10.1021/bi00078a005 (alanine hot spots; dDG >3 kcal/mol loss)"
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
    hotspot_labels = ", ".join(
        f"{PARENT[position - 1]}{position}{mutant}"
        for position, mutant in sorted(ALANINE_HOTSPOT_MUTATIONS)
    )
    hotspot_citation = short_citation(
        literature[sorted(ALANINE_HOTSPOT_MUTATIONS)[0]]
    )

    for index, (mode, literature_specs, consensus_specs) in enumerate(DESIGN_SPECS, 1):
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
        if set(positions) & ALANINE_HOTSPOT_POSITIONS:
            raise ValueError(f"published alanine hot spot included: {positions}")
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
        filler_list = ", ".join(
            f"{PARENT[position - 1]}{position}{mutant}"
            for position, mutant in consensus_specs
        )
        if mode == "d102w":
            if literature_specs != (AFFINITY_IMPROVING_MUTATION,):
                raise ValueError(f"bad D102W specification: {literature_specs}")
            inclusion_citation = short_citation(
                literature[AFFINITY_IMPROVING_MUTATION]
            )
            rationale = (
                f"{mutation_list} - D102W published affinity improvement "
                f"[{inclusion_citation}]; both-model consensus fillers: "
                f"{filler_list}; published alanine-scan hot spots excluded "
                f"[{hotspot_labels}: {hotspot_citation}]; max 2 per CDR and "
                "all mutated CAs >=10 A apart"
            )
        elif mode == "hotspot_avoid":
            if literature_specs:
                raise ValueError(
                    f"avoidance design has literature mutation: {literature_specs}"
                )
            rationale = (
                f"{mutation_list} - consensus mutations, published hot spots "
                f"excluded [{hotspot_labels}: {hotspot_citation}]; both models "
                "prefer every substitution; max 2 per CDR and all mutated CAs "
                ">=10 A apart"
            )
        else:
            raise ValueError(f"unknown literature design mode: {mode}")
        rows.append(
            {
                "name": name,
                "sequence": sequence,
                "rationale": rationale,
                "tool": "literature",
            }
        )

    candidates_path = REPO / "results/candidates.csv"
    with candidates_path.open(newline="") as handle:
        existing = list(csv.DictReader(handle))
    replacements = {row["name"]: row for row in rows}
    existing_literature_names = {
        row["name"] for row in existing if row["tool"] == "literature"
    }
    if existing_literature_names != set(replacements):
        raise SystemExit(
            "expected exactly literature_01 through literature_04, found "
            f"{sorted(existing_literature_names)}"
        )
    rebuilt = [replacements.get(row["name"], row) for row in existing]

    with candidates_path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["name", "sequence", "rationale", "tool"]
        )
        writer.writeheader()
        writer.writerows(rebuilt)

    controls_path = REPO / "results/controls.csv"
    with controls_path.open(newline="") as handle:
        antidesigns = [
            row for row in csv.DictReader(handle) if row["name"] == "ctrl_antidesign"
        ]
    if len(antidesigns) != 1:
        raise ValueError(f"expected one ctrl_antidesign, found {len(antidesigns)}")
    wave_rows = [
        {"name": row["name"], "sequence": row["sequence"]} for row in rows
    ]
    wave_rows.append(
        {
            "name": antidesigns[0]["name"],
            "sequence": antidesigns[0]["sequence"],
        }
    )
    wave_path = REPO / "results/wave_05.csv"
    with wave_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["name", "sequence"])
        writer.writeheader()
        writer.writerows(wave_rows)

    for row in rows:
        print(f"{row['name']}: {row['rationale']}")
    print(f"rebuilt {len(rows)} literature designs; candidates total {len(rebuilt)}")
    print(f"wrote {wave_path.relative_to(REPO)} with {len(wave_rows)} rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
