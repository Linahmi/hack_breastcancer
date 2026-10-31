"""Build an eight-design first wave under the current spreading rule.

The published wave_01.csv predates the CA-distance matrix and is retained as
legacy work. This builder now enforces the >=10 A CA rule and therefore is not
expected to reproduce that legacy file.
"""

import csv
import itertools
import random
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from Bio.Align import PairwiseAligner, substitution_matrices
from anarci.germlines import all_germlines


REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO))

from scripts.parent import PARENT, VH, VL  # noqa: E402
from scripts.validate import check, mutations  # noqa: E402


CANON = "ACDEFGHIKLMNPQRSTVWY"
AA_INDEX = {aa: i for i, aa in enumerate(CANON)}
HARD_EXCLUSIONS = set(range(121, 136)) | {22, 96, 158, 223}
MIN_CA_DISTANCE_A = 10.0
CA_DISTANCES = np.load(REPO / "structure/ca_distances.npy")
CDR_RANGES = {
    "H1": (31, 35),
    "H2": (50, 66),
    "H3": (99, 109),
    "L1": (159, 169),
    "L2": (185, 191),
    "L3": (224, 232),
}


@dataclass(frozen=True)
class MutationChoice:
    position: int
    mutant: str
    quality: float


@dataclass(frozen=True)
class Design:
    sequence: str
    choices: tuple[MutationChoice, ...]
    quality: float
    spread: int


def cdr_for(position: int) -> str | None:
    for name, (start, end) in CDR_RANGES.items():
        if start <= position <= end:
            return name
    return None


def region_label(position: int) -> str:
    loop = cdr_for(position)
    if loop:
        return loop
    return "VH-framework" if position <= VH[1] else "VL-framework"


def rank_normalised(scores: np.ndarray) -> np.ndarray:
    """Best score gets 1, worst gets 0; ranks are within one position."""
    order = np.argsort(-scores, kind="mergesort")
    ranks = np.empty(len(scores), dtype=float)
    ranks[order] = np.arange(1, len(scores) + 1)
    return (len(scores) - ranks) / (len(scores) - 1)


def build_choice_table(
    esm_logp: np.ndarray, ablang_logp: np.ndarray
) -> tuple[dict[int, list[MutationChoice]], set[int]]:
    choices: dict[int, list[MutationChoice]] = {}
    position_consensus: set[int] = set()
    scored_positions = set(range(1, 121)) | set(range(136, 243))
    for position in sorted(scored_positions - HARD_EXCLUSIONS):
        wt_index = AA_INDEX[PARENT[position - 1]]
        mutant_indices = [i for i in range(len(CANON)) if i != wt_index]
        esm_scores = esm_logp[position - 1, mutant_indices] - esm_logp[position - 1, wt_index]
        ablang_scores = (
            ablang_logp[position - 1, mutant_indices]
            - ablang_logp[position - 1, wt_index]
        )
        if np.any(esm_scores > 0) and np.any(ablang_scores > 0):
            position_consensus.add(position)

        esm_rank = rank_normalised(esm_scores)
        ablang_rank = rank_normalised(ablang_scores)
        shared = []
        for k, aa_index in enumerate(mutant_indices):
            if esm_scores[k] > 0 and ablang_scores[k] > 0:
                shared.append(
                    MutationChoice(
                        position,
                        CANON[aa_index],
                        float((esm_rank[k] + ablang_rank[k]) / 2),
                    )
                )
        if shared:
            choices[position] = sorted(
                shared, key=lambda choice: (-choice.quality, choice.mutant)
            )
    return choices, position_consensus


def germline_map() -> dict[int, str]:
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5

    references = (
        (
            PARENT[VH[0] - 1 : VH[1]],
            all_germlines["V"]["H"]["human"]["IGHV3-66*01"].replace("-", ""),
            VH[0],
        ),
        (
            PARENT[VL[0] - 1 : VL[1]],
            all_germlines["V"]["K"]["human"]["IGKV1-39*01"].replace("-", ""),
            VL[0],
        ),
    )
    mapped: dict[int, str] = {}
    for parent_domain, germline, parent_offset in references:
        alignment = aligner.align(parent_domain, germline)[0]
        parent_blocks, germline_blocks = alignment.aligned
        for (parent_start, parent_end), (germ_start, germ_end) in zip(
            parent_blocks, germline_blocks
        ):
            if parent_end - parent_start != germ_end - germ_start:
                raise ValueError("aligned block lengths differ")
            for delta in range(parent_end - parent_start):
                mapped[parent_offset + parent_start + delta] = germline[
                    germ_start + delta
                ]
    return mapped


def respects_cdr_cap(positions: tuple[int, ...]) -> bool:
    counts = Counter(cdr_for(position) for position in positions if cdr_for(position))
    return all(count <= 2 for count in counts.values())


def respects_ca_distance_floor(positions: tuple[int, ...]) -> bool:
    """Require every mutated-position CA pair to be at least 10 A apart."""
    for left, right in itertools.combinations(positions, 2):
        distance = CA_DISTANCES[left - 1, right - 1]
        if not np.isfinite(distance) or distance < MIN_CA_DISTANCE_A:
            return False
    return True


def respects_spreading_rule(positions: tuple[int, ...]) -> bool:
    return respects_cdr_cap(positions) and respects_ca_distance_floor(positions)


def sequence_for(choices: tuple[MutationChoice, ...]) -> str:
    sequence = list(PARENT)
    for choice in choices:
        sequence[choice.position - 1] = choice.mutant
    return "".join(sequence)


def positional_spread(positions: tuple[int, ...]) -> int:
    return sum(
        abs(right - left)
        for left, right in itertools.combinations(sorted(positions), 2)
    )


def make_design(position_set: tuple[int, ...], choices_by_position) -> Design | None:
    if not respects_spreading_rule(position_set):
        return None
    selected = tuple(choices_by_position[position][0] for position in sorted(position_set))
    sequence = sequence_for(selected)
    if check("candidate", sequence):
        return None
    return Design(
        sequence=sequence,
        choices=selected,
        quality=sum(choice.quality for choice in selected),
        spread=positional_spread(position_set),
    )


def enumerate_arm(position_pool, choices_by_position, sizes=(6, 5)) -> list[Design]:
    designs = []
    seen = set()
    for size in sizes:
        for position_set in itertools.combinations(sorted(position_pool), size):
            design = make_design(position_set, choices_by_position)
            if design and design.sequence not in seen:
                seen.add(design.sequence)
                designs.append(design)
    return designs


def sample_aggressive(position_pool, choices_by_position) -> list[Design]:
    rng = random.Random(20261031)
    positions = sorted(position_pool)
    designs = []
    seen_sets = set()
    for _ in range(50000):
        position_set = tuple(sorted(rng.sample(positions, 9)))
        if position_set in seen_sets:
            continue
        seen_sets.add(position_set)
        design = make_design(position_set, choices_by_position)
        if design:
            designs.append(design)
    return designs


def hamming(left: str, right: str) -> int:
    return sum(a != b for a, b in zip(left, right))


def choose_two(candidates: list[Design], already_selected: list[Design], spread_first=False):
    if len(candidates) < 2:
        raise ValueError(f"only {len(candidates)} valid internal candidates")
    if spread_first:
        first = max(candidates, key=lambda d: (d.spread, d.quality))
    else:
        first = max(candidates, key=lambda d: (d.quality, d.spread))
    reference = already_selected + [first]
    remaining = [candidate for candidate in candidates if candidate.sequence != first.sequence]
    second = max(
        remaining,
        key=lambda d: (
            min(hamming(d.sequence, prior.sequence) for prior in reference),
            d.spread if spread_first else d.quality,
            d.quality,
        ),
    )
    return first, second


def mutation_text(design: Design) -> str:
    return ";".join(mutations(design.sequence))


def rationale(arm: str, design: Design) -> str:
    labels = []
    for choice in design.choices:
        label = region_label(choice.position)
        if label not in labels:
            labels.append(label)
    spread = "/".join(labels)
    clauses = {
        "consensus_green": (
            "GREEN consensus, both models prefer, max 2 per CDR and all mutated "
            f"CAs >=10 A apart, spread across {spread}"
        ),
        "contact_test": (
            "RED contact test, both models prefer, max 2 per CDR and all mutated "
            f"CAs >=10 A apart, spread across {spread}"
        ),
        "germline_revert": (
            "germline reversion toward IGHV3-66/IGKV1-39, both models prefer, "
            f"max 2 per CDR and all mutated CAs >=10 A apart, spread across {spread}"
        ),
        "aggressive": (
            f"aggressive consensus, both models prefer, {len(design.choices)} mutations, "
            "max 2 per CDR and all mutated CAs >=10 A apart, "
            f"maximally spread across {spread}"
        ),
    }
    return f"{mutation_text(design)} - {clauses[arm]}"


def main() -> int:
    candidates_path = REPO / "results/candidates.csv"
    existing = list(csv.DictReader(candidates_path.open(newline="")))
    if existing:
        raise SystemExit(
            f"refusing to append wave 01: results/candidates.csv already has {len(existing)} rows"
        )

    esm_logp = np.load(REPO / "src/lina/esm2_logprobs.npy")
    ablang_logp = np.load(REPO / "src/lina/ablang2_logprobs.npy")
    choices, consensus_positions = build_choice_table(esm_logp, ablang_logp)
    if len(consensus_positions) != 35:
        raise ValueError(f"expected 35 consensus positions, found {len(consensus_positions)}")

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
        germline_aa = germline.get(position)
        if germline_aa is None or germline_aa == PARENT[position - 1]:
            continue
        matching = [choice for choice in position_choices if choice.mutant == germline_aa]
        if matching:
            germline_choices[position] = matching

    print(f"consensus positions: {len(consensus_positions)}")
    print(f"positions with at least one shared qualifying residue: {len(choices)}")
    print(f"GREEN shared-residue positions: {len(green_positions)}")
    print(f"RED shared-residue positions: {len(red_positions)}")
    print(f"germline shared-residue positions: {len(germline_choices)}")
    print(
        "germline candidates: "
        + ",".join(
            f"{PARENT[p-1]}{p}{germline_choices[p][0].mutant}"
            for p in sorted(germline_choices)
        )
    )

    arm_candidates = {
        "consensus_green": enumerate_arm(green_positions, choices),
        "contact_test": enumerate_arm(red_positions, choices),
        "germline_revert": enumerate_arm(germline_choices, germline_choices),
        "aggressive": sample_aggressive(set(choices), choices),
    }
    for arm, designs in arm_candidates.items():
        print(f"{arm} valid internal candidates: {len(designs)}")

    selected_by_arm = {}
    selected_so_far = []
    for arm in ("consensus_green", "contact_test", "germline_revert", "aggressive"):
        selected = choose_two(
            arm_candidates[arm], selected_so_far, spread_first=(arm == "aggressive")
        )
        selected_by_arm[arm] = selected
        selected_so_far.extend(selected)

    # Interleave arms so the first wave remains balanced even if read top-down.
    ordered = []
    for replicate in range(2):
        for arm in ("consensus_green", "contact_test", "germline_revert", "aggressive"):
            ordered.append((arm, replicate + 1, selected_by_arm[arm][replicate]))

    full_rows = []
    wave_rows = []
    for arm, replicate, design in ordered:
        name = f"w01_{arm}_{replicate:02d}"
        errors = check(name, design.sequence)
        if errors:
            raise ValueError(f"{name} failed validation before write: {errors}")
        wave_rows.append({"name": name, "sequence": design.sequence})
        full_rows.append(
            {
                "name": name,
                "sequence": design.sequence,
                "rationale": rationale(arm, design),
                "tool": arm,
            }
        )

    wave_path = REPO / "results/wave_01.csv"
    with wave_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["name", "sequence"])
        writer.writeheader()
        writer.writerows(wave_rows)

    with candidates_path.open("a", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["name", "sequence", "rationale", "tool"]
        )
        writer.writerows(full_rows)

    print(f"wrote {wave_path.relative_to(REPO)} with {len(wave_rows)} rows")
    print(f"appended {len(full_rows)} rows to {candidates_path.relative_to(REPO)}")
    for row in full_rows:
        print(f"{row['name']}: {row['rationale']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
