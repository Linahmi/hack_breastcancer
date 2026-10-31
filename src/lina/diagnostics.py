"""Numerical diagnostics for the ESM2 and AbLang2 masked-marginal scans.

This script reports diagnostics only. It does not select substitutions or
assemble candidates. Shannon entropy is calculated in bits after renormalising
each position over the 20 canonical amino acids.
"""

import csv
import sys
from pathlib import Path

import numpy as np


REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO))

from scripts.parent import PARENT, VH, VL  # noqa: E402


CANON = "ACDEFGHIKLMNPQRSTVWY"
AA_INDEX = {aa: i for i, aa in enumerate(CANON)}
STRUCTURAL_CYS = (22, 96, 158, 223)
CDR_RANGES = (
    (31, 35),
    (50, 66),
    (99, 109),
    (159, 169),
    (185, 191),
    (224, 232),
)
CDR_POSITIONS = {
    position
    for start, end in CDR_RANGES
    for position in range(start, end + 1)
}
DOMAIN_POSITIONS = set(range(VH[0], VH[1] + 1)) | set(range(VL[0], VL[1] + 1))
FRAMEWORK_POSITIONS = DOMAIN_POSITIONS - CDR_POSITIONS


def canonical_probabilities(logp: np.ndarray) -> np.ndarray:
    """Renormalise stored log probabilities over the 20 canonical residues."""
    shifted = logp - np.max(logp, axis=1, keepdims=True)
    weights = np.exp(shifted)
    return weights / weights.sum(axis=1, keepdims=True)


def entropy_bits(probabilities: np.ndarray) -> np.ndarray:
    return -(probabilities * np.log2(probabilities)).sum(axis=1)


def non_wt_outranks_wt(logp: np.ndarray, position: int) -> bool:
    wt_index = AA_INDEX[PARENT[position - 1]]
    competitors = np.delete(logp[position - 1], wt_index)
    return bool(np.max(competitors) > logp[position - 1, wt_index])


def average_ranks(values: np.ndarray) -> np.ndarray:
    """One-based average ranks, including exact-tie handling."""
    order = np.argsort(values, kind="mergesort")
    sorted_values = values[order]
    ranks = np.empty(len(values), dtype=float)
    start = 0
    while start < len(values):
        end = start + 1
        while end < len(values) and sorted_values[end] == sorted_values[start]:
            end += 1
        ranks[order[start:end]] = (start + 1 + end) / 2.0
        start = end
    return ranks


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    return float(np.corrcoef(average_ranks(x), average_ranks(y))[0, 1])


def score_pairs(
    esm_logp: np.ndarray,
    ablang_logp: np.ndarray,
    positions: set[int],
    include_wt: bool,
) -> tuple[np.ndarray, np.ndarray]:
    esm_scores = []
    ablang_scores = []
    for position in sorted(positions):
        wt_index = AA_INDEX[PARENT[position - 1]]
        for aa_index in range(len(CANON)):
            if not include_wt and aa_index == wt_index:
                continue
            esm_scores.append(
                esm_logp[position - 1, aa_index]
                - esm_logp[position - 1, wt_index]
            )
            ablang_scores.append(
                ablang_logp[position - 1, aa_index]
                - ablang_logp[position - 1, wt_index]
            )
    return np.asarray(esm_scores), np.asarray(ablang_scores)


def distribution(values: np.ndarray) -> dict[str, float]:
    quantiles = np.quantile(values, [0.0, 0.05, 0.25, 0.5, 0.75, 0.95, 1.0])
    return {
        "n": int(len(values)),
        "min": float(quantiles[0]),
        "p05": float(quantiles[1]),
        "p25": float(quantiles[2]),
        "median": float(quantiles[3]),
        "mean": float(np.mean(values)),
        "p75": float(quantiles[4]),
        "p95": float(quantiles[5]),
        "max": float(quantiles[6]),
    }


def main() -> int:
    esm_logp = np.load(REPO / "src/lina/esm2_logprobs.npy")
    ablang_logp = np.load(REPO / "src/lina/ablang2_logprobs.npy")
    if esm_logp.shape != (len(PARENT), len(CANON)):
        raise ValueError(f"unexpected ESM2 shape: {esm_logp.shape}")
    if ablang_logp.shape != (len(PARENT), len(CANON)):
        raise ValueError(f"unexpected AbLang2 shape: {ablang_logp.shape}")

    esm_prob = canonical_probabilities(esm_logp)
    # AbLang2's linker rows are zero placeholders and must never be used.
    ablang_prob = canonical_probabilities(ablang_logp)
    esm_entropy = entropy_bits(esm_prob)
    ablang_entropy = entropy_bits(ablang_prob)

    entropy_path = REPO / "src/lina/entropy.csv"
    with entropy_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(
            ["position", "wt_aa", "esm2_entropy", "ablang2_entropy"]
        )
        for position, wt_aa in enumerate(PARENT, start=1):
            ablang_value = (
                f"{ablang_entropy[position - 1]:.6f}"
                if position in DOMAIN_POSITIONS
                else ""
            )
            writer.writerow(
                [position, wt_aa, f"{esm_entropy[position - 1]:.6f}", ablang_value]
            )

    print("A. STRUCTURAL CYSTEINE RANKINGS")
    print("position model top_aa top_raw_p C_rank C_raw_p")
    for position in STRUCTURAL_CYS:
        c_index = AA_INDEX["C"]
        for model_name, logp in (("ESM2", esm_logp), ("AbLang2", ablang_logp)):
            row = logp[position - 1]
            top_index = int(np.argmax(row))
            c_rank = 1 + int(np.sum(row > row[c_index]))
            print(
                f"{position:8d} {model_name:7s} {CANON[top_index]:6s} "
                f"{np.exp(row[top_index]):.8f} {c_rank:6d} "
                f"{np.exp(row[c_index]):.8f}"
            )

    cdr_indices = np.asarray([p - 1 for p in sorted(CDR_POSITIONS)])
    framework_indices = np.asarray([p - 1 for p in sorted(FRAMEWORK_POSITIONS)])
    print("\nB. SHANNON ENTROPY (BITS; 20-WAY CANONICAL RENORMALISATION)")
    print("model group n mean_entropy")
    for model_name, values in (("ESM2", esm_entropy), ("AbLang2", ablang_entropy)):
        cdr_mean = float(np.mean(values[cdr_indices]))
        framework_mean = float(np.mean(values[framework_indices]))
        print(f"{model_name:7s} CDR       {len(cdr_indices):3d} {cdr_mean:.6f}")
        print(
            f"{model_name:7s} framework {len(framework_indices):3d} "
            f"{framework_mean:.6f}"
        )
        print(f"{model_name:7s} CDR_gt_framework {cdr_mean > framework_mean}")

    esm_flags = {
        position: non_wt_outranks_wt(esm_logp, position)
        for position in range(1, len(PARENT) + 1)
    }
    ablang_flags = {
        position: non_wt_outranks_wt(ablang_logp, position)
        for position in DOMAIN_POSITIONS
    }

    with (REPO / "structure/red_green.csv").open(newline="") as handle:
        structure_rows = list(csv.DictReader(handle))
    class_by_position = {int(row["position"]): row["class"] for row in structure_rows}

    print("\nC. POSITIONS WITH A NON-WT RESIDUE ABOVE WT")
    print(f"ESM2 {sum(esm_flags.values())}/{len(esm_flags)}")
    print(f"AbLang2 {sum(ablang_flags.values())}/{len(ablang_flags)}")
    print(
        "both "
        f"{sum(esm_flags[p] and ablang_flags[p] for p in DOMAIN_POSITIONS)}"
        f"/{len(DOMAIN_POSITIONS)}"
    )
    print("class n_parent n_ablang esm2_yes ablang2_yes both_yes")
    for class_name in ("RED", "GREEN", "UNLISTED"):
        parent_positions = [
            p
            for p in range(1, len(PARENT) + 1)
            if class_by_position.get(p, "UNLISTED") == class_name
        ]
        ablang_positions = [p for p in parent_positions if p in DOMAIN_POSITIONS]
        print(
            f"{class_name:8s} {len(parent_positions):8d} {len(ablang_positions):8d} "
            f"{sum(esm_flags[p] for p in parent_positions):8d} "
            f"{sum(ablang_flags[p] for p in ablang_positions):11d} "
            f"{sum(esm_flags[p] and ablang_flags[p] for p in ablang_positions):8d}"
        )
    print("common-position combination cross-tab")
    print("class neither esm2_only ablang2_only both total")
    for class_name in ("RED", "GREEN", "UNLISTED"):
        positions = [
            p
            for p in DOMAIN_POSITIONS
            if class_by_position.get(p, "UNLISTED") == class_name
        ]
        combinations = {"neither": 0, "esm2_only": 0, "ablang2_only": 0, "both": 0}
        for position in positions:
            esm_flag = esm_flags[position]
            ablang_flag = ablang_flags[position]
            if esm_flag and ablang_flag:
                combinations["both"] += 1
            elif esm_flag:
                combinations["esm2_only"] += 1
            elif ablang_flag:
                combinations["ablang2_only"] += 1
            else:
                combinations["neither"] += 1
        print(
            f"{class_name:8s} {combinations['neither']:7d} "
            f"{combinations['esm2_only']:9d} {combinations['ablang2_only']:12d} "
            f"{combinations['both']:4d} {len(positions):5d}"
        )

    print("\nD. SPEARMAN CORRELATION OF MUTATION SCORES")
    print("scope rows rho")
    for scope_name, positions in (
        ("overall_nonWT", DOMAIN_POSITIONS),
        ("CDR_nonWT", CDR_POSITIONS),
    ):
        esm_scores, ablang_scores = score_pairs(
            esm_logp, ablang_logp, positions, include_wt=False
        )
        print(f"{scope_name:18s} {len(esm_scores):4d} {spearman(esm_scores, ablang_scores):.6f}")
    for scope_name, positions in (
        ("overall_withWT", DOMAIN_POSITIONS),
        ("CDR_withWT", CDR_POSITIONS),
    ):
        esm_scores, ablang_scores = score_pairs(
            esm_logp, ablang_logp, positions, include_wt=True
        )
        print(f"{scope_name:18s} {len(esm_scores):4d} {spearman(esm_scores, ablang_scores):.6f}")

    print("\nE. RAW P(WT) DISTRIBUTION")
    print("model n min p05 p25 median mean p75 p95 max")
    esm_pwt = np.asarray(
        [np.exp(esm_logp[p - 1, AA_INDEX[PARENT[p - 1]]]) for p in range(1, len(PARENT) + 1)]
    )
    ablang_pwt = np.asarray(
        [np.exp(ablang_logp[p - 1, AA_INDEX[PARENT[p - 1]]]) for p in sorted(DOMAIN_POSITIONS)]
    )
    for model_name, values in (("ESM2", esm_pwt), ("AbLang2", ablang_pwt)):
        stats = distribution(values)
        print(
            f"{model_name:7s} {stats['n']:3d} {stats['min']:.8f} {stats['p05']:.8f} "
            f"{stats['p25']:.8f} {stats['median']:.8f} {stats['mean']:.8f} "
            f"{stats['p75']:.8f} {stats['p95']:.8f} {stats['max']:.8f}"
        )

    print(f"\nF. WROTE {entropy_path.relative_to(REPO)} ({len(PARENT)} rows)")
    print("AbLang2 entropy is blank for linker positions 121-135.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
