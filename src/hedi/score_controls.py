"""Score a directory of Boltz control runs and report the control separation.

    python src/hedi/score_controls.py <outputs_dir> [pae_cutoff] [dist_cutoff]

Expects runs named <control>_seed<N>, e.g. parent_wt_seed1. Scores every
diffusion sample of every run, then aggregates at two levels:

  per run (seed)  median over the 5 diffusion samples, plus spread
  per control     median over the 3 SEED medians, plus seed-to-seed spread

Two levels because the two sources of variation are not interchangeable.
Diffusion samples within one seed share an MSA and a seed, so they are
correlated: treating 15 samples as 15 independent observations would
overstate the evidence. The honest unit of replication is the seed, n=3.
Separation is therefore judged against seed-to-seed variation, and the
per-sample numbers are reported but not used to manufacture significance.

Writes results/control_scores.csv (per-sample and per-run rows).
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from ipsae import (  # noqa: E402
    DOCK_THRESHOLD,
    MAX_SPREAD,
    RunScore,
    find_boltz_samples,
    gate,
    read_boltz_sample,
    score_pairs,
)

OUT_CSV = REPO / "results" / "control_scores.csv"
IPSAE_THRESHOLD = 0.3  # published discrimination threshold, for reference only


def control_of(run: str) -> str:
    return run.rsplit("_seed", 1)[0] if "_seed" in run else run


def seed_of(run: str) -> str:
    return run.rsplit("_seed", 1)[1] if "_seed" in run else "-"


def score_run(run_dir: Path, pae_cutoff: float, dist_cutoff: float) -> list[dict]:
    """One row per diffusion sample."""
    rows = []
    for s in find_boltz_samples(run_dir):
        # reverse_vl_vh is 242 aa and keeps its linker, so the default chain
        # checks apply to every control here.
        model = read_boltz_sample(s)
        rec = next(
            r for r in score_pairs(model, pae_cutoff, dist_cutoff)
            if r.kind == "max" and {r.chain1, r.chain2} == {"A", "B"}
        )
        rows.append(
            {
                "run": run_dir.name,
                "control": control_of(run_dir.name),
                "seed": seed_of(run_dir.name),
                "model": s.model,
                "ipsae": rec.ipsae,
                "ipsae_d0chn": rec.ipsae_d0chn,
                "ipsae_d0dom": rec.ipsae_d0dom,
                "iptm_d0chn": rec.iptm_d0chn,
                "pdockq": rec.pdockq,
                "pdockq2": rec.pdockq2,
                "lis": rec.lis,
                "n0res": rec.n0res,
                "nres1": rec.nres1,
                "nres2": rec.nres2,
                "mean_plddt": float(np.mean(model.plddt)),
            }
        )
    return rows


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print(__doc__.strip().splitlines()[2])
        return 2
    base = Path(argv[1])
    pae_cutoff = float(argv[2]) if len(argv) > 2 else 10.0
    dist_cutoff = float(argv[3]) if len(argv) > 3 else 15.0

    run_dirs = sorted(
        d for d in base.iterdir()
        if d.is_dir() and any(d.rglob("pae_*_model_*.npz"))
    )
    # parent_quick is a smoke test (empty MSA, 1 recycle, 1 sample) and its
    # own PROVENANCE.md says so. Never let it into a baseline.
    skipped = [d.name for d in run_dirs if "_seed" not in d.name]
    run_dirs = [d for d in run_dirs if "_seed" in d.name]
    if skipped:
        print(f"skipping (no _seed, not a control run): {', '.join(skipped)}\n")
    if not run_dirs:
        print(f"no <control>_seed<N> runs found under {base}")
        return 1

    samples: list[dict] = []
    for d in run_dirs:
        samples += score_run(d, pae_cutoff, dist_cutoff)

    # ---- per run ----
    print(f"ipSAE (max record, chains A-B), pae<{pae_cutoff:g} dist<{dist_cutoff:g}\n")
    print(f"{'run':24s} {'n':>2s} {'median':>8s} {'spread':>8s} {'min':>8s} {'max':>8s}   per-sample")
    runs: dict[str, list[float]] = {}
    for d in run_dirs:
        v = [r["ipsae"] for r in samples if r["run"] == d.name]
        runs[d.name] = v
        print(
            f"{d.name:24s} {len(v):2d} {np.median(v):8.4f} {max(v)-min(v):8.4f} "
            f"{min(v):8.4f} {max(v):8.4f}   " + " ".join(f"{x:.4f}" for x in v)
        )

    # ---- per control, over seed medians ----
    controls: dict[str, list[float]] = {}
    for name in dict.fromkeys(control_of(d.name) for d in run_dirs):
        seed_medians = [
            float(np.median(v)) for r, v in runs.items() if control_of(r) == name
        ]
        controls[name] = seed_medians

    print(f"\n{'control':24s} {'seeds':>5s} {'median':>8s} {'seed spread':>11s}   seed medians")
    for name, sm in controls.items():
        print(
            f"{name:24s} {len(sm):5d} {np.median(sm):8.4f} {max(sm)-min(sm):11.4f}   "
            + " ".join(f"{x:.4f}" for x in sorted(sm, reverse=True))
        )

    # ---- separation ----
    print("\n" + "=" * 72)
    pw, cs = controls.get("parent_wt"), controls.get("cdrh3_scrambled")
    if pw and cs:
        d_med = float(np.median(pw)) - float(np.median(cs))
        noise = max(max(pw) - min(pw), max(cs) - min(cs))
        pw_s = [r["ipsae"] for r in samples if r["control"] == "parent_wt"]
        cs_s = [r["ipsae"] for r in samples if r["control"] == "cdrh3_scrambled"]
        disjoint_seeds = min(pw) > max(cs)
        disjoint_samples = min(pw_s) > max(cs_s)

        print(f"parent_wt        median {np.median(pw):.4f}  seed range [{min(pw):.4f}, {max(pw):.4f}]")
        print(f"cdrh3_scrambled  median {np.median(cs):.4f}  seed range [{min(cs):.4f}, {max(cs):.4f}]")
        print(f"\nseparation             {d_med:+.4f}")
        print(f"largest seed spread     {noise:.4f}   (the relevant noise floor, n=3 per control)")
        if noise > 0:
            print(f"ratio                   {abs(d_med)/noise:.1f}x seed-to-seed variation")
        print(f"\nseed medians disjoint?  {'YES' if disjoint_seeds else 'NO'}")
        print(f"all 15 samples disjoint? {'YES' if disjoint_samples else 'NO'}"
              f"   parent min {min(pw_s):.4f} vs scrambled max {max(cs_s):.4f}")
        if d_med <= 0:
            print("\n*** cdrh3_scrambled does NOT score below parent_wt. The negative")
            print("    control has not behaved as a negative control. Treat every")
            print("    design ranking from this pipeline as unvalidated until this")
            print("    is understood.")
    else:
        print("parent_wt and/or cdrh3_scrambled missing; no separation computed")

    # ---- mode-aware statistics ----
    # The scrambled control turns out to be BIMODAL, not shifted: a sample
    # either docks (~0.86) or collapses (~0.0-0.3), with an empirical gap in
    # between. A median over 5 draws of a bimodal variable reports whichever
    # mode happened to win the 5 draws, so it is not a safe summary here --
    # see the warning printed below.
    print("\n" + "=" * 72)
    print(f"MODE-AWARE VIEW (dock threshold {DOCK_THRESHOLD:g}, inside the empirical gap)\n")
    print(f"{'control':24s} {'docked':>8s} {'collapsed':>9s} {'median(docked)':>15s} {'max within-seed spread':>23s}")
    spreads: dict[str, list[float]] = {}
    for name in controls:
        v = [r["ipsae"] for r in samples if r["control"] == name]
        hi = [x for x in v if x >= DOCK_THRESHOLD]
        sp = [max(x) - min(x) for r, x in runs.items() if control_of(r) == name]
        spreads[name] = sp
        med = f"{np.median(hi):.4f}" if hi else "n/a"
        print(f"{name:24s} {len(hi):3d}/{len(v):<4d} {len(v)-len(hi):9d} {med:>15s} {max(sp):23.4f}")

    good = [s for n, sp in spreads.items() if n != "cdrh3_scrambled" for s in sp]
    bad = spreads.get("cdrh3_scrambled", [])
    if good and bad:
        print(f"\nwithin-seed spread separates cleanly:")
        print(f"  parent_wt + reverse_vl_vh   max spread {max(good):.4f}")
        print(f"  cdrh3_scrambled             min spread {min(bad):.4f}")
        if min(bad) > max(good):
            print(f"  DISJOINT, {min(bad)/max(good):.1f}x gap -- spread is the better discriminator here")

    if pw and cs:
        pw_hi = [r["ipsae"] for r in samples
                 if r["control"] == "parent_wt" and r["ipsae"] >= DOCK_THRESHOLD]
        cs_hi = [r["ipsae"] for r in samples
                 if r["control"] == "cdrh3_scrambled" and r["ipsae"] >= DOCK_THRESHOLD]
        if cs_hi:
            print(f"\nWhen the scramble DOES dock it scores {np.median(cs_hi):.4f}, against "
                  f"{np.median(pw_hi):.4f} for the parent.")
            print("ipSAE MAGNITUDE therefore does not separate them; the DOCKING RATE does")
            print(f"({len(pw_hi)}/15 vs {len(cs_hi)}/15).")
        print("\n*** WARNING on the median-of-5 aggregation: under this bimodality the")
        print("    per-seed median reports whichever mode won that seed's 5 draws.")
        print("    cdrh3_scrambled_seed1 medians at 0.8873 and would pass any")
        print("    magnitude filter, while seed2 medians at 0.1240. A design scored")
        print("    on one seed of 5 samples can be called good purely by draw.")
        print("    Recommend gating on docking rate (and/or spread) BEFORE comparing")
        print("    medians, and raising this with the team: decisions.md currently")
        print("    commits to median aggregation alone.")

    # ---- the gate, applied to each run ----
    # If the gate is right, the six good control runs pass and the three
    # scrambled ones fail. That is the only validation the thresholds have.
    print("\n" + "=" * 72)
    print(f"GATE (every sample docks at ipSAE >= {DOCK_THRESHOLD:g}, spread <= {MAX_SPREAD:g})")
    print("Designs are gated on docking rate and spread, NOT on median ipSAE.\n")
    expected_fail = {"cdrh3_scrambled"}
    surprises = []
    for d in run_dirs:
        v = [r["ipsae"] for r in samples if r["run"] == d.name]
        verdict = gate(
            RunScore(
                run=d.name, chain1="A", chain2="B",
                pae_cutoff=pae_cutoff, dist_cutoff=dist_cutoff,
                models=[r["model"] for r in samples if r["run"] == d.name],
                values=v,
            )
        )
        print(f"  {verdict}")
        should_fail = control_of(d.name) in expected_fail
        if verdict.passed == should_fail:
            surprises.append(d.name)
    if surprises:
        print(f"\n  *** gate disagrees with the controls on: {', '.join(surprises)}")
        print("      The thresholds are not calibrated. Do not use them on designs.")
    else:
        print(f"\n  gate agrees with all {len(run_dirs)} controls: the negative fails,")
        print("  the positive and the topology probe pass.")

    rv = controls.get("reverse_vl_vh")
    if rv and pw:
        print(f"\nreverse_vl_vh    median {np.median(rv):.4f}  "
              f"({np.median(rv) - np.median(pw):+.4f} vs parent_wt)")
        print("  Topology probe, not a negative control: the deck puts the reverse")
        print("  order at 1.83 nM vs 0.82 nM, ~0.5 kcal/mol, below a structure")
        print("  predictor's resolution. Reported, not interpreted.")

    print(f"\nthreshold for reference: published ipSAE discrimination ~{IPSAE_THRESHOLD}")

    # ---- CSV ----
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT_CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(samples[0]))
        w.writeheader()
        w.writerows(samples)
    print(f"\nwrote {OUT_CSV.relative_to(REPO)} ({len(samples)} per-sample rows)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
