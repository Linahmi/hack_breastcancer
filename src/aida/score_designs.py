"""Score Boltz-2 cofold runs with ipSAE and apply the stability filter.

    python src/aida/score_designs.py outputs/<run> [outputs/<run> ...]
    python src/aida/score_designs.py --all
    python src/aida/score_designs.py --all --csv results/scored.csv

Uses src/hedi/ipsae.py for the ipSAE computation. ipTM is read from Boltz's
own confidence_*.json, NOT from ipsae.py -- its iptm_af returns -1.0
deliberately, because pair_chains_iptm is keyed by undocumented integer
indices and guessing the chain mapping would be worse than reporting
nothing.

WHY WE FILTER ON STABILITY, NOT SCORE
-------------------------------------
Broken designs do not score low. They score UNSTABLE. Per diffusion
sample a broken design either docks (~0.87) or collapses (<0.3), with an
empirical gap between roughly 0.42 and 0.75 and nothing in between. The
median hides this completely: the scrambled-CDR-H3 control has median
ipSAE 0.8773, inside the parent's own seed range of 0.8615-0.9218. A
median-only pipeline passes a design whose CDR-H3 has been shuffled.

So the reject rule runs BEFORE the median, on the 5 individual samples,
and the median is reported as a secondary number only. Never collapse a
design to its median.

THRESHOLDS ARE PROVISIONAL
--------------------------
Fitted on n=3 controls (parent_wt, reverse_vl_vh, cdrh3_scrambled).
Margins are clean but the sample size is three. Hedi re-fits at 50
designs. Do not treat these as settled.
"""

import csv
import glob
import json
import statistics as st
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent.parent


def _resolve_ipsae():
    """Make src/hedi/ipsae.py importable from any branch.

    It lives on Hedi's branch, so a checkout of `aida` has only
    src/hedi/.gitkeep. Rather than copy his file into this branch, extract
    it from the ref into a temp dir. Tries, in order: a local
    src/hedi/ipsae.py, then origin/hedi, then any ref that has it.
    """
    local = REPO / "src" / "hedi"
    if (local / "ipsae.py").exists():
        sys.path.insert(0, str(local))
        return f"{local / 'ipsae.py'} (local checkout)"

    import subprocess
    import tempfile
    refs = ["origin/hedi", "hedi"] + [
        r for r in subprocess.run(
            ["git", "-C", str(REPO), "for-each-ref", "--format=%(refname:short)",
             "refs/remotes", "refs/heads"],
            capture_output=True, text=True).stdout.split()
    ]
    seen = set()
    for ref in refs:
        if ref in seen:
            continue
        seen.add(ref)
        got = subprocess.run(
            ["git", "-C", str(REPO), "show", f"{ref}:src/hedi/ipsae.py"],
            capture_output=True, text=True)
        if got.returncode == 0 and got.stdout:
            d = Path(tempfile.mkdtemp(prefix="ipsae_"))
            (d / "ipsae.py").write_text(got.stdout)
            sys.path.insert(0, str(d))
            return f"{ref}:src/hedi/ipsae.py (extracted to {d})"
    raise SystemExit(
        "FAIL: cannot find src/hedi/ipsae.py locally or on any ref.\n"
        "      Fetch Hedi's branch:  git fetch origin hedi"
    )


IPSAE_SOURCE = _resolve_ipsae()

PAE_CUTOFF = 10.0
DIST_CUTOFF = 15.0

# --- Gate, adopted from src/hedi/ipsae.py on branch hedi (commit 60e2178) ---
# DOCK_THRESHOLD = 0.5 : a sample has docked if its ipSAE is at least this.
# MAX_SPREAD     = 0.30: reject a run whose samples disagree by more.
# This replaces the earlier 0.20-spread / 0.65-floor rule, which is gone.
DOCK_FLOOR = 0.50
MAX_SAMPLE_SPREAD = 0.30
PROVISIONAL = (
    "PROVISIONAL: controls are extreme by construction. If real designs "
    "cluster near 0.5 ipSAE or 0.30 spread, refit against the actual "
    "distribution rather than applying these blindly."
)

FIELDS = [
    "run", "verdict", "reason", "n_samples", "docking_rate",
    "max_sample_spread", "min_sample", "ipsae_median", "iptm_median",
    "ipsae_samples", "iptm_samples", "collapse_diag", "thresholds",
]


def iptm_samples(run_dir: Path):
    """ipTM per diffusion sample, from Boltz's confidence json."""
    run = run_dir.name
    out = []
    pat = run_dir / f"boltz_results_{run}" / "predictions" / run
    for p in sorted(pat.glob("confidence_*_model_*.json"),
                    key=lambda q: int(q.stem.rsplit("_", 1)[1])):
        out.append(json.load(open(p))["iptm"])
    return out


def plddt_means(run_dir: Path):
    """Mean pLDDT per diffusion sample, rescaled to 0-100."""
    import numpy as np
    run = run_dir.name
    pat = run_dir / f"boltz_results_{run}" / "predictions" / run
    out = []
    for q in sorted(pat.glob("plddt_*_model_*.npz"),
                    key=lambda x: int(x.stem.rsplit("_", 1)[1])):
        a = np.load(q)["plddt"]
        out.append(float(a.mean()) * 100.0)
    return out


def ipsae_samples(run_dir: Path):
    """ipSAE per diffusion sample, plus the collapse diagnostic.

    Returns (values, diagnostics). The diagnostic per sample carries n0res
    and the HER2-side interface residue count, which both go to ~0 on a
    collapsed sample while mean pLDDT holds in the low 80s -- the structure
    stays confidently folded, it just stops touching the target. That is why
    pLDDT cannot be used to detect collapse and the docking rate can.
    """
    import ipsae
    _, every = ipsae.score_boltz_run(run_dir, PAE_CUTOFF, DIST_CUTOFF)
    vals, diag = [], []
    for recs in every:
        mx = [r for r in recs if r.kind == "max"]
        if not mx:
            raise RuntimeError(f"{run_dir.name}: no 'max' record in a sample")
        r = mx[0]
        vals.append(r.ipsae)
        # asym B->A carries the HER2-side interface residue count
        ba = [q for q in recs if q.kind == "asym" and q.chain1 == "B"]
        diag.append({
            "ipsae": r.ipsae,
            "n0res": getattr(r, "n0res", None),
            "nres_her2": getattr(ba[0], "nres1", None) if ba else None,
            "dist1": getattr(r, "dist1", None),
        })
    return vals, diag


def judge(samples):
    """Hedi's gate: EVERY sample must dock, and spread must be within bounds.

    The median does not gate. It only ranks designs that already passed.
    """
    n = len(samples)
    spread = max(samples) - min(samples)
    lo = min(samples)
    docked = sum(1 for v in samples if v >= DOCK_FLOOR)
    stats = {
        "n_samples": n,
        "docking_rate": f"{docked}/{n}",
        "max_sample_spread": round(spread, 4),
        "min_sample": round(lo, 4),
        "ipsae_median": round(st.median(samples), 4),
    }
    reasons = []
    if docked < n:
        failed = [f"{v:.4f}" for v in samples if v < DOCK_FLOOR]
        reasons.append(
            f"only {docked}/{n} samples docked at ipSAE >= {DOCK_FLOOR:g} "
            f"(collapsed: {', '.join(failed)})"
        )
    if spread > MAX_SAMPLE_SPREAD:
        reasons.append(f"sample-spread {spread:.4f} > {MAX_SAMPLE_SPREAD:g}")
    verdict = "REJECT" if reasons else "PASS"
    return verdict, "; ".join(reasons) or "all samples docked, spread in bounds", stats


def score(run_dir: Path):
    sam, diag = ipsae_samples(run_dir)
    pl = plddt_means(run_dir)
    for i, d in enumerate(diag):
        d["plddt"] = pl[i] if i < len(pl) else float("nan")
    it = iptm_samples(run_dir)
    verdict, reason, stats = judge(sam)
    row = {
        "run": run_dir.name, "verdict": verdict, "reason": reason,
        "iptm_median": round(st.median(it), 4) if it else "",
        "ipsae_samples": " ".join(f"{v:.4f}" for v in sam),
        "iptm_samples": " ".join(f"{v:.4f}" for v in it),
        "thresholds": PROVISIONAL,
        "collapse_diag": " | ".join(
            f"m{i}:ipsae={d['ipsae']:.3f},n0res={d['n0res']},"
            f"her2_iface={d['nres_her2']},plddt={d['plddt']:.1f}"
            for i, d in enumerate(diag)),
        **stats,
    }
    return row


def main(argv):
    args = [a for a in argv[1:] if not a.startswith("--")]
    out_csv = None
    if "--csv" in argv:
        out_csv = argv[argv.index("--csv") + 1]
        args = [a for a in args if a != out_csv]

    if "--all" in argv:
        runs = sorted(
            Path(p) for p in glob.glob(str(REPO / "outputs" / "*"))
            if Path(p).is_dir()
            and not Path(p).name.startswith(("_", "parent_quick"))
            and list(Path(p).glob(f"boltz_results_{Path(p).name}"))
        )
    else:
        runs = [Path(a) for a in args]
    if not runs:
        print(__doc__.strip().splitlines()[1])
        return 2

    print(f"ipSAE via {IPSAE_SOURCE}\n  pae_cutoff={PAE_CUTOFF} "
          f"dist_cutoff={DIST_CUTOFF}")
    print(f"GATE: every one of the samples must dock at ipSAE >= {DOCK_FLOOR:g}"
          f" AND spread <= {MAX_SAMPLE_SPREAD:g}")
    print("The median does NOT gate -- it only ranks designs that passed.")
    print(f"[{PROVISIONAL}]\n")

    rows = []
    for d in runs:
        try:
            rows.append(score(d))
        except Exception as e:
            print(f"  ERROR {d.name}: {type(e).__name__}: {e}")
            rows.append({"run": d.name, "verdict": "ERROR",
                         "reason": f"{type(e).__name__}: {e}",
                         "thresholds": PROVISIONAL})

    hdr = (f"{'run':24} {'verdict':8} {'dock':>6} {'spread':>8} "
           f"{'min':>7} {'median':>8} {'ipTM':>7}")
    print(hdr)
    print("-" * len(hdr))
    for r in rows:
        if r["verdict"] == "ERROR":
            print(f"{r['run']:24} {'ERROR':8} {r['reason'][:48]}")
            continue
        print(f"{r['run']:24} {r['verdict']:8} {r['docking_rate']:>6} "
              f"{r['max_sample_spread']:8.4f} {r['min_sample']:7.4f} "
              f"{r['ipsae_median']:8.4f} {r['iptm_median']:7.4f}")
    print("\nper-sample ipSAE (all 5, never collapsed):")
    for r in rows:
        if r["verdict"] != "ERROR":
            print(f"  {r['run']:24} {r['ipsae_samples']}")
            if r["verdict"] == "REJECT":
                print(f"  {'':24} ^ {r['reason']}")
    print("\ncollapse diagnostic (n0res and HER2 interface residues go to ~0")
    print("on a collapsed sample while pLDDT holds -- folded but not touching):")
    for r in rows:
        if r["verdict"] != "ERROR":
            print(f"  {r['run']}")
            print(f"    {r['collapse_diag']}")

    npass = sum(1 for r in rows if r["verdict"] == "PASS")
    print(f"\n{npass}/{len(rows)} PASS   "
          f"{sum(1 for r in rows if r['verdict']=='REJECT')} REJECT   "
          f"{sum(1 for r in rows if r['verdict']=='ERROR')} ERROR")

    if out_csv:
        Path(out_csv).parent.mkdir(parents=True, exist_ok=True)
        with open(out_csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=FIELDS)
            w.writeheader()
            for r in rows:
                w.writerow({k: r.get(k, "") for k in FIELDS})
        print(f"wrote {out_csv}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
