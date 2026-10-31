"""Batch Boltz-2 cofold runner: scFv candidates against HER2 domain IV.

Reads a CSV with columns `name,sequence`. For each row it builds the cofold
YAML (chain A = scFv, chain B = HER2 domain IV 511-652), runs Boltz-2 on
Modal, and writes the outputs under the agreed team path:

    outputs/<run>/boltz_results_<run>/predictions/<run>/pae_<run>_model_0.npz
    outputs/<run>/boltz_results_<run>/predictions/<run>/plddt_<run>_model_0.npz
    outputs/<run>/boltz_results_<run>/predictions/<run>/pde_<run>_model_0.npz
    outputs/<run>/boltz_results_<run>/predictions/<run>/confidence_<run>_model_0.json
    outputs/<run>/boltz_results_<run>/predictions/<run>/<run>_model_0.cif
    outputs/<run>/run_meta.json

PATH CONTRACT -- VERIFIED AGAINST A REAL RUN (boltz 2.2.1):
Boltz inserts a `boltz_results_<stem>` level that the upstream docs omit.
The team's circulating path `outputs/<name>/predictions/<name>/...` is
missing that level and matches nothing on disk. This runner preserves
Boltz's real layout rather than flattening it, so a committed example run
and live pipeline output have the same shape.

where <run> is the CSV `name`, plus `_seed<N>` when more than one seed is
requested. Nothing in <run> may contain the substring "pae".

Resumable: a row whose pae_*.npz already exists is skipped. Per-row
wall-clock is logged. A row that fails is reported loudly and recorded as
FAILED in outputs/_manifest.csv; it is never silently dropped, and the
process exits non-zero if any row failed.

NOTE ON SCOPE: this script runs a shortlist through cofold. Per the
organiser constraint, cofold is a FILTER, not an exploration tool --
screen with sequence models first and bring 50-100 candidates (or at
minimum the final 15) here.

Usage:
    # audit the plan without spending GPU time
    modal run run_cofold.py --csv shortlist.csv --dry-run

    # real run
    modal run run_cofold.py --csv shortlist.csv --seeds 1 --concurrency 8
"""

import csv
import json
import os
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from modal_boltz2 import app, boltz  # noqa: E402

HERE = Path(__file__).resolve().parent
DOMAIN_IV_FASTA = HERE / "her2_domainIV.fasta"

# --write_full_pae is MANDATORY: without it Boltz writes no pae_*.npz and the
# whole ipSAE half of the pipeline has nothing to read.
BASE_PARAMS = (
    "--use_msa_server "
    "--write_full_pae "
    "--recycling_steps 10 "
    "--diffusion_samples 5 "
    "--step_scale 1.0"
)

ALLOWED_AAS = set("ACDEFGHIKLMNPQRSTVWY")
MAX_LEN = 260
LINKER = "GGGGSGGGGSGGGGS"

# UniProt P04626 disulfides fully contained in domain IV 511-652.
HER2_SS_UNIPROT = [(511, 520), (515, 528), (531, 540), (544, 560),
                   (563, 576), (567, 584), (587, 596), (600, 623),
                   (626, 634), (630, 642)]
HER2_LO = 511


def load_her2():
    """Return the HER2 domain IV sequence and its local-index disulfide pairs."""
    if not DOMAIN_IV_FASTA.exists():
        raise SystemExit(f"FAIL: missing {DOMAIN_IV_FASTA}")
    seq = "".join(l.strip() for l in DOMAIN_IV_FASTA.read_text().splitlines()
                  if not l.startswith(">"))
    pairs = [(a - HER2_LO + 1, b - HER2_LO + 1) for a, b in HER2_SS_UNIPROT]
    for la, lb in pairs:
        if not (seq[la - 1] == "C" and seq[lb - 1] == "C"):
            raise SystemExit(
                f"FAIL: domain IV construct has no Cys at local {la}/{lb}; "
                f"her2_domainIV.fasta does not match the expected 511-652 window"
            )
    return seq, pairs


def check_run_name(run: str):
    """Reject names that would break the output contract."""
    if not run:
        raise ValueError("empty name")
    if "pae" in run.lower():
        raise ValueError(
            f"name {run!r} contains 'pae'; it would collide with pae_*.npz globbing"
        )
    bad = set(run) - set(
        "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789._-"
    )
    if bad:
        raise ValueError(f"name {run!r} has unsafe characters: {sorted(bad)}")


def scfv_disulfides(seq: str):
    """Pair scFv cysteines in order. Requires exactly 4 (VH C22-C96, VL C158-C223).

    Positions are found dynamically rather than hardcoded, because designs may
    use the 18 aa length headroom and shift them. Exactly 4 also enforces the
    'no new cysteines' rule.
    """
    cys = [i + 1 for i, a in enumerate(seq) if a == "C"]
    if len(cys) != 4:
        raise ValueError(
            f"{len(cys)} Cys at {cys}, expected exactly 4 "
            f"(conserved VH + VL pairs); refusing to guess the pairing"
        )
    return [(cys[0], cys[1]), (cys[2], cys[3])]


def validate_row(name: str, seq: str):
    """Validate one CSV row. Raises ValueError with a specific reason."""
    check_run_name(name)
    seq = seq.strip().upper()
    if not seq:
        raise ValueError("empty sequence")
    bad = sorted(set(seq) - ALLOWED_AAS)
    if bad:
        raise ValueError(f"non-standard residues {bad}")
    if len(seq) > MAX_LEN:
        raise ValueError(f"length {len(seq)} > {MAX_LEN} aa limit")
    n = seq.count(LINKER)
    if n != 1:
        raise ValueError(f"(GGGGS)3 linker appears {n}x, expected exactly 1")
    return seq, scfv_disulfides(seq)


def build_yaml(scfv: str, a_ss, her2: str, b_ss) -> str:
    """Emit the cofold YAML by hand: deterministic chain order, no pyyaml dep."""
    out = ["version: 1\n", "sequences:\n"]
    out += ["  - protein:\n", "      id: A\n", f"      sequence: {scfv}\n"]
    out += ["  - protein:\n", "      id: B\n", f"      sequence: {her2}\n"]
    out += ["constraints:\n"]
    for ch, pairs in (("A", a_ss), ("B", b_ss)):
        for x, y in pairs:
            out += ["  - bond:\n",
                    f"      atom1: [{ch}, {x}, SG]\n",
                    f"      atom2: [{ch}, {y}, SG]\n"]
    return "".join(out)


def expected_npz(out_root: Path, run: str) -> Path:
    """The pae npz path Boltz 2.2.1 actually produces (verified on a real run)."""
    return (out_root / run / f"boltz_results_{run}" / "predictions" / run
            / f"pae_{run}_model_0.npz")


def one(run: str, yaml_str: str, params: str, out_root: Path, meta: dict):
    """Run a single cofold and write its outputs. Returns (run, seconds, nfiles)."""
    t0 = time.monotonic()
    outputs = boltz.remote(yaml_str, params_str=params, name=run)
    elapsed = time.monotonic() - t0

    dest = out_root / run
    written = []
    for rel, content in outputs:
        target = dest / Path(rel)
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(content)
        written.append(str(Path(rel)))

    npz = expected_npz(out_root, run)
    if not npz.exists():
        raise RuntimeError(
            f"{run}: run returned {len(written)} files but no {npz.name}. "
            f"Is --write_full_pae in params? got files: {written[:12]}"
        )

    meta.update(seconds=round(elapsed, 1), n_files=len(written),
               files=sorted(written), pae_npz=str(npz.relative_to(out_root)))
    (dest / "run_meta.json").write_text(json.dumps(meta, indent=2) + "\n")
    return run, elapsed, len(written)


@app.local_entrypoint()
def main(
    csv_path: str,
    seeds: str = "1",
    out_dir: str = "outputs",
    concurrency: int = 8,
    params_str: str | None = None,
    dry_run: bool = False,
    force: bool = False,
):
    """Run every row of `csv_path` (columns name,sequence) through Boltz-2 cofold."""
    from concurrent.futures import ThreadPoolExecutor, as_completed

    her2, b_ss = load_her2()
    base = params_str or BASE_PARAMS
    if "--write_full_pae" not in base:
        raise SystemExit(
            "FAIL: --write_full_pae missing from params. Without it Boltz writes "
            "no pae_*.npz and ipSAE has nothing to read. Refusing to run."
        )

    seed_list = [int(s) for s in seeds.split(",") if s.strip()]
    if not seed_list:
        raise SystemExit("FAIL: no seeds given")
    suffix = len(seed_list) > 1

    src = Path(csv_path)
    if not src.exists():
        raise SystemExit(f"FAIL: no such CSV {src}")
    rows = list(csv.DictReader(src.open()))
    if not rows:
        raise SystemExit(f"FAIL: {src} has no data rows")
    cols = set(rows[0].keys())
    for need in ("name", "sequence"):
        if need not in cols:
            raise SystemExit(f"FAIL: {src} missing column {need!r}; has {sorted(cols)}")

    out_root = Path(out_dir)
    out_root.mkdir(parents=True, exist_ok=True)

    # ---- validate everything up front, so a bad row fails before GPU spend ----
    plan, rejected, seen = [], [], set()
    for i, r in enumerate(rows, start=2):
        nm = (r["name"] or "").strip()
        try:
            seq, a_ss = validate_row(nm, r["sequence"] or "")
        except ValueError as e:
            rejected.append((i, nm, str(e)))
            continue
        if nm in seen:
            rejected.append((i, nm, "duplicate name"))
            continue
        seen.add(nm)
        for s in seed_list:
            run = f"{nm}_seed{s}" if suffix else nm
            check_run_name(run)
            plan.append({
                "run": run, "name": nm, "seed": s, "seq": seq,
                "yaml": build_yaml(seq, a_ss, her2, b_ss),
                "params": f"{base} --seed {s}",
                "a_ss": a_ss,
            })

    if rejected:
        print(f"\n!! {len(rejected)} REJECTED ROW(S) -- fix these, they are not skipped silently:")
        for ln, nm, why in rejected:
            print(f"   line {ln}: name={nm!r}: {why}")

    done = [p for p in plan if expected_npz(out_root, p["run"]).exists() and not force]
    todo = [p for p in plan if p not in done]

    print(f"\nCSV {src}  rows={len(rows)}  valid={len(seen)}  seeds={seed_list}")
    print(f"chain A scFv from CSV | chain B HER2 domain IV 511-652 ({len(her2)} aa)")
    print(f"constraints per job: A:2 + B:{len(b_ss)} = {2+len(b_ss)} disulfides")
    print(f"params: {base} --seed <N>")
    print(f"out contract: {out_root}/<run>/predictions/<run>/pae_<run>_model_0.npz")
    print(f"jobs: {len(plan)}  already done (resume): {len(done)}  to run: {len(todo)}")
    for p in done:
        print(f"   SKIP {p['run']} (pae npz present)")

    if dry_run:
        print("\n--dry-run: nothing submitted. First job's YAML:\n")
        if plan:
            print(plan[0]["yaml"])
        if rejected:
            raise SystemExit(1)
        return

    if not todo:
        print("\nnothing to do.")
        return

    results, failures = [], []
    t0 = time.monotonic()
    with ThreadPoolExecutor(max_workers=max(1, concurrency)) as ex:
        futs = {}
        for p in todo:
            meta = {
                "run": p["run"], "name": p["name"], "seed": p["seed"],
                "scfv_len": len(p["seq"]), "scfv_sequence": p["seq"],
                "scfv_disulfides": p["a_ss"],
                "her2_construct": "P04626 511-652 (domain IV, 10 disulfides)",
                "her2_len": len(her2),
                "params": p["params"],
                "gpu": os.environ.get("GPU", "L40S"),
                "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            }
            futs[ex.submit(one, p["run"], p["yaml"], p["params"], out_root, meta)] = p
        for fut in as_completed(futs):
            p = futs[fut]
            try:
                run, secs, n = fut.result()
                print(f"  OK   {run}: {n} files in {secs:.1f}s")
                results.append((run, secs, n))
            except Exception as e:
                print(f"  FAIL {p['run']}: {type(e).__name__}: {e}")
                failures.append((p["run"], f"{type(e).__name__}: {e}"))
    wall = time.monotonic() - t0

    man = out_root / "_manifest.csv"
    new = not man.exists()
    with man.open("a", newline="") as fh:
        w = csv.writer(fh)
        if new:
            w.writerow(["run", "status", "seconds", "n_files", "detail"])
        for run, secs, n in results:
            w.writerow([run, "OK", f"{secs:.1f}", n, ""])
        for run, why in failures:
            w.writerow([run, "FAILED", "", "", why])

    print(f"\n{len(results)} ok, {len(failures)} failed, "
          f"{len(done)} resumed. wall-clock {wall:.1f}s")
    if results:
        times = sorted(s for _, s, _ in results)
        print(f"per-job seconds: min {times[0]:.1f} / "
              f"median {times[len(times)//2]:.1f} / max {times[-1]:.1f}")
    print(f"manifest: {man}")

    if failures or rejected:
        raise SystemExit(
            f"EXIT 1: {len(failures)} failed job(s), {len(rejected)} rejected row(s)"
        )
