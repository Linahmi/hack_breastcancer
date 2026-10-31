"""Parity test: src/hedi/ipsae.py vs the Dunbrack lab reference script.

Runs the reference on synthetic fixtures and diffs every reported quantity
against ours -- both the per-residue arrays (the _byres file: n0res, d0res,
ipTM_pae, ipSAE_d0chn, ipSAE_d0dom, ipSAE for every residue of every ordered
chain pair) and the summary rows (asym and max).

The per-residue comparison is the one that matters: a bug in the row-wise d0
or in the mean-over-columns can easily leave the final max-over-rows scalar
looking plausible. Checking ~1300 intermediate values per fixture does not.

The reference is NOT vendored into this repo -- point IPSAE_REF at a local
copy (MIT licensed, https://github.com/DunbrackLab/IPSAE):

    IPSAE_REF=/path/to/ipsae.py python -I src/hedi/test_parity.py

Without it the test skips rather than silently passing.
"""

from __future__ import annotations

import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import fixtures  # noqa: E402
from ipsae import calc_d0, calc_d0_array, read_af2, score_pairs  # noqa: E402

# The reference prints these to 4 decimals in the byres file and 6 in the
# summary, so exact equality is not available; parity is to the printed width.
BYRES_TOL = 5e-5
SUMMARY_TOL = 5e-7
D0_TOL = 5e-4  # d0 columns print to 3 decimals

CUTOFFS = [(10.0, 15.0), (15.0, 15.0), (5.0, 10.0)]


class Failure(Exception):
    pass


def _run_reference(ref: Path, pae: Path, pdb: Path, pae_cut: float, dist_cut: float):
    """Run the reference and return (summary_rows, byres_rows) as dict lists."""
    proc = subprocess.run(
        [sys.executable, "-I", str(ref), str(pae), str(pdb), str(pae_cut), str(dist_cut)],
        capture_output=True,
        text=True,
        cwd=str(ref.parent.parent),
    )
    if proc.returncode != 0:
        raise Failure(f"reference exited {proc.returncode}: {proc.stderr[-800:]}")

    stem = pdb.with_suffix("")
    # The reference zero-pads a cutoff below 10, so 5 -> "05".
    tag = "_" + "_".join(f"{int(c):02d}" for c in (pae_cut, dist_cut))
    summary_path = Path(f"{stem}{tag}.txt")
    byres_path = Path(f"{stem}{tag}_byres.txt")
    if not summary_path.exists():
        raise Failure(f"reference wrote no summary at {summary_path}")

    def _table(path: Path) -> list[dict]:
        lines = [l for l in path.read_text().splitlines() if l.strip()]
        header = lines[0].split()
        rows = []
        for line in lines[1:]:
            parts = line.split()
            if len(parts) < len(header):
                continue
            rows.append(dict(zip(header, parts)))
        return rows

    return _table(summary_path), _table(byres_path)


def _check(label: str, mine, theirs, tol: float, problems: list[str]) -> None:
    if isinstance(theirs, str):
        if str(mine).strip() != theirs.strip():
            problems.append(f"{label}: ours={mine!r} ref={theirs!r}")
        return
    if not np.isclose(mine, theirs, atol=tol, rtol=0.0):
        problems.append(f"{label}: ours={mine!r} ref={theirs!r} (tol {tol})")


def check_d0_helpers(problems: list[str]) -> None:
    """The d0 discontinuity we deliberately preserve.

    calc_d0 takes the cube-root branch only for L > 27, so it returns exactly
    1.0 at L=27 and ~1.039 just above. calc_d0_array clamps L at 26 instead,
    so the two disagree on 26 < L <= 27. If someone "tidies" either function
    these assertions are what catches it.
    """
    if calc_d0(27) != 1.0:
        problems.append(f"calc_d0(27)={calc_d0(27)}, expected exactly 1.0")
    if not np.isclose(calc_d0(27.001), 1.0391, atol=1e-3):
        problems.append(f"calc_d0(27.001)={calc_d0(27.001)}, expected ~1.0391")
    if not np.isclose(calc_d0(28), 1.1160, atol=1e-3):
        problems.append(f"calc_d0(28)={calc_d0(28)}, expected ~1.1160")
    # array version: clamped below 26, so 0 and 26 agree
    arr = calc_d0_array([0, 10, 26, 27, 28, 40])
    if not np.isclose(arr[0], arr[2]):
        problems.append("calc_d0_array: L=0 and L=26 should both be clamped")
    if not np.isclose(arr[2], 1.0, atol=1e-6):
        problems.append(f"calc_d0_array(26)={arr[2]}, expected 1.0 (post-2026-01-03 fix)")
    if not np.isclose(arr[5], 1.8260, atol=1e-3):
        problems.append(f"calc_d0_array(40)={arr[5]}, expected ~1.8260")
    # the deliberate divergence between the two
    if np.isclose(calc_d0(27), float(calc_d0_array([27])[0]), atol=1e-6):
        problems.append(
            "calc_d0(27) and calc_d0_array(27) now agree -- the reference's "
            "clamp discontinuity was removed, parity with it is broken"
        )


def check_fixture(ref: Path, manifest: dict, pae_cut: float, dist_cut: float) -> list[str]:
    problems: list[str] = []
    pae_path, pdb_path = Path(manifest["json"]), Path(manifest["pdb"])
    tag = f"{manifest['name']} @ pae<{pae_cut:g} dist<{dist_cut:g}"

    summary, byres = _run_reference(ref, pae_path, pdb_path, pae_cut, dist_cut)
    model = read_af2(pae_path, pdb_path)
    records = score_pairs(model, pae_cut, dist_cut)

    # ---- per-residue arrays ----
    by_pair: dict[tuple[str, str], list[dict]] = {}
    for row in byres:
        by_pair.setdefault((row["AlignChn"], row["ScoredChain"]), []).append(row)

    n_values = 0
    for rec in records:
        if rec.kind != "asym":
            continue
        rows = by_pair.get((rec.chain1, rec.chain2))
        if not rows:
            problems.append(f"{tag}: reference has no byres rows for {rec.pair}")
            continue
        for row in rows:
            i = int(row["i"]) - 1  # reference prints 1-based
            where = f"{tag} {rec.pair} res{row['i']}"
            _check(f"{where} n0res", int(rec.byres["n0res"][i]), int(row["n0res"]), 0, problems)
            _check(f"{where} d0res", rec.byres["d0res"][i], float(row["d0res"]), D0_TOL, problems)
            _check(f"{where} n0chn", rec.n0chn, int(row["n0chn"]), 0, problems)
            _check(f"{where} n0dom", rec.n0dom, int(row["n0dom"]), 0, problems)
            _check(f"{where} d0chn", rec.d0chn, float(row["d0chn"]), D0_TOL, problems)
            _check(f"{where} d0dom", rec.d0dom, float(row["d0dom"]), D0_TOL, problems)
            _check(f"{where} plddt", model.plddt[i], float(row["AlignRespLDDT"]), 5e-3, problems)
            _check(f"{where} iptm", rec.byres["iptm_d0chn"][i], float(row["ipTM_pae"]), BYRES_TOL, problems)
            _check(f"{where} ipsae_d0chn", rec.byres["ipsae_d0chn"][i], float(row["ipSAE_d0chn"]), BYRES_TOL, problems)
            _check(f"{where} ipsae_d0dom", rec.byres["ipsae_d0dom"][i], float(row["ipSAE_d0dom"]), BYRES_TOL, problems)
            _check(f"{where} ipsae", rec.byres["ipsae"][i], float(row["ipSAE"]), BYRES_TOL, problems)
            n_values += 11

    # ---- summary rows ----
    want = {
        "ipSAE": ("ipsae", SUMMARY_TOL),
        "ipSAE_d0chn": ("ipsae_d0chn", SUMMARY_TOL),
        "ipSAE_d0dom": ("ipsae_d0dom", SUMMARY_TOL),
        "ipTM_d0chn": ("iptm_d0chn", SUMMARY_TOL),
        "pDockQ": ("pdockq", 5e-5),
        "pDockQ2": ("pdockq2", 5e-5),
        "LIS": ("lis", 5e-5),
        "n0res": ("n0res", 0),
        "n0chn": ("n0chn", 0),
        "n0dom": ("n0dom", 0),
        "d0res": ("d0res", 5e-3),
        "d0chn": ("d0chn", 5e-3),
        "d0dom": ("d0dom", 5e-3),
        "nres1": ("nres1", 0),
        "nres2": ("nres2", 0),
        "dist1": ("dist1", 0),
        "dist2": ("dist2", 0),
    }
    mine = {(r.chain1, r.chain2, r.kind): r for r in records}
    seen = set()
    for row in summary:
        key = (row["Chn1"], row["Chn2"], row["Type"])
        # The reference emits its max row under one ordering only.
        rec = mine.get(key) or mine.get((row["Chn2"], row["Chn1"], row["Type"]))
        if rec is None:
            problems.append(f"{tag}: we produced no record for {key}")
            continue
        seen.add(key)
        for col, (attr, tol) in want.items():
            _check(f"{tag} {key} {col}", getattr(rec, attr), float(row[col]), tol, problems)
            n_values += 1

    if not seen:
        problems.append(f"{tag}: no summary rows compared -- parity proved nothing")
    print(f"  {tag}: {n_values} values compared, {len(problems)} mismatch(es)")
    return problems


def main() -> int:
    ref_env = os.environ.get("IPSAE_REF")
    if not ref_env:
        print("SKIP: set IPSAE_REF to a local copy of the reference ipsae.py")
        print("      https://github.com/DunbrackLab/IPSAE")
        return 0
    ref = Path(ref_env)
    if not ref.is_file():
        print(f"FAIL: IPSAE_REF={ref} is not a file")
        return 1

    problems: list[str] = []
    print("d0 helpers:")
    check_d0_helpers(problems)
    print(f"  {len(problems)} mismatch(es)")

    with tempfile.TemporaryDirectory(prefix="ipsae_parity_") as tmp:
        manifests = fixtures.build_all(Path(tmp))
        print(f"fixtures in {tmp}")
        for m in manifests:
            for pae_cut, dist_cut in CUTOFFS:
                problems += check_fixture(ref, m, pae_cut, dist_cut)

    print()
    if problems:
        print(f"PARITY FAILED: {len(problems)} mismatch(es)")
        for p in problems[:40]:
            print(f"  {p}")
        if len(problems) > 40:
            print(f"  ... and {len(problems) - 40} more")
        return 1
    print("PARITY OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
