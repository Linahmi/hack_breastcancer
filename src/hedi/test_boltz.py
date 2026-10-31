"""Tests for the Boltz 2.2.1 reader.

Array keys are no longer stubbed -- they were confirmed against a real 2.2.1
tree (outputs/parent_quick on the aida branch): 'pae' in pae_*.npz and
'plddt' in plddt_*.npz, both indexed by residue for a protein-only complex.

Covered here:
  * the real layout, which has a boltz_results_<run> level the upstream docs
    omit, AND the docs' layout without it -- the recursive glob handles both
  * --diffusion_samples N: every sample found and scored, not just model_0
  * median aggregation, and that it resists a single collapsed sample
  * a run with no PAE at all (--write_full_pae not passed) fails loudly
  * a renamed array key fails naming the keys that ARE present
  * token-vs-residue PAE shape mismatch
  * chain identity from the structure, never from confidence_*.json

Also re-scores a real Boltz run, which is REQUIRED, not optional:

    IPSAE_BOLTZ_RUN=/path/to/outputs/parent_wt_seed1 python -I src/hedi/test_boltz.py

The synthetic trees are written by this same file, so alone they only prove
the reader is self-consistent. They cannot catch a wrong array key or a real
layout change. Without IPSAE_BOLTZ_RUN this suite FAILS.
"""

from __future__ import annotations

import os
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import ipsae  # noqa: E402
from ipsae import (  # noqa: E402
    Model,
    assert_chain_layout,
    chain_sequence,
    find_boltz_samples,
    read_boltz_sample,
    read_cif_atoms,
    score_boltz_run,
)
from scripts.parent import PARENT  # noqa: E402

AA1TO3 = {v: k for k, v in ipsae.AA3TO1.items()}

# Stand-in for HER2 domain IV (the real trim is P04626 511-652, 142 aa).
# Its identity is never asserted; only that chain B exists and is non-empty.
TARGET_STUB = (
    "CHQLCARGHCWGPGPTQCVNCSQFLRGQECVEECRVLQGLPREYVNARHCLPCHPECQPQ"
    "NGSVTCFGPEADQCVACAHYKDPPFCVARCPSGVKPDLSYMPIWKFPDEEGACQPCPINC"
    "THSCVDLDDKGCPAEQRASPLT"
)


def _cif(chains: dict[str, str]) -> str:
    """Minimal valid mmCIF. Column order is deliberately not the reference's,
    to prove the parser reads fields by declared name rather than position."""
    head = [
        "data_test", "#", "loop_",
        "_atom_site.group_PDB", "_atom_site.id", "_atom_site.label_atom_id",
        "_atom_site.label_comp_id", "_atom_site.label_asym_id",
        "_atom_site.auth_asym_id", "_atom_site.label_seq_id",
        "_atom_site.auth_seq_id", "_atom_site.Cartn_x", "_atom_site.Cartn_y",
        "_atom_site.Cartn_z", "_atom_site.B_iso_or_equiv",
    ]
    rows, serial = [], 1
    for ci, (chain, seq) in enumerate(chains.items()):
        for i, aa in enumerate(seq):
            rn = AA1TO3.get(aa, "ALA")
            x, y, z = 3.8 * i, 9.0 * ci, 0.0
            rows.append(
                f"ATOM {serial} CA {rn} {chain} {chain} {i+1} {i+1} {x:.3f} {y:.3f} {z:.3f} 90.00"
            )
            serial += 1
            if rn != "GLY":
                rows.append(
                    f"ATOM {serial} CB {rn} {chain} {chain} {i+1} {i+1} "
                    f"{x:.3f} {y+1.4:.3f} {z+0.5:.3f} 90.00"
                )
                serial += 1
    return "\n".join(head + rows) + "\n#\n"


def _tree(
    root: Path,
    stem: str,
    chains: dict[str, str],
    *,
    models: list[int] | None = None,
    boltz_results_level: bool = True,
    inter_pae: list[float] | None = None,
    pae_key: str = "pae",
    pae_size: int | None = None,
    write_pae: bool = True,
) -> Path:
    """Write a Boltz-shaped output tree.

    boltz_results_level mirrors the REAL 2.2.1 layout (an extra
    boltz_results_<stem> directory); False gives the layout the upstream docs
    describe. Both must resolve.
    """
    models = models or [0]
    base = root / f"boltz_results_{stem}" if boltz_results_level else root
    d = base / "predictions" / stem
    d.mkdir(parents=True, exist_ok=True)

    n = sum(len(s) for s in chains.values())
    sizes = [len(s) for s in chains.values()]
    for k, m in enumerate(models):
        (d / f"{stem}_model_{m}.cif").write_text(_cif(chains))
        if write_pae:
            size = pae_size or n
            pae = np.full((size, size), 1.0, dtype=np.float32)
            if inter_pae is not None and size == n:
                # constant inter-chain PAE, distinct per sample
                a = sizes[0]
                pae[:a, a:] = inter_pae[k]
                pae[a:, :a] = inter_pae[k]
            np.savez(d / f"pae_{stem}_model_{m}.npz", **{pae_key: pae})
        np.savez(
            d / f"plddt_{stem}_model_{m}.npz",
            plddt=np.full(n, 0.9, dtype=np.float32),
        )
        (d / f"confidence_{stem}_model_{m}.json").write_text(
            # Present on purpose: nothing may read chain identity out of this.
            '{"pair_chains_iptm": {"0": {"1": 0.77}, "1": {"0": 0.77}}}'
        )
    return d


def _expect(fn, exc, needle: str, label: str, problems: list[str]) -> None:
    try:
        fn()
    except exc as e:
        if needle.lower() not in str(e).lower():
            problems.append(f"{label}: {exc.__name__} message lacks {needle!r}: {e}")
        return
    except Exception as e:  # noqa: BLE001
        problems.append(f"{label}: raised {type(e).__name__} not {exc.__name__}: {e}")
        return
    problems.append(f"{label}: did not raise {exc.__name__}")


def main() -> int:
    problems: list[str] = []
    chains = {"A": PARENT, "B": TARGET_STUB}
    n = len(PARENT) + len(TARGET_STUB)

    with tempfile.TemporaryDirectory(prefix="boltz_test_") as tmp:
        root = Path(tmp)

        # ---- 1. the REAL layout, with the boltz_results_ level ----
        real = root / "run_real"
        d = _tree(real, "parent_quick", chains, boltz_results_level=True)
        s = find_boltz_samples(real)
        if len(s) != 1 or s[0].pae_npz != d / "pae_parent_quick_model_0.npz":
            problems.append(f"real layout resolved {[str(x.pae_npz) for x in s]}")
        if s and s[0].stem != "parent_quick":
            problems.append(f"stem was {s[0].stem!r}")

        # ---- 2. the docs' layout, without that level, must also resolve ----
        docs = root / "run_docs"
        _tree(docs, "design_a", chains, boltz_results_level=False)
        if len(find_boltz_samples(docs)) != 1:
            problems.append("layout without boltz_results_ level did not resolve")

        # ---- 3. --diffusion_samples 5 -> all five, in model order ----
        multi = root / "run_multi"
        _tree(multi, "d5", chains, models=[0, 1, 2, 3, 4])
        got = find_boltz_samples(multi)
        if [x.model for x in got] != [0, 1, 2, 3, 4]:
            problems.append(f"diffusion samples resolved as {[x.model for x in got]}")
        if len({x.structure for x in got}) != 5:
            problems.append("samples do not have distinct structure files")

        # ---- 4. no PAE anywhere: the --write_full_pae failure mode ----
        nopae = root / "run_nopae"
        _tree(nopae, "nopae", chains, write_pae=False)
        _expect(
            lambda: find_boltz_samples(nopae),
            FileNotFoundError, "write_full_pae", "run without --write_full_pae", problems,
        )

        # ---- 5. several prediction directories under one out_dir ----
        two = root / "run_two"
        _tree(two, "design_001", chains)
        _tree(two, "design_002", chains)
        _expect(
            lambda: find_boltz_samples(two),
            ValueError, "expected one run", "two prediction dirs", problems,
        )

        # ---- 6. real keys, end to end, no monkeypatching ----
        m = read_boltz_sample(find_boltz_samples(real)[0])
        if m.numres != n:
            problems.append(f"numres {m.numres} != {n}")
        if chain_sequence(m, "A") != PARENT:
            problems.append("chain A does not round-trip the parent sequence")
        if not np.isclose(m.plddt.max(), 90.0):
            problems.append(f"plddt not rescaled 0-1 -> 0-100 (max {m.plddt.max()})")
        if m.iptm_af != -1.0:
            problems.append("iptm_af must stay -1.0, never read from confidence json")

        # ---- 7. a renamed array key must name what IS present ----
        renamed = root / "run_renamed"
        _tree(renamed, "renamed", chains, pae_key="predicted_aligned_error")
        _expect(
            lambda: read_boltz_sample(find_boltz_samples(renamed)[0]),
            KeyError, "predicted_aligned_error", "renamed PAE key", problems,
        )

        # ---- 8. PAE sized over tokens, not residues ----
        tok = root / "run_tokens"
        _tree(tok, "tok", chains, pae_size=n + 9)
        _expect(
            lambda: read_boltz_sample(find_boltz_samples(tok)[0]),
            ValueError, "token", "token-sized PAE", problems,
        )

        # ---- 9. chain checks, from the structure only ----
        sw = root / "run_swapped"
        _tree(sw, "swapped", {"A": TARGET_STUB, "B": PARENT})
        _expect(
            lambda: read_boltz_sample(find_boltz_samples(sw)[0]),
            ValueError, "swapped", "swapped chains", problems,
        )
        three = root / "run_three"
        _tree(three, "three", {"A": PARENT, "B": TARGET_STUB, "C": TARGET_STUB})
        _expect(
            lambda: read_boltz_sample(find_boltz_samples(three)[0]),
            ValueError, "exactly 2 chains", "three chains", problems,
        )
        short = root / "run_short"
        _tree(short, "short", {"A": PARENT[:-6], "B": TARGET_STUB})
        _expect(
            lambda: read_boltz_sample(find_boltz_samples(short)[0]),
            ValueError, "236 aa", "unexpected scFv length", problems,
        )
        read_boltz_sample(find_boltz_samples(short)[0], scfv_len=None)  # opt out

        # ---- 10. aggregation: median over samples, and its spread ----
        # Four good samples and one collapsed one. The median must ignore the
        # outlier; the mean would not. That is the whole reason for median.
        agg = root / "run_agg"
        _tree(agg, "agg", chains, models=[0, 1, 2, 3, 4],
              inter_pae=[2.0, 2.5, 3.0, 3.5, 25.0])
        run, every = score_boltz_run(agg, 10.0, 15.0)
        if run.n_samples != 5:
            problems.append(f"aggregated {run.n_samples} samples, expected 5")
        if run.models != [0, 1, 2, 3, 4]:
            problems.append(f"models {run.models}")
        if len(every) != 5:
            problems.append("per-sample records not all retained")
        if not np.isclose(run.median, float(np.median(run.values))):
            problems.append("median property disagrees with its own values")
        if not np.isclose(run.spread, max(run.values) - min(run.values)):
            problems.append("spread is not max-min")
        mean = float(np.mean(run.values))
        if run.median <= mean:
            problems.append(
                f"median {run.median:.4f} should sit ABOVE the mean {mean:.4f} "
                "when one sample collapses -- that is why we use it"
            )
        if run.minimum != min(run.values) or run.maximum != max(run.values):
            problems.append("min/max properties disagree with values")

    # ---- 11. regression against a real run (REQUIRED) ----
    # Not optional. The synthetic trees above are written by this file, so on
    # their own they only prove the reader is self-consistent -- they cannot
    # catch a wrong array key, a layout Boltz actually emits, or a shape
    # convention change. Only a real run can. Passing without it would be a
    # green result that validated nothing real.
    real_run = os.environ.get("IPSAE_BOLTZ_RUN")
    if not real_run:
        problems.append(
            "IPSAE_BOLTZ_RUN is not set, so the reader was never run against "
            "real Boltz output. The synthetic fixtures are written by this "
            "same file and cannot catch a wrong array key or a real layout "
            "change. Point it at a Boltz run directory, e.g. "
            "outputs/parent_wt_seed1."
        )
    else:
        run, _ = score_boltz_run(real_run, 10.0, 15.0)
        print(f"real run {run.run}: n={run.n_samples} median={run.median:.6f} "
              f"spread={run.spread:.6f} [{run.per_sample()}]")
        if not 0.0 <= run.median <= 1.0:
            problems.append(f"real run median {run.median} out of range")
        if run.n_samples < 1:
            problems.append("real run produced no samples")

    print(f"boltz reader: {len(problems)} problem(s)")
    for p in problems:
        print(f"  {p}")
    if problems:
        print("\nBOLTZ TESTS FAILED")
        return 1
    print("BOLTZ TESTS OK")
    return 0


if __name__ == "__main__":
    sys.exit(main())
