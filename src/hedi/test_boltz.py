"""Tests for the Boltz 2.2.1 reader -- everything except the array keys.

`_boltz_arrays` is still stubbed (we have no real output yet), so these tests
monkeypatch it to prove that every step AROUND it already works and keeps
working once it is filled in:

  * the PAE file is resolved by glob, with exactly-one-match enforced
  * the Boltz-1 boltz_results_*/ layout fails with a message saying so
  * mmCIF parsing by declared column name
  * chain identity derived from the structure, not confidence_*.json
  * the --write_full_pae failure mode (no PAE matrix) fails loudly
  * a token-vs-residue PAE shape mismatch fails loudly

    python -I src/hedi/test_boltz.py
"""

from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import ipsae  # noqa: E402
from ipsae import (  # noqa: E402
    BoltzArrays,
    Model,
    assert_chain_layout,
    chain_sequence,
    find_boltz_outputs,
    read_boltz,
    read_cif_atoms,
)
from scripts.parent import PARENT  # noqa: E402

AA1TO3 = {v: k for k, v in ipsae.AA3TO1.items()}

# Stand-in for HER2 domain IV. Its identity is never asserted (the real trim
# is a teammate's job); only that chain B exists and is non-empty.
TARGET_STUB = ("TQVCTGTDMKLRLPASPETHLDMLRHLYQGCQVVQGNLELTYLPTNASLSFLQDIQEVQGYVLIAHNQVRQVPLQRLRIVRGTQLFEDNYALAVLDNGDPLNNTTPVTGASPG")


def _cif(chains: dict[str, str]) -> str:
    """Minimal but valid mmCIF with an _atom_site loop (CA, plus CB for non-GLY).

    Column order here is deliberately NOT the reference's, to prove the parser
    reads fields by declared name rather than by position.
    """
    head = [
        "data_test",
        "#",
        "loop_",
        "_atom_site.group_PDB",
        "_atom_site.id",
        "_atom_site.label_atom_id",
        "_atom_site.label_comp_id",
        "_atom_site.label_asym_id",
        "_atom_site.auth_asym_id",
        "_atom_site.label_seq_id",
        "_atom_site.auth_seq_id",
        "_atom_site.Cartn_x",
        "_atom_site.Cartn_y",
        "_atom_site.Cartn_z",
        "_atom_site.B_iso_or_equiv",
    ]
    rows, serial = [], 1
    for ci, (chain, seq) in enumerate(chains.items()):
        for i, aa in enumerate(seq):
            rn = AA1TO3.get(aa, "ALA")
            x, y, z = 3.8 * i, 9.0 * ci, 0.0
            rows.append(
                f"ATOM {serial} CA {rn} {chain} {chain} {i + 1} {i + 1} "
                f"{x:.3f} {y:.3f} {z:.3f} 90.00"
            )
            serial += 1
            if rn != "GLY":
                rows.append(
                    f"ATOM {serial} CB {rn} {chain} {chain} {i + 1} {i + 1} "
                    f"{x:.3f} {y + 1.4:.3f} {z + 0.5:.3f} 90.00"
                )
                serial += 1
    return "\n".join(head + rows) + "\n#\n"


def _tree(root: Path, stem: str, chains: dict[str, str], model: int = 0) -> Path:
    """Write a Boltz 2.2.1-shaped output tree: <root>/predictions/<stem>/."""
    d = root / "predictions" / stem
    d.mkdir(parents=True, exist_ok=True)
    (d / f"{stem}_model_{model}.cif").write_text(_cif(chains))
    np.savez(d / f"pae_{stem}_model_{model}.npz", placeholder=np.zeros((2, 2)))
    np.savez(d / f"plddt_{stem}_model_{model}.npz", placeholder=np.zeros(2))
    (d / f"confidence_{stem}_model_{model}.json").write_text(
        # Present on purpose: nothing may read chain identity out of this.
        '{"pair_chains_iptm": {"0": {"1": 0.77}, "1": {"0": 0.77}}}'
    )
    return d


def _expect(fn, exc, needle: str, label: str, problems: list[str]) -> None:
    try:
        fn()
    except exc as e:
        if needle.lower() not in str(e).lower():
            problems.append(f"{label}: raised {exc.__name__} but message lacks {needle!r}: {e}")
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

        # ---- 1. unique resolution ----
        d = _tree(root / "run1", "parent_her2_dIV", chains)
        paths = find_boltz_outputs(root / "run1")
        if paths.pae_npz != d / "pae_parent_her2_dIV_model_0.npz":
            problems.append(f"resolved wrong PAE file: {paths.pae_npz}")
        if paths.structure != d / "parent_her2_dIV_model_0.cif":
            problems.append(f"resolved wrong structure: {paths.structure}")
        for name in ("plddt_npz", "confidence_json"):
            if getattr(paths, name) is None:
                problems.append(f"sibling {name} not resolved")
        if paths.stem != "parent_her2_dIV":
            problems.append(f"stem was {paths.stem!r}")

        # ---- 2. zero matches, and the stale Boltz-1 layout ----
        stale = root / "run_stale"
        (stale / "boltz_results_x").mkdir(parents=True)
        np.savez(stale / "boltz_results_x" / "pae_x_model_0.npz", pae=np.zeros((2, 2)))
        _expect(
            lambda: find_boltz_outputs(stale),
            FileNotFoundError, "boltz 2.2.1", "stale boltz_results layout", problems,
        )

        # ---- 3. more than one match ----
        multi = root / "run_multi"
        _tree(multi, "design_001", chains)
        _tree(multi, "design_002", chains)
        _expect(
            lambda: find_boltz_outputs(multi),
            ValueError, "expected exactly 1", "two predictions", problems,
        )

        # ---- 4. mmCIF parsing ----
        atoms = read_cif_atoms(paths.structure)
        ca = [a for a in atoms if a["atom_name"] == "CA"]
        if len(ca) != n:
            problems.append(f"CIF gave {len(ca)} CA atoms, expected {n}")
        if {a["chain"] for a in atoms} != {"A", "B"}:
            problems.append(f"CIF chains were {sorted({a['chain'] for a in atoms})}")

        # ---- 5. the stub itself must refuse, pointing at the fix ----
        _expect(
            lambda: read_boltz(root / "run1"),
            NotImplementedError, "describe_boltz_npz", "stubbed _boltz_arrays", problems,
        )

        # ---- 6. guards around the stub, via monkeypatch ----
        real = ipsae._boltz_arrays
        try:
            # 6a. no PAE at all -- the --write_full_pae failure mode
            ipsae._boltz_arrays = lambda p: BoltzArrays(pae=None, plddt=None)
            _expect(
                lambda: read_boltz(root / "run1"),
                ValueError, "--write_full_pae", "missing PAE matrix", problems,
            )

            # 6b. empty array must not pass as "present"
            ipsae._boltz_arrays = lambda p: BoltzArrays(pae=np.array([]), plddt=None)
            _expect(
                lambda: read_boltz(root / "run1"),
                ValueError, "--write_full_pae", "empty PAE matrix", problems,
            )

            # 6c. token-length PAE (ligand/PTM tokens) must not be scored
            ipsae._boltz_arrays = lambda p: BoltzArrays(
                pae=np.zeros((n + 7, n + 7)), plddt=None
            )
            _expect(
                lambda: read_boltz(root / "run1"),
                ValueError, "token", "token-vs-residue PAE shape", problems,
            )

            # 6d. plddt length mismatch
            ipsae._boltz_arrays = lambda p: BoltzArrays(
                pae=np.full((n, n), 5.0), plddt=np.zeros(n - 3)
            )
            _expect(
                lambda: read_boltz(root / "run1"),
                ValueError, "plddt", "plddt length mismatch", problems,
            )

            # From here on the trees have different residue counts, so the
            # patch must size the PAE from the structure -- otherwise the
            # shape guard fires first and masks the check under test.
            def _sized(with_plddt: bool = False):
                def f(p):
                    total = len([a for a in read_cif_atoms(p.structure) if a["atom_name"] == "CA"])
                    rng = np.random.default_rng(0)
                    return BoltzArrays(
                        pae=np.full((total, total), 5.0),
                        plddt=rng.uniform(0.3, 0.95, size=total) if with_plddt else None,
                    )
                return f

            # 6e. the happy path, with 0-1 plddt rescaled to 0-100
            ipsae._boltz_arrays = _sized(with_plddt=True)
            m = read_boltz(root / "run1")
            if m.numres != n:
                problems.append(f"happy path numres {m.numres} != {n}")
            if m.plddt.max() <= 1.0:
                problems.append("plddt was not rescaled from 0-1 to 0-100")
            if m.iptm_af != -1.0:
                problems.append("iptm_af should stay -1.0 (never read from confidence json)")
            if chain_sequence(m, "A") != PARENT:
                problems.append("chain A sequence does not round-trip the parent")

            # 6f. chains swapped -> must be caught, and say they look swapped
            sw = root / "run_swapped"
            _tree(sw, "swapped", {"A": TARGET_STUB, "B": PARENT})
            _expect(
                lambda: read_boltz(sw),
                ValueError, "swapped", "swapped chains", problems,
            )

            # 6g. a third chain -> refuse (untrimmed target, or two copies)
            three = root / "run_three"
            _tree(three, "three", {"A": PARENT, "B": TARGET_STUB, "C": TARGET_STUB})
            _expect(
                lambda: read_boltz(three),
                ValueError, "exactly 2 chains", "three chains", problems,
            )

            # 6h. a design of non-parent length is allowed only explicitly
            short = root / "run_short"
            _tree(short, "short", {"A": PARENT[:-6], "B": TARGET_STUB})
            _expect(
                lambda: read_boltz(short),
                ValueError, "236 aa", "unexpected scFv length", problems,
            )
            m2 = read_boltz(short, scfv_len=None)  # opt out
            if m2.numres != len(PARENT) - 6 + len(TARGET_STUB):
                problems.append("scfv_len=None did not accept the shorter design")
        finally:
            ipsae._boltz_arrays = real

        # ---- 7. assert_chain_layout rejects a chain with no linker ----
        m3 = Model(
            pae=np.zeros((n, n)),
            chains=np.array(["A"] * len(PARENT) + ["B"] * len(TARGET_STUB)),
            resnums=np.array(list(range(1, len(PARENT) + 1)) + list(range(1, len(TARGET_STUB) + 1))),
            resnames=[AA1TO3.get(c, "ALA") for c in PARENT.replace("G", "A") + TARGET_STUB],
            plddt=np.zeros(n),
            distances=np.zeros((n, n)),
        )
        _expect(
            lambda: assert_chain_layout(m3),
            ValueError, "no (GGGGS)x3 linker", "linker-less chain A", problems,
        )

    print(f"boltz reader: {len(problems)} problem(s)")
    for p in problems:
        print(f"  {p}")
    if problems:
        print("\nBOLTZ TESTS FAILED")
        return 1
    print("BOLTZ TESTS OK (array keys still stubbed, by design)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
