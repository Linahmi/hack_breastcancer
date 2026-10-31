"""Deterministic synthetic AF2-style fixtures for ipSAE parity testing.

We cannot parity-test against the Dunbrack reference without an input both
it and our implementation can read. Rather than wait on a real Boltz-2 run,
we synthesise PDB + PAE-JSON pairs we fully control, so we can deliberately
hit the edge cases that a real prediction might not:

  * rows whose interface is small enough that d0res clamps to 1.0, and rows
    large enough that the cube-root formula is actually evaluated
  * at least one row with NO inter-chain PAE below the cutoff (the 0.0 branch)
  * an asymmetric PAE matrix, so A->B and B->A genuinely differ
  * a glycine (no CB), which the reference scores using its CA instead
  * both AF2 JSON spellings: 'pae' and 'predicted_aligned_error'

Fixture sizes are chosen so n0res spans the L=26/27 clamp boundary in
calc_d0 / calc_d0_array. A 15-residue partner chain would clamp every row
to d0res=1.0 and the test would pass while checking nothing.

    python src/hedi/fixtures.py <outdir>
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

# 19 of the 20, deliberately excluding GLY so we can place glycines by hand.
AA3 = [
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "HIS", "ILE", "LEU",
    "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
]

# PAE values are capped at 31.75 A in real AF2 output.
PAE_MAX = 31.75


def _resnames(n: int, gly_at: tuple[int, ...], seed: int) -> list[str]:
    """Residue names for one chain, with GLY forced at the given 0-based slots."""
    rng = np.random.default_rng(seed)
    names = [AA3[i] for i in rng.integers(0, len(AA3), size=n)]
    for g in gly_at:
        if 0 <= g < n:
            names[g] = "GLY"
    return names


def _pae(
    sizes: dict[str, int],
    interface: dict[tuple[str, str], tuple[range, range]],
    dead_rows: dict[str, tuple[int, ...]],
    seed: int,
) -> np.ndarray:
    """Build an asymmetric PAE matrix with planted low-PAE interface patches.

    Intra-chain blocks get low PAE (a confidently folded chain). Inter-chain
    blocks are high by default; for each (c1, c2) in `interface` the given
    row/column window is driven low so those rows have many passing columns.
    Rows listed in `dead_rows` are forced high across every other chain, to
    exercise the "no valid pairs -> 0.0" branch.
    """
    rng = np.random.default_rng(seed)
    order = list(sizes)
    total = sum(sizes.values())

    # 0-based start offset of each chain in the matrix
    start, off = {}, 0
    for c in order:
        start[c] = off
        off += sizes[c]

    # Default: high PAE everywhere, then carve out structure.
    pae = rng.uniform(18.0, PAE_MAX, size=(total, total))

    # Intra-chain: low, rising with sequence separation.
    for c in order:
        s, n = start[c], sizes[c]
        idx = np.arange(n)
        sep = np.abs(idx[:, None] - idx[None, :]).astype(float)
        block = 0.5 + 0.22 * sep + rng.uniform(0.0, 0.9, size=(n, n))
        pae[s:s + n, s:s + n] = np.clip(block, 0.3, PAE_MAX)

    # Inter-chain interface patches: low PAE, decaying away from the patch centre.
    for (c1, c2), (rows, cols) in interface.items():
        for i in rows:
            for j in cols:
                r, k = start[c1] + i, start[c2] + j
                d = abs((i - rows.start) - (j - cols.start))
                pae[r, k] = np.clip(1.1 + 0.55 * d + rng.uniform(0.0, 0.7), 0.3, PAE_MAX)
                # Asymmetric on purpose: the transpose gets its own, larger value
                # so ipSAE(A->B) != ipSAE(B->A) and the max/asym split is tested.
                pae[k, r] = np.clip(pae[r, k] + rng.uniform(1.5, 4.5), 0.3, PAE_MAX)

    # Dead rows: nothing inter-chain passes a cutoff of 10 or 15.
    for c, rows in dead_rows.items():
        for i in rows:
            r = start[c] + i
            for c2 in order:
                if c2 == c:
                    continue
                s2, n2 = start[c2], sizes[c2]
                pae[r, s2:s2 + n2] = rng.uniform(26.0, PAE_MAX, size=n2)

    np.fill_diagonal(pae, 0.3)
    return np.round(pae, 2)


def _coords(sizes: dict[str, int], seed: int) -> dict[str, np.ndarray]:
    """CA coordinates: each chain an extended strand, chains stacked in y.

    Chain spacing is ~9 A so that CB-CB contacts land near the hard-coded
    8.0 A pDockQ cutoff -- close enough that some pairs pass and some do not,
    which is what we want the distance-dependent columns to be tested on.
    """
    rng = np.random.default_rng(seed)
    out = {}
    for k, (c, n) in enumerate(sizes.items()):
        xyz = np.zeros((n, 3))
        xyz[:, 0] = 3.8 * np.arange(n)
        xyz[:, 1] = 9.0 * k
        xyz[:, 2] = rng.uniform(-0.6, 0.6, size=n)
        out[c] = np.round(xyz, 3)
    return out


def _pdb(
    sizes: dict[str, int],
    names: dict[str, list[str]],
    ca: dict[str, np.ndarray],
    plddt: np.ndarray,
) -> str:
    """Write CA (+CB for non-GLY) atoms in fixed-column PDB format.

    The reference indexes its CA list and its CB list in parallel, so every
    residue must contribute exactly one CB entry. Glycine contributes none,
    and the reference falls back to its CA -- emitting a CB for GLY here
    would silently misalign the two lists.
    """
    lines, serial, flat = [], 1, 0
    for c in sizes:
        for i in range(sizes[c]):
            rn = names[c][i]
            x, y, z = ca[c][i]
            b = plddt[flat]
            lines.append(
                f"ATOM  {serial:>5d}  CA  {rn:>3s} {c}{i + 1:>4d}    "
                f"{x:>8.3f}{y:>8.3f}{z:>8.3f}  1.00{b:>6.2f}           C"
            )
            serial += 1
            if rn != "GLY":
                # CB nudged toward the neighbouring chain, so interface CB-CB
                # distances are shorter than the CA-CA spacing.
                sign = 1.0 if list(sizes).index(c) == 0 else -1.0
                lines.append(
                    f"ATOM  {serial:>5d}  CB  {rn:>3s} {c}{i + 1:>4d}    "
                    f"{x:>8.3f}{y + sign * 1.4:>8.3f}{z + 0.5:>8.3f}  1.00{b:>6.2f}           C"
                )
                serial += 1
            flat += 1
    lines.append("TER")
    lines.append("END")
    return "\n".join(lines) + "\n"


def build(name: str, outdir: Path, *, pae_key: str = "pae") -> dict:
    """Write <name>.pdb and <name>.json into outdir. Returns a manifest."""
    if name == "pair":
        # A=30, B=40: n0res reaches 40, so d0res spans the clamp boundary.
        sizes = {"A": 30, "B": 40}
        interface = {("A", "B"): (range(4, 26), range(6, 38))}
        dead = {"A": (0, 29), "B": (0,)}
        gly = {"A": (7, 18), "B": (3, 22, 39)}
        seed = 20260131
    elif name == "trimer":
        # Three chains: checks pair enumeration and that chain C never leaks
        # into the A/B interface mask.
        sizes = {"A": 25, "B": 30, "C": 20}
        interface = {
            ("A", "B"): (range(3, 20), range(5, 28)),
            ("B", "C"): (range(2, 18), range(1, 16)),
        }
        dead = {"A": (24,), "C": (19,)}
        gly = {"A": (5,), "B": (11, 29), "C": (0, 8)}
        seed = 777
    else:
        raise ValueError(f"unknown fixture {name!r}")

    total = sum(sizes.values())
    rng = np.random.default_rng(seed + 1)
    plddt = np.round(rng.uniform(32.0, 98.0, size=total), 2)
    names = {c: _resnames(n, gly.get(c, ()), seed + i) for i, (c, n) in enumerate(sizes.items())}
    pae = _pae(sizes, interface, dead, seed)
    ca = _coords(sizes, seed)

    outdir.mkdir(parents=True, exist_ok=True)
    pdb_path = outdir / f"{name}.pdb"
    json_path = outdir / f"{name}.json"
    pdb_path.write_text(_pdb(sizes, names, ca, plddt))
    json_path.write_text(
        json.dumps(
            {
                pae_key: pae.tolist(),
                "plddt": plddt.tolist(),
                # Reference reads these straight through as ipTM_af / pTM.
                "iptm": 0.6421,
                "ptm": 0.7788,
            }
        )
    )
    return {
        "name": name,
        "pdb": str(pdb_path),
        "json": str(json_path),
        "sizes": sizes,
        "numres": total,
        "pae_key": pae_key,
    }


FIXTURES = {"pair": {"pae_key": "pae"}, "trimer": {"pae_key": "predicted_aligned_error"}}


def build_all(outdir: Path) -> list[dict]:
    return [build(n, outdir, **kw) for n, kw in FIXTURES.items()]


if __name__ == "__main__":
    out = Path(sys.argv[1] if len(sys.argv) > 1 else "fixtures")
    for m in build_all(out):
        print(f"{m['name']:8s} {m['numres']:4d} res  chains={m['sizes']}  key={m['pae_key']}")
        print(f"         {m['pdb']}")
        print(f"         {m['json']}")
