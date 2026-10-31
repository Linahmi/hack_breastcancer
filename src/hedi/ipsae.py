"""ipSAE scoring for scFv / HER2-domain-IV cofolds.

A reimplementation of the Dunbrack lab reference (DunbrackLab/IPSAE, ipsae.py
version 4, 2026-01-03; https://doi.org/10.1101/2025.02.10.637595), kept
numerically identical to it. `test_parity.py` diffs every reported quantity
against the reference script column-for-column; if you change anything in
`score_pairs` or the d0/ptm helpers, that test is what tells you whether you
broke parity.

Why reimplement rather than shell out to the reference:
  * we need the per-residue arrays as objects, not as a text file to re-parse
  * we score batches of designs in one process
  * the reference has no library entry point; it runs off sys.argv and writes
    files next to its input

What this is FOR: triage. ipSAE separates a broken interface from a plausible
one (the paper's discrimination threshold is ~0.3). It does NOT resolve small
differences between two good designs -- a 0.02 gap is noise, not affinity.
Use it to throw designs out, then pick a diverse shortlist on other grounds.

Note that ipSAE values are much lower than ipTM values for the same model;
do not carry an ipTM-style cutoff (0.8 and such) over to ipSAE.

Scoring is asymmetric: ipSAE(A->B) != ipSAE(B->A). For ranking designs use
the symmetric `max` record (`Record.kind == "max"`), which is what the
reference reports on its "max" row.

    python src/hedi/ipsae.py <pae.json> <model.pdb> <pae_cutoff> <dist_cutoff>
"""

from __future__ import annotations

import json
import math
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

# The 20 standard residues, as the reference spells them.
RESIDUE_SET = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
}

# Hard-coded in the reference; not a tunable.
PDOCKQ_CONTACT_CUTOFF = 8.0
LIS_CUTOFF = 12.0
PDOCKQ2_D0 = 10.0


# --------------------------------------------------------------------------
# The TM-score kernel. Both exactly as the reference defines them.
# --------------------------------------------------------------------------

def ptm_func(x, d0):
    """1 / (1 + (x/d0)^2) -- the per-pair TM-score transform of a PAE value."""
    return 1.0 / (1 + (x / d0) ** 2.0)


def calc_d0(L, pair_type: str = "protein") -> float:
    """Scalar d0. Used for d0chn and d0dom.

    Mirrors the reference exactly, including the discontinuity: this returns
    1.0 for L <= 27 and jumps to ~1.039 just above it, because the cube-root
    branch is only taken for L > 27. calc_d0_array does NOT have this
    discontinuity -- it clamps L at 26 instead. The two therefore disagree
    in 26 < L <= 27. That is the reference's behaviour and we keep it, so
    our numbers match; do not "fix" it here.
    """
    L = float(L)
    min_value = 2.0 if pair_type == "nucleic_acid" else 1.0
    if L > 27:
        d0 = 1.24 * (L - 15) ** (1.0 / 3.0) - 1.8
    else:
        d0 = 1.0
    return max(min_value, d0)


def calc_d0_array(L, pair_type: str = "protein") -> np.ndarray:
    """Vectorised d0 for the per-residue d0res. Clamps L at 26, not 27."""
    L = np.array(L, dtype=float)
    L = np.maximum(26, L)
    min_value = 2.0 if pair_type == "nucleic_acid" else 1.0
    return np.maximum(min_value, 1.24 * (L - 15) ** (1.0 / 3.0) - 1.8)


# --------------------------------------------------------------------------
# Model: the parsed structure + confidence arrays the scorer needs.
# --------------------------------------------------------------------------

@dataclass
class Model:
    """Everything scoring needs, decoupled from the file format it came from.

    chains[i] is the chain id of residue i; resnums[i] its author residue
    number; pae is (n, n) and asymmetric; distances is the (n, n) CB-CB
    matrix (CA used for glycine, which has no CB).
    """

    pae: np.ndarray
    chains: np.ndarray
    resnums: np.ndarray
    resnames: list[str]
    plddt: np.ndarray
    distances: np.ndarray
    iptm_af: float = -1.0
    ptm_af: float = -1.0
    name: str = ""

    @property
    def numres(self) -> int:
        return len(self.chains)

    @property
    def unique_chains(self) -> list[str]:
        return sorted(set(self.chains.tolist()))

    def label(self, i: int) -> str:
        """Reference-style residue label, e.g. 'ARG   A     17'."""
        return f"{self.resnames[i]:3}   {self.chains[i]:3} {self.resnums[i]:4}"

    def validate(self) -> None:
        n = self.numres
        for nm, arr in (("pae", self.pae), ("distances", self.distances)):
            if arr.shape != (n, n):
                raise ValueError(f"{nm} is {arr.shape}, expected ({n}, {n})")
        if len(self.plddt) != n:
            raise ValueError(f"plddt has {len(self.plddt)} entries, expected {n}")
        if len(self.unique_chains) < 2:
            raise ValueError(
                f"need >= 2 chains to score an interface, found {self.unique_chains}"
            )


@dataclass
class Record:
    """One scored chain pair. `kind` is 'asym' (ordered) or 'max' (symmetric)."""

    chain1: str
    chain2: str
    kind: str
    pae_cutoff: float
    dist_cutoff: float

    ipsae: float = 0.0           # the headline score (d0res)
    ipsae_d0chn: float = 0.0
    ipsae_d0dom: float = 0.0
    iptm_d0chn: float = 0.0
    iptm_af: float = -1.0
    pdockq: float = 0.0
    pdockq2: float = 0.0
    lis: float = 0.0

    n0res: int = 0
    n0chn: int = 0
    n0dom: int = 0
    d0res: float = 0.0
    d0chn: float = 0.0
    d0dom: float = 0.0
    nres1: int = 0
    nres2: int = 0
    dist1: int = 0
    dist2: int = 0

    asymres: str = "None"
    byres: dict[str, np.ndarray] = field(default_factory=dict, repr=False)

    @property
    def pair(self) -> str:
        return f"{self.chain1}-{self.chain2}"


def _pair_type(model: Model, c1: str, c2: str) -> str:
    """Protein-only. Nucleic-acid pairs change d0's floor to 2.0; we have no
    test coverage for that path, so refuse rather than return a wrong number.
    """
    for c in (c1, c2):
        bad = {
            model.resnames[i]
            for i in range(model.numres)
            if model.chains[i] == c and model.resnames[i] not in RESIDUE_SET
        }
        if bad:
            raise NotImplementedError(
                f"chain {c} has non-standard residues {sorted(bad)}; this "
                "implementation is protein-only (see _pair_type)"
            )
    return "protein"


# --------------------------------------------------------------------------
# Scoring
# --------------------------------------------------------------------------

def score_pairs(
    model: Model, pae_cutoff: float = 10.0, dist_cutoff: float = 15.0
) -> list[Record]:
    """Score every ordered chain pair, then add the symmetric max records.

    Returns asym records first (one per ordered pair), then one max record
    per unordered pair -- the same rows, in the same order, as the
    reference's summary file.
    """
    model.validate()
    chains, pae = model.chains, model.pae
    dist, plddt = model.distances, model.plddt
    numres = model.numres
    uniq = model.unique_chains

    asym: dict[tuple[str, str], Record] = {}

    for c1 in uniq:
        for c2 in uniq:
            if c1 == c2:
                continue
            ptype = _pair_type(model, c1, c2)
            rec = Record(c1, c2, "asym", pae_cutoff, dist_cutoff, iptm_af=model.iptm_af)

            is1, is2 = chains == c1, chains == c2

            # --- d0chn: both chain lengths, no PAE cutoff anywhere in it ---
            rec.n0chn = int(is1.sum() + is2.sum())
            rec.d0chn = calc_d0(rec.n0chn, ptype)
            ptm_d0chn = ptm_func(pae, rec.d0chn)

            # --- the interface mask. Strict '<', and direction-specific,
            # because PAE is asymmetric. ---
            valid = np.outer(is1, is2) & (pae < pae_cutoff)

            # --- d0dom: unique residues on each side with >=1 passing partner.
            # The reference counts author residue NUMBERS, not indices. ---
            rows_hit = valid.any(axis=1)
            cols_hit = valid.any(axis=0)
            res1 = {int(model.resnums[i]) for i in np.where(rows_hit)[0]}
            res2 = {int(model.resnums[j]) for j in np.where(cols_hit)[0]}
            rec.nres1, rec.nres2 = len(res1), len(res2)
            rec.n0dom = len(res1) + len(res2)
            rec.d0dom = calc_d0(rec.n0dom, ptype)
            ptm_d0dom = ptm_func(pae, rec.d0dom)

            # --- d0res: per-row, from that row's count of passing columns ---
            n0res_byres = valid.sum(axis=1)
            d0res_byres = calc_d0_array(n0res_byres, ptype)

            iptm_byres = np.zeros(numres)
            ipsae_chn_byres = np.zeros(numres)
            ipsae_dom_byres = np.zeros(numres)
            ipsae_res_byres = np.zeros(numres)

            for i in np.where(is1)[0]:
                vp = valid[i]
                # ipTM: every chain2 column, no cutoff -- this is the
                # dilution that ipSAE exists to avoid.
                iptm_byres[i] = ptm_d0chn[i, is2].mean() if is2.any() else 0.0
                if vp.any():
                    ipsae_chn_byres[i] = ptm_d0chn[i, vp].mean()
                    ipsae_dom_byres[i] = ptm_d0dom[i, vp].mean()
                    ipsae_res_byres[i] = ptm_func(pae[i], d0res_byres[i])[vp].mean()

            # --- distance-gated tallies. Reported only; they do not enter
            # ipSAE. Note '<' here but '<=' for the pDockQ contact cutoff. ---
            dvalid = np.outer(is1, is2) & (pae < pae_cutoff) & (dist < dist_cutoff)
            rec.dist1 = len({int(model.resnums[i]) for i in np.where(dvalid.any(axis=1))[0]})
            rec.dist2 = len({int(model.resnums[j]) for j in np.where(dvalid.any(axis=0))[0]})

            # --- aggregate: max over chain1 rows.
            # argmax runs over the FULL array, as in the reference. Non-chain1
            # rows are left at 0.0, so this is the max over chain1 -- except
            # when every chain1 row scores 0.0, where argmax returns index 0
            # and the reported n0res/d0res/asymres may belong to another
            # chain. Faithful to the reference; such a pair has no interface
            # and its score is 0.0 regardless.
            mi = int(np.argmax(ipsae_res_byres))
            rec.ipsae = float(ipsae_res_byres[mi])
            rec.n0res = int(n0res_byres[mi])
            rec.d0res = float(d0res_byres[mi])
            rec.asymres = model.label(mi)
            rec.ipsae_d0chn = float(ipsae_chn_byres[np.argmax(ipsae_chn_byres)])
            rec.ipsae_d0dom = float(ipsae_dom_byres[np.argmax(ipsae_dom_byres)])
            rec.iptm_d0chn = float(iptm_byres[np.argmax(iptm_byres)])

            # --- pDockQ / pDockQ2: 8 A CB contacts, '<=' ---
            contact_res: set[int] = set()
            npairs = 0
            ptm_sum = 0.0
            for i in np.where(is1)[0]:
                cvalid = is2 & (dist[i] <= PDOCKQ_CONTACT_CUTOFF)
                npairs += int(cvalid.sum())
                if cvalid.any():
                    contact_res.add(int(i))
                    contact_res.update(int(j) for j in np.where(cvalid)[0])
                    ptm_sum += ptm_func(pae[i][cvalid], PDOCKQ2_D0).sum()
            if npairs > 0:
                mean_plddt = plddt[list(contact_res)].mean()
                x = mean_plddt * math.log10(npairs)
                rec.pdockq = 0.724 / (1 + math.exp(-0.052 * (x - 152.611))) + 0.018
                x2 = mean_plddt * (ptm_sum / npairs)
                rec.pdockq2 = 1.31 / (1 + math.exp(-0.075 * (x2 - 84.733))) + 0.005

            # --- LIS: mean of (12-PAE)/12 over the whole cross block, PAE<12 ---
            sel = pae[np.outer(is1, is2)]
            low = sel[sel < LIS_CUTOFF]
            rec.lis = float(((LIS_CUTOFF - low) / LIS_CUTOFF).mean()) if low.size else 0.0

            rec.byres = {
                "n0res": n0res_byres,
                "d0res": d0res_byres,
                "iptm_d0chn": iptm_byres,
                "ipsae_d0chn": ipsae_chn_byres,
                "ipsae_d0dom": ipsae_dom_byres,
                "ipsae": ipsae_res_byres,
            }
            asym[(c1, c2)] = rec

    return list(asym.values()) + _max_records(asym, uniq)


def _max_records(
    asym: dict[tuple[str, str], Record], uniq: list[str]
) -> list[Record]:
    """Symmetric records: each metric's better direction, taken independently.

    Different metrics may win in different directions for the same pair --
    the reference resolves each one separately, and so do we. LIS is the
    exception: it is AVERAGED over the two directions, not maximised.
    """
    out = []
    for a in uniq:
        for b in uniq:
            if a >= b:
                continue
            fwd, rev = asym[(a, b)], asym[(b, a)]
            m = Record(a, b, "max", fwd.pae_cutoff, fwd.dist_cutoff, iptm_af=fwd.iptm_af)

            # ipSAE (d0res) carries n0res/d0res/asymres from the winner.
            win = fwd if fwd.ipsae >= rev.ipsae else rev
            m.ipsae, m.n0res, m.d0res = win.ipsae, win.n0res, win.d0res
            m.asymres = win.asymres

            # d0dom likewise carries its own winner's n0dom/d0dom.
            dwin = fwd if fwd.ipsae_d0dom >= rev.ipsae_d0dom else rev
            m.ipsae_d0dom, m.n0dom, m.d0dom = dwin.ipsae_d0dom, dwin.n0dom, dwin.d0dom

            m.ipsae_d0chn = max(fwd.ipsae_d0chn, rev.ipsae_d0chn)
            m.iptm_d0chn = max(fwd.iptm_d0chn, rev.iptm_d0chn)
            m.pdockq = max(fwd.pdockq, rev.pdockq)
            m.pdockq2 = max(fwd.pdockq2, rev.pdockq2)
            m.lis = (fwd.lis + rev.lis) / 2.0  # averaged, not maximised

            m.n0chn, m.d0chn = fwd.n0chn, fwd.d0chn  # direction-independent
            m.nres1, m.nres2 = fwd.nres1, fwd.nres2
            m.dist1, m.dist2 = fwd.dist1, fwd.dist2
            out.append(m)
    return out


# --------------------------------------------------------------------------
# AF2 reader: fixed-column PDB + PAE JSON.
# --------------------------------------------------------------------------

def _parse_pdb_atom(line: str) -> dict | None:
    """Fixed-column PDB ATOM/HETATM record. Columns, not whitespace splits."""
    resname = line[17:20].strip()
    if resname == "LIG":
        return None
    return {
        "atom_num": int(line[6:11]),
        "atom_name": line[12:16].strip(),
        "resname": resname,
        "chain": line[21].strip(),
        "resnum": int(line[22:26]),
        "xyz": (float(line[30:38]), float(line[38:46]), float(line[46:54])),
    }


def read_pdb(pdb_path: str | Path) -> dict:
    """Pull the CA list (scoring order) and the CB list (distances) from a PDB.

    These two lists are indexed in PARALLEL, so every residue must contribute
    exactly one CB entry -- glycine has no CB and falls back to its CA. A
    residue with a missing CB (truncated side chain, alternate location) would
    shift the CB list relative to the CA list and silently corrupt every
    distance-derived number. The reference does not check for this; we do.
    """
    ca, cb = [], []
    for line in Path(pdb_path).read_text().splitlines():
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue
        atom = _parse_pdb_atom(line)
        if atom is None:
            continue
        name, resname = atom["atom_name"], atom["resname"]
        if name == "CA" or "C1" in name:
            ca.append(atom)
        if name == "CB" or "C3" in name or (resname == "GLY" and name == "CA"):
            cb.append(atom)

    if len(ca) != len(cb):
        raise ValueError(
            f"{Path(pdb_path).name}: {len(ca)} CA atoms but {len(cb)} CB atoms. "
            "These are indexed in parallel, so distances would be misaligned. "
            "Check for residues with missing side chains."
        )
    for i, (a, b) in enumerate(zip(ca, cb)):
        if a["chain"] != b["chain"] or a["resnum"] != b["resnum"]:
            raise ValueError(
                f"{Path(pdb_path).name}: CA/CB lists diverge at index {i} "
                f"({a['chain']}{a['resnum']} vs {b['chain']}{b['resnum']})"
            )
    if not ca:
        raise ValueError(f"{Path(pdb_path).name}: no CA atoms found")
    return {"ca": ca, "cb": cb}


def read_af2(pae_path: str | Path, pdb_path: str | Path) -> Model:
    """AF2-style model: .pdb structure + .json confidences."""
    parsed = read_pdb(pdb_path)
    ca, cb = parsed["ca"], parsed["cb"]
    n = len(ca)

    data = json.loads(Path(pae_path).read_text())
    # Key precedence matches the reference: 'pae' wins over the long spelling.
    if "pae" in data:
        pae = np.array(data["pae"], dtype=float)
    elif "predicted_aligned_error" in data:
        pae = np.array(data["predicted_aligned_error"], dtype=float)
    else:
        raise ValueError(
            f"{Path(pae_path).name}: no 'pae' or 'predicted_aligned_error' key "
            f"(found {sorted(data)[:8]})"
        )

    plddt = (
        np.array(data["plddt"], dtype=float) if "plddt" in data else np.zeros(n)
    )
    coords = np.array([a["xyz"] for a in cb], dtype=float)
    distances = np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=-1)

    model = Model(
        pae=pae,
        chains=np.array([a["chain"] for a in ca]),
        resnums=np.array([a["resnum"] for a in ca]),
        resnames=[a["resname"] for a in ca],
        plddt=plddt,
        distances=distances,
        iptm_af=float(data.get("iptm", -1.0)),
        ptm_af=float(data.get("ptm", -1.0)),
        name=Path(pdb_path).stem,
    )
    model.validate()
    return model


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------

COLUMNS = (
    "chain1", "chain2", "kind", "ipsae", "ipsae_d0chn", "ipsae_d0dom",
    "iptm_af", "iptm_d0chn", "pdockq", "pdockq2", "lis",
    "n0res", "n0chn", "n0dom", "d0res", "d0chn", "d0dom",
    "nres1", "nres2", "dist1", "dist2",
)


def main(argv: list[str]) -> int:
    if len(argv) < 3:
        print(__doc__.strip().splitlines()[-1])
        return 2
    pae_path, pdb_path = argv[1], argv[2]
    pae_cutoff = float(argv[3]) if len(argv) > 3 else 10.0
    dist_cutoff = float(argv[4]) if len(argv) > 4 else 15.0

    model = read_af2(pae_path, pdb_path)
    records = score_pairs(model, pae_cutoff, dist_cutoff)

    print(f"# {model.name}  pae_cutoff={pae_cutoff}  dist_cutoff={dist_cutoff}")
    print("\t".join(COLUMNS))
    for r in records:
        row = []
        for c in COLUMNS:
            v = getattr(r, c)
            row.append(f"{v:.6f}" if isinstance(v, float) else str(v))
        print("\t".join(row))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
