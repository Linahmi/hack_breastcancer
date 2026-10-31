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
import re
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

            # Interface residue counts are per CHAIN, and each chain has two
            # direction-specific counts: chain A appears as nres1 of A->B and
            # as nres2 of B->A. The max row reports the larger of the two for
            # each chain -- NOT the forward direction's pair.
            #
            # Easy to get wrong and easy to miss: it only shows up when the
            # two directions disagree, which a symmetric test fixture will
            # never produce. The real parent_quick cofold is what caught it
            # (ref 242/124/97/75 vs 226/124/97/72 from taking fwd's pair).
            m.nres1 = max(fwd.nres1, rev.nres2)
            m.nres2 = max(fwd.nres2, rev.nres1)
            m.dist1 = max(fwd.dist1, rev.dist2)
            m.dist2 = max(fwd.dist2, rev.dist1)
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
# Boltz 2.2.1 reader.
#
# Layout, confirmed against a real 2.2.1 tree (outputs/parent_quick on the
# aida branch) -- the upstream docs are wrong about this. Boltz inserts a
# boltz_results_<run> level that the docs omit:
#
#   <out_dir>/boltz_results_<run>/predictions/<run>/pae_<run>_model_<N>.npz
#
# The path that circulated in the team, <out_dir>/predictions/<run>/, is
# missing that level and matches nothing on disk. We therefore glob
# RECURSIVELY for */predictions/*/pae_*, which matches with or without the
# boltz_results_ level and so survives either layout.
#
# --diffusion_samples N produces model_0 .. model_(N-1), each with its own
# pae/plddt/confidence/structure. One sample is not a result: we score every
# sample and aggregate (see score_boltz_run).
#
# Sibling files are located by globbing inside the prediction directory, not
# by string-substituting the PAE filename. The reference does the latter
# (pae -> plddt, pae -> confidence), which rewrites the FIRST occurrence of
# "pae" anywhere in the path -- including in a parent directory name -- and so
# silently resolves to a file that does not exist, or worse, the wrong one.
# --------------------------------------------------------------------------

BOLTZ_PAE_GLOB = "**/predictions/*/pae_*_model_*.npz"
BOLTZ_MODEL_RE = re.compile(r"_model_(\d+)\.npz$")
PARENT_SCFV_LEN = 242  # the 242 aa trastuzumab scFv parent

AA3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}


@dataclass
class BoltzOutputs:
    """Resolved paths for ONE diffusion sample. Optional siblings may be None."""

    out_dir: Path
    pred_dir: Path
    stem: str
    model: int
    pae_npz: Path
    structure: Path
    plddt_npz: Path | None = None
    confidence_json: Path | None = None


def find_boltz_samples(out_dir: str | Path) -> list[BoltzOutputs]:
    """Resolve every diffusion sample of one Boltz run, in model order.

    At least one sample is required, not exactly one: --diffusion_samples 5
    yields model_0..model_4 and all of them are scored. What must be unique
    is the PREDICTION DIRECTORY -- several directories under one out_dir
    means out_dir points above several runs, and silently scoring a mixture
    of them would be worse than failing.
    """
    out_dir = Path(out_dir)
    if not out_dir.is_dir():
        raise NotADirectoryError(f"out_dir {out_dir} is not a directory")

    matches = sorted(out_dir.glob(BOLTZ_PAE_GLOB))
    if not matches:
        present = sorted(p.name for p in out_dir.rglob("*") if p.is_file())[:25]
        raise FileNotFoundError(
            f"no PAE npz matching {BOLTZ_PAE_GLOB!r} under {out_dir}. "
            f"Files present: {present or '(none)'}. Expected "
            "<out_dir>/boltz_results_<run>/predictions/<run>/pae_<run>_model_<N>.npz "
            "(Boltz inserts the boltz_results_ level; the upstream docs omit it). "
            "A run without --write_full_pae produces no PAE at all."
        )

    by_dir: dict[Path, list[Path]] = {}
    for m in matches:
        by_dir.setdefault(m.parent, []).append(m)
    if len(by_dir) > 1:
        raise ValueError(
            f"found predictions in {len(by_dir)} directories under {out_dir}, "
            f"expected one run: {[str(d) for d in sorted(by_dir)]}. "
            "Point out_dir at a single run."
        )

    pred_dir, paes = next(iter(by_dir.items()))
    stem = pred_dir.name

    samples: list[BoltzOutputs] = []
    for pae in paes:
        m = BOLTZ_MODEL_RE.search(pae.name)
        if m is None:
            raise ValueError(f"cannot read a model index from {pae.name}")
        idx = int(m.group(1))

        struct = sorted(pred_dir.glob(f"*_model_{idx}.cif")) or sorted(
            pred_dir.glob(f"*_model_{idx}.pdb")
        )
        # Exclude the npz siblings that also end in _model_<idx>
        struct = [p for p in struct if p.suffix in (".cif", ".pdb")]
        if len(struct) != 1:
            raise ValueError(
                f"expected exactly 1 structure for model {idx} in {pred_dir}, "
                f"found {[p.name for p in struct]}"
            )

        def _optional(pat: str) -> Path | None:
            hits = sorted(pred_dir.glob(pat))
            return hits[0] if len(hits) == 1 else None

        samples.append(
            BoltzOutputs(
                out_dir=out_dir,
                pred_dir=pred_dir,
                stem=stem,
                model=idx,
                pae_npz=pae,
                structure=struct[0],
                plddt_npz=_optional(f"plddt_*_model_{idx}.npz"),
                confidence_json=_optional(f"confidence_*_model_{idx}.json"),
            )
        )

    samples.sort(key=lambda s: s.model)
    return samples


def describe_boltz_npz(path: str | Path) -> str:
    """Report an npz's keys, shapes and dtypes.

    Run this first when the real output lands -- its output is exactly what
    `_boltz_arrays` below needs to be filled in with.

        python -c "import sys; sys.path.insert(0,'src/hedi'); \
import ipsae; print(ipsae.describe_boltz_npz(sys.argv[1]))" <file.npz>
    """
    with np.load(path, allow_pickle=False) as z:
        lines = [f"{path}:"]
        for k in z.files:
            a = z[k]
            lines.append(f"  {k!r}: shape={a.shape} dtype={a.dtype}")
        if not z.files:
            lines.append("  (no arrays)")
    return "\n".join(lines)


@dataclass
class BoltzArrays:
    """Raw arrays pulled out of Boltz's npz files, before any masking."""

    pae: np.ndarray
    plddt: np.ndarray | None


BOLTZ_PAE_KEY = "pae"
BOLTZ_PLDDT_KEY = "plddt"


def _boltz_arrays(paths: BoltzOutputs) -> BoltzArrays:
    """THE ONLY place Boltz array keys are read.

    Verified against boltz 2.2.1 output (outputs/parent_quick, aida branch):

        pae_<run>_model_<N>.npz    -> 'pae',   (384, 384) float32
        plddt_<run>_model_<N>.npz  -> 'plddt', (384,)     float32, 0-1 scale

    384 is exactly the residue count (242 aa scFv + 142 aa HER2 domain IV),
    so for a protein-only complex Boltz's PAE is indexed by RESIDUE and no
    token mask is needed. That is asserted rather than assumed: read_boltz
    compares the PAE shape against the parsed structure and refuses a
    mismatch. A complex containing a ligand or a modified residue would
    introduce extra tokens and trip that guard -- which is the correct
    outcome, since the mask logic to handle it is not written.

    pLDDT came back on a 0-1 scale here; read_boltz rescales to 0-100 when
    the max is <= 1.0, because Boltz has shipped both conventions.

    If a key is missing we raise with the keys that ARE present, rather than
    letting a KeyError surface with no context.
    """
    def _get(path: Path, key: str, what: str) -> np.ndarray:
        with np.load(path, allow_pickle=False) as z:
            if key not in z.files:
                raise KeyError(
                    f"{path.name} has no {key!r} array (holds {list(z.files)}). "
                    f"If boltz renamed the {what} key, fix it in _boltz_arrays "
                    "-- it is the only place these keys are read."
                )
            return np.array(z[key])

    pae = _get(paths.pae_npz, BOLTZ_PAE_KEY, "PAE")
    plddt = None
    if paths.plddt_npz is not None:
        plddt = _get(paths.plddt_npz, BOLTZ_PLDDT_KEY, "pLDDT")
    return BoltzArrays(pae=pae, plddt=plddt)


def _cif_tokens(line: str) -> list[str]:
    """Split an mmCIF data line, honouring single/double quoted values.

    Atom names like C1' and O5' are quoted in mmCIF, so a plain .split()
    would mangle them.
    """
    out, cur, quote = [], "", None
    for ch in line:
        if quote:
            if ch == quote:
                quote = None
            else:
                cur += ch
        elif ch in "'\"":
            quote = ch
        elif ch.isspace():
            if cur:
                out.append(cur)
                cur = ""
        else:
            cur += ch
    if cur:
        out.append(cur)
    return out


def read_cif_atoms(path: str | Path) -> list[dict]:
    """Parse the _atom_site loop of an mmCIF file.

    mmCIF declares its own column names in the loop_ header, so this reads
    them by name rather than by position -- no column-order assumptions, and
    nothing here depends on Boltz's particular field ordering.
    """
    if str(path).endswith(".pdb"):
        raise ValueError("read_cif_atoms called on a .pdb file; use read_pdb")
    text = Path(path).read_text()

    fields: list[str] = []
    atoms: list[dict] = []
    in_loop = False

    for raw in text.splitlines():
        line = raw.strip()
        if line.startswith("_atom_site."):
            fields.append(line.split(".", 1)[1].split()[0])
            in_loop = True
            continue
        if in_loop and (line.startswith("ATOM") or line.startswith("HETATM")):
            parts = _cif_tokens(line)
            if len(parts) < len(fields):
                continue
            row = dict(zip(fields, parts))

            def pick(*names: str) -> str | None:
                for n in names:
                    v = row.get(n)
                    if v is not None and v not in (".", "?"):
                        return v
                return None

            resname = pick("label_comp_id", "auth_comp_id")
            if resname is None or resname == "LIG":
                continue
            chain = pick("auth_asym_id", "label_asym_id")
            resnum = pick("auth_seq_id", "label_seq_id")
            if chain is None or resnum is None:
                continue
            atoms.append(
                {
                    "atom_name": pick("label_atom_id", "auth_atom_id") or "",
                    "resname": resname,
                    "chain": chain,
                    "resnum": int(resnum),
                    "xyz": (
                        float(row["Cartn_x"]),
                        float(row["Cartn_y"]),
                        float(row["Cartn_z"]),
                    ),
                }
            )
        elif in_loop and line and not line.startswith(("#", "ATOM", "HETATM")):
            if atoms:  # loop finished
                break

    if not atoms:
        raise ValueError(f"{Path(path).name}: no _atom_site ATOM/HETATM rows parsed")
    return atoms


def _ca_cb(atoms: list[dict], source: str) -> tuple[list[dict], list[dict]]:
    """Split an atom list into the parallel CA and CB lists, and check them."""
    ca, cb = [], []
    for atom in atoms:
        name, resname = atom["atom_name"], atom["resname"]
        if name == "CA" or "C1" in name:
            ca.append(atom)
        if name == "CB" or "C3" in name or (resname == "GLY" and name == "CA"):
            cb.append(atom)
    if len(ca) != len(cb):
        raise ValueError(
            f"{source}: {len(ca)} CA atoms but {len(cb)} CB atoms. These are "
            "indexed in parallel, so distances would be misaligned."
        )
    for i, (a, b) in enumerate(zip(ca, cb)):
        if a["chain"] != b["chain"] or a["resnum"] != b["resnum"]:
            raise ValueError(
                f"{source}: CA/CB lists diverge at index {i} "
                f"({a['chain']}{a['resnum']} vs {b['chain']}{b['resnum']})"
            )
    if not ca:
        raise ValueError(f"{source}: no CA atoms found")
    return ca, cb


def chain_sequence(model: Model, chain: str) -> str:
    """One-letter sequence of one chain, in structure order. X for unknown."""
    return "".join(
        AA3TO1.get(model.resnames[i], "X")
        for i in range(model.numres)
        if model.chains[i] == chain
    )


def assert_chain_layout(
    model: Model,
    scfv_chain: str = "A",
    target_chain: str = "B",
    scfv_len: int | None = PARENT_SCFV_LEN,
) -> None:
    """Check chain A is the scFv and chain B the target, from the STRUCTURE.

    Chain identity is derived from the chain definitions in the output, never
    from confidence_*.json: `pair_chains_iptm` there is keyed by stringified
    integers whose mapping onto the YAML chain IDs is undocumented, so using
    it to decide which chain is the binder is a silent-wrong-answer risk.

    The scFv is identified positively by its (GGGGS)x3 linker, which is the
    one feature of our construct that cannot be coincidental. Length is
    checked too, but a design is allowed to differ from the parent's 242 --
    pass scfv_len=None to skip that, and scripts/validate.py is what actually
    enforces sequence rules.
    """
    present = model.unique_chains
    missing = [c for c in (scfv_chain, target_chain) if c not in present]
    if missing:
        raise ValueError(
            f"expected chains {scfv_chain} (scFv) and {target_chain} (target), "
            f"missing {missing}; structure has {present}"
        )
    if len(present) > 2:
        raise ValueError(
            f"expected exactly 2 chains, found {present}. Trim the target to "
            "HER2 domain IV and cofold a single scFv against it."
        )

    seq = chain_sequence(model, scfv_chain)
    linker = "GGGGS" * 3
    if linker not in seq:
        other = chain_sequence(model, target_chain)
        hint = (
            f" Chain {target_chain} DOES contain it -- the chains look swapped."
            if linker in other
            else ""
        )
        raise ValueError(
            f"chain {scfv_chain} has no (GGGGS)x3 linker, so it is not the "
            f"scFv (length {len(seq)}).{hint}"
        )
    if scfv_len is not None and len(seq) != scfv_len:
        raise ValueError(
            f"chain {scfv_chain} is {len(seq)} aa, expected {scfv_len}. "
            "Pass scfv_len=None if this design is deliberately a different "
            "length, then check it with scripts/validate.py."
        )
    if not chain_sequence(model, target_chain):
        raise ValueError(f"chain {target_chain} (target) has no residues")


def read_boltz_sample(
    paths: BoltzOutputs,
    *,
    scfv_chain: str = "A",
    target_chain: str = "B",
    scfv_len: int | None = PARENT_SCFV_LEN,
    check_chains: bool = True,
) -> Model:
    """Read ONE diffusion sample into a Model."""
    if str(paths.structure).endswith(".cif"):
        atoms = read_cif_atoms(paths.structure)
    else:
        parsed = read_pdb(paths.structure)
        atoms = None
        ca, cb = parsed["ca"], parsed["cb"]
    if atoms is not None:
        ca, cb = _ca_cb(atoms, Path(paths.structure).name)
    n = len(ca)

    arrays = _boltz_arrays(paths)

    # --- PAE must exist. A Boltz run without --write_full_pae produces the
    # other confidence files but no PAE matrix, and ipSAE is undefined without
    # it. This must never degrade to a default or a zero matrix. ---
    if arrays.pae is None or getattr(arrays.pae, "size", 0) == 0:
        raise ValueError(
            f"{paths.pae_npz} contains no PAE matrix. ipSAE cannot be computed "
            "without one. This is the --write_full_pae failure mode: re-run "
            "Boltz with --write_full_pae."
        )
    pae = np.asarray(arrays.pae, dtype=float)
    if pae.ndim != 2 or pae.shape[0] != pae.shape[1]:
        raise ValueError(f"PAE from {paths.pae_npz} has shape {pae.shape}, expected square 2-D")
    if pae.shape != (n, n):
        raise ValueError(
            f"PAE is {pae.shape} but the structure has {n} residues. The PAE is "
            "probably over TOKENS rather than residues, which needs a token "
            "mask selecting CA rows/columns -- see the TODO in _boltz_arrays."
        )

    if arrays.plddt is None:
        plddt = np.zeros(n)
    else:
        plddt = np.asarray(arrays.plddt, dtype=float)
        if len(plddt) != n:
            raise ValueError(
                f"plddt has {len(plddt)} entries but the structure has {n} residues"
            )
        if plddt.size and plddt.max() <= 1.0:
            plddt = plddt * 100.0  # Boltz has shipped both 0-1 and 0-100

    coords = np.array([a["xyz"] for a in cb], dtype=float)
    distances = np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=-1)

    m = Model(
        pae=pae,
        chains=np.array([a["chain"] for a in ca]),
        resnums=np.array([a["resnum"] for a in ca]),
        resnames=[a["resname"] for a in ca],
        plddt=plddt,
        distances=distances,
        # iptm_af is deliberately left at -1.0: it would come from
        # confidence_*.json's pair_chains_iptm, whose integer-keyed chain
        # mapping is undocumented. ipSAE does not need it.
        name=f"{paths.stem}_model_{paths.model}",
    )
    m.validate()
    if check_chains:
        assert_chain_layout(
            m, scfv_chain=scfv_chain, target_chain=target_chain, scfv_len=scfv_len
        )
    return m


def read_boltz(out_dir: str | Path, model: int = 0, **kw) -> Model:
    """Read a single named diffusion sample of a Boltz run."""
    samples = find_boltz_samples(out_dir)
    for s in samples:
        if s.model == model:
            return read_boltz_sample(s, **kw)
    raise ValueError(
        f"model {model} not found in {out_dir}; available: "
        f"{[s.model for s in samples]}"
    )


@dataclass
class RunScore:
    """ipSAE across every diffusion sample of one run, for one chain pair.

    Aggregated by MEDIAN, not mean: with 5 samples a single collapsed pose
    would drag a mean down and misrepresent the run, while the median
    reports what the sampler typically produced.

    `values` keeps every per-sample score and `spread` (max - min) is
    retained deliberately: a tight ensemble is a more trustworthy number
    than a scattered one, so spread is available as a tie-breaker between
    designs whose medians are too close to separate on their own.
    """

    run: str
    chain1: str
    chain2: str
    pae_cutoff: float
    dist_cutoff: float
    models: list[int]
    values: list[float]
    records: list[Record] = field(default_factory=list, repr=False)

    @property
    def median(self) -> float:
        return float(np.median(self.values))

    @property
    def spread(self) -> float:
        return float(max(self.values) - min(self.values)) if self.values else 0.0

    @property
    def minimum(self) -> float:
        return float(min(self.values)) if self.values else 0.0

    @property
    def maximum(self) -> float:
        return float(max(self.values)) if self.values else 0.0

    @property
    def stdev(self) -> float:
        return float(np.std(self.values, ddof=1)) if len(self.values) > 1 else 0.0

    @property
    def n_samples(self) -> int:
        return len(self.values)

    def per_sample(self) -> str:
        return ", ".join(
            f"m{m}={v:.4f}" for m, v in zip(self.models, self.values)
        )


def score_boltz_run(
    out_dir: str | Path,
    pae_cutoff: float = 10.0,
    dist_cutoff: float = 15.0,
    *,
    chain1: str = "A",
    chain2: str = "B",
    **kw,
) -> tuple[RunScore, list[list[Record]]]:
    """Score every diffusion sample of a run and aggregate by median.

    Returns (RunScore for the chain1-chain2 `max` record, per-sample full
    record lists) so nothing is thrown away -- the individual values, the
    other chain-pair records and all the secondary metrics stay available.
    """
    samples = find_boltz_samples(out_dir)
    models, values, chosen, every = [], [], [], []

    for s in samples:
        model = read_boltz_sample(s, **kw)
        records = score_pairs(model, pae_cutoff, dist_cutoff)
        every.append(records)
        hit = next(
            (
                r
                for r in records
                if r.kind == "max"
                and {r.chain1, r.chain2} == {chain1, chain2}
            ),
            None,
        )
        if hit is None:
            raise ValueError(
                f"no max record for {chain1}-{chain2} in {s.pae_npz.name}; "
                f"chains present: {model.unique_chains}"
            )
        models.append(s.model)
        values.append(hit.ipsae)
        chosen.append(hit)

    return (
        RunScore(
            run=samples[0].stem,
            chain1=chain1,
            chain2=chain2,
            pae_cutoff=pae_cutoff,
            dist_cutoff=dist_cutoff,
            models=models,
            values=values,
            records=chosen,
        ),
        every,
    )


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------

COLUMNS = (
    "chain1", "chain2", "kind", "ipsae", "ipsae_d0chn", "ipsae_d0dom",
    "iptm_af", "iptm_d0chn", "pdockq", "pdockq2", "lis",
    "n0res", "n0chn", "n0dom", "d0res", "d0chn", "d0dom",
    "nres1", "nres2", "dist1", "dist2",
)


def _main_boltz(out_dir: str, pae_cutoff: float, dist_cutoff: float) -> int:
    """Score a whole Boltz run: every diffusion sample, aggregated by median."""
    run, every = score_boltz_run(out_dir, pae_cutoff, dist_cutoff)
    print(f"# run {run.run}  pae_cutoff={pae_cutoff}  dist_cutoff={dist_cutoff}")
    print(f"# {run.n_samples} diffusion sample(s): {run.per_sample()}")
    print()
    print("\t".join(("run", "chain1", "chain2", "n", "median", "min", "max", "spread", "stdev")))
    print(
        "\t".join(
            (
                run.run, run.chain1, run.chain2, str(run.n_samples),
                f"{run.median:.6f}", f"{run.minimum:.6f}", f"{run.maximum:.6f}",
                f"{run.spread:.6f}", f"{run.stdev:.6f}",
            )
        )
    )
    print()
    print("# per sample, all metrics")
    print("\t".join(("model",) + COLUMNS))
    for m, records in zip(run.models, every):
        for r in records:
            row = [str(m)]
            for c in COLUMNS:
                v = getattr(r, c)
                row.append(f"{v:.6f}" if isinstance(v, float) else str(v))
            print("\t".join(row))
    return 0


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print(__doc__.strip().splitlines()[-1])
        return 2

    # A directory means a Boltz run; two files mean an AF2-style pair.
    if Path(argv[1]).is_dir():
        pae_cutoff = float(argv[2]) if len(argv) > 2 else 10.0
        dist_cutoff = float(argv[3]) if len(argv) > 3 else 15.0
        return _main_boltz(argv[1], pae_cutoff, dist_cutoff)

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
