"""Submission gate for candidate scFv designs.

    python src/aida/validate_candidates.py results/candidates.csv

Takes a CSV with at least `name,sequence`. Prints one line per broken rule
per row, then a PASS/FAIL summary. Exits non-zero if any row fails.

NOTE: this is deliberately a SEPARATE file from scripts/validate.py, which
is Lina's (skeleton commit a452d0a). Hers checks format/developability;
this one adds the submission-gate rules (name pattern, duplicates,
conserved-cysteine identity, indel-tolerant mutation counting, CDR-scoped
NG). Run both, or fold these rules into hers once she agrees.

Rules enforced
--------------
R1  >= 5 mutations vs parent (indel-tolerant, counted programmatically)
R2  length <= 260 aa
R3  single chain in VH-(GGGGS)3-VL order
R4  (GGGGS)3 linker present exactly once and unmutated
R5  no conserved cysteine mutated (VH C22/C96, VL C158/C223)
R6  no unpaired cysteine introduced (total must stay 4)
R7  no new N-glycosylation sequon N-X-S/T, X != P
R8  no new NG motif created inside a CDR
R9  name matches ^her2_\\d{3}$ and contains no "pae"
R10 no duplicate names, no duplicate sequences
"""

import csv
import re
import sys
from difflib import SequenceMatcher

PARENT = (
    "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFT"
    "ISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSSGGGGSGGGGSGGGGSDIQ"
    "MTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFT"
    "LTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK"
)
LINKER = "GGGGSGGGGSGGGGS"
VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")
MIN_MUTATIONS = 5
MAX_LENGTH = 260
NAME_RE = re.compile(r"^her2_\d{3}$")

_li = PARENT.find(LINKER)
VH, VL = PARENT[:_li], PARENT[_li + len(LINKER):]
CONSERVED_CYS = [i + 1 for i, a in enumerate(PARENT) if a == "C"]  # 22,96,158,223

# CDRs located by exact substring in the parent, so the positions are derived
# rather than hand-typed. 1-based inclusive ranges.
_CDR_SEQS = {
    "H1": "GFNIKDTYIH", "H2": "RIYPTNGYTRYADSVKG", "H3": "WGGDGFYAMDY",
    "L1": "RASQDVNTAVA", "L2": "SASFLYS", "L3": "QQHYTTPPT",
}
CDRS = {}
for _n, _s in _CDR_SEQS.items():
    _i = PARENT.find(_s)
    if _i < 0:
        raise SystemExit(f"internal error: CDR-{_n} {_s} not found in parent")
    CDRS[_n] = (_i + 1, _i + len(_s))


def identity(a: str, b: str) -> float:
    """Fraction of aligned positions that match (0..1)."""
    if not a or not b:
        return 0.0
    return SequenceMatcher(None, a, b, autojunk=False).ratio()


def parent_index_map(cand: str):
    """Map 1-based parent position -> 1-based candidate position for aligned
    (equal or replaced) blocks. Deleted parent positions map to None."""
    sm = SequenceMatcher(None, PARENT, cand, autojunk=False)
    m = {}
    for tag, i1, i2, j1, j2 in sm.get_opcodes():
        if tag in ("equal", "replace"):
            for k in range(min(i2 - i1, j2 - j1)):
                m[i1 + k + 1] = j1 + k + 1
    return m


def count_mutations(cand: str) -> int:
    """Differences vs parent, tolerant of insertions/deletions."""
    sm = SequenceMatcher(None, PARENT, cand, autojunk=False)
    n = 0
    for tag, i1, i2, j1, j2 in sm.get_opcodes():
        if tag == "replace":
            n += max(i2 - i1, j2 - j1)
        elif tag in ("delete", "insert"):
            n += max(i2 - i1, j2 - j1)
    return n


def sequons(s: str):
    """1-based start positions of N-X-S/T, X != P."""
    return {m.start() + 1 for m in re.finditer(r"(?=N[^P][ST])", s)}


def ng_sites(s: str):
    return {m.start() + 1 for m in re.finditer(r"(?=NG)", s)}


def in_any_cdr(pos: int) -> str | None:
    for n, (lo, hi) in CDRS.items():
        if lo <= pos <= hi:
            return n
    return None


def check(name: str, seq: str) -> list[str]:
    """Return a list of broken-rule messages. Empty list == row passes."""
    errs = []
    seq = (seq or "").strip().upper()

    # R9 name
    if not NAME_RE.match(name or ""):
        errs.append(f"R9 name {name!r} does not match ^her2_\\d{{3}}$")
    if "pae" in (name or "").lower():
        errs.append(f"R9 name {name!r} contains 'pae'")

    if not seq:
        errs.append("R0 empty sequence")
        return errs
    bad = sorted(set(seq) - VALID_AA)
    if bad:
        errs.append(f"R0 invalid residues {bad}")
        return errs

    # R2 length
    if len(seq) > MAX_LENGTH:
        errs.append(f"R2 length {len(seq)} > {MAX_LENGTH}")

    # R4 linker present exactly once and unmutated
    hits = [m.start() for m in re.finditer(re.escape(LINKER), seq)]
    if len(hits) != 1:
        errs.append(
            f"R4 (GGGGS)3 linker found {len(hits)}x, expected exactly 1 "
            f"(a mutated linker will not match)"
        )
    else:
        # R3 chain order: segment before linker must look like VH, not VL
        pre, post = seq[:hits[0]], seq[hits[0] + len(LINKER):]
        id_pre_vh, id_pre_vl = identity(pre, VH), identity(pre, VL)
        if id_pre_vh < id_pre_vl:
            errs.append(
                f"R3 chain order looks reversed: segment before the linker "
                f"matches VL better than VH ({id_pre_vl:.2f} vs {id_pre_vh:.2f})"
            )
        elif identity(post, VL) < identity(post, VH):
            errs.append("R3 chain order wrong: segment after linker matches VH")

    # R1 mutations
    n_mut = count_mutations(seq)
    if n_mut < MIN_MUTATIONS:
        errs.append(f"R1 only {n_mut} mutations vs parent, need >= {MIN_MUTATIONS}")

    # R5 conserved cysteines still cysteine
    imap = parent_index_map(seq)
    for cp in CONSERVED_CYS:
        q = imap.get(cp)
        if q is None:
            errs.append(f"R5 conserved Cys{cp} deleted")
        elif seq[q - 1] != "C":
            errs.append(f"R5 conserved Cys{cp} mutated to {seq[q - 1]}")

    # R6 cysteine pairing
    nc = seq.count("C")
    if nc != PARENT.count("C"):
        errs.append(
            f"R6 cysteine count {nc} != parent {PARENT.count('C')} "
            f"({'unpaired cysteine introduced' if nc % 2 else 'pairing changed'})"
        )

    # R7 new N-glycosylation sequons
    new_seq = len(sequons(seq)) - len(sequons(PARENT))
    if new_seq > 0:
        errs.append(f"R7 {new_seq} new N-glycosylation sequon(s) N-X-S/T")

    # R8 new NG inside a CDR.
    # Map each CDR's parent range into candidate coordinates, rather than
    # mapping individual NG positions -- a 2-residue substitution can be
    # emitted by SequenceMatcher as delete+insert, leaving those positions
    # absent from imap and silently skipping the check.
    n_parent_cdr_ng = len([p for p in ng_sites(PARENT) if in_any_cdr(p)])
    n_cand_cdr_ng = 0
    cand_cdr = {}
    for cname, (lo, hi) in CDRS.items():
        mapped = [imap[p] for p in range(lo, hi + 1) if p in imap]
        if mapped:
            clo, chi = min(mapped), max(mapped)
        else:
            sub = PARENT[lo - 1:hi]
            j = seq.find(sub)
            if j < 0:
                continue
            clo, chi = j + 1, j + len(sub)
        cand_cdr[cname] = (clo, chi)
        for p in ng_sites(seq):
            if clo <= p <= chi:
                n_cand_cdr_ng += 1
    new_cdr_ng = n_cand_cdr_ng - n_parent_cdr_ng
    if new_cdr_ng > 0:
        where = [f"CDR-{c}" for c, (lo, hi) in cand_cdr.items()
                 if any(lo <= p <= hi for p in ng_sites(seq))]
        errs.append(
            f"R8 {new_cdr_ng} new NG motif(s) inside a CDR "
            f"(NG present in {', '.join(sorted(where))})"
        )

    return errs


def main(path: str) -> int:
    try:
        rows = list(csv.DictReader(open(path)))
    except OSError as e:
        print(f"cannot read {path}: {e}")
        return 1
    if not rows:
        print(f"{path}: no data rows")
        return 1
    if "sequence" not in rows[0]:
        print(f"{path}: missing 'sequence' column; has {sorted(rows[0])}")
        return 1

    seen_name, seen_seq = {}, {}
    failed = 0
    for i, row in enumerate(rows, start=2):
        name = (row.get("name") or "").strip()
        seq = (row.get("sequence") or "").strip().upper()
        errs = check(name, seq)

        # R10 duplicates
        if name in seen_name:
            errs.append(f"R10 duplicate name, first seen on line {seen_name[name]}")
        else:
            seen_name[name] = i
        if seq and seq in seen_seq:
            errs.append(
                f"R10 duplicate sequence, identical to {seen_seq[seq][0]!r} "
                f"on line {seen_seq[seq][1]}"
            )
        else:
            seen_seq[seq] = (name, i)

        if errs:
            failed += 1
            print(f"FAIL  line {i}  {name or '<no name>'}")
            for e in errs:
                print(f"        {e}")

    total = len(rows)
    print(f"\n{total - failed}/{total} passed")
    print("PASS" if not failed else f"FAIL ({failed} row(s) broke at least one rule)")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else "results/candidates.csv"))
