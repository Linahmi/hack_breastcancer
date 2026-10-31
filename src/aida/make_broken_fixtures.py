"""Build a CSV of deliberately broken candidates, one per rule.

Usage: python3 -I mkbroken.py <out.csv>
"""
import csv
import re
import sys

PARENT = (
    "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFT"
    "ISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSSGGGGSGGGGSGGGGSDIQ"
    "MTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFT"
    "LTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK"
)
LINKER = "GGGGSGGGGSGGGGS"
li = PARENT.find(LINKER)
VH, VL = PARENT[:li], PARENT[li + len(LINKER):]


def mut(seq, changes):
    s = list(seq)
    for pos, aa in changes:
        s[pos - 1] = aa
    return "".join(s)


def sequons(s):
    return len(re.findall(r"(?=N[^P][ST])", s))


def ngs(s):
    return len(re.findall(r"(?=NG)", s))


# 5 safe framework mutations, chosen programmatically: outside CDRs and the
# linker, not a cysteine, actually changes the residue, and adds no sequon/NG.
CDR_SEQS = ["GFNIKDTYIH", "RIYPTNGYTRYADSVKG", "WGGDGFYAMDY",
            "RASQDVNTAVA", "SASFLYS", "QQHYTTPPT"]
cdr_pos = set()
for c in CDR_SEQS:
    i = PARENT.find(c)
    assert i >= 0, c
    cdr_pos |= set(range(i + 1, i + len(c) + 1))
link_pos = set(range(li + 1, li + len(LINKER) + 1))

SAFE = []
for p in range(1, len(PARENT) + 1):
    if len(SAFE) == 5:
        break
    if p in cdr_pos or p in link_pos or PARENT[p - 1] in "CA":
        continue
    cand = mut(PARENT, SAFE + [(p, "A")])
    if sequons(cand) > sequons(PARENT) or ngs(cand) > ngs(PARENT):
        continue
    SAFE.append((p, "A"))
assert len(SAFE) == 5, SAFE
print(f"safe framework mutations: {[(p, PARENT[p-1] + str(p) + 'A') for p, _ in SAFE]}")
good = mut(PARENT, SAFE)
assert len(good) == 242
assert sequons(good) <= sequons(PARENT), "test seq accidentally adds sequon"
assert ngs(good) <= ngs(PARENT), "test seq accidentally adds NG"
assert good.count("C") == 4
assert sum(1 for a, b in zip(good, PARENT) if a != b) == 5

four = mut(PARENT, SAFE[:4])                      # R1: only 4 mutations
too_long = good + "AAAAAAAAAAAAAAAAAAA"           # R2: 261 aa
bad_linker = good.replace(LINKER, "GGGGSGGAGSGGGGS")  # R4: linker mutated
cys_mut = mut(good, [(22, "S")])                  # R5 + R6: conserved Cys gone
reverse = VL + LINKER + VH                        # R3: reversed order
# R7: new N-glyc sequon. Put N-A-T into framework.
glyc = mut(good, [(60, "N"), (61, "A"), (62, "T")])
# R8: new NG inside CDR-H3 (WGGDGFYAMDY at 99-109); make pos 103-104 = NG
h3 = PARENT.find("WGGDGFYAMDY") + 1
ng_cdr = mut(good, [(h3 + 4, "N"), (h3 + 5, "G")])

# her2_009: destroy conserved Cys22 BUT add a cysteine elsewhere, so the
# total count stays 4 and a count-only check cannot see it.
_swap = None
for p in range(1, len(PARENT) + 1):
    if p in cdr_pos or p in link_pos:
        continue
    if PARENT[p - 1] == "S":
        _swap = p
        break
assert _swap
cys_swap = mut(good, [(22, "S"), (_swap, "C")])
assert cys_swap.count("C") == 4, cys_swap.count("C")

# her2_010: a LEGITIMATE design using the 18 aa headroom the README
# encourages -- 5 substitutions plus a 6 aa insertion, 248 aa total.
_h3end = PARENT.find("WGGDGFYAMDY") + len("WGGDGFYAMDY")
headroom = good[:_h3end] + "AAAAAA" + good[_h3end:]
assert len(headroom) == 248, len(headroom)
assert headroom.count("C") == 4
assert sequons(headroom) <= sequons(PARENT)
assert ngs(headroom) <= ngs(PARENT)

rows = [
    ("her2_001", good,       "should PASS"),
    ("her2_009", cys_swap,   "conserved Cys22 destroyed, count still 4"),
    ("her2_010", headroom,   "LEGITIMATE 248 aa design using headroom"),
    ("her2_002", four,       "R1 only 4 mutations"),
    ("her2_003", too_long,   "R2 261 aa"),
    ("her2_004", bad_linker, "R4 linker mutated"),
    ("her2_005", cys_mut,    "R5/R6 conserved Cys22 mutated"),
    ("her2_006", reverse,    "R3 reversed VL-linker-VH"),
    ("her2_007", glyc,       "R7 new N-glycosylation sequon"),
    ("her2_008", ng_cdr,     "R8 new NG inside CDR-H3"),
    ("design_9", good,       "R9 bad name + R10 dup sequence"),
    ("her2_pae1", good,      "R9 name contains 'pae' + bad pattern"),
    ("her2_001", four,       "R10 duplicate name her2_001"),
]

with open(sys.argv[1], "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["name", "sequence", "note"])
    for n, s, note in rows:
        w.writerow([n, s, note])

print(f"wrote {sys.argv[1]} with {len(rows)} rows")
for n, s, note in rows:
    print(f"  {n:12} len={len(s):3}  {note}")
