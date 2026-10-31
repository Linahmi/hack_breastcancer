# What I need from Ghita, and what I can do without it

## The file: `results/positions.csv`

Columns, exactly as the README fixes them:

```
position,wt_aa,region,imgt,category,min_distance_A,delta_sasa
```

What I actually consume:

| Column | Need | Why |
|---|---|---|
| `position` | **required** | must be **1-based on the 242 aa parent**, not PDB numbering, not IMGT |
| `wt_aa` | **required** | my integrity check, see below |
| `category` | **required** | exactly `contact` / `second_shell` / `free` — three values. This is what forms the arms |
| `min_distance_A` | nice to have | ranks *within* the contact set |
| `delta_sasa` | nice to have | ranks *within* second shell |
| `region` | not needed | I derive it from `scripts.parent.region()` |
| `imgt` | not needed | hers; I can't produce it anyway (no `hmmscan` on this machine) |

If she's short on time, **`position` + `category` alone unblocks me.** Even a list of contact positions in chat is enough to start the contact arm — the CSV can follow.

## The one thing most likely to go wrong

Numbering. 1N8Z is the trastuzumab **Fab** bound to HER2: VH and VL are separate chains with their own PDB numbering, and there is no linker. Our parent is a single 242 aa chain where VL starts at 136. Any of three slips — PDB numbering, IMGT numbering, or VL numbered 1–107 instead of 136–242 — produces a file that looks perfectly valid and is silently wrong by a constant offset.

So on load I will assert, for every row:

```python
assert row["wt_aa"] == PARENT[int(row["position"]) - 1]
```

and reject the whole file if any row fails, rather than build candidates on a shifted frame. A VL-offset error would show up as a near-total mismatch; a PDB-numbering error as scattered mismatches. **This is the single highest-risk integration point in my half of the pipeline** — worth 30 seconds of agreement with her before she exports, not after.

Two follow-ups for her, cheap now and expensive later:
- Does `position` for VL residues run 136–242?
- Does `category` use all three values, or did second shell collapse into contact?

## What's blocked and what isn't

| Arm (`tool`) | Needs positions.csv? |
|---|---|
| `plm_consensus` | **no** — pure sequence, runs on the two matrices alone |
| `format` | **no** — linker placement and the 18 aa headroom are from the rules |
| `literature` | **no** — known 4D5 / trastuzumab affinity-maturation positions |
| `contact` | **yes** — that's the definition of the arm |
| `second_shell` | **yes** — and specifically needs the three-way `category`, not a red/green split |

So roughly half the candidate pool does not depend on her at all. The scan (deliverables 1 and 2) is entirely independent — her file is only needed at assembly, deliverable 3.

## If her file is late

1. **Proceed with the unblocked arms** and assemble those candidates first. No waiting.
2. **Do not substitute CDR positions for contact positions.** I can get Kabat/Chothia CDR boundaries from sequence alone, but CDRs are not contacts — plenty of CDR residues point away from the interface, and some framework residues touch it. If I use a sequence-derived set as a stand-in, it ships tagged `plm_consensus`, never `contact`, so her arm is never silently faked.
3. **I won't recompute 1N8Z contacts myself.** It duplicates her work, and the repo rule keeps me in `src/lina/` — I don't write `results/positions.csv` under any circumstances.
4. If it's still missing near the 14:00 milestone, the honest version is a pool that's short on the structural arms, with that stated in the pitch — not a pool padded with a guess wearing her label.
