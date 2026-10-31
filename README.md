# HERoes — Adaptyv track

Breast Cancer Awareness Hackathon, Zurich, 31 October 2026.
Affinity maturation of the trastuzumab scFv against the HER2 ectodomain.

**Team:** Aida (Boltz-2 cofold) · Hedi (ipSAE + ranking) · Lina (PLM scans) · Ghita (1N8Z contacts)

## The task

Parent: 242 aa scFv, VH–(GGGGS)₃–VL, **KD = 0.82 nM** by BLI. That's the number to beat.
Target: HER2 ectodomain, UniProt P04626, Thr23–Thr652. Structure: PDB 1N8Z.

Every submitted design must have **≥5 mutations**, be a **single chain in VH–(GGGGS)₃–VL order**,
and be **≤260 aa**. The parent is 242, so there are **18 aa of headroom** — worth using on at
least one design, since most teams will treat 242 as fixed.

Deliverable: **15 sequences**, each with an ipSAE score and a rationale, plus a 2-minute pitch.
Judged today on the ipSAE of the pool and the pitch. Real KD comes back in ~3 weeks.

## Numbering — read this before you write any code

Positions are **1-based on the 242 aa parent scFv**, never PDB numbering.
Mutations are written `D31A`, wild-type letter, position, mutant letter.

```
VH       1–120    (120 aa)
linker   121–135  (15 aa, fixed)
VL       136–242  (107 aa)
```

Do not hardcode these. Import them:

```python
from scripts.parent import PARENT, VH, LINKER, VL, region
```

`scripts/parent.py` reads the sequence from `data/parent.fasta` and finds the linker itself,
so there is one source of truth and no retyping. Run `python scripts/parent.py` to print them.

## Setup

```bash
python scripts/fetch_inputs.py      # 1N8Z + HER2 ectodomain, both gitignored
python scripts/parent.py            # sanity check
```

## Shared files — do not rename these columns

| file | owner | columns |
|---|---|---|
| `results/positions.csv` | Ghita | `position,wt_aa,region,imgt,category,min_distance_A,delta_sasa` |
| `results/substitutions.csv` | Lina | `position,wt_aa,mut_aa,esm2_score,ablang2_score,germline_flag` |
| `results/candidates.csv` | Lina | `name,sequence,rationale,tool` |
| `results/scored.csv` | Aida/Hedi | candidates + `ipsae` |

`category` in positions.csv is one of `contact` / `second_shell` / `free` — three categories,
not two. Second shell (within 8 Å of the interface, not in direct contact) is where a lot of
affinity maturation actually happens, and a red/green split throws it away.

`rationale` should open with the explicit mutation list so nobody has to re-derive it:

```
D31A;Y33S;N55K;Y57Q;W99V — contact set, both models agree, position tolerant in alanine scan
```

`tool` carries which arm the design came from: `plm_consensus`, `contact`, `second_shell`,
`abmpnn`, `literature`, `format`.

## Validate before you ship

```bash
python scripts/validate.py results/candidates.csv
```

Checks ≥5 mutations, ≤260 aa, linker present once and centred, no new cysteines, no new
N-glycosylation sequons, no new NG or DG motifs. Expression is reported publicly on
ProteinBase, so a design that doesn't express is a visible zero. Run it again on the final 15.

## Controls — owner: Hedi, verdict needed by 13:00

Before trusting any ipSAE ranking, push three sequences through the pipeline:

1. **parent** — 0.82 nM, measured
2. **reverse order** VL–(GGGGS)₃–VH — 1.83 nM, measured (per the Adaptyv deck)
3. **scrambled CDR-H3** — should be clearly worse than both

If ipSAE doesn't separate 1 from 2, it can't resolve a known two-fold difference and it won't
resolve our variants either. If it doesn't rank 3 last, it isn't measuring binding at all.
Either way we want to know at 13:00, not 16:00. Whatever the answer, it goes in the pitch —
testing your own metric is worth more than ranking noise.

## Working rules

- One folder per person under `src/`. Stay in yours.
- Each stage reads a file and writes a file. Nothing passed in memory, so any stage can be
  rerun alone when something breaks.
- Branch per person. Don't push to main, don't merge without asking.
- **16:00 is the freeze.** No sequence changes after that. Pitch gets built 16:00–16:45.

## Links

- Adaptyv protein-design-skills plugin (24 skills) — installed
- `huggingface.co/datasets/Anthropic/claude-protein-binder-design`, see
  `prompts/single_target/EGFR.md` for the campaign prompt structure
- adaptyvbio.com/blog/po101 and po102 — optimising on a small test budget
