# structure/ — mutable vs non-mutable positions (Ghita)

Structural constraints for the trastuzumab scFv, derived from **PDB 1N8Z**
(HER2 ectodomain + trastuzumab Fab, 2.52 Å), plus a literature arm of
mutations with published experimental evidence.

All positions are **1-based on the 242 aa parent** (`scripts/parent.py` is the
source of truth), never PDB numbering.

## Files

| File | What it is |
|---|---|
| `red_green.csv` | 60 rows: 20 RED (contact, conserve) + 40 GREEN (Kabat CDR, mutable), `position,aa,class,min_dist_A`. |
| `contacts.py` | Produces `red_green.csv` from the structure. numpy only. |
| `ca_distances.npy` | 242×242 CA–CA distance matrix in Å, indexed by parent position − 1. |
| `ca_distances.py` | Produces `ca_distances.npy`. |
| `literature_mutations.csv` | 5 mutations with published experimental measurements, translated to parent numbering. |
| `check_literature.py` | Asserts every `wt_aa` matches the parent; cross-references RED/GREEN. |
| `ghita_structure_report.md` | Method + per-CDR breakdown (written against the original water-inclusive run). |
| `her2_domainIV.pdb` | HER2 domain IV (489–607), chain C, the Boltz-2 cofolding target. |
| `view_heroes.cxc` | ChimeraX script: HER2 surface, RED in red, GREEN in lime. |

## Numbering

**The brief's boundaries were off by one.** True layout:

- **VH 1–120**
- **linker 121–135** (GGGGS ×3)
- **VL 136–242**

Verified residue by residue: `contacts.py` aligns each chain to the parent and
asserts identity at every mapped position, so a numbering slip fails loudly
instead of producing a plausible wrong answer. In 1N8Z, chain B = VH, chain A =
VL, chain C = HER2.

Two numbering facts worth knowing before trusting any translation:

- **1N8Z numbers VH sequentially, not Kabat.** Chain B's author numbering is
  identical to the parent index for all 120 positions, and there is no `100a`
  insertion code anywhere in the file.
- **Chain A is Kabat.** VL author numbering + 135 = parent position
  (Kabat L30 → parent 165 = Asn, Kabat L91 → parent 226 = His).

So for CDR-H3, where Kabat inserts (`100a`–`100c`), the translation cannot be
read off the structure. It was derived from the Kabat insertion convention and
then **validated by residue identity** — see the literature arm below.

`view_heroes.cxc` uses the structure's own chain/residue numbering, not parent
numbering, since that is what ChimeraX reads from 1N8Z.

## How RED and GREEN were defined

- **RED** — any scFv position with a heavy atom within **4.5 Å** of any HER2
  heavy atom. Direct paratope, do not mutate.
- **GREEN (40)** — Kabat CDR position *not* in contact. Safe to scan.
- Framework positions outside the CDRs are neither; leave them alone.

Kabat CDRs in parent numbering: H1 31–35, H2 50–66, H3 99–109, L1 159–169,
L2 185–191, L3 224–232.

GREEN contains **no cysteines** (22, 96, 158, 223) and **no linker positions**,
asserted in code, so it can be fed to a PLM scan as-is.

### Water correction — what changed and why

The **first** committed `red_green.csv` listed 21 RED. That run counted
crystallographic waters of chain C as HER2 atoms. Waters in a 2.52 Å structure
are solvent, not epitope, and including them inflates the paratope:

| Position | First CSV | Nearest HER2 **protein** atom | What was at 3–4 Å |
|---|---|---|---|
| 58 T | 3.96 | 4.11 (Glu558 OE2) | HOH 1016 |
| 163 D | 2.47 | 4.50 (Glu598 OE1) | HOH 1045 |
| **201 R** | **3.52** | **5.45** (Glu598 OE1) | HOH 1026 |

**Position 201 had no HER2 protein contact under 4.5 Å at all** — it was RED
because of one water molecule. It was also the only framework position in the
RED list, which was the visible tell.

The CSV in this directory is the **corrected, water-free** version: **20 RED**.
58 and 163 stay RED at their true protein distances; 201 is dropped. The
ChimeraX script and the per-CDR report were updated to match.

**GREEN is unchanged (40)**, because 201 is framework and never was a GREEN
candidate — the design whitelist was never affected by this.

`contacts.py` defaults to excluding waters. `--include-waters` reproduces the
original 21-RED CSV exactly, so the correction is auditable rather than
asserted. Commit `f01f504` holds the pre-correction file.

### Cutoff robustness

`python3 structure/contacts.py 1N8Z.cif --sweep`

Water-free, RED grows 14 → 21 across 4.0–5.0 Å; the **4.0 Å core is 14
residues**. (The "17-residue core" in `ghita_structure_report.md` is the
water-inclusive number — 17 at 4.0 Å with waters counted. Prefer 14.)
Position 184 enters only at ≥4.75 Å.

## Literature arm

`literature_mutations.csv` — mutations with **published experimental**
measurements, translated onto parent numbering. Columns:
`position,wt_aa,mut_aa,reported_effect,source_numbering,source`.

Every row is checked with `python3 structure/check_literature.py`, which
asserts `wt_aa == PARENT[position-1]`. This matters because a wrong
Kabat→parent offset is silent: it still yields a plausible row. All 5 rows
pass, and the residue-identity match is what validates the CDR-H3 insertion
translation (Kabat H95 → W at parent 99, H98 → D at 102, H100a → Y at 105).

The result is an independent check on the structural work:

- **All four published alanine-scan hotspots land on RED positions** —
  Kabat H50, H95, H100a, L91 → parent 50, 99, 105, 226. An experiment from
  1993 agrees with the 2003 crystal structure's contact map.
- **The one published affinity-*improving* mutation lands on GREEN** —
  Kabat H98 → parent 102, Asp→Trp, ~3-fold better than parental hu4D5.

Sources: Kelley & O'Connell, *Biochemistry* 1993;32(27):6828-35
(doi:10.1021/bi00078a005) for the alanine scan; Gerstner, Carter & Lowman,
*J Mol Biol* 2002;321(5):851-62 (doi:10.1016/s0022-2836(02)00677-0) for the
phage-display plasticity scan and the Asp→Trp gain.

**Not included, deliberately.** Mason et al., *Nat Biomed Eng*
2021;5:600-612 (doi:10.1038/s41551-021-00699-9) is a deep mutational scan of
trastuzumab CDR-H3 (~36k classified variants), and Chinery et al. 2024 derive
a >524k-variant CDR-H3 affinity set from it. Both are class-labelled
(high/medium/low binder), not per-mutation ΔΔG, and the per-variant data was
not pulled here. They are the right place to mine CDR-H3 if the literature arm
needs more rows — but rows were not invented from them. Computational-only
predictions were excluded; this arm is measurement-grounded by definition.

ANARCI was not used: it is not installed in any env on this machine
(nor is Biopython), which is why every script here is numpy-only.

## Rerunning

```bash
curl -O https://files.rcsb.org/download/1N8Z.cif
```

```bash
python3 structure/contacts.py 1N8Z.cif > structure/red_green.csv
```

```bash
python3 structure/ca_distances.py 1N8Z.cif structure/ca_distances.npy
```

```bash
python3 structure/check_literature.py
```

`1N8Z.cif` is not committed — it is a public download, not our data.

## Using the distance matrix

```python
import numpy as np
D = np.load("structure/ca_distances.npy")   # (242, 242), Angstrom
D[21, 95]          # positions 22 and 96 -> 6.47 A (VH disulfide)
np.nanmin(D[32])   # closest CA to position 33
```

Rows and columns 120–134 (positions 121–135, the linker) are **NaN**: the
linker has no coordinates because 1N8Z is a Fab, not an scFv. The array is kept
at 242×242 so callers index by parent position without an offset. Use
`np.nan*` functions. Verified on both disulfides: C22–C96 = 6.47 Å,
C158–C223 = 6.63 Å.

## Downstream

`red_green.csv` is the position whitelist for the mutation/scanning step. The
shared submission file stays `results/candidates.csv`
(`name,sequence,rationale,tool`, pipeline adds `ipsae`) and is not owned here.
