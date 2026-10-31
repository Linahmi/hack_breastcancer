# structure/ — mutable vs non-mutable positions (Ghita)

Structural constraints for the trastuzumab scFv, derived from **PDB 1N8Z**
(HER2 ectodomain + trastuzumab Fab, 2.52 Å). These files tell the design side
which positions are free to mutate and which must be conserved, and give the
cofolding side a trimmed HER2 target.

## Files

| File | What it is |
|---|---|
| `red_green.csv` | The deliverable. 61 rows: 21 RED + 40 GREEN, with `position,aa,class,min_dist_A`. |
| `ghita_structure_report.md` | Method + per-CDR breakdown of the RED/GREEN split, with per-residue minimum distances. |
| `her2_domainIV.pdb` | HER2 domain IV (residues 489–607), chain C, extracted from 1N8Z — the binding target for Boltz-2 cofolding. |
| `contacts.py` | The analysis script that produces `red_green.csv` from the 1N8Z structure. |
| `view_heroes.cxc` | ChimeraX script: HER2 surface, RED contacts in red, GREEN mutable positions in lime. |

## Numbering

**The brief's boundaries were off by one.** The true scFv layout on our parent
sequence (1–242) is:

- **VH 1–120**
- **linker 121–135**
- **VL 136–242**

All positions in `red_green.csv` use this parent numbering, verified residue by
residue against `data/parent.fasta` — every mapped PDB residue matches.
`view_heroes.cxc` uses the structure's own chain/residue numbering instead
(chain A = VL, chain B = VH, chain C = HER2), since that is what ChimeraX reads
from 1N8Z.

## How RED and GREEN were defined

- **RED (21)** — any scFv residue with an atom within **4.5 Å** of any HER2 atom.
  These form the proven paratope; do not mutate them.
- **GREEN (40)** — CDR positions (**Kabat** definition) that are *not* in contact.
  These are the ones that are safe to scan.
- Framework positions outside the CDRs are neither: leave them alone.

The 4.5 Å cutoff is robust — a core of 17 residues stays RED across the whole
4.0–5.0 Å range, so the split does not hinge on the exact threshold.

GREEN was additionally checked to contain **no cysteines** (22, 96, 158, 223)
and **no linker positions** (121–135), so it can be fed to a PLM scan as-is
without extra filtering.

## Rerunning the contact analysis

Needs `biopython`. Fetch the structure once, then run the script with the CIF
as its argument:

```bash
# from the repo root
curl -O https://files.rcsb.org/download/1N8Z.cif
python3 structure/contacts.py 1N8Z.cif > structure/red_green.csv
```

Output is `red_green.csv` as above. `1N8Z.cif` itself is not committed — it is a
public download, not our data.

## Downstream

`red_green.csv` is the position whitelist for the mutation/scanning step; the
shared submission file stays `results/candidates.csv`
(`name,sequence,rationale,tool`, pipeline adds `ipsae`) and is not owned here.
