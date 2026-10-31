# rationale.md

One section per teammate. Append only; do not edit anyone else's section.

## Aida — cofold (Boltz-2) and the stability filter

### 1. The HER2 construct was wrong, and is now UniProt 511–652

The original target file, `inputs/her2_541_620_from_1n8z.pdb`, was named for
PDB residues 541–620 but **contained 541–607**. Converted to UniProt
numbering (PDB + 22) that is 563–629, and it had three defects:

- a **10-residue internal gap** at UniProt 603–612, a loop that is
  disordered and unresolved in 1N8Z chain C;
- **Cys626 left unpaired** — the 626–634 disulfide is cut, because 634 is
  outside the window;
- 4 intact disulfides where the domain has 10.

Replaced with `src/aida/her2_domainIV.fasta`, **UniProt P04626 511–652**,
142 aa. Chosen by auditing every candidate window against three criteria:

| window | aa | intact SS | cut SS | free Cys | contacts missing |
|---|---|---|---|---|---|
| 563–642 (old YAML) | 80 | 6 | 0 | 0 | 0 |
| 563–629 (old PDB) | 67 | 4 | **1** | **1** | 0 |
| 501–629 (cryo-EM def.) | 129 | 8 | **2** | **2** | 0 |
| 509–604 (DARPin def.) | 96 | 7 | **1** | **1** | **8** |
| **511–652 (chosen)** | **142** | **10** | **0** | **0** | **0** |

511–652 is the only window that is disulfide-closed *and* contains every
HER2 residue within 5 Å of the trastuzumab Fab in 1N8Z. Those contacts
span UniProt **579–627** (19 residues, computed by heavy-atom distance,
not taken from literature), leaving 68 aa of N-terminal and 25 aa of
C-terminal margin.

The end point matters: the C-terminus **must** reach 642, because the
630–642 disulfide closes there. The commonly cited "domain IV ends at
N629" definition cuts two disulfides and leaves two free cysteines, which
is why we did not use it.

Disulfide constraints are **derived per sequence**, never copied. The
reverse-order control's cysteines sit at 23/88/144/218 rather than
22/96/158/223, and the runner finds them by cysteine order, so a
length-changing design shifts its own pairing correctly.

### 2. Multi-seed baseline — and a control that failed

Parent, 3 seeds × 5 diffusion samples, identical settings
(`--use_msa_server --write_full_pae --recycling_steps 10
--diffusion_samples 5 --step_scale 1.0`):

| control | ipSAE (max) | ipTM | seed-spread (ipSAE) |
|---|---|---|---|
| `parent_wt` | **0.9078** | 0.9321 | 0.0603 |
| `reverse_vl_vh` | 0.8842 | 0.9191 | 0.0579 |
| `cdrh3_scrambled` | 0.8002 | 0.8810 | 0.3357 |

**Control 2 fails.** The parent-vs-reverse-order gap is **+0.0236** ipSAE
against a parent seed-spread of **0.0603** — the gap is 0.39× the noise.
The measured difference between these two is real and known: 0.82 nM
versus 1.83 nM, roughly two-fold. Neither metric resolves it.

This is not a metric artefact. We checked ipTM (+0.0129, 0.34× noise),
ipSAE `max` (+0.0236, 0.39×), and ipSAE asymmetric d0res (+0.0372,
0.34×). Three metrics, same answer: the reverse-order gap sits at about a
third of seed noise. Both are dominated by the same seed-to-seed variance
in the underlying Boltz prediction, so swapping the score does not help.

**ipSAE is noisier than ipTM** on this system — seed-spread 0.0603 vs
0.0377, a factor of 1.6 — so it buys discrimination on broken designs
(below) at the cost of some resolution between good ones.

Control 3 passes on both metrics, and decisively once stability is used
rather than score.

### 3. Broken designs are unstable, not low-scoring

The finding that changed our filter. Per-sample ipSAE, never collapsed:

```
parent_wt_seed1          0.9102 0.9165 0.9078 0.8875 0.8885   5/5 docked
parent_wt_seed2          0.9333 0.9428 0.9218 0.8127 0.9177   5/5
parent_wt_seed3          0.9050 0.8268 0.8794 0.8615 0.8513   5/5
reverse_vl_vh_seed1      0.8818 0.9055 0.8842 0.8975 0.7901   5/5
reverse_vl_vh_seed2      0.8985 0.8807 0.8444 0.7779 0.8355   5/5
reverse_vl_vh_seed3      0.9156 0.8870 0.9023 0.9204 0.8596   5/5

cdrh3_scrambled_seed1    0.9276 0.8773 0.9272 0.7153 0.5795   4/5 docked
cdrh3_scrambled_seed2    0.8223 0.8002 0.8290 0.5440 0.0404   3/5
cdrh3_scrambled_seed3    0.8984 0.5416 0.6764 0.0698 0.0000   2/5
```

A broken design does not score low. Each diffusion sample either **docks**
(~0.87) or **collapses** (<0.3), with an empirically empty band between
roughly 0.42 and 0.75. The scrambled control's *best* samples (0.9276,
0.8984) are indistinguishable from the parent's. What distinguishes it is
that only 2–4 of 5 samples dock at all.

**The median hides this completely.** `cdrh3_scrambled_seed1` has median
ipSAE **0.8773**, inside the parent's own seed range of 0.8615–0.9218. A
median-only pipeline passes a design whose CDR-H3 has been shuffled. Had
we run one seed and reported the median, we would have concluded that a
scrambled CDR-H3 binds as well as trastuzumab.

### 4. Why we filter on stability, not on score

The gate runs **before** the median, on the 5 individual samples
(adopted from `src/hedi/ipsae.py`, commit 60e2178):

- **every one of the 5 samples must dock**, ipSAE ≥ **0.50**
- and **spread across the 5 samples ≤ 0.30**

The median does **not** gate. It only ranks designs that already passed.
Applied to the controls: **6/9 PASS, 3/9 REJECT** — every scrambled seed
rejected, every parent and reverse seed kept.

The collapse is visible in the interface geometry, not in confidence:

```
                      ipSAE   n0res   HER2 iface   mean pLDDT
parent, typical       0.910     242          141         85.5
scrambled, docked     0.898     242          141         83.0
scrambled, partial    0.070      43           81         82.7
scrambled, collapsed  0.000       0            0         80.6
```

`n0res` goes 242 → 0 and the HER2 interface goes 142 → 0, while **mean
pLDDT holds at 80–85 throughout**. The model is still confidently folded;
it has simply stopped touching the target. That is why pLDDT cannot
detect collapse and the docking rate can.

**These thresholds are PROVISIONAL.** The controls are extreme by
construction — a shuffled CDR-H3 is a far harsher negative than anything
a PLM scan will produce. If real designs cluster near 0.5 ipSAE or 0.30
spread, refit against the actual distribution rather than applying these
blindly. We report, per design, all 5 individual samples, the docking
rate, the max sample-spread, and the median as a secondary number only —
a design is never collapsed to its median.

One calibration caveat: `cdrh3_scrambled_seed1` docks 5/5 under the 0.50
floor and is rejected on spread alone, 0.3481 against a limit of 0.30 —
a margin of 0.048. The gate catches it, but not with room to spare.

Scoring is `src/aida/score_designs.py`, which uses `src/hedi/ipsae.py` at
`pae_cutoff=10`, `dist_cutoff=15`. ipTM comes from Boltz's own
`confidence_*.json`, not from ipsae.py: its `iptm_af` returns −1.0 by
design, because `pair_chains_iptm` is keyed by undocumented integer
indices and guessing the chain mapping would be worse than reporting
nothing.

### 5. Throughput

Measured over 26 real cofolds at production settings, 384 tokens, L40S:
median **207 s**, p75 **243 s**, range 73–253 s. The variance is MSA-server
latency, not compute.

Concurrency: **8 containers ran with zero queueing** — wall-clock equalled
the slowest single job (253.0 s for 8 jobs; 248.9 s for 6; 142.0 s for 9).
8 is a measured floor, not a ceiling; nothing queued, so the real limit is
higher and untested.

The 5 samples the filter needs come from **one** cofold, not from three
seeds, so the stability filter costs nothing extra. At 8-way concurrency
and 1 cofold per design that is **~270 designs in the 143 minutes before
the freeze**; at 3 seeds per design, ~90. Screen at 1 seed, spend seeds on
the final 15.

### Open issue at time of writing

Which CDR-H3 scramble is canonical is **unresolved**. The outputs on disk
and every number in section 3 use a 13-residue permutation of positions
97–109 (`SRWGGDGFYAMDY` → `RMADDWGFYGSYG`). An 11-residue permutation of
99–109 (`WGGDGFYAMDY` → `GFDMAWYGYDG`) was also prepared. Both preserve
whole-sequence composition. The thresholds in section 4 were fitted on the
97–109 version and will need refitting if we switch — which they need
anyway at n=3.
