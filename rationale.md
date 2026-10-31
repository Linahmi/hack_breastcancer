# rationale.md

One section per teammate. Append only; do not edit anyone else's section.

## Aida — cofold (Boltz-2) and the stability gate

### The construct was wrong before it was right

The target file I inherited was named `her2_541_620_from_1n8z.pdb`. It did
not contain 541–620. It contained PDB 541–607 — UniProt 563–629 — with two
defects that matter for cofolding:

- a **10-residue internal gap at UniProt 603–612**, a loop that is
  disordered and unresolved in 1N8Z chain C;
- **Cys626 unpaired**, because its partner Cys634 falls outside the window.
  Four intact disulfides where domain IV has ten.

Replaced with **UniProt P04626 511–652**, 142 aa. I audited five candidate
windows against three criteria — disulfides closed, no free cysteines, and
every Fab contact residue present:

| window | aa | intact SS | cut SS | free Cys | contacts missing |
|---|---|---|---|---|---|
| 563–642 (old YAML) | 80 | 6 | 0 | 0 | 0 |
| 563–629 (old PDB) | 67 | 4 | 1 | 1 | 0 |
| 501–629 (cryo-EM definition) | 129 | 8 | 2 | 2 | 0 |
| 509–604 (DARPin construct) | 96 | 7 | 1 | 1 | 8 |
| **511–652 (chosen)** | **142** | **10** | **0** | **0** | **0** |

The contacts are computed, not cited: 19 HER2 residues within 5 Å of the
Fab heavy atoms in 1N8Z, spanning UniProt 579–627. 511–652 contains all of
them with 68 aa of N-terminal and 25 aa of C-terminal margin.

The C-terminus must reach 642 because the 630–642 disulfide closes there.
The widely cited "domain IV ends at N629" definition cuts two disulfides
and leaves two free cysteines. That is why we did not use it.

### We measured our own noise before we scored anything

Before ranking a single design we ran the parent three times with three
seeds, 5 diffusion samples each — 15 predictions of one sequence — to find
out how much of any score difference is just the model talking to itself.

Parent ipSAE, per-seed medians: 0.9078, 0.9218, 0.8615. Noise, by three
defensible definitions: seed-spread of medians **0.0603**, mean
within-seed sample spread **0.0791**, stdev of per-seed medians **0.0316**.

That number is the floor under every comparison that follows.

### The reverse-order control failed. Reporting it as a negative result.

VL–(GGGGS)₃–VH is measured at 1.83 nM against the parent's 0.82 nM, a
**2.2× difference in KD**. If the pipeline cannot separate those two, it
cannot resolve anything finer.

```
parent_wt       median ipSAE 0.9078
reverse_vl_vh   median ipSAE 0.8842
gap                          0.0236
```

The gap is **0.0236 against a noise floor of 0.0603–0.0791** — between
0.30× and 0.39× the noise, depending on definition. Under every
definition the gap is smaller than the noise. ipSAE cannot resolve a
two-fold KD difference on this system.

We checked whether this was an artefact of the metric. It is not. ipTM
gives +0.0129 against 0.0377 noise (0.34×); ipSAE asymmetric d0res gives
+0.0372 against 0.1081 (0.34×). Three metrics, the same answer, because
all three are dominated by the same seed-to-seed variance in the
underlying prediction. Changing the score does not buy resolution.

So we do not claim our rankings resolve two-fold affinity differences.
They do not.

### The finding: broken designs score unstable, not low

A broken design does not come back with a low score. It comes back
**bimodal**. Each diffusion sample either docks at roughly 0.87 or
collapses below 0.3, with an empirically empty band between about 0.42 and
0.75.

The same sequence, two different seeds:

```
cdrh3_scrambled_strong_seed1   median 0.8857   4/5 docked
cdrh3_scrambled_strong_seed2   median 0.1240   2/5 docked
```

Identical input. Medians differing by 0.76. A single-seed run of that
sequence would have reported either "excellent" or "dead" depending
entirely on the seed.

The collapse is a real structural event, not a scoring artefact:

```
                      ipSAE   n0res   HER2 iface   mean pLDDT
parent, typical       0.910     242          141         85.5
scrambled, docked     0.898     242          141         83.0
scrambled, partial    0.070      43           81         82.7
scrambled, collapsed  0.000       0            0         80.6
```

`n0res` goes 242 → 0 and the HER2 interface goes 142 → 0, while **mean
pLDDT holds between 80 and 85 throughout**. The model is still
confidently folded. It has simply stopped touching the target. pLDDT
cannot see this; confidence metrics are blind to it.

### Therefore: gate on docking rate and spread, never on the median

```
every one of the 5 samples must dock at ipSAE >= 0.50
AND spread across the 5 samples <= 0.30
```

The median does not gate. It only ranks designs that already passed.

Demonstrated on a real design, not a constructed control.
**`w01_consensus_green_01`** has median ipSAE **0.6553** — mid-pack, the
sort of number a ranking would quietly carry forward. Its weakest sample
is **0.0180** and only 3 of 5 docked. A median-only pipeline ships it.

The thresholds were fitted on three controls and we expected them to need
refitting. They transferred:

```
PASS   spread max  0.2293
REJECT spread min  0.3229
```

Nothing occupies the gap between. 34 of 58 designs pass.

### Two negatives, and we measured which was harsher rather than assuming

We built two scrambled-CDR-H3 controls, both preserving whole-sequence
composition, differing only in window:

```
                  docking rate (3 seeds)    min spread
97-109, 13 aa     5/5, 4/5, 3/5             0.348
99-109, 11 aa     4/5, 2/5, 1/5             0.703
```

The 11-residue scramble is decisively the harsher negative. We had
predicted the opposite — more perturbed residues ought to mean a more
broken protein — and we were wrong. Shuffling S97–R98 moved the
permutation into a region that preserves enough paratope to keep docking.
**More mutations did not mean a harsher control.** The weaker negative is
what sets the gate's real margin, so the stronger one anchors it.

### ESM2-only positions: stated carefully

Our design pool required ESM2 and AbLang2 to agree. Lina built an arm from
the ~40 positions **ESM2 flagged and AbLang2 did not** — positions the
pool had never touched.

```
esm2_only          8/8  PASS  (first batch, 5 mutations each)
consensus_green    8/10 PASS  (80%)
```

This supports **"requiring both models to agree was over-cautious"**. It
does **not** support "ESM2-only is better". Fisher's exact on 8/8 vs 8/10
gives **p ≈ 0.47** — entirely consistent with chance. The finding is the
**absence of a penalty**: positions we had been discarding on principle
performed no worse than the ones we kept. That is a reason to widen the
pool, not a claim of superiority.

### The confound, and the check we ran on it

The first 8 `esm2_only` designs all carried exactly 5 mutations, the
minimum, while other arms ranged 5–11. An ESM2-only effect and a
low-perturbation effect are not separable by design in that comparison.
We bounded it instead of hand-waving:

```
r(mutations, PASS)    = -0.03
r(mutations, spread)  = +0.14
PASS mean 6.38 mutations   vs   REJECT mean 6.50
```

Mutation count does not predict passing. The decisive comparison is at
matched perturbation:

```
arm                mean mutations   pass rate
esm2_only                5.00         100%  (first 8)
literature               5.00          25%
consensus_green          5.60          80%
aggressive               9.38          75%
```

`literature` sits at **exactly the same 5.00 mutations** as the first
`esm2_only` batch and passes at **25%**. Same perturbation, a 75-point
gap. Low mutation count is not what confers passing. Lina then built
`esm2_only` designs at 7–9 mutations from the same fixed position set,
turning the confound into a dose-response within one arm. Those four were
still folding when this was written; the dose-response result is not in
this section.

### Literature at 25% is a result, not a gap in our data

The `literature` arm draws on a published alanine scan. It is the
**worst-performing arm we scored**. That is a finding worth stating rather
than an arm to quietly drop: mutations reported in the literature for this
target transferred least well into this scFv format. Format matters, and
published hotspot data does not arrive format-free.

### Cross-validation

The ipSAE numbers were computed twice, independently — Hedi's
implementation and my pipeline reading the same outputs — and they agree
**digit for digit on every seed and every spread**. The scramble
provenance was also settled against the structures rather than against
either of our notes: 15 CIF files and 3 npz agree on which sequence was
folded.

### Throughput

26 cofolds measured at production settings, 384 tokens, L40S: median
**207 s**, p75 **243 s**, range 73–253 s. The variance is MSA-server
latency, not compute.

8 containers ran with **zero queueing** — wall-clock equalled the slowest
single job (253.0 s for 8 jobs, 248.9 s for 6, 142.0 s for 9). 8 is a
measured floor, not a ceiling.

The 5 samples the gate needs come from **one** cofold, not three seeds, so
stability filtering costs nothing extra. At 8-way concurrency that is
roughly 270 designs in a 143-minute window. We used 58.
