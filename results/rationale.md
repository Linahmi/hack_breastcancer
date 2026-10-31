# Design rationale

## 1. Why two models

ESM2 650M for general protein constraints, AbLang2 for antibody-specific priors. The choice is measured on this parent, not assumed from training sets. Masked wild-type recovery, over the 223 positions both models score, was 64.6% for ESM2 and 80.7% for AbLang2, so the antibody-trained model has the sharper priors here. Both rank cysteine first at all four structural positions, 22, 96, 158 and 223, above p = 0.98, so both priors are calibrated on structural constraint. CDR entropy exceeds framework entropy for both models, 2.83 against 1.25 bits for ESM2 and 2.51 against 1.49 for AbLang2, which is the expected hypervariability signature. Spearman correlation between the two scores is 0.65 across 4,237 non-wild-type substitutions. Enough of the signal is model-specific that requiring agreement is a real filter: 79 positions pass ESM2, 43 pass AbLang2, 35 pass both.

## 2. The finding

Requiring both models to agree that some non-wild-type residue outranks the wild type, the filter fires on 15 of 20 HER2-contact residues, 75%, against 8 of 167 unlisted positions, 4.8%. That is roughly a sixteen-fold enrichment, Fisher's exact p around 1e-12, concentrated exactly where the antibody was engineered. The contact set is Ghita's, from `structure/red_green.csv` at commit `f0ee860`, which excludes water-mediated contacts. An earlier revision of that file counted position 201 as a contact and gives 16 of 21 against 7 of 166 instead; the figures here are the water-excluded ones. Trastuzumab is a murine antibody humanised onto a human framework and then affinity matured, so its contact residues are deliberately atypical. For an already-optimised therapeutic, model surprise is ambiguous: the residue may be suboptimal, or engineered on purpose, and sequence alone cannot separate the two. The strongest signal our scans produced is also the one we trust least.

## 3. How we handled it

We could not resolve the ambiguity today, so we built it into the pool as a test. The `contact_test` arm deliberately mutates those residues; the measurements in three weeks answer the question. Seven predictions were registered in `results/predictions.md` before any design was scored, in commits `e771c4c` at 14:22 and `8b22bc1` at 14:35. The early returns are consistent with the central one: `contact_test` is passing at 50% against 80% for `consensus_green`. At ten designs per arm that is directionally consistent and underpowered. Fisher's exact on 5 of 10 against 8 of 10 gives p = 0.35, which is not evidence.

## 4. The pool and the controls

58 designs, six arms, 70 unique mutated positions. Each arm answers a different question: `consensus_green`, whether safe CDR-adjacent consensus mutations help; `contact_test`, engineered or suboptimal; `germline_revert`, whether AbLang2's germline pull helps or hurts; `aggressive`, whether mutation count helps or epistasis bites; `literature`, whether published mutations reproduce here; `esm2_only`, whether requiring consensus was justified or merely conservative.

Two comparisons were designed rather than assembled:

- Four matched D102W pairs, each differing from a `consensus_green` design at position 102 alone, every background docking 5 of 5. One pair then failed the stability gate itself, so three usable head-to-heads remain. Our models and the published record converged on 102 independently: both models put serine top with no structural input, 0.35 for ESM2 and 0.37 for AbLang2 against 0.04 for wild-type aspartate, while Gerstner's phage-displayed Fab library selection found 102 the most sequence-plastic site and D102W about three-fold affinity-improving.
- `esm2_only` spans 5 to 9 mutations across a fixed 40-position set, so mutation count varies with position choice held constant.

The negative control is a scrambled CDR-H3 with composition preserved, run at settings identical to the parent. ipSAE magnitude does not separate it: when the scramble docks it scores 0.873 against the parent's 0.905. It fails bimodally rather than low. Seven of its fifteen samples docked against fifteen of fifteen for parent and reverse-order, and the collapses are physical rather than a scoring artefact, with interface residues falling from 242 to 0 while mean pLDDT holds at 82 to 85. Both chains fold and only the interface fails. What separates the control is docking rate and intra-run spread, not magnitude, so we use ipSAE as a filter and not a ranking. We require all five diffusion samples to dock above 0.5 with spread at most 0.30. A median-only filter would have kept `w01_consensus_green_01`, which docked 3 times in 5 at a median of 0.6553. ipSAE also does not resolve the 2.2-fold parent against reverse-order difference, 0.82 against 1.83 nM measured: +0.0236 against a parent seed spread of 0.0603. That is roughly 0.5 kcal/mol, below what a structure predictor resolves, so we report it rather than interpret it. The thresholds remain provisional, calibrated on controls extreme by construction, and the scramble now sits only 0.048 above the spread gate.

## 5. What we do not know

- The pool is exploratory, not comparative. Two to three working designs per arm cannot separate a 50% hit rate from 25%.
- In the `aggressive` arm, mutation count and position choice are confounded. The arm needed nine mutations per design, and only 20 non-contact positions passed the consensus filter, so it ranged outside them across 28 positions, including hot spots 50, 99 and 226. `esm2_only` is the only arm that avoids this.
- Neither model has ever seen HER2. They score how typical a residue is, which tracks stability and foldability more than binding affinity.
- The second-shell region is unexplored. Only 8 of 167 unlisted positions passed the consensus filter, so we had no design space there. That emptiness is a consequence of our own filter choice, not a fact about the molecule.

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
PASS   spread max  0.2608
REJECT spread min  0.3229
```

Nothing occupies the gap between, but the margin is thin: **0.0621**.
43 of 58 designs pass. One passing design, `lit_hotspot_excluded_01` at
0.2608, sits 0.06 from a rejection — close enough that the threshold is
now adjudicating rather than comfortably separating, and it should be
refit against this distribution rather than the controls it came from.

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
esm2_only         12/12 PASS  (100%)
consensus_green    8/10 PASS  (80%)
```

This supports **"requiring both models to agree was over-cautious"**. It
does **not** support "ESM2-only is better". Fisher's exact on 12/12 vs
8/10 gives **p = 0.195** — consistent with chance. The finding is the
**absence of a penalty**: positions we had been discarding on principle
performed no worse than the ones we kept. That is a reason to widen the
pool, not a claim of superiority.

The arm is also the most *stable* of the six — median sample spread
**0.0956** against 0.1091 for `consensus_green` and 0.3479 for
`contact_test`. Stability, not score, is where the difference shows.

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

At the time of that check `literature` sat at **exactly the same 5.00
mutations** as the first `esm2_only` batch and passed at **25%**. Same
perturbation, a 75-point gap. Low mutation count is not what confers
passing. (With six further literature designs the arm finished at 60%,
mean 5.30 mutations — still below `esm2_only` at matched perturbation.) Lina then built
`esm2_only` designs at 7–9 mutations from the same fixed position set,
turning the confound into a dose-response within one arm.

**The dose-response settles it.** Within that fixed position set, varying
mutation count from 5 to 9 changed nothing about pass rate, and the
higher-count designs were *tighter*, not looser:

```
n_mut    n   PASS   median spread
    5    8      8          0.1020
    7    1      1          0.0581
    8    1      1          0.1160
    9    2      2          0.0619

5 mutations    median spread 0.1020
7-9 mutations  median spread 0.0619   <- tighter
```

Position choice, not perturbation size, is what mattered. The
low-perturbation explanation for the `esm2_only` result is ruled out.

### The literature arm is two hypotheses, and they must not be blended

Our first reading — "published mutations transfer poorly" — was wrong,
because the arm is heterogeneous. Of the 10 designs labelled `literature`,
only 6 carry an actual published mutation; the other 4 are
hotspot-excluded designs, which are consensus designs with a constraint.
Those test different questions and have to be reported separately.

**Does a published affinity-improving mutation reproduce? D102W.**

```
D102W            4/6 PASS (67%)   median spread 0.1792
  lit_d102w_pair_04   PASS  5/5  spread 0.0535  median 0.8812
  lit_d102w_pair_02   PASS  5/5  spread 0.0737  median 0.9240
  lit_d102w_pair_03   PASS  5/5  spread 0.1727  median 0.7751
  literature_02       PASS  5/5  spread 0.1857  median 0.8534
  lit_d102w_pair_01   REJECT 5/5 spread 0.3530  median 0.7387
  literature_01       REJECT 3/5 spread 0.4773  median 0.8285
```

Yes, mostly. D102W reproduces on 4 of 6 matched backgrounds. Against
D102S — the residue the PLMs chose at the same position — it is
**indistinguishable**: D102W 4/6 versus D102S 8/11, Fisher two-sided
**p = 1.000**. Two independent routes to the same position, picked by a
published alanine scan and by two protein language models, give
statistically identical results. We do not claim a winner; n is too small.

**Do hotspot-excluded designs work? A different question.**

```
hotspot_excluded 2/4 PASS (50%)   median spread 0.4390
  lit_hotspot_excluded_02  PASS   5/5  spread 0.0345  median 0.9176
  lit_hotspot_excluded_01  PASS   5/5  spread 0.2608  median 0.8490
  literature_04            REJECT 4/5  spread 0.6173  median 0.5776
  literature_03            REJECT 2/5  spread 0.8654  median 0.4955
```

Half pass, and the failures are severe — median spreads of 0.62 and 0.87.
Excluding the published hotspot positions and designing around them is the
weaker strategy here, which is the opposite conclusion from the D102W
result and is exactly why the blended 60% figure was uninformative.

One construction caveat we should fix rather than paper over:
`lit_hotspot_excluded_01` carries **D102S**, so it mutates the very
hotspot its label says it excludes. It is also the design that tripped our
spread watchpoint at 0.2608. Treat it as mislabelled, not as evidence
either way.

### Cross-validation

The ipSAE numbers were computed twice, independently — Hedi's
implementation and my pipeline reading the same outputs — and they agree
**digit for digit on every seed and every spread**. The scramble
provenance was also settled against the structures rather than against
either of our notes: 15 CIF files and 3 npz agree on which sequence was
folded.

### Throughput

58 designs plus 21 control runs were cofolded. Timing from 26 measured
cofolds at production settings, 384 tokens, L40S: median
**207 s**, p75 **243 s**, range 73–253 s. The variance is MSA-server
latency, not compute.

8 containers ran with **zero queueing** — wall-clock equalled the slowest
single job (253.0 s for 8 jobs, 248.9 s for 6, 142.0 s for 9). 8 is a
measured floor, not a ceiling.

The 5 samples the gate needs come from **one** cofold, not three seeds, so
stability filtering costs nothing extra. At 8-way concurrency that is
roughly 270 designs in a 143-minute window. We used 58.
