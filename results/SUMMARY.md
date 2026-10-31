# HERoes — closing summary (Aida, cofold and scoring)

Breast Cancer Awareness Hackathon, Zurich, 31 October 2026.
Affinity maturation of the trastuzumab scFv against HER2 domain IV.

Submission: `results/final_15.csv`. Scores: `results/scored.csv`.
Full reasoning: `results/rationale.md`. Decisions: `decisions.md`.

---

## The final 15

| name | arm | mutations | dock | median ipSAE | spread |
|---|---|---|---|---|---|
| `consensus_green_03` | consensus_green | T54G;D102S;A106Y;V164I;A169S;Y190Q | 5/5 | 0.9265 | 0.0247 |
| `w01_consensus_green_02` | consensus_green | T32Y;T54G;M107L;V164I;A169S | 5/5 | 0.9246 | 0.0721 |
| `contact_test_07` | contact_test | F104D;T166Y;A167Y;R201G;H226Y;T229Y | 5/5 | 0.9077 | 0.0577 |
| `contact_test_10` | contact_test | R59Y;T166Y;A167Y;F188S;T228S;T229Y | 5/5 | 0.8785 | 0.0674 |
| `contact_test_09` | contact_test | R59Y;F104D;N165S;A167Y;S185G;T229Y | 5/5 | 0.9194 | 0.0873 |
| `w01_germline_revert_02` | germline_revert | N28T;N55G;R59Y;S97A;V164I;A167Y | 5/5 | 0.9369 | 0.0638 |
| `germline_revert_06` | germline_revert | N28T;K30S;T32N;R59Y;N165S;T166S | 5/5 | 0.8839 | 0.0870 |
| `germline_revert_03` | germline_revert | T32N;A72R;S97A;N165S;V168L;H226S | 5/5 | 0.8788 | 0.0924 |
| `aggressive_08` | aggressive | T32Y;T54G;N55G;A72R;A79V;D163G;A167Y;S185G;H226Y | 5/5 | 0.9336 | 0.0789 |
| `w01_aggressive_02` | aggressive | Y57D;A72R;A79V;S97A;D102S;F104D;D163G;A169S;T229Y | 5/5 | 0.8757 | 0.1019 |
| `lit_d102w_pair_04` | literature | T32Y;T54G;N55G;D102W;V164I;Y190Q | 5/5 | 0.8812 | 0.0535 |
| `lit_d102w_pair_02` | literature | N55G;D102W;M107L;A169S;Y190Q | 5/5 | 0.9240 | 0.0737 |
| `esm2_only_07` | esm2_only | G100S;Q162S;A186S;L213V;Q235S | 5/5 | 0.9054 | 0.0574 |
| `esm2_only_06` | esm2_only | D31S;F68A;S120G;M139L;P215A | 5/5 | 0.8948 | 0.0706 |
| `esm2_only_08` | esm2_only | F27Y;T91S;S144A;V150P;P230A | 5/5 | 0.9097 | 0.0897 |

All 15 dock 5/5. 15 unique sequences, all 242 aa, 5–9 mutations each.

## Pool statistics

```
spread    min 0.0247   median 0.0721   max 0.1019
margin to the 0.30 gate                  0.1981
ipSAE     min 0.8757   median 0.9077   max 0.9369
docking   15/15 designs at 5/5 samples
composition  consensus_green 2, contact_test 3, germline_revert 3,
             aggressive 2, literature 2, esm2_only 3
```

## The six arms, all 58 designs

```
arm                PASS/n    rate
esm2_only          12/12     100%
consensus_green     8/10      80%
germline_revert     6/ 8      75%
aggressive          6/ 8      75%
literature          6/10      60%
contact_test        5/10      50%
TOTAL              43/58      74%
```

`contact_test` keeps 3 of its 5 passing designs despite the worst pass
rate. Deliberate: dropping a whole strategy concentrates risk on the arms
that happened to score well. Its weakness is stated, not hidden.

## Key findings, in order

### 1. Broken designs score unstable, not low — and the gate that follows

A broken design does not return a low score. It returns a **bimodal** one.
Each diffusion sample either docks at ~0.87 or collapses below 0.3, with an
empirically empty band between ~0.42 and ~0.75. The same sequence, two seeds:

```
cdrh3_scrambled_strong_seed1   median 0.8857   4/5 docked
cdrh3_scrambled_strong_seed2   median 0.1240   2/5 docked
```

Identical input, medians differing by 0.76.

The collapse is structural, not a scoring artefact: `n0res` goes 242 → 0
and the HER2 interface 142 → 0, while **mean pLDDT holds at 80–85**. The
model stays confidently folded and simply stops touching the target, so
confidence metrics cannot see it.

Therefore the gate is **every one of 5 samples must dock at ipSAE ≥ 0.50,
and spread ≤ 0.30**. The median never gates; it only ranks what passed.
Demonstrated on a real design: `w01_consensus_green_01` has median 0.6553
— mid-pack — with a weakest sample of 0.0180 and 3/5 docked. A median-only
pipeline ships it.

### 2. ESM2-only positions: the absence of a penalty

Our pool required ESM2 and AbLang2 to agree. An arm built from the ~40
positions ESM2 flagged and AbLang2 did not scored **12/12 PASS** across
5–9 mutations, with the tightest median spread of any arm (0.0956).

This supports **“requiring both models to agree was over-cautious.”** It
does **not** support “ESM2-only is better”: Fisher’s exact against
`consensus_green` 8/10 gives **p = 0.195**. The finding is the absence of
a penalty on positions we had been discarding on principle.

The low-perturbation confound was checked, not waved away. Within the
fixed position set, 5-mutation designs had median spread 0.1020 and 7–9
mutation designs 0.0619 — tighter at higher count, pass rate unchanged.
Across all 58, r(mutations, PASS) = −0.03.

### 3. The literature arm gave two opposite answers

The arm is two hypotheses and its blended 60% was uninformative:

```
D102W, a published affinity mutation   4/6 PASS   median spread 0.1792
hotspot-excluded designs               2/4 PASS   median spread 0.4390
```

D102W reproduces on 4 of 6 matched backgrounds. Designing *around* the
published hotspots is the weaker strategy, and its clean failures are
severe (spreads 0.6173 and 0.8654). So “published mutations transfer
poorly” is not supported — if anything the reverse.

### 4. The reverse-order control failed

VL–(GGGGS)₃–VH is measured at 1.83 nM against the parent’s 0.82 nM, a 2.2×
KD difference. Our gap was **0.0236 against a noise floor of 0.0603–0.0791**
— between 0.30× and 0.39× the noise. Checked against ipTM (+0.0129 vs
0.0377) and ipSAE asym (+0.0372 vs 0.1081): same answer under every
definition. **ipSAE cannot resolve a two-fold KD difference on this system,**
and we do not claim our rankings do.

### 5. The construct was wrong, and corrected

The inherited target file was named 541–620 but contained PDB 541–607
(UniProt 563–629), with a 10-residue internal gap at UniProt 603–612
(disordered in 1N8Z) and **Cys626 unpaired**. Replaced by UniProt P04626
**511–652**: 10 intact disulfides, zero free cysteines, and all 19
Fab-contact residues (579–627) with 68/25 aa margin. The C-terminus must
reach 642 because the 630–642 disulfide closes there — which is why the
common “domain IV ends at N629” definition was not used.

## What we do NOT claim

- **No superiority claim for `esm2_only`.** 12/12 vs 8/10 is p = 0.195.
  The claim is the absence of a penalty, nothing more.
- **No conclusion on D102W vs D102S.** 4/6 against 8/11, Fisher **p = 1.000**.
  Two independent routes to the same position, statistically identical. n
  is too small to call a winner and we do not.
- **No claim to resolve two-fold affinity differences.** The reverse-order
  control says we cannot.
- **Thresholds were never moved after seeing the data.** 0.50/0.30 was
  fitted on the controls and kept. Tightening once we could see where the
  designs fell would have been fitting the threshold to the answer we
  wanted. The margin is stated rather than hidden: the closest passing
  design sat 0.09 from a rejection, and the only borderline point was a
  mislabelled design (`lit_hotspot_excluded_01`, which carries D102S and
  so mutates the hotspot its name excludes). It is not in the 15.
- **The gate detects broken predictions, not weak binders.** See below.

## The limitation we found last, and it matters

`ctrl_antidesign` — five mutations chosen as a deliberate negative
(`I34D;G66M;A160W;S191M;Q224P`) — **PASSED the gate 3/3**:

```
ctrl_antidesign_seed1   PASS   5/5   spread 0.0228   median 0.8897
ctrl_antidesign_seed2   PASS   5/5   spread 0.1190   median 0.8759
ctrl_antidesign_seed3   PASS   5/5   spread 0.1413   median 0.8544
```

Its tightest seed, spread **0.0228**, is tighter than every one of our
final 15 and tighter than the parent itself.

The gate measures whether the predicted complex forms **stably**. It does
not measure affinity, and it does not flag a design that folds and docks
confidently while being functionally dead. Stability is a necessary
condition we can test, not a sufficient one. Anyone reading our 15 should
know that our filter removes broken predictions and nothing else.

## Scale

58 designs and 24 control runs cofolded, each 5 diffusion samples.
Median 207 s per cofold on an L40S, 8 containers concurrent with zero
queueing. The 5 samples the gate needs come from one cofold, not three
seeds, so stability filtering cost nothing extra.

ipSAE was computed by two independent implementations agreeing digit for
digit, and Hedi’s was separately validated against the Dunbrack reference
(10,170 synthetic values plus a real cofold, zero deviation) and
mutation-tested (9 injected bugs, 9 caught).

