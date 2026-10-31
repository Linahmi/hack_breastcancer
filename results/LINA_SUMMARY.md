# Lina — reading guide

Everything of mine Ghita needs for the report and the pitch, in reading order.
Compiled 31 October 2026. Numbers below were recomputed from the raw matrices
for this document; see **Verification** at the end for the three places where
the prose and the recomputation disagree.

## Read in this order

1. **`src/lina/rationale.md`** — the argument in prose, no numbers. Why model
   surprise is not evidence of improvable binding in a humanised,
   affinity-matured therapeutic, and why structural context has to come before
   the scan is actionable. This is the reasoning the rest of the files quantify.

2. **`results/rationale.md`** — the same argument with the measurements, in five
   sections: why two models, the finding, how we handled it, the pool and the
   controls, what we do not know. If Ghita reads one file, this is it.

3. **`results/predictions.md`** — seven predictions registered *before* any
   design was scored, locked at **2026-10-31 14:21:44 CET**, in commits
   `e771c4c` (14:22) and `8b22bc1` (14:35). Predictions 6 and 7 were added
   after the consensus-position analysis but still before any design scores
   existed. This is what makes the three-week readout a test rather than a
   story.

4. **`results/candidates.csv`** — the 58 designs with per-design rationale and
   arm. The metadata source of record; the wave CSVs carry only
   `name,sequence` because that is the GPU input contract.

5. **`results/controls.csv`** and **`results/geometry_audit.csv`** — the four
   controls, and per-design CA-distance covariates. The audit is a covariate
   table, not a pass/fail verdict.

6. **`results/LINA_HANDOFF.md`** — which wave files are ready and what Aida and
   Hedi each run. Operational; skip for the report, needed for the pitch demo.

7. **`src/lina/approaches.md`** and **`src/lina/dependencies.md`** — the
   decision reasoning, written before the work was done. Options considered
   with their trade-offs and failure modes, and the integration contract with
   Ghita's `positions.csv`. These are pre-decision snapshots and are superseded
   in places (see below); they are in the report because they show the choices
   were made deliberately rather than by default.

8. **`src/lina/scan_esm2.py`**, **`scan_ablang2.py`**, **`diagnostics.py`**, the
   `build_*.py` scripts, and the `*_masked_marginals.csv` / `*_logprobs.npy`
   matrices — the scan itself and everything reproducible from it.

## The two models, and what validates them

ESM2 650M for general protein constraints, AbLang2 for antibody-specific
priors. Both were scored by masked marginals: mask the parent residue, compare
the model's probability for each substitution against the wild-type residue.
The choice of models is justified on this parent, not assumed from training
sets:

- **Wild-type recovery.** On the 223 positions both models score, ESM2 ranks
  the parent residue top at 64.6% and AbLang2 at 80.7%. The antibody-trained
  model has the sharper priors here.
- **Cysteine.** Both models rank cysteine first at all four structural
  positions — 22, 96, 158, 223 — at p ≥ 0.98 (ESM2 0.992–0.999, AbLang2
  0.982–0.995). Both sets of priors are calibrated on structural constraint.
- **Entropy.** CDR entropy exceeds framework entropy for both models: 2.83
  against 1.25 bits for ESM2, 2.51 against 1.49 for AbLang2, over 60 CDR and
  167 framework positions. That is the expected hypervariability signature.
- **Independence.** Spearman between the two score columns is 0.65 across 4,237
  non-wild-type substitutions. The models overlap, but enough signal is
  model-specific that requiring agreement is a real filter: 79 positions pass
  ESM2, 43 pass AbLang2, 35 pass both.

## The finding: 75% against 4.8%

Requiring both models to agree that some non-wild-type residue outranks the
wild type, the filter fires on **15 of 20 HER2-contact residues (75%)** against
**8 of 167 unlisted positions (4.8%)** — roughly a sixteen-fold enrichment at
Fisher's exact p ≈ 1.2e-12, concentrated exactly where the antibody was
engineered.

The contact set is Ghita's, from `structure/red_green.csv` at commit `f0ee860`,
*"exclude crystallographic waters from the contact definition"*. **That commit
is on `origin/ghita-structure` and is not yet merged**; the copy currently on
`main` is the earlier `f01f504`, which counts position 201 as a contact through
a water-mediated contact and yields 16 of 21 against 7 of 166 instead. The
figures above are the water-excluded ones, because excluding waters is correct.
Anyone recomputing this from `main` before that merge will get the other
number — see **Verification**.

What it means: trastuzumab's contact residues are deliberately atypical,
murine loops grafted onto a human framework and then affinity matured. So model
surprise at those positions is ambiguous — the residue may be suboptimal, or it
may be engineered on purpose, and sequence alone cannot separate the two. The
strongest signal the scans produced is also the one we trust least. This is the
central result and the reason the pool is built as a test rather than a
shortlist.

## The six arms

58 designs, 70 unique mutated positions. Each arm answers one question:

| Arm | n | Question |
|---|---|---|
| `consensus_green` | 10 | Do safe CDR-adjacent consensus mutations help? |
| `contact_test` | 10 | Engineered on purpose, or genuinely suboptimal? |
| `germline_revert` | 8 | Does AbLang2's germline pull help or hurt? |
| `aggressive` | 8 | Does mutation count help, or does epistasis bite? |
| `literature` | 10 | Do published mutations reproduce here? |
| `esm2_only` | 12 | Was requiring consensus justified, or merely conservative? |

The literature arm splits by hypothesis: 4 matched D102W pairs, 4 plain
literature designs, 2 with the published alanine-scan hot spots excluded.

Two comparisons were designed rather than assembled. The **D102W pairs** differ
from a `consensus_green` design at position 102 alone, so the backgrounds are
identical — models and the published record converged on position 102
independently and disagree on the residue (both models put serine top with no
structural input, 0.35 ESM2 and 0.37 AbLang2 against 0.04 for wild-type
aspartate; Gerstner's phage-display selection found D102W about three-fold
improving). One pair later failed the stability gate, leaving three usable
head-to-heads. **`esm2_only`** spans 5 to 9 mutations across a fixed 40-position
set, so mutation count varies with position choice held constant — the only arm
that is not confounded that way.

## The controls, and why ipSAE is a filter and not a ranking

Four controls: parent, reverse-order scFv, scrambled CDR-H3 with composition
preserved, and `ctrl_antidesign` built from the mutations both models most
dislike. The anti-design is the signal check — if it binds as well as our
designs, our scores carried nothing.

**ipSAE magnitude does not discriminate.** When the scrambled-CDR-H3 control
docks, it scores 0.873 against the parent's 0.905. It fails *bimodally*, not
low: 7 of 15 samples docked against 15 of 15 for parent and reverse-order, and
the collapses are physical rather than a scoring artefact — interface residues
drop from 242 to 0 while mean pLDDT holds at 82–85. Both chains fold; only the
interface fails. ipSAE also does not resolve the known 2.2-fold parent against
reverse-order affinity difference (0.82 against 1.83 nM measured): +0.0236
against a parent seed spread of 0.0603, about 0.5 kcal/mol, below what a
structure predictor resolves.

So the discriminating quantities are **docking rate and intra-run spread**, not
score magnitude. The gate is all five diffusion samples docking above 0.5 with
spread at most 0.30. A median-only filter would have passed
`w01_consensus_green_01`, which docked 3 times in 5 at a median of 0.6553.
The thresholds are provisional — calibrated on controls that are extreme by
construction, and the scramble sits only 0.048 above the spread gate.

## Honest limitations

- **The pool is exploratory, not comparative.** Two to three working designs per
  arm cannot separate a 50% hit rate from 25%.
- **Neither model has ever seen HER2.** They score how typical a residue is,
  which tracks stability and foldability more than binding affinity. Nothing
  here is an affinity prediction.
- **The `aggressive` arm is confounded.** Mutation count and position choice
  move together: it needed nine mutations per design and only 20 non-contact
  positions passed the consensus filter, so it ranged across 28 positions
  including the hot spots at 50, 99 and 226.
- **The second shell is unexplored.** Only 8 of 167 unlisted positions passed
  the consensus filter, so there was no design space there. That emptiness is a
  consequence of our filter choice, not a fact about the molecule — and
  Ghita's three-way file shows it is not even true: see the gradient below.
- **The early returns are not evidence.** `contact_test` passing at 50% against
  80% for `consensus_green` is directionally consistent with prediction 5, but
  Fisher's exact on 5 of 10 against 8 of 10 gives p = 0.35.
- **Those early returns are not in this repository.** There is no
  `results/scored.csv` or any other design-score file here; the 50%/80% figures
  came from Aida directly, around 14:45. They should be sourced to her run
  output before they go in the report.

## Verification

Every figure in this file and in `results/rationale.md` was recomputed from
`*_masked_marginals.csv`, `*_logprobs.npy`, `entropy.csv` and
`structure/red_green.csv` on 31 October. The entropy figures, the cysteine
probabilities, the Spearman coefficient and the 35 positions passing both
models reproduced unchanged, as did the enrichment. Three figures did not, and
`results/rationale.md` has been corrected to match this file. The superseded
values are recorded here because they appear in the branch history and in
anything drafted before 31 October 15:30:

- **Recovery, previously 66.1% and 81.1%.** Those numbers gave each model its
  own denominator (242 and 227) and counted the excluded cysteine and linker
  positions as recovered, which made the two models not directly comparable.
  On the 223 positions both models actually score it is 64.6% and 80.7%. The
  conclusion — AbLang2 has the sharper priors — is unchanged.
- **Positions passing ESM2, previously 82, now 79.** A three-position
  difference in how the pass set was counted; the recomputation excludes the
  wild-type self-row, which scores exactly zero and is not a passing
  substitution. 43 for AbLang2 and 35 for both were unaffected.
### The structure file, and one retracted correction

`structure/red_green.csv` is Ghita's and has two revisions. `f01f504` (13:27)
is what sits on `main`; `f0ee860` (13:58), on the unmerged
`origin/ghita-structure`, excludes water-mediated contacts and drops position
201 from the contact set. **Hers is the later and correct one, and this file
uses it.**

Between 15:24 and 15:30 on 31 October the enrichment in `rationale.md` was
changed to 16 of 21 against 7 of 166 — computed from `main`'s copy without
checking for a newer revision on her branch. That was wrong and has been
reverted. The original 15 of 20 against 8 of 167 was correct all along, and so
was the 167 figure, which is the unlisted count under her classification rather
than a framework count. Recorded here because the wrong numbers are in the
history between `55a4bde` and the revert.

**Anyone recomputing the enrichment must take `red_green.csv` from
`origin/ghita-structure`, not from `main`,** until that branch is merged.

### The gradient, once her three-way file lands

`results/positions.csv` on `main` is an empty header row. The full 242-row
three-way file is on `origin/ghita-structure`, and it gives a strictly stronger
result than the binary split — the pass rate falls monotonically with distance
from the interface:

| Category | n | Passing both models |
|---|---|---|
| `contact` | 20 | 15 (75.0%) |
| `second_shell` | 21 | 8 (38.1%) |
| `free` | 186 | 12 (6.5%) |

This is worth building the pitch around instead of contact-against-everything.
It also retires the second-shell limitation above: there are 8 passing
positions there, so the design space exists and the red/green split was hiding
it. `src/lina/dependencies.md` called this risk in advance — the second-shell
arm needs the three-way `category`, not a red/green split — which is why that
arm does not exist in the pool.
