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
