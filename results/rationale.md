# Design rationale

## 1. Why two models

ESM2 650M for general protein constraints, AbLang2 for antibody-specific priors. The choice is measured on this parent, not assumed from training sets. Masked wild-type recovery was 66.1% for ESM2 (160 of 242) and 81.1% for AbLang2 (184 of 227 scored positions), so the antibody-trained model has the sharper priors here. Both rank cysteine first at all four structural positions, 22, 96, 158 and 223, at probability above 0.98, so both sets of priors are calibrated on structural constraint. CDR entropy exceeds framework entropy for both models: 2.83 against 1.25 bits for ESM2, 2.51 against 1.49 for AbLang2. That is the expected hypervariability signature. Spearman correlation between the two score columns is 0.65 across 4,313 non-wild-type substitutions. The models overlap, but enough of the signal is model-specific that requiring agreement is a real filter: 82 positions pass ESM2, 43 pass AbLang2, 35 pass both.

## 2. The finding

Requiring both models to agree that some non-wild-type residue outranks the wild type, the filter fires on 15 of 20 HER2-contact residues, 75%, against 8 of 167 unlisted positions, 4.8%. That is roughly a sixteen-fold enrichment, Fisher's exact p around 1e-12, concentrated exactly where the antibody was engineered. Trastuzumab is a murine antibody humanised onto a human framework and then affinity matured, so its contact residues are deliberately atypical. For an already-optimised therapeutic, model surprise is ambiguous: the residue may be suboptimal, or engineered on purpose, and sequence alone cannot separate the two. The strongest signal our scans produced is also the one we trust least.

## 3. How we handled it

We could not resolve the ambiguity today, so we built it into the pool as a test. The `contact_test` arm deliberately mutates those residues; the measurements in three weeks answer the question. Seven predictions were registered in `results/predictions.md` before any design was scored, in commits `e771c4c` at 14:22 and `8b22bc1` at 14:35. The early returns are consistent with the central prediction: `contact_test` is passing at 50% against 80% for `consensus_green`. At ten designs per arm that is directionally consistent and underpowered. Fisher's exact on 5 of 10 against 8 of 10 gives p = 0.35, which is not evidence.

## 4. The pool and the controls

58 designs, six arms, 70 unique mutated positions. Each arm answers a different question: `consensus_green`, whether safe CDR-adjacent consensus mutations help; `contact_test`, engineered or suboptimal; `germline_revert`, whether AbLang2's germline pull helps or hurts; `aggressive`, whether mutation count helps or epistasis bites; `literature`, whether published mutations reproduce here; `esm2_only`, whether requiring consensus was justified or merely conservative.

Two comparisons were designed rather than assembled:

- Four matched D102W pairs, each differing from a `consensus_green` design at position 102 alone, with every background docking 5 of 5. Our models and the published record converged on 102 independently. Both models put serine top there with no structural input, at 0.35 for ESM2 and 0.37 for AbLang2 against 0.04 for wild-type aspartate, while Gerstner's phage-displayed Fab library selection found 102 the most sequence-plastic site and D102W about three-fold affinity-improving. Same position, different residue, identical backgrounds.
- `esm2_only` spans 5 to 9 mutations across a fixed 40-position set, so mutation count varies with position choice held constant.

For controls, `ctrl_antidesign` is built from the mutations both models most dislike: if it binds as well as our designs, our scores carried no signal. ipSAE does not resolve the known 2.2-fold affinity difference between parent and reverse-order scFv, 0.82 against 1.83 nM measured, separating them by +0.0129 against 0.0377 of seed noise. It does cleanly separate a broken design at about -0.28. We therefore use it as a filter, not a ranking. Broken designs also score unstably rather than low: one had a median of 0.6553 but docked only 3 times in 5, and a median filter would have passed it. We filter on docking rate and seed variance across five seeds instead.

## 5. What we do not know

- The pool is exploratory, not comparative. Two to three working designs per arm cannot separate a 50% hit rate from 25%.
- In the `aggressive` arm, mutation count and position choice are confounded. The arm needed nine mutations per design, and only 20 non-contact positions passed the consensus filter, so it ranged outside them across 28 positions, including the published hot spots at 50, 99 and 226. `esm2_only` is the only arm that avoids this.
- Neither model has ever seen HER2. They score how typical a residue is, which tracks stability and foldability more than binding affinity.
- The second-shell region is unexplored. Only 8 of 167 unlisted positions passed the consensus filter, so we had no design space there. That emptiness is a consequence of our own filter choice, not a fact about the molecule.
