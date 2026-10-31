# Design rationale

## 1. Why two models

I used ESM2 650M for general protein constraints and AbLang2 for antibody-specific priors. That choice is measured on this parent, not assumed from their training sets. Masked wild-type recovery was 66.1% for ESM2 and 81.1% for AbLang2, so the antibody-trained model has sharper priors here. Both rank cysteine first at all four structural positions, 22, 96, 158 and 223, with probability above 0.98. CDR entropy also exceeds framework entropy for both; for ESM2 it is 2.83 versus 1.25 bits. These checks do not demonstrate affinity prediction, but they show that the scores reflect antibody sequence constraints.

## 2. The main finding

I asked where both models agree that a non-wild-type residue outranks the parent. The filter fires on 16 of 21 HER2-contact residues, or 76%, but only 7 of 166 unlisted positions, or 4%: a nineteen-fold enrichment where the antibody was engineered. Trastuzumab was a murine antibody humanised onto a human framework and then affinity matured, so its contact residues are deliberately atypical. In an already optimised therapeutic, model surprise can mean either suboptimal or engineered for binding. Sequence alone cannot separate those explanations. The strongest model signal was therefore also unsafe to trust without testing.

## 3. How we handled the ambiguity

We could not resolve that ambiguity today, so I made it part of the experiment. The `contact_test` arm deliberately mutates surprising contact residues rather than excluding them or calling them improvements. Measurements in three weeks can show whether the models found sequence liabilities or proposed undoing affinity maturation. We also registered seven predictions in `results/predictions.md` before any design was scored, in commits `e771c4c` at 14:22 and `8b22bc1` at 14:35. This does not make them right; it prevents us rewriting the reasoning after seeing the answer.

## 4. The pool

The pool contains 48 designs across six arms and covers 70 mutated positions. `consensus_green` asks whether lower-risk CDR-adjacent consensus mutations help. `contact_test` asks whether flagged contacts are suboptimal or engineered. `germline_revert` tests whether AbLang2's germline pull helps or removes binding adaptations. `aggressive` compares eight to ten mutations with five or six and tests the cost of epistasis. `literature` asks whether published effects reproduce here. `esm2_only` tests whether two-model agreement was justified or merely conservative. The `ctrl_antidesign` combines mutations both models strongly dislike; if it binds well, their scores carried little signal. Position 102 is a direct comparison. With no structural input, the models selected D102S in seven consensus designs. Gerstner's phage-display scan identified D102W. They converged on the position but chose different residues, so measurement decides between them.

## 5. What we do not know

The pool is exploratory, not comparative. Two or three effective designs per arm cannot distinguish a 50% hit rate from 25%. In the aggressive arm, mutation count and position choice are confounded: eight to ten mutations from only 30 consensus positions forced it outside that set. If it performs better, we cannot tell which factor caused it. Neither model has seen HER2. They score residue typicality, which tracks foldability and stability more directly than binding affinity. ipSAE also cannot resolve the known 2.2-fold affinity difference between parent and reverse-order scFv: +0.0129 against 0.0377 seed noise. It does separate a broken design by about -0.28. We therefore use ipSAE as a failure filter, not a ranking of plausible binders.
