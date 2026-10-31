# Pre-score predictions

**Prediction lock timestamp:** 2026-10-31 14:21:44 CET (UTC+01:00, Europe/Berlin)

These predictions were recorded before any design scores existed. At the lock
timestamp, no `results/scored.csv` or other design-score output was present in
the repository.

1. **`aggressive_04` is the design most likely to fail.** It carries W99D, a
   Trp-to-Asp substitution at a published hot spot, together with R50Y in a
   nine-mutation design. Kelley and O'Connell identified W99 as a hot spot by
   alanine scanning, with dDG greater than 3 kcal/mol for W99A. Replacing Trp
   with a charged Asp is therefore the strongest individual failure signal in
   the current pool.

2. **The four R50Y designs are at elevated risk, but are not predicted to fail
   outright.** The affected designs are `w01_aggressive_01`,
   `w01_contact_test_02`, `aggressive_04`, and `contact_test_05`. The published
   scan measured R50A, which removes the side chain entirely, whereas R50Y
   retains a large aromatic side chain. If these designs survive, the result
   will distinguish a requirement for side-chain bulk from a requirement for
   positive charge at this hot spot.

3. **`ctrl_antidesign` will show high seed-to-seed variance rather than a
   uniformly low score.** This follows Aida's bimodal failure model: individual
   seeds may alternate between plausible and failed docking modes, producing a
   broad or bimodal score distribution.

4. **Designs with minimum pairwise mutation distance under 5 A will have a
   lower docking rate than designs over 10 A.** The pre-score groups contain 17
   designs under 5 A and 4 designs over 10 A. Distance is treated as a
   covariate, not as a retrospective filter.

5. **If `contact_test` fails disproportionately relative to
   `consensus_green`, the language models were flagging engineered binding
   residues rather than genuine sequence suboptimality.** That result would
   support the central hypothesis of this approach: model surprise becomes
   useful only after structural context distinguishes intentional interface
   engineering from sequence liability.

Predictions 6 and 7 were added after the consensus-arm position analysis but
still before any design was scored.

6. **If `consensus_green` failures cluster on the seven designs carrying
   D102S, the position is the problem. If failures spread evenly across all
   ten designs, the first failure was stochastic.** The first failure was
   `w01_consensus_green_01`: docking rate 3/5, minimum seed 0.0180, median
   0.6553. Aida confirms this with all ten consensus designs around 14:45.

7. **D102S versus D102W is a direct head-to-head at the same position.** D102S
   appears in seven designs and was chosen independently by both language
   models. D102W appears in two designs and comes from the published Gerstner
   alanine-scan study. Models and literature converged on the position and
   disagree on the residue; the measurements settle which was right.

## Literature source

Kelley RF, O'Connell MP. Thermodynamic analysis of an antibody functional
epitope. *Biochemistry*. 1993;32(27):6828-6835.
doi:10.1021/bi00078a005.
