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

## Literature source

Kelley RF, O'Connell MP. Thermodynamic analysis of an antibody functional
epitope. *Biochemistry*. 1993;32(27):6828-6835.
doi:10.1021/bi00078a005.
