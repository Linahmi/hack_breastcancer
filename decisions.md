# HERoes — decisions log

One block per teammate. Append only; do not edit anyone else's entry.

## Ghita - structure constraints (1N8Z, 4.5 A)
- RED (conserve, direct HER2 contact), 20 positions: 33,50,52,57,58,59,99,103,104,105,163,165,166,167,185,188,226,227,228,229
- GREEN (CDR, mutable): 31,32,34,35,51,53,54,55,56,60,61,62,63,64,65,66,100,101,102,106,107,108,109,159,160,161,162,164,168,169,186,187,189,190,191,224,225,230,231,232
- tool: structure/contacts.py (numpy only) on PDB 1N8Z; chains classified by sequence; numbering asserted residue by residue vs parent.
- reason: preserve the proven paratope (0.82 nM), let PLM scans mutate only non-contact CDR positions (lower risk to affinity).
- correction appliquee: les eaux cristallographiques de la chaine C etaient comptees comme atomes HER2. Corrige -> RED passe de 21 a 20 (la position 201 n avait que ce contact medie par une eau, 5.45 A au plus proche atome proteique) ; GREEN reste 40, inchange, la whitelist de design n a jamais ete affectee. contacts.py --include-waters reproduit l ancien CSV ; commit f01f504 conserve la version pre-correction.
- litterature: 4 hotspots d alanine scan (Kelley 1993) tombent tous sur des positions RED ; la mutation ameliorante publiee D->W (Gerstner 2002, ~3x) tombe sur la position GREEN 102.

## Aida - cofold, gate thresholds, and pool composition

- **Construct.** `inputs/her2_541_620_from_1n8z.pdb` rejected: it is named
  541-620 but contains PDB 541-607 (UniProt 563-629), with a 10-residue
  internal gap at UniProt 603-612 (disordered in 1N8Z chain C) and Cys626
  unpaired because Cys634 falls outside the window. Replaced by
  `src/aida/her2_domainIV.fasta`, UniProt P04626 511-652: 10 intact
  disulfides, zero free cysteines, all 19 Fab-contact residues
  (UniProt 579-627) present with 68/25 aa margin. The C-terminus must
  reach 642 because the 630-642 disulfide closes there, which is why the
  common "domain IV ends at N629" definition was not used.

- **Gate thresholds: KEPT at 0.50 / 0.30. No refit.** Every sample must
  dock at ipSAE >= 0.50 and sample-spread must be <= 0.30. Hedi's
  reasoning, recorded as given: tightening the threshold now that we can
  see where the designs fall would be fitting the threshold to the answer
  we want. The band is empty (closest PASS 0.2608, closest REJECT
  0.3530) and the criterion still separates - it is simply tighter than
  it was on the controls. We state the margin explicitly rather than
  hide it: **the closest passing design sits 0.09 from a rejection.**

- **The threshold question is settled, not worked around.** The only
  borderline point was `lit_hotspot_excluded_01` at spread 0.2608, and
  that design is mislabelled - it carries D102S and so mutates the
  published hotspot its name says it excludes. It is excluded from the
  final pool. No threshold was moved to accommodate anything.

- **The median never gates.** It ranks designs that already passed.
  Demonstrated on a real design: `w01_consensus_green_01` has median
  ipSAE 0.6553 - mid-pack - with a weakest sample of 0.0180 and only 3/5
  docked. A median-only pipeline ships it.

- **`contact_test` keeps all 3 slots despite being the weakest arm.**
  50% pass rate (5/10), and its median design is rejected at 0.3479
  spread. This is a deliberate diversity choice, not an oversight. We
  take its three tightest designs - 0.0577, 0.0674, 0.0873 - because
  dropping a whole strategy from the pool concentrates risk on the arms
  that happened to score well. The arm's weakness is stated, not hidden.

- **Literature arm split by hypothesis.** The blended 60% figure mixed
  two questions and must not be pitched as "published mutations transfer
  poorly". D102W (published mutation): 4/6 PASS, median spread 0.1792.
  Hotspot-excluded (consensus designs with a constraint): 2/4 PASS,
  median spread 0.4390, and effectively n=3 once the mislabelled design
  is removed. Both literature slots test D102W, one coherent hypothesis.
  Against D102S at the same position: 4/6 vs 8/11, Fisher p = 1.000 -
  indistinguishable, n too small to call a winner.

- **esm2_only: absence of a penalty, not superiority.** 12/12 PASS
  across 5-9 mutations within a fixed position set. Fisher vs
  consensus_green 8/10 gives p = 0.195 - not significant. The
  low-perturbation confound is ruled out by the dose-response: 5-mutation
  designs have median spread 0.1020, 7-9 mutation designs 0.0619
  (tighter, not looser), so position choice rather than perturbation size
  is what mattered.

- **Reverse-order control FAILED, reported as a negative result.** Gap
  parent vs reverse 0.0236 against a noise floor of 0.0603-0.0791
  depending on definition - between 0.30x and 0.39x the noise. ipSAE
  cannot resolve the known 2.2x KD difference (0.82 vs 1.83 nM). Checked
  against ipTM (+0.0129 vs 0.0377) and ipSAE asym (+0.0372 vs 0.1081):
  same answer. We do not claim our rankings resolve two-fold affinity.
