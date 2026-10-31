# Submission — results/final_15.csv

Parent: 242 aa scFv, VH-(GGGGS)3-VL, KD 0.82 nM by BLI.
Target: HER2 domain IV, UniProt P04626 residues 511-652 (142 aa).
Each design was cofolded against that target with Boltz-2 (v2.2.1),
5 diffusion samples per design, identical settings throughout.
ipsae_median_of_5 is interface quality (higher better; parent = 0.908).
ipsae_spread_across_5 is prediction CONSISTENCY (lower better).
All 15 passed our stability gate: all 5 samples docked, spread <= 0.30.
The lab needs only design_name and sequence; the rest is supporting data.
Full method: results/SUMMARY.md and results/rationale.md.

## Columns

| column | what it is |
|---|---|
| `design_name` | Our identifier for the design. The prefix is the design strategy used. |
| `sequence` | Full scFv amino acid sequence, 242 aa, VH-(GGGGS)3-VL in that order. |
| `n_mutations` | Number of substitutions versus the parent trastuzumab scFv. |
| `mutations` | The substitutions, wild-type/position/mutant, 1-based on the 242 aa parent. Semicolon separated. |
| `sequence_length_aa` | Length in amino acids. All designs are 242 aa, the same as the parent. |
| `design_strategy` | Which of six design approaches produced this candidate. |
| `rationale` | Why this design was selected, opening with its mutation list. |
| `ipsae_median_of_5` | Median ipSAE interface score over 5 structure predictions. Higher is better; 0 to 1. The parent control scores 0.908 under the same settings. |
| `ipsae_all_5_samples` | The 5 individual ipSAE values, one per structure prediction, space separated. Reported in full rather than collapsed to the median. |
| `ipsae_spread_across_5` | Highest minus lowest of those 5 ipSAE values. A measure of prediction CONSISTENCY, not quality: low means all 5 predictions agreed. We reject designs above 0.30 because broken designs predict inconsistently rather than scoring low. All 15 here are at or below 0.102. |
| `samples_docked_of_5` | How many of the 5 predictions placed the scFv on the target at all (ipSAE >= 0.50). All 15 designs are 5 of 5. A design where some predictions collapse is rejected even if its median looks healthy. |
| `iptm_median_of_5` | Median ipTM from Boltz-2, for reference. A second interface confidence measure, 0 to 1, independent of ipSAE. |
| `passed_stability_gate` | yes for all 15: all 5 predictions docked AND spread at or below 0.30. |

## The one thing worth knowing about our scoring

A broken design does not score low. It scores *inconsistently* — some
of the 5 structure predictions place the binder correctly and others
collapse entirely. A scrambled-CDR-H3 negative control had a median
ipSAE of 0.877, inside the range of our real designs, while only 3 of
its 5 predictions docked. So we filter on consistency and docking
rate, never on the median, and we report all 5 values per design so
you can see the variation yourself.

Caveat we would rather state than have found: our filter tests whether
the predicted complex forms stably. It does not predict affinity. A
deliberate anti-design control passed it.

