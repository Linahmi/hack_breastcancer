# Lina handoff for cofold and ipSAE

Validated on 31 October 2026. `results/candidates.csv` is the metadata source
for rationale and arm; the wave CSVs are the GPU inputs and intentionally carry
only `name,sequence`.

## Ready inputs

- `wave_01.csv` through `wave_04.csv`: 8 rows each, legacy designs produced
  before the CA-distance matrix arrived. All 32 are sequence-valid, but every
  row has at least one mutated-position CA pair below 10 A.
- `wave_05.csv`: 4 literature-arm designs. All four pass the current spreading
  rule (max 2 mutations per CDR and every mutated-position CA pair >=10 A).
- `controls.csv`: parent, reverse-order, scrambled-CDR-H3, and five-mutation
  anti-design controls. These are controls, not submission candidates.
- `geometry_audit.csv`: the authoritative per-design geometry verdict. Join on
  `name`; do not infer compliance from the wave number.

The four literature designs use exact published mutations. Four of the five
literature rows are loss-of-affinity alanine scans; only D102W was reported as
affinity-improving. Each candidate rationale records the source and effect.

## Aida

The wave and control CSVs satisfy `run_cofold.py`'s `name,sequence` contract.
Audit before spending GPU time:

```bash
modal run src/aida/run_cofold.py --csv results/wave_05.csv --dry-run
```

Then remove `--dry-run` and use the desired seed/concurrency settings. Waves
01-04 are useful as legacy exploration but are not compliant with the current
CA-spacing rule.

## Hedi

Score each single-run directory produced by Aida, for example:

```bash
python src/hedi/ipsae.py outputs/literature_01 10 15
```

Use the symmetric `max` ipSAE record for ranking, then join results back to
`candidates.csv` and `geometry_audit.csv` by `name`.
