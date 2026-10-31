# Lina handoff for cofold and ipSAE

Validated on 31 October 2026. `results/candidates.csv` is the metadata source
for rationale and arm; the wave CSVs are the GPU inputs and intentionally carry
only `name,sequence`.

## Ready inputs

- `wave_01.csv` through `wave_04.csv`: 8 rows each. All 32 are sequence-valid.
  CA distance is a post-hoc covariate, not an exclusion rule.
- `wave_05.csv`: 4 rebuilt literature-arm designs plus `ctrl_antidesign`. The
  four candidates were prospectively constructed with max 2 mutations per CDR
  and every mutated-position CA pair >=10 A; the anti-design is a negative
  control.
- `controls.csv`: parent, reverse-order, scrambled-CDR-H3, and five-mutation
  anti-design controls. These are controls, not submission candidates.
- `wave_06.csv`: the four previously undispatched consensus/contact designs.
- `wave_07.csv`: 8 ESM2-only designs using 40 positions absent from the prior
  pool; all obey the prospective CDR and CA-separation constraints.
- `geometry_audit.csv`: per-design CA-distance covariates for correlation with
  docking rate. It is not a compliance or filtering verdict.

Two literature designs use D102W, the only affinity-improving published result,
with both-model consensus fillers. Two use consensus mutations while excluding
the published alanine-scan hot spots. Each candidate rationale cites the source
for the D102W inclusion and/or the hot-spot exclusions.

## Aida

The wave and control CSVs satisfy `run_cofold.py`'s `name,sequence` contract.
Audit before spending GPU time:

```bash
modal run src/aida/run_cofold.py --csv results/wave_07.csv --dry-run
```

Then remove `--dry-run` and use the desired seed/concurrency settings. Analyze
minimum CA distance as a covariate across all waves; do not use 10 A as a
retrospective filter.

## Hedi

Score each single-run directory produced by Aida, for example:

```bash
python src/hedi/ipsae.py outputs/literature_01 10 15
```

Use the symmetric `max` ipSAE record for ranking, then join results back to
`candidates.csv` and `geometry_audit.csv` by `name`.
