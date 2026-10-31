# parent_quick — structural smoke test (NOT a scientific result)

Force-added to git (`outputs/` is gitignored) for one reason only: to give
Hedi a **real** Boltz-2 output tree to write the ipSAE reader against.

## Do not use these numbers

This run was deliberately configured for speed, not quality:

| setting | value | normal pipeline value |
|---|---|---|
| MSA | `msa: empty` (single-sequence mode) | `--use_msa_server` |
| `--recycling_steps` | 1 | 10 |
| `--diffusion_samples` | 1 | 5 |

Boltz itself warns: *"will run these in single sequence mode … predictions
will be suboptimal without an MSA."* The `iptm` of 0.843 in
`confidence_parent_quick_model_0.json` is **not** a usable parent baseline.
The real 3-seed baseline is a separate run.

## What it is

Parent trastuzumab scFv cofolded against HER2 domain IV.

- chain A = parent scFv, 242 aa, PAE rows/cols **0–241**
- chain B = HER2 domain IV, UniProt P04626 **511–652**, 142 aa, PAE rows/cols **242–383**
- 12 disulfide bond constraints applied (2 scFv + 10 HER2); confirmed by
  `processed/structures/parent_quick.npz` → `bonds` has shape `(12,)`
- GPU L40S, 22 s prediction wall-clock, seed 1, boltz 2.2.1
- input: `src/aida/inputs/parent_smoketest.yaml`

## Path layout — read this

Boltz inserts a `boltz_results_<stem>` level that the upstream docs omit.
The path circulating in the team, `outputs/<name>/predictions/<name>/…`,
is **missing that level** and matches nothing on disk. The real path is:

```
outputs/<run>/boltz_results_<run>/predictions/<run>/pae_<run>_model_0.npz
```

`src/aida/run_cofold.py` preserves this real layout, so this example and
live pipeline output have the same shape.
