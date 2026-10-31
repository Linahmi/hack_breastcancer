# Controls — force-added for Hedi's ipSAE reader

`outputs/` is gitignored; these nine runs were added with `git add -f`.

## What is here

Three controls × seeds 1, 2, 3:

| run prefix | sequence |
|---|---|
| `parent_wt` | parent trastuzumab scFv, 242 aa |
| `reverse_vl_vh` | VL-(GGGGS)3-VH, control only, **not submittable** |
| `cdrh3_scrambled` | CDR-H3 `WGGDGFYAMDY` → `GFDMAWYGYDG`, same composition |

Path pattern, `N = 0..4`:

```
outputs/<run>/boltz_results_<run>/predictions/<run>/pae_<run>_model_N.npz
outputs/<run>/boltz_results_<run>/predictions/<run>/plddt_<run>_model_N.npz
outputs/<run>/boltz_results_<run>/predictions/<run>/confidence_<run>_model_N.json
outputs/<run>/boltz_results_<run>/predictions/<run>/<run>_model_N.cif
outputs/<run>/boltz_results_<run>/processed/structures/<run>.npz   <- chain table
outputs/<run>/boltz_results_<run>/processed/constraints/<run>.npz
outputs/<run>/run_meta.json                                        <- settings
```

## npz array keys

One array per file, key equals the metric name:

| file | key | shape | dtype |
|---|---|---|---|
| `pae_*.npz` | `pae` | (384, 384) | float32 |
| `plddt_*.npz` | `plddt` | (384,) | float32 |

**`plddt` is 0–1, not 0–100.** Scale it if ipsae.py expects 0–100.

## Chain layout for slicing the 384×384 PAE

```
rows/cols   0-241  = chain A, scFv, 242 aa
rows/cols 242-383  = chain B, HER2 domain IV (UniProt 511-652), 142 aa
```

Machine-readable in `processed/structures/<run>.npz`, key `chains`:
`[('A', ..., res_start=0, res_count=242), ('B', ..., res_start=242, res_count=142)]`

## Settings — identical across all nine runs

```
--use_msa_server --write_full_pae --recycling_steps 10 --diffusion_samples 5 --step_scale 1.0
```

construct UniProt P04626 511–652, 142 aa, 10 disulfides; scFv chain 242 aa.

Disulfide constraints are derived per sequence, never copied:

```
parent_wt        [[22, 96], [158, 223]]
cdrh3_scrambled  [[22, 96], [158, 223]]
reverse_vl_vh    [[23, 88], [144, 218]]
```

## What is NOT here, and why

Full trees are 191 MB, over the ~50 MB budget. Excluded:

- `msa/` — 14 MB per run, 126 MB total. Regenerate with `--use_msa_server`.
- `pde_*.npz` — 2.5 MB per run, not used by ipSAE.
- `processed/mols`, `processed/msa`, `lightning_logs`.

Kept: 207 files, 33.6 MB. CIFs are included because ipSAE needs CB
coordinates as well as the PAE matrix.

## Do not confuse with `parent_quick`

`outputs/parent_quick/` is a structural smoke test — empty MSA, 1 recycle,
1 diffusion sample. **Not a baseline.** See its own `PROVENANCE.md`. The
valid parent baseline is `parent_wt_seed{1,2,3}`.
