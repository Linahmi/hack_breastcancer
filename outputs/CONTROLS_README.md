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

---

## Stability filter — PROVISIONAL (fitted on n=3)

Broken designs are **unstable, not low-scoring**. Each diffusion sample
either docks (~0.87) or collapses (<0.3); the empirical band between
~0.42 and ~0.75 is empty.

Reject rule, applied **BEFORE** the median, on the 5 individual samples:

- sample-spread > **0.20**
- or any single sample < **0.65**

```
                     good (parent + reverse)    broken (scrambled)
max sample-spread              0.1301                 0.3481
worst single sample            0.7779                 0.0000
```

Applied to these controls: 6/9 PASS, 3/9 REJECT.

**Never collapse a design to its median.** `cdrh3_scrambled_seed1` has
median ipSAE 0.8773, inside the parent's own range (0.8615-0.9218), and a
median-only filter passes it. Report all 5 samples, the docking rate, the
max sample-spread, and the median as secondary only.

**PROVISIONAL: n=3 control sequences. Hedi re-fits at 50 designs.**

Scoring: `python3 src/aida/score_designs.py --all --csv results/scored.csv`

### Unresolved: which CDR-H3 scramble is canonical

The `cdrh3_scrambled_*` outputs here, and the thresholds above, use a
13-residue permutation of positions **97-109**
(`SRWGGDGFYAMDY` -> `RMADDWGFYGSYG`), from `src/hedi/inputs/`.

An 11-residue permutation of **99-109**
(`WGGDGFYAMDY` -> `GFDMAWYGYDG`) was also prepared and was pushed in
commit dac6c33, then replaced in 5327aaa. Both preserve whole-sequence
composition. If we revert to the 99-109 version these thresholds must be
refitted.
