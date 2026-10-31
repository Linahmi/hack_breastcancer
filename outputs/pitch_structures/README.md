# pitch_structures — two cofold models for the ChimeraX figure

`outputs/` is gitignored, so these two files are force-added (`git add -f`)
to make them visible without shipping the 34 MB of run data.

## The two files

| file | design | run | seed | model | ipSAE of THIS sample |
|---|---|---|---|---|---|
| `consensus_green_03_seed1_model_2.cif` | `consensus_green_03` | `outputs/consensus_green_03` | 1 | **2** | **0.9265** |
| `w01_germline_revert_02_seed1_model_1.cif` | `w01_germline_revert_02` | `outputs/w01_germline_revert_02` | 1 | **1** | **0.9369** |

## These are REPRESENTATIVE, not cherry-picked

Each design was cofolded once with 5 diffusion samples. The sample chosen
is the one whose ipSAE **equals the design's reported median** — so the
figure matches the number quoted in `results/final_15.csv`, not a
flattering outlier.

```
consensus_green_03      m0=0.9360  m1=0.9230  m2=0.9265  m3=0.9360  m4=0.9113
                        median 0.9265  -> model_2 chosen
                        (best sample was model_0 at 0.9360, deliberately NOT used)

w01_germline_revert_02  m0=0.9454  m1=0.9369  m2=0.9141  m3=0.8872  m4=0.9510
                        median 0.9369  -> model_1 chosen
                        (best sample was model_4 at 0.9510, deliberately NOT used)
```

Note neither is `model_0`. Boltz writes `model_0` first but does not rank
samples by ipSAE, so `model_0` is not the median one in either case.

Both designs pass the gate at 5/5 samples docked, spreads 0.0247 and
0.0638 against a 0.30 limit.

## Chains — colour the right one

Verified by parsing the files themselves, not assumed:

```
chain A = scFv              label_seq_id 1-242   (242 aa)
chain B = HER2 domain IV    label_seq_id 1-142   (142 aa)
```

Chain B is **UniProt P04626 511-652** renumbered 1-142 in the file. To
convert a chain B residue number to UniProt numbering, add 510. To PDB
1N8Z chain C numbering, add 488.

The trastuzumab epitope sits in chain B around UniProt 579-627, which is
**chain B residues 69-117** in these files — that is the interface worth
showing.

The scFv is a single chain: VH 1-120, (GGGGS)3 linker 121-135, VL 136-242.

## Provenance

Production settings, identical for both:

```
--use_msa_server --write_full_pae --recycling_steps 10
--diffusion_samples 5 --step_scale 1.0 --seed 1
```

Boltz 2.2.1, L40S. 12 disulfide bond constraints applied per structure
(2 scFv, 10 HER2 domain IV).

Full scores for both designs, including all 5 samples, are in
`results/scored.csv`. The mutation lists are in `results/final_15.csv`.
