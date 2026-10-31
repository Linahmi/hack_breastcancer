"""Runs the HER2/scFv cofold across several seeds in parallel on Modal.

Reuses the lean Boltz-2 app/function from modal_boltz2.py and fans the seeds
out with .starmap, so all seeds occupy separate containers concurrently.

The entrypoint is named `fanout`, not `main`: modal_boltz2 already registers a
`main` local entrypoint on this same app, and two entrypoints sharing a name
is a hard InvalidError raised at import time. The name must therefore be
given explicitly after `::`.

Usage:
    modal run run_seeds.py::fanout --input-yaml inputs/her2_scfv_cofold.yaml --seeds 1,2,3

NOTE: every seed here reaches boltz() without a `name`, so all of them use the
default stem "fixed" and land on identical paths -- the files are written to
separate out_dir/seed<N>/ directories, so they do not clobber, but each tree
is internally named boltz_results_fixed/predictions/fixed/pae_fixed_*.npz.
For named, resumable, per-row output use run_cofold.py instead.
"""

from pathlib import Path

from modal_boltz2 import app, boltz

BASE_PARAMS = "--use_msa_server --write_full_pae --recycling_steps 10 --diffusion_samples 5 --step_scale 1.0"


@app.local_entrypoint()
def fanout(
    input_yaml: str = "inputs/her2_scfv_cofold.yaml",
    seeds: str = "1,2,3",
    out_dir: str = "./out/boltz",
):
    import time

    input_str = open(input_yaml).read()
    seed_list = [int(s) for s in seeds.split(",")]

    args = [(input_str, f"{BASE_PARAMS} --seed {s}") for s in seed_list]

    t0 = time.monotonic()
    results = list(boltz.starmap(args))
    elapsed = time.monotonic() - t0

    for seed, outputs in zip(seed_list, results):
        dest = Path(out_dir) / f"seed{seed}"
        for out_file, out_content in outputs:
            target = dest / Path(out_file)
            target.parent.mkdir(parents=True, exist_ok=True)
            target.write_bytes(out_content)
        print(f"seed {seed}: {len(outputs)} files -> {dest}")

    print(f"\n{len(seed_list)} seeds wall-clock: {elapsed:.1f}s")
