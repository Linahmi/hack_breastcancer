"""ESM2 650M masked-marginals over the full parent: 242 positions x 20 residues.

    python src/lina/scan_esm2.py

Writes src/lina/esm2_masked_marginals.csv incrementally (full matrix, nothing
filtered) and src/lina/esm2_logprobs.npy at successful completion (242 x 20
log-probabilities, canonical AA order). The CSV is flushed and synced every
10 positions so an interrupted run retains its completed checkpoints.

Frame guard: every parent position is checked token-by-token against PARENT
before it is masked, so a frame shift aborts the run instead of producing a
plausible-looking matrix shifted by a constant.
"""

import csv
import os
import sys
import time
from pathlib import Path

import numpy as np
import torch

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO))
import esm  # noqa: E402

from scripts.parent import PARENT, region  # noqa: E402

CANON = "ACDEFGHIKLMNPQRSTVWY"
LINKER_POS = set(range(121, 136))
CYS_POS = {22, 96, 158, 223}
EXCLUDED = LINKER_POS | CYS_POS
BATCH = 8
CHECKPOINT = 10
OUT_DIR = Path(__file__).resolve().parent


def main() -> int:
    t_start = time.time()
    print(f"parent {len(PARENT)} aa")

    model, alphabet = esm.pretrained.esm2_t33_650M_UR50D()
    model = model.eval()
    bc = alphabet.get_batch_converter()
    _, _, toks = bc([("parent", PARENT)])
    print(f"tokens {tuple(toks.shape)}  (cls at 0, parent position p at token p)")

    # --- frame guard, exact. Position p (1-based) must sit at token index p. ---
    for p in range(1, len(PARENT) + 1):
        tok_letter = alphabet.get_tok(int(toks[0, p]))
        if tok_letter != PARENT[p - 1]:
            raise SystemExit(
                f"FRAME SHIFT: parent position {p} is '{PARENT[p-1]}' but token {p} "
                f"decodes to '{tok_letter}'. Aborting before any scoring."
            )
    print(f"frame guard OK: all {len(PARENT)} positions decode to the right residue")

    aa_idx = torch.tensor([alphabet.get_idx(a) for a in CANON])
    logp = np.zeros((len(PARENT), len(CANON)), dtype=np.float32)

    out_csv = OUT_DIR / "esm2_masked_marginals.csv"
    positions = list(range(1, len(PARENT) + 1))
    active_batch = BATCH
    done = 0

    with open(out_csv, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow([
            "position", "wt_aa", "mut_aa", "region", "excluded",
            "esm2_logp_mut", "esm2_logp_wt", "esm2_score",
        ])
        fh.flush()
        os.fsync(fh.fileno())

        # Checkpoint in groups of 10 positions. Within each checkpoint, score
        # masked copies in batches of 8 (then 2). If the first batched forward
        # pass runs out of memory, retry that position and the rest sequentially.
        for checkpoint_start in range(0, len(positions), CHECKPOINT):
            checkpoint_positions = positions[
                checkpoint_start : checkpoint_start + CHECKPOINT
            ]
            offset = 0
            while offset < len(checkpoint_positions):
                chunk = checkpoint_positions[offset : offset + active_batch]
                batch = toks.repeat(len(chunk), 1)
                for row, p in enumerate(chunk):
                    batch[row, p] = alphabet.mask_idx
                try:
                    with torch.no_grad():
                        out = model(batch)["logits"]
                except RuntimeError as exc:
                    if active_batch == 1 or "out of memory" not in str(exc).lower():
                        raise
                    print(
                        f"batch size {active_batch} ran out of memory; "
                        "falling back to sequential scoring",
                        flush=True,
                    )
                    active_batch = 1
                    del batch
                    if torch.cuda.is_available():
                        torch.cuda.empty_cache()
                    continue

                for row, p in enumerate(chunk):
                    lp = torch.log_softmax(out[row, p], dim=-1)
                    logp[p - 1] = lp[aa_idx].cpu().numpy()
                offset += len(chunk)
                del batch, out

            # Append this completed checkpoint to disk before continuing.
            for p in checkpoint_positions:
                wt = PARENT[p - 1]
                lp_wt = float(logp[p - 1, CANON.index(wt)])
                for k, mut in enumerate(CANON):
                    w.writerow([
                        p, wt, mut, region(p), int(p in EXCLUDED),
                        f"{logp[p-1, k]:.6f}", f"{lp_wt:.6f}",
                        f"{logp[p-1, k] - lp_wt:.6f}",
                    ])
            fh.flush()
            os.fsync(fh.fileno())

            done += len(checkpoint_positions)
            elapsed = time.time() - t_start
            eta = elapsed / done * (len(positions) - done)
            print(
                f"  {done:3}/{len(positions)} positions  "
                f"{elapsed/60:.1f} min elapsed  eta {eta/60:.1f} min  "
                f"CSV checkpoint synced",
                flush=True,
            )

    np.save(OUT_DIR / "esm2_logprobs.npy", logp)

    # Report (do not assert) masked argmax WT recovery: informative, but on an
    # affinity-matured antibody it is legitimately well under 90%.
    recovered = sum(
        1 for p in range(1, len(PARENT) + 1)
        if CANON[int(np.argmax(logp[p - 1]))] == PARENT[p - 1]
    )
    print(f"\nmasked argmax recovers WT at {recovered}/{len(PARENT)} = "
          f"{recovered/len(PARENT):.1%} of positions (reported, not asserted)")

    n_rows = len(PARENT) * len(CANON)
    print(f"wrote {out_csv.name}  ({n_rows} rows = {len(PARENT)} x {len(CANON)})")
    print(f"excluded positions marked: {sorted(EXCLUDED)[:6]}... "
          f"({len(EXCLUDED)} positions, {len(EXCLUDED)*len(CANON)} rows)")
    print(f"total {(time.time()-t_start)/60:.1f} min")
    return 0


if __name__ == "__main__":
    sys.exit(main())
