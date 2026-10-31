"""AbLang2 paired masked-marginals over the parent VH/VL: 227 positions x 20 residues.

    python src/lina/scan_ablang2.py

Writes src/lina/ablang2_masked_marginals.csv and src/lina/ablang2_logprobs.npy.

Masking is AbLang2's own `_predict_logits_with_step_masking` (encodings.py:28-47):
genuine masked-marginals - one token masked at a time, prediction read at the
masked position. The softmax is applied here explicitly, because the raw method
returns logits (mode="likelihood" returns those unnormalised; p(wt) came back >1
during verification).

The public `mode=...`/`align=True` wrapper is deliberately bypassed: align=True
routes through ANARCI numbering, which needs hmmscan and dies on this machine.
The internal path takes the same formatted "<VH>|<VL>" string and skips numbering.

Frame guard, two layers:
  1. exact - every mapped column must decode to the right parent residue;
  2. calibrated - the UNMASKED (wt-marginals) pass must recover WT at >90%,
     the quantity measured at 91.2% during environment verification. The
     masked recovery rate is reported but not asserted: masking legitimately
     lowers it, so a threshold there would be a quality check wearing a
     frame-guard's clothes.

The linker is absent by construction here - AbLang2 scores VH and VL only - so
linker positions 121-135 simply have no AbLang2 row.
"""

import csv
import sys
import time
from pathlib import Path

import numpy as np
import torch

REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(REPO))
import ablang2  # noqa: E402
from ablang2.pretrained import format_seq_input  # noqa: E402

from scripts.parent import PARENT, VH, VL, region  # noqa: E402

CANON = "ACDEFGHIKLMNPQRSTVWY"
CYS_POS = {22, 96, 158, 223}
LINKER_POS = set(range(121, 136))
EXCLUDED = CYS_POS | LINKER_POS
SPECIALS = set("<>|*X-")
OUT_DIR = Path(__file__).resolve().parent


def main() -> int:
    t0 = time.time()
    vh = PARENT[VH[0] - 1 : VH[1]]
    vl = PARENT[VL[0] - 1 : VL[1]]
    print(f"VH {len(vh)} aa, VL {len(vl)} aa, {len(vh)+len(vl)} scored positions")

    model = ablang2.pretrained(
        model_to_use="ablang2-paired", random_init=False, ncpu=1, device="cpu"
    )
    tok = model.tokenizer
    inv = {v: k for k, v in tok.aa_to_token.items()}

    # --- token layout, and position -> column map ---
    seqs, chain = format_seq_input([[vh, vl]], fragmented=False)
    print(f"formatted input: {len(seqs[0])} chars, chain {chain}")
    ids = np.asarray(tok(seqs, pad=True, w_extra_tkns=False, device="cpu"))[0]
    toks = [inv[int(i)] for i in ids]
    print(f"token layout: {len(toks)} columns, specials at "
          f"{[(i, t) for i, t in enumerate(toks) if t in SPECIALS]}")

    pos2col, k = {}, 0
    for i, t in enumerate(toks):
        if t in SPECIALS:
            continue
        pos = (VH[0] + k) if k < len(vh) else (VL[0] + k - len(vh))
        pos2col[pos] = i          # verified offset 0 against this layout
        k += 1

    # guard 1: exact frame check
    for p, col in pos2col.items():
        if toks[col] != PARENT[p - 1]:
            raise SystemExit(
                f"FRAME SHIFT: parent position {p} is '{PARENT[p-1]}' but column "
                f"{col} decodes to '{toks[col]}'. Aborting before any scoring."
            )
    print(f"frame guard OK: all {len(pos2col)} VH+VL positions map correctly")

    # guard 2: calibrated, on the unmasked pass (the 91.2% figure)
    with torch.no_grad():
        unmasked = torch.softmax(model._predict_logits(seqs), dim=-1).cpu().numpy()[0]
    rec = sum(
        1 for p, col in pos2col.items()
        if inv[int(np.argmax(unmasked[col]))] == PARENT[p - 1]
    )
    frac = rec / len(pos2col)
    print(f"unmasked WT recovery {rec}/{len(pos2col)} = {frac:.1%}")
    if frac < 0.90:
        raise SystemExit(
            f"GUARD FAILED: unmasked WT recovery {frac:.1%} < 90%. Verified at "
            f"91.2%; a drop this large means the frame or the model changed."
        )

    # --- masked-marginals ---
    print(f"running stepwise-masked pass ({len(toks)} masked forward passes)...",
          flush=True)
    with torch.no_grad():
        logits = model._predict_logits_with_step_masking(seqs)
    probs = torch.softmax(logits, dim=-1).cpu().numpy()[0]   # softmax applied here
    print(f"masked probability matrix {probs.shape} in {time.time()-t0:.0f}s")

    masked_rec = sum(
        1 for p, col in pos2col.items()
        if inv[int(np.argmax(probs[col]))] == PARENT[p - 1]
    )
    print(f"masked argmax recovers WT at {masked_rec}/{len(pos2col)} = "
          f"{masked_rec/len(pos2col):.1%} (reported, not asserted)")

    # renormalise over the 20 canonical residues, then store log-probabilities
    aa_cols = [tok.aa_to_token[a] for a in CANON]
    logp = np.zeros((len(PARENT), len(CANON)), dtype=np.float32)
    have = np.zeros(len(PARENT), dtype=bool)
    for p, col in pos2col.items():
        sub = probs[col][aa_cols]
        sub = sub / sub.sum()
        logp[p - 1] = np.log(np.clip(sub, 1e-12, None))
        have[p - 1] = True
    np.save(OUT_DIR / "ablang2_logprobs.npy", logp)

    out_csv = OUT_DIR / "ablang2_masked_marginals.csv"
    n = 0
    with open(out_csv, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow([
            "position", "wt_aa", "mut_aa", "region", "excluded",
            "ablang2_logp_mut", "ablang2_logp_wt", "ablang2_score",
        ])
        for p in range(1, len(PARENT) + 1):
            if not have[p - 1]:
                continue                      # linker: AbLang2 never sees it
            wt = PARENT[p - 1]
            lp_wt = float(logp[p - 1, CANON.index(wt)])
            for j, mut in enumerate(CANON):
                w.writerow([
                    p, wt, mut, region(p), int(p in EXCLUDED),
                    f"{logp[p-1, j]:.6f}", f"{lp_wt:.6f}",
                    f"{logp[p-1, j] - lp_wt:.6f}",
                ])
                n += 1
    print(f"wrote {out_csv.name}  ({n} rows = {int(have.sum())} x {len(CANON)})")
    print(f"total {time.time()-t0:.0f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
