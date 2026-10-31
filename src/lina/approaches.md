# PLM scan — decisions to make

Options and trade-offs only. Nothing chosen, nothing run beyond the environment check.

## Verified environment (31 Oct, before any scan)

| Thing | Status |
|---|---|
| `python scripts/parent.py` | `VH 1-120`, `linker 121-135`, `VL 136-242`, 242 aa, 18 aa headroom |
| ESM2 650M (`fair-esm` 2.0.0) | loads on CPU, 2.6 GB cached; `logits (1, 244, 33)` = 242 + cls + eos; **7.3 s** per forward pass |
| AbLang2 paired (`ablang2` 0.2.1) | loads on CPU; `likelihood (1, 232, 26)`; 0.3 s per forward pass |
| AbLang2 position mapping | **proven**, not assumed. Layout is `< VH > \| < VL >`, specials at 0/121/122/123/231. Offset 0 recovers WT at 91.2%; every other offset in ±4 gives ~7%. All 227 VH+VL positions round-trip to the right parent letter. |
| AbLang2 output scale | raw **logits**, not probabilities (p(wt) came back >1). Needs an explicit softmax. |
| Germline reference | offline, no download needed — ANARCI ships the IMGT DB as Python: 249 human IGHV, 82 IGKV, 86 IGLV, IMGT-gapped |
| ANARCI *numbering* | **unavailable on this machine.** It shells out to `hmmscan`; HMMER is not installed and `pyhmmer` is POSIX-only. The HMM profiles ship but the binary does not. Germline assignment therefore has to be alignment-based (biopython `PairwiseAligner`, verified working). The `imgt` column is Ghita's anyway. |

Alignment top hits for the parent: VH → `IGHV3-66` (3-66\*01/\*02/\*04 tied), VL → `IGKV1-39`/`IGKV1D-7-1`. Consistent with trastuzumab being humanised onto a human VH3 / Vκ1 framework. Note the ties — see (e).

Cost anchor: a full masked-marginals scan is 242 forward passes ≈ **30 min** CPU unbatched, less when batched. AbLang2 is ~25× cheaper per pass.

---

## a. Scoring scheme

| Option | Cost | What it assumes / how it fails |
|---|---|---|
| **masked-marginals** — mask position *i*, read log P(aa) at *i* | 242 passes, ~30 min | The ESM-paper default for variant effect. No WT leakage. Still scores every mutation in the *parent* background, so it says nothing about combinations. |
| **wt-marginals** — one unmasked pass, read logits at every position | 1 pass, 7 s | The WT token is visible at its own position, so the model largely copies it: p(WT) inflates and mutation effects compress. Correlates with masked-marginals but ranks worse on DMS benchmarks. Cheap enough to compute as a cross-check either way. |
| **full pseudo-likelihood per variant** | infeasible for singles (242×19 variants × 242 passes); ~36k passes for 150 candidates | The only one that scores a *sequence* rather than a substitution. Out of reach for the matrix; reachable for a shortlist. |
| **masked passes at mutated positions only, in the mutant background** | ~5–10 passes per variant | Middle option: catches sequence-model epistasis between the chosen mutations at ~1 min per variant. Does not model structure. |

Note for AbLang2: its `mode="likelihood"` is a single unmasked pass, so **it is wt-marginals by construction**. Matching ESM2's masked-marginals means masking AbLang2's input by hand (`*` token) — a separate 227-pass run, cheap at 0.3 s/pass. If we don't, the two columns are not the same kind of number, which feeds straight into (b).

## b. Combining ESM2 and AbLang2

They are on different scales, and AbLang2 is far more peaked (antibody-specialised).

| Option | Keeps | Loses / fails when |
|---|---|---|
| **rank-average (Borda)** | scale-free, outlier-robust | discards magnitude — a position where one model is emphatic and the other indifferent looks like mild agreement |
| **z-score sum** | magnitude | assumes comparable distribution shape, which these two do not have. Sub-choice: *global* z keeps "this position is conserved at all", *per-position* z removes it and scores only "which residue, given this position" |
| **intersection of top-k** | precision | recall, and *k* is arbitrary. May not yield enough distinct mutations to build 150 candidates |
| **reciprocal rank fusion** (Σ 1/(60+rank)) | agreement near the top | still rank-only; the constant is a free parameter |
| **don't combine** — one arm per model + an agreement arm | disagreement stays visible as its own category | nothing, except that "consensus" then isn't a single number. The `tool` column is already built for arms |

This choice is coupled to (e): **any consensus that includes AbLang2 inherits the germline pull.** Agreement between a general protein LM and an antibody repertoire LM is partly just "this looks like a normal antibody", which for an affinity-matured therapeutic is the wrong objective.

## c. Which mutations qualify

| Option | Rows | Bias it introduces |
|---|---|---|
| **all 242 × 19** | 4598 | none; keeps the indifference signal the brief asks for. Linker rows are meaningless but harmless |
| **only where some non-WT outranks WT** | few hundred | **this filter preferentially proposes undoing the engineering.** It fires where the parent residue is unusual — which in a humanised, affinity-matured scFv is precisely the CDR positions that were deliberately engineered |
| **margin threshold** (log-ratio > +0.5) | fewer | same bias, sharper |
| **exclude linker 121–135 always** | −285 rows | free; the rules fix the linker and the validator checks it |

Orthogonal: the filter can live *in* `substitutions.csv` (fewer rows) or *downstream* at assembly (all rows, filter later). The column spec has no "qualifies" field, which argues for keeping all rows and filtering at assembly — but that's a call, not a constraint.

## d. Singles → 5+ mutation variants (the epistasis problem)

| Option | Assumes | How it fails |
|---|---|---|
| **stack the top-N singles** | additivity | the brief's explicit warning. Breaks on same-position/adjacent pairs, charge clusters in one CDR, and any correlated pair. Produces one variant per N, so it also wastes the budget |
| **sample many combinations from a qualified pool** | that the downstream scorer can separate them | this is a bet on ipSAE — and whether ipSAE resolves anything is exactly what Hedi's 13:00 control test decides. If it can't split parent from reverse-order, broad sampling is rolling dice |
| **spread across regions, cap per CDR** (≤2/CDR, ≥1 in VL) | dispersed mutations are more independent | true for packing and charge clustering, false for long-range electrostatics and for VH/VL interface pairs |
| **contact + second-shell pairing** (needs Ghita) | positions.csv numbering matches ours | gives every variant a mechanism instead of a score. Blocked on her file — see `dependencies.md` |
| **build, then rescore each combination in its own background** | sequence-model epistasis is the epistasis that matters | ~1 min/variant, catches the worst additivity failures; blind to structural epistasis |

On the 8–10 mutation gamble arm: more mutations monotonically raises the chance of losing binding *and* of losing expression, and expression is publicly reported. Worth having as **variance**, not as expectation — a couple of designs, not a tranche.

## e. Germline flag

**Reference**

- ANARCI's bundled IMGT human V genes — offline, deterministic, already verified. Covers V only.
- IMGT/V-QUEST or OAS over the network — authoritative, slow, no real gain over the bundled IMGT for V-gene assignment.
- Granularity: single best allele / top-k set / whole family. The top hits came back **tied** (`IGHV3-66*01/*02/*04` all at 377), so "the single best allele" is an arbitrary tie-break. Family level (`IGHV3-66`, or all `IGHV3`) is robust but blurs.

**The CDR3 hole:** CDR-H3 is the V-D-J junction and has **no germline V counterpart**. Those positions cannot be flagged either way, so the column needs an explicit third value — not a silent `no`. Given CDR-H3 is where a lot of the affinity lives, this is not an edge case.

**Defining "moves toward germline"**

- **strict identity** — mutant letter == germline letter at the aligned position, and parent differs. Unambiguous, binary.
- **similarity** — mutant is closer to germline than parent is, by BLOSUM. Catches conservative drift; fuzzier, needs a threshold.
- **per-variant count** — germline reversions per sequence rather than per mutation. Different object: flags candidates, not substitutions. Possibly both.

**Alignment caveat:** IMGT germlines are IMGT-*gapped*, and we cannot produce IMGT numbering here (no `hmmscan`). So either align ungapped with BLOSUM62 and accept indel misregistration risk in CDR1/CDR2, or recover registration from the germline's own gap pattern. The first is verified working; the second is more faithful.

**Suggested value set** (the brief says flagged mutations become their own category, so the column probably can't be a bool): `toward_germline` / `away_from_germline` / `no_germline_equivalent` / `already_germline`.

**Why this matters enough to be deliverable 2:** trastuzumab's CDRs carry affinity-driven, deliberately non-germline residues from the murine 4D5 parent. AbLang2 scores germline-like residues as more natural, so ESM2 × AbLang2 consensus will propose reverting them. The flag is what separates *"the model likes this because it's natural"* from *"the model likes this because it's better"* — and those point in opposite directions here.
