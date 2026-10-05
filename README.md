# Confined accelerated-Collatz words — companion computations

Exact-arithmetic computations, pre-registrations and reports supporting

> E. De Jesús, *From Sturmian Capacity to Entropy-Deficit Survival: Exact Enumeration and Large
> Deviations in Confined Accelerated Collatz Words*, Zenodo, 2026.
> [doi:10.5281/zenodo.23167068](https://doi.org/10.5281/zenodo.23167068)

and its companion

> E. De Jesús, *A Sturmian Complexity Clock in Accelerated Collatz Words: Beatty Thresholds,
> Cumulative Multiplicity, and Bulk Resolution Rates*, Zenodo, 2026.
> [doi:10.5281/zenodo.22096151](https://doi.org/10.5281/zenodo.22096151)

**Nothing here bears on the Collatz conjecture or on any `O(log m)` occupation bound.** Every
confinement, prune, return and tie decision is an exact integer comparison; floating point appears
only in descriptive summaries, and is labelled where it does.

---

## Status of the work after v2 (October 2026)

| finding | status |
|---|---|
| The checkpoint decision tree of the clock note is a staircase-constrained lattice-path count; both horizon operators (`T0` identity, `T1` cumulative) are exact for every window | **PROVED** — Appendix C of the capacity note, v2 |
| Scalar reduction `Pi_h = T1^{K(h)} Pi_{h0}` on its base range `j <= h0` | **PROVED**; does *not* extend above `h0` |
| A reported failure of window-universality ("16 of 56", gaps 5, 7, 12) | **WITHDRAWN** — indexing artifact, reconstructed exactly |
| Global shape does not depend on the Diophantine type of the slope | **VERIFIED** on held-out data against a pre-registration |
| Unique maximum of the weighted profile for `r >= 20`, slopes `θ >= 0.1384` | **CONJECTURE**; threshold frozen in advance, 0 exceptions |
| Peak location `r - 4` | specific to `θ = log2 3 - 1` and the reference window; corrects the convention-dependent "`K - 1`" |
| Ties caused by near-balance | **REJECTED** |
| Fixed-modulus 2-adic exclusion of ties | **ROUTE CLOSED** |

Details, with effect sizes and permutation p-values, in the reports below.

## Layout

```
REPORT.md                      the operator study: reduction, T0/T1 proved, artifact diagnosis
PREDICTION.md                  its pre-registration (committed before the holdout ran)
code/                          exact tree, DP, chain count, sweeps, held-out tests
data/ out/                     tables and raw logs

shape/REPORT.md                slope-independence; the narrowed uniqueness conjecture (§7, §7b)
shape/REPORT_balance.md        why ties are not a balance phenomenon
shape/REPORT_2adic.md          why no fixed-modulus 2-adic argument excludes ties
shape/PREDICTION.md            pre-registration for the slope study
shape/PREDICTION_balance.md    pre-registration for the balance study
shape/NEXT_REVISION_NOTE.md    draft wording for a future revision of Remark C.16
shape/code/ shape/out/ shape/figs/

ERRATA_clock_note.md           errata for the clock note, keyed to its own claim numbers
ZENODO_DESCRIPTION_UPDATE.md   draft text for the Zenodo record
```

The LaTeX of the capacity note itself is not mirrored here; the published version is the Zenodo
record linked above. `shape/NEXT_REVISION_NOTE.md` holds draft wording for its next revision.

## Reproducing

Python 3.12, `numpy`, `matplotlib` (figures only). No other dependencies; everything else is the
standard library with exact integer arithmetic.

```bash
python3 code/stage0_validate.py        # cross-validate tree vs DP vs chain count
python3 code/stage3_holdout.py         # the operator law, held out
python3 shape/code/stage0_shape.py     # reproduce the published W0 numbers; local law, all slopes
python3 shape/code/stage4_holdout.py   # slope-independence, held out
python3 shape/code/stage_b4_2adic.py   # the 2-adic test
```

Scripts regenerate their own CSVs; the large ones are not tracked — `shape/data/README.md` gives
the run order. Every number quoted in a report is preserved in `out/` and `shape/out/`, so claims
can be checked without rerunning anything. Scripts resolve their own paths and may be run from any
working directory.

## Licence

MIT for the code (`LICENSE`); CC BY 4.0 for the reports, notes and figures (`LICENSE-DOCS.md`).

## Method notes

* **Pre-registration.** Each study froze its predictions in a `PREDICTION*.md` committed *before*
  the held-out computation. Where a frozen prediction failed it is reported as failed — see
  `shape/REPORT.md` §7b (a tie threshold of `K >= 17` fitted on exploration broke on holdout and
  was corrected to `K >= 19`), and `shape/REPORT_balance.md` §4 (two of my own numeric brackets
  missed while the hypothesis under test was still rejected).
* **Exactness.** `floor(n*theta)` is computed by bracketing between successive continued-fraction
  convergents, which *proves* the value; for `theta = log2 3 - 1` it is cross-checked against
  `(3^n).bit_length() - 1 - n`.
* **Labels.** Claims are tagged PROVED / VERIFIED / OBSERVED / CONJECTURE throughout, and each
  report ends with a claim audit.
