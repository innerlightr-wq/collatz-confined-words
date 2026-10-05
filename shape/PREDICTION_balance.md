# PREDICTION_balance.md — frozen before the HOLDOUT run

Frozen at the timestamp of the git commit adding this file. HOLDOUT = 8 new random-CF slopes
(seed 20261007, none sharing a CF with any earlier slope), `theta` spread over `[0.15, 0.9]`,
windows `(2,1), (4,2), (7,1)`, horizons `h <= 400`. Nothing below has been computed on them.

## Definitions (fixed in Stage 1, not retuned)

`w(j) = Pi_h(j) 2^{-j}`; `rho(j) = w(j+1)/w(j) = Pi(j+1)/(2 Pi(j))` as an exact Fraction;
peak `j*` = argmax `w` (smallest index if tied); balance margin
`m = min(|log rho(j*-1)|, |log rho(j*)|)`; top-irregularity `D` = number of staircase drops in
`j in [j*-5, j*+5]`; `k_fit = argmin_k TV(w, w*_k)` over `k` within 20 of `round(mean)`,
`S4` = that TV; `Delta = Pi(j*+1) - 2 Pi(j*)` (exact integer; a tie is `Delta = 0`).

## The mechanism's predictions

**P1 — all ties are adjacent.** Predicted: **0 non-adjacent ties** on HOLDOUT. Note the Stage-1
wording correction: tied sets are adjacent *runs*, of length 2 (78%) or 3 (22%), not only pairs.

**P2 — margin `m` increases with top-irregularity `D`.** The mechanism requires a positive rank
correlation after `theta` and `r` are held fixed (pooled `D` is ~`11 theta` and would merely
restate `theta`). EXPLORE gave median within-(slope, r-band) Spearman `-0.018` with 44% of cells
positive. **Predicted on HOLDOUT: P2 FAILS** — median within-cell Spearman in `[-0.10, +0.10]`
and fraction of positive cells in `[0.35, 0.65]`, i.e. indistinguishable from no relation.

**P3 — at matched `r`, median `m` and median `S4` both increase with `theta`.** EXPLORE, `r` band
10–20: `m` rose 0.0009 → 0.0152 (17x) and `S4` 0.0226 → 0.0786 (3.5x) from the lowest to the
highest `theta` band. **Predicted on HOLDOUT: P3 HOLDS** — in the `r` band 10–20, `m` rises by a
factor `>= 3` and `S4` by `>= 1.5` from the lowest to the highest `theta` band, with
Spearman(`theta`, `m`) `> 0` at `p < 0.05` (permutation).

**P4 — ties occur only where the ideal balance point is within one index of `j*`.**
**As stated this is vacuous**: in EXPLORE `|j* - (k_fit - 2)| <= 1` for **100.0%** of rows, tied
and untied alike, so it cannot discriminate. Predicted on HOLDOUT: `>= 99%` for both classes.
**Sharpened P4'**, which can discriminate: `P(j* = k_fit - 2 exactly | tied)` exceeds
`P(j* = k_fit - 2 exactly | untied)`. EXPLORE: 0.730 vs 0.458. **Predicted: the tied rate exceeds
the untied rate by at least 0.15.**

## The competing hypothesis (ties as small-integer coincidence)

A tie requires `Pi(j*+1) = 2 Pi(j*)` **exactly**. Stage 1 found the relative margin `m` *shrinks*
with `r` (0.47 at `r<5` to 0.0011 at `r>=100`: the profile becomes steadily more balanced) while
the absolute integer gap `|Delta|` *explodes* (median 3 to 3e66). Ties track `|Delta|`, not `m`.

**P5 — ties are an arithmetic small-number effect.** Predicted on HOLDOUT:
* **no tie at any row with `Pi(j*) >= 10^10`** (EXPLORE: the largest `Pi(j*)` carrying a tie was
  300,540,195, 9 digits; untied rows reach 101 digits);
* **tie rate exactly 0 for every `r >= 20`**;
* conditioning on the **10% most balanced rows** of each `r` band, the tie rate falls to **0 for
  every `r` band with `r >= 20`**, despite those rows having `m <= ~4e-4` — i.e. maximal balance
  with no ties. (EXPLORE: that conditional tie rate was 1.00, 1.00, 0.82, 0.15, 0.00, 0.00, 0.00,
  0.00 across `r` bands 1-5, 5-10, 10-15, 15-20, 20-30, 30-50, 50-100, 100-400.)
* median `|Delta|` increases monotonically across `r` bands.

## Verdict rule (frozen)

* **SUPPORTED** if P1–P4 hold *and* P5's decisive clause fails, i.e. ties do occur at large `r`
  where balance is tightest.
* **REJECTED** if P2 fails and P5 holds: balance improves with `r` while ties vanish, so balance
  is not what selects a tie.
* **PARTIAL** otherwise — stating exactly which parts hold.

A single tie at `r >= 20`, or at `Pi(j*) >= 10^10`, would falsify P5 and reopen the mechanism.
