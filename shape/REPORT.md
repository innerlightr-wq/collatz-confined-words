# VERDICT: **UNIVERSAL** — the global shape is a smooth function of the slope's *value*, not of its Diophantine type

Across 37 slopes spanning badly-approximable quadratics, `beta = log2 3 - 1`, `e-2`, `pi-3`,
Liouville-like slopes with partial quotients 50 and 500, **rational** slopes, and 20 random
continued fractions, every global shape statistic of the Haar-weighted staircase profile is
reproduced by a smooth function of `(theta, r)` plus a small local-phase term; after `theta` is in
the model, continued-fraction features add **1.0–1.3%** of out-of-sample RMSE, no CF feature
correlates with per-slope residuals (all `p >= 0.17`), and the one resonance effect that reached
`p = 0.007` on held-out data is an `h`-confound that vanishes under stratification (`p = 0.17`).
The sharpest available control is decisive: at matched `r`, the golden slope (hardest to
approximate) and the rationals `55/89` and `144/233` agree in `Var/r` to **five decimal places**,
and a Liouville-like slope `0.0026` away in `theta` differs by only `0.004` — exactly its
smooth-`theta` prediction.
The uniqueness of the Haar mode — the open half of Remark C.16 — holds for **every** slope tested
once `K >= 19` (0 ties in 20,327 rows), so over the tested slope range `theta in [0.1384, 0.8791]`
it is a statement about monotone Beatty staircases waiting to be proved, with no arithmetic input
(§7 states it with its scope; slopes below the tested range are not claimed); but the *location*
`mode = K(h) - 1` is **not** universal, it is specific to `beta` with window `(M,E) = (2,1)` and
the reference `h0 = 5`.

So: **STOP for the Diophantine hypothesis; GO for proving Remark C.16 in general.**

No claim here concerns the Collatz conjecture or any `O(log m)` bound. For `theta != log2 3 - 1`
these are staircase combinatorics with a different slope, not accelerated-Collatz windows.

---

## 1. Verdict justification (3 sentences)

1. With `theta` itself in the baseline, the continued-fraction block raises held-out `R^2` by
   0.002–0.018 and cuts RMSE by 1.0–1.3%, while the same block without `theta` reached `R^2`
   0.43–0.55 purely by proxying for `theta` — the Diophantine "signal" in the naive analysis was
   the slope's value in disguise.
2. At matched `r`, one cubic in `theta` (4 parameters) fits 12–15 slopes with residual sd
   0.0011–0.0092, the fit *tightening* as `r` grows, and rational slopes sit on the same curve as
   badly-approximable ones.
3. The only statistic that is slope-independent outright is the uniqueness of the mode, which holds
   for every slope at `K >= 19`; everything else varies smoothly with `theta`, which is not a
   Diophantine property.

## 2. Stage 0 — validation

| check | result |
|---|---|
| exact `floor(n*theta)` for every slope | two successive convergents give identical floors for all `n <= 460`, which **proves** the common value is `floor(n*theta)` since `n*theta` lies strictly between them. Asserted at construction for all 37 slopes |
| independent cross-check for `beta` | CF-based `G(n)` equals the bit-length formula `(3^n).bit_length()-1-n` for all `n < 450` |
| reproduce the W0 numbers (`theta = beta`, `(M,E)=(2,1)`) | `mode = K-1` untied at `h = 40,60,80,100`; `mean-K = +0.0682, +0.0116, -0.0018, -0.0061` (note: 0.068, 0.012, -0.002, -0.006); `Var/K = 1.6590, 1.7474, 1.8025, 1.8426` (note: 1.659, 1.747, 1.803, 1.843) — **exact agreement** |
| local law `T0`/`T1` for every slope | 37 slopes x 5 windows x 119 horizons: identity on 11,612/11,612 `eta=0` steps, `T1` on 10,318/10,318 `eta=1` steps, residual `-1` on 85/85 degenerate steps. **Zero exceptions** |

The last row is the "local rules are slope-neutral" half of the hypothesis, and it is confirmed
exactly: the transition rules carry no arithmetic of `theta` whatsoever (as proved in Appendix C
of the September note for an arbitrary monotone staircase).

## 3. Slope families x windows x statistics

Windows: `(M,E) in {(2,1),(3,1),(5,2),(8,1),(13,3)}`; horizons `h <= 150` (EXPLORE) and
`h <= 400` (HOLDOUT). Statistics, restricted to `K >= 20` unless stated:

* `S1 = mode - K(h)`, `S1r = mode - r` (convention-free: `K` carries an arbitrary reference `h0`)
* `S2r = mean - r`, normalized `U_loc = S2r/sqrt(r)`
* `S3r = Var/r = U_var`
* `S4m` = total variation to the **mean-matched** `NB(round(mean), 1/2)`; `U_tv = S4m*sqrt(r)`
* `S5` = Haar mass above `h0 = 5`

`S5 = 1.000000` to six places in essentially every row: the Haar mass lives entirely in the region
*above* the reference horizon, i.e. exactly the region Corollary C.12 does **not** control. That
is Remark C.16's point, measured.

Per-slope means of `U_var` (EXPLORE, `K>=20`), ordered by `theta` — the dependence is monotone and
smooth, and families interleave rather than separate:

| theta | slope | family | U_var |
|---|---|---|---|
| 0.192123 | random#03 | random | +1.9619 |
| 0.302776 | [0;3,3,3,…] | bounded | +1.9260 |
| 0.414214 | sqrt2-1 | bounded | +1.8899 |
| 0.500499 | liouville a_3=500 | liouville | +1.8165 |
| 0.563009 | random#13 | random | +1.7469 |
| 0.600781 | liouville a_5=50 | liouville | +1.6890 |
| 0.615385 | rational 8/13 | rational | +1.6743 |
| 0.615500 | liouville a_7=50 | liouville | +1.6743 |
| 0.617647 | rational 21/34 | rational | +1.6713 |
| 0.718282 | e-2 | named | +1.4607 |

Cubic in `theta` through the 18 EXPLORE slopes: `R^2 = 0.9993` (`U_var`), `0.9964` (`U_loc`);
residual sd 0.0034 and 0.0059. Full table: `out/stage2b.log`, data in `data/explore_shape2.csv`.

## 4. Baselines, EXPLORE and HOLDOUT

`B3t` = smooth in `(r, theta)`; `B1t` = `B3t` + local phase (`frac((M+h)theta)`, `frac(M theta)`,
`r mod 2`, `r mod 3`); `B2t` = `B1t` + star discrepancy; `CFt` = `B2t` + 5 CF features
(log max `a_k`, log sum `a_k`, CF depth, distance of `M+h` and of `h` to the nearest `q_k`).

**EXPLORE, in-sample `R^2`:**

| target | B3t | B1t | B2t | CFt |
|---|---|---|---|---|
| `U_loc` | 0.6590 | 0.6968 | 0.6971 | 0.7043 |
| `U_var` | 0.9060 | 0.9123 | 0.9124 | 0.9146 |
| `U_tv` | 0.7635 | 0.7772 | 0.7781 | 0.7958 |

Without `theta` in the baseline the same CF block reached `R^2` 0.4295 / 0.5481 / 0.5464 against a
`B3` of 0.0010 / 0.0012 / 0.0806 — **that entire apparent Diophantine signal was `theta` leaking
through the CF features** (CF depth and `sum a_k` both grow with `theta`). Correlation of `theta`
with per-slope residuals of the no-`theta` model: `r = -0.911, -0.920, +0.696`, all `p < 0.0005`.

**HOLDOUT (19 slopes, 16,036 rows, 2,375 at `h = 160..400`), out-of-sample:**

| target | B1t RMSE | CFt RMSE | change |
|---|---|---|---|
| `U_var` | 0.0538 | 0.0532 | **+1.0%** |
| `U_loc` | 0.0845 | 0.0834 | **+1.3%** |

CF features vs per-slope held-out residuals of `B1t` (4000-permutation):
`U_var`: log max `a_k` `r=+0.333 p=0.167`, log sum `a_k` `r=+0.191 p=0.416`, discrepancy
`r=-0.053 p=0.921`. `U_loc`: `+0.263 p=0.266`, `+0.139 p=0.559`, `+0.031 p=0.839`. Nothing.

**Matched-`r` universality test** (one cubic in `theta` per `r`-band; 15 slopes, 4 parameters,
11 degrees of freedom):

| r-band | `U_var` resid sd | `U_loc` resid sd |
|---|---|---|
| 20–30 | 0.0087 | 0.0092 |
| 30–40 | 0.0057 | 0.0068 |
| 40–55 | 0.0038 | 0.0063 |
| 55–75 | 0.0029 | 0.0047 |
| 75–105 | 0.0023 | 0.0055 |
| 105–150 | 0.0011 | 0.0023 |
| 150–260 | 0.0012 | 0.0025 |

The fit tightens monotonically with `r`: whatever slope-dependence is not captured by `theta`
shrinks as the profile grows.

## 5. Calibration and the two things the holdout flagged

**Convergent-coincidence null.** For these slopes, 0.8–6.8% of integers in `1..400` lie within
tolerance 1–2 of a convergent denominator, so "near a convergent" is not automatic here; the
resonance test below is scored against a random-denominator null of matched count and magnitude.

**Shuffled-CF control.** 16 slopes built by shuffling an EXPLORE slope's partial-quotient multiset
(same multiset, different order, hence different `theta` and different Diophantine type) deviate
from the EXPLORE `theta`-curve by mean `+0.0071`, sd `0.0219` (`U_var`) — i.e. they move *along*
the curve, they do not leave it.

**(a) The frozen tie threshold was an EXPLORE artifact.** P1 claimed no ties for `K >= 17`; the
holdout produced 28, all at `K = 17, 18`, all in window `(13,3)`, all from the two smallest-`theta`
slopes (`pi-3`, `random#19`). **P1 as frozen is FALSIFIED.** The corrected statement survives on
both halves: *no ties for `K >= 19`* — 0 of 8,885 EXPLORE and 0 of 11,442 HOLDOUT rows.

**(b) A resonance test reached `p = 0.007` and then dissolved.** On HOLDOUT, `|residual|` for rows
with `min_k |h - q_k| <= 2` was 0.0841 vs 0.0729 elsewhere. But near-convergent rows sit at mean
`h = 85.7` against `h = 138.1` for the rest, and the residual is negatively correlated with `h`
(`r = -0.323`). Within `h`-bins the difference **reverses** to `-0.0030`, and an `h`-stratified
permutation test gives `p = 0.1694`. Accounting: 4 resonance tests were run (2 statistics x 2
halves) and the EXPLORE version of this very test had the opposite sign (`p = 0.43`). It is an
`h`-confound, not a Diophantine resonance.

**P2 verdict: HOLDS.** In-range held-out slopes: `U_var` residual sd 0.0133, max 0.0460
(predicted `<= 0.02 / <= 0.08`); `U_loc` sd 0.0076, max 0.0469 (predicted `<= 0.03 / <= 0.08`).
Restricted to `h <= 150` — the range the curve was fitted on — the residuals collapse to
**mean `+0.0006`, sd `0.0043`, max `0.0084`**. The systematic positive bias in the full set is a
finite-size effect in `r`: the extended horizons reach `r <= 260` where the curve was fitted at
`r <= 108`, and the shift is shared by all slopes (common shift `-0.135`, between-slope spread
`0.045`).

**P4 verdict: HOLDS**, and is the single most decisive table in the study — `U_var` at matched `r`:

| r-band | golden (bounded) | rational 144/233 | rational 55/89 | liouville a_7=500 | spread |
|---|---|---|---|---|---|
| 20–40 | 1.50256 | 1.50256 | 1.50256 | 1.50775 | 0.0052 |
| 40–70 | 1.68362 | 1.68362 | 1.68362 | 1.68741 | 0.0038 |
| 70–110 | 1.78615 | 1.78615 | 1.78615 | 1.78991 | 0.0038 |
| 110–260 | 1.89744 | 1.89744 | 1.89742 | 1.89892 | 0.0015 |

Caveat stated plainly: golden and the two rationals have `theta` agreeing to `5.6e-5`, so their
staircases coincide over most of the tested range — their agreement to five decimals partly
reflects that, and is a consistency check rather than an independent test. The informative column
is `liouville a_7=500`, whose `theta` differs by `2.7e-3`: it differs by `0.004`, which is what the
smooth `theta`-curve predicts for that gap. Diophantine type (rational vs badly-approximable vs
Liouville-like) produces no separation at any `r`.

**P3 verdict: HOLDS** (table in §4). **P5 (extended horizons) HOLDS** for the tie statement
(0 ties in 2,375 extended rows with `K >= 17`) and for the slope-agreement statement at matched `r`;
the frozen *curve* does not extrapolate to `r > 108`, a finite-size failure common to all slopes.

## 6. Claim audit

| # | claim | label |
|---|---|---|
| 1 | `floor(n*theta)` computed here is exact for every slope and every `n <= 460` | **PROVED** (bracketing by successive convergents) + asserted in code |
| 2 | The local transition rules (`T0` identity, `T1` cumulative) carry no arithmetic of `theta` | **PROVED** (Appendix C, September note) + VERIFIED on 22,015 transitions over 37 slopes |
| 3 | The W0 shape numbers of the note are reproduced exactly | **VERIFIED** |
| 4 | The Haar mode is unique for every slope tested once `K >= 19`, equivalently `r >= 20` | **VERIFIED** (0 ties in 20,327 rows of the main study, both halves, `h <= 400`; 0 ties in 16,366 rows with `r >= 20` in the confirmatory re-test, `h <= 600`); **CONJECTURE** over the tested slope range `theta in [0.1384, 0.8791]` (§7), not claimed below it |
| 5 | `mode = K(h) - 1` | **NOT universal** — specific to `beta`, window `(2,1)`, reference `h0 = 5`. There it holds for 150/162 rows to `h = 400`, and `mode - r = -4` exactly for all 136 rows with `K >= 20`. Across held-out slopes `mode - K` ranges over `{-6,…,-1}` |
| 6 | `Var/r` and `(mean-r)/sqrt(r)` are smooth functions of `theta` at fixed `r` | **VERIFIED** (matched-`r` cubic, 11 df, residual sd 0.001–0.009, tightening with `r`) |
| 7 | Diophantine type (CF features, discrepancy, convergent resonance) adds nothing | **VERIFIED** on held-out data: +1.0%/+1.3% RMSE, all correlations `p >= 0.17`, resonance `p = 0.17` stratified |
| 8 | The apparent CF signal in the naive model was `theta` leaking through CF features | **VERIFIED** (`R^2` 0.43–0.55 collapses to +0.002–0.018 once `theta` is included) |
| 9 | P1 as frozen (`K >= 17`) | **FALSIFIED** on holdout; corrected to `K >= 19` |
| 10 | The `p = 0.007` resonance | **ARTIFACT** (`h`-confound; stratified `p = 0.169`, opposite sign on EXPLORE) |
| 11 | Anything about Collatz or `O(log m)` | **NOT CLAIMED** |

## 7. Next step — the statement to prove

The Diophantine hypothesis is closed for these statistics, and what is left is a clean
combinatorial problem with no arithmetic in it. Stated over the slope range actually tested
(see Scope below), it is the open half of Remark C.16:

> **Conjecture (unique Haar-weighted maximum, tested range).** [OBSERVED / CONJECTURE — not
> proved.] Let `theta_min = 0.138380349690391...` be the smallest slope tested. For every
> `theta` in `[theta_min, 1)`, rational or irrational, every window `(M,E)` among those tested,
> and the staircase
> `ell_0 = r`, `ell_j = min(ell_0, floor((M+h)theta) - floor((M+j)theta))`,
> let `Q_L(j)` be the number of chains `x_0 >= ... >= x_j >= 0` with `x_i <= ell_i - 1`.
> Then once `r >= 20` the sequence `w(j) = Q_L(j) 2^{-j}` has a **unique** maximum.
> In particular this covers `theta = beta = log2(3) - 1 = 0.5849625...`, the Collatz case.

**Scope.** The evidence is 46,911 rows over **47 distinct slopes** with
`theta in [0.1383803, 0.8790732]`, **8 windows** `(M,E) in {(2,1), (3,1), (5,2), (8,1), (13,3)}`
(main study) and `{(4,2), (7,1), (11,2)}` (confirmatory re-test), `r` from 1 to 527 and horizons
`h` from 1 to 600. Of these, **36,863 rows have `r >= 20` and none has a tie**. The threshold
`r >= 20` was frozen in this section before the confirmatory run of §7b and held there with
**0 ties** in 16,366 rows. The largest `r` at which any tie was observed is **17**, in the main
study (slopes `random#19` and `pi-3`); in the confirmatory run alone the largest was **12**. Ties
concentrate at the smallest slopes tested — the two slopes attaining a tie at `r = 17` are the two
smallest, `theta = 0.1384` and `theta = 0.1416` — so the threshold may well grow as `theta -> 0`,
and **slopes below `theta_min` are not claimed and are left open**. The *location* of the maximum
is **not** claimed: it varies with slope and window (`mode - r` takes values in `{-6,…,-3}`). For
`theta = beta` and window `(2,1)` the observed location is `mode = r - 4` (all 136 rows with
`K >= 20`, `h <= 400`), which is the `mode = K(h) - 1` of the September note under that note's
`h0 = 5` convention.

Two further observations a proof should explain, both VERIFIED but unproved:

* the ratio `w(j+1)/w(j) = Q_L(j+1)/(2 Q_L(j))` appears to cross 1 exactly once — **log-concavity
  of `j -> Q_L(j)`** remains the natural lemma to target, a statement about staircase-bounded chain
  counts with no arithmetic in it. Note that log-concavity gives *unimodality* but not uniqueness:
  a separate argument is needed to exclude a two-way tie `w(j) = w(j+1)` at the summit, which is
  exactly the behaviour observed at small `r`;
* the location of the maximum is *not* slope-free (see Scope), so a proof of uniqueness should not
  attempt to pin the location to a constant.

Log-concavity of `Q_L`, plus the tie-exclusion argument, would give the single-mode statement for
the September note's Appendix C profile as the special case `theta = log2 3 - 1`, closing
Remark C.16 over this range.

## 7b. Confirmatory re-test of the uniqueness threshold (new data, thresholds fixed in advance)

Run after §7 was written and committed, on data that played no part in setting the threshold:
**10 new random-CF slopes** (seed 20261006; none shares a continued fraction with the earlier 37),
**3 new windows** `(M,E) in {(4,2),(7,1),(11,2)}`, horizons `h = 1..600`. 17,920 rows, `r` reaching
527. Exactness of `floor(n*theta)` asserted to `n <= 620` for every slope.

Profiles were advanced with the proved operators (`eta = 0` identity, `eta = 1` `T1`) rather than
recomputed from scratch, and spot-checked against a direct chain-count recomputation at 210
horizons: **0 failures** — which incidentally re-verifies the local law on 10 new slopes and 3 new
windows.

**Ties, reported in `r`:**

| window `(M,E)` | rows | `r` range | #ties | largest `r` with a tie | tie-free from |
|---|---|---|---|---|---|
| (4,2) | 5,967 | 1–526 | 98 | **9** | `r >= 10` |
| (7,1) | 5,989 | 1–527 | 42 | **6** | `r >= 7` |
| (11,2) | 5,964 | 1–527 | 120 | **12** | `r >= 13` |
| **all** | **17,920** | 1–527 | 260 | **12** | `r >= 13` |

The same ties in `K`: largest `K` with a tie is 9, 5, 13 by window; **13** overall.

**Thresholds under test (set before this run, not adjusted after it):**

| threshold | rows | ties | verdict |
|---|---|---|---|
| `r >= 20` (the conjecture of §7) | 16,366 | **0** | **HOLDS** |
| `K >= 19` (observed in the earlier study) | 16,353 | **0** | **HOLDS** |

Tie frequency pooled over the three windows: 23.1% at `r in [1,5)`, 41.3% at `[5,10)`, 3.9% at
`[10,15)`, and **0% at every `r >= 15`** across 15,673 rows. The largest tie seen anywhere in this
run is at `r = 12`, so the pre-set threshold `r >= 20` clears this run's data by 7 in `r`; across
all three runs the largest tie is at `r = 17` (main study), so the margin over the whole study is
3. The threshold is left where it was.

Ties concentrate in the small-`theta` slopes — `new-random#07` (`theta = 0.1848`, 86 ties, largest
at `r = 12`) and `new-random#00` (`theta = 0.1751`, 79 ties) — against 3–6 ties each for the slopes
with `theta > 0.55`. This matches the earlier holdout, where the only ties above `K = 17` came from
`pi-3` and `random#19`, the two smallest-`theta` slopes. The largest tie `r` also varies with the
window (6, 9, 12 for the three tested), so the threshold is not a single clean constant across
windows at this resolution; `r >= 20` covers all of them.

*(Labels: the tie counts and thresholds here are VERIFIED on this data; the conjecture of §7
remains CONJECTURE.)*

## 8. Files

```
shape/code/slopes.py              CF slopes, exact floor(n*theta), split, Diophantine features
shape/code/shape_core.py          staircase, chain counts, Haar statistics (exact)
shape/code/stage0_shape.py        0.1-0.2  W0 reproduction; local law over all 37 slopes
shape/code/stage1_sweep.py        1.1-1.5  EXPLORE sweep
shape/code/stage1b_analysis.py    1.6-1.9  convention-free statistics, K-collapse, tie anatomy
shape/code/stage2_baselines.py    2.1-2.5  B3/B1/B2/CF without theta  (the confound)
shape/code/stage2b_theta.py       2.6-2.9  the same with theta  (the resolution)
shape/code/stage3_calibration.py  3.1-3.3  convergent null, shuffled-CF control, resonance
shape/code/freeze_predictions.py  fits and freezes the theta-curve
shape/code/stage4_holdout.py      4.1-4.6  frozen P1-P5 on HOLDOUT + extended horizons
shape/code/stage5_followup.py     5.1-5.5  tie threshold; the p=0.007 resonance dissected
shape/code/stage6_finitesize.py   6.1-6.3  finite-size law in r
shape/code/stage7_matched_r.py    7.1-7.2  matched-r universality; the four-slope control
```
Raw outputs `shape/out/*.log`, data `shape/data/*.csv`, pre-registration `shape/PREDICTION.md`
(committed 2026-10-05 12:48:49 -0400, before any holdout computation).
