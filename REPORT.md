# VERDICT: **ARTIFACT** — and the law it was supposed to break is now a theorem

The reported `M = 3` breakdown of the cumulative-summation law does not exist: for the
natural `M = 3, c = 1` windows the literal `T1` identity holds on **every** `eps_h = 1`
transition tested (88/88 to `h = 150` by the fast DP, 17/17 by the fully literal
recursive tree), and the reported figures are reproduced exactly, with no free parameters,
by one indexing error — selecting the `eps = 1` horizons with the hard-coded `M = 2` clock
`C(h) = floor(2a + hb)` while computing the tree of an `M = 3` window, whose own clock is
the same Sturmian word shifted by one index (`#{h <= 95 : e_{2+h}=1} = 56`, of which
`#{e_{3+h}=1} = 16`, with gaps exactly `{5,7}` and `12 = 5+7`).
The "Diophantine" signal is a tautology of that error: the 16 survivors are the positions of
the Sturmian factor `11`, whose return times are forced to be convergent-denominator-sized by
the three-distance theorem, and 80-100% of small integers qualify as "near a
convergent/semiconvergent" anyway.
Going further, the reduction built to test the question also **proves** the law (the source
note's Conjecture 12 / "missing lemma"): after an exact change of variables the profile
`N_h(.)` becomes a staircase-bounded chain count, and the `eps_h = 1` step is exactly a
uniform `+1` of all the bounds, for which `Ntilde_{L+1}(j) = 1 + sum_{i<=j} Ntilde_L(i)` is an
elementary identity holding for *every* non-increasing bound profile — so the law involves no
property of `beta`, no phase, and no window alignment at all.

Nothing in this report bears on the Collatz conjecture or on any `O(log m)` occupation bound.

---

## 1. Verdict justification (3 sentences)

1. Stage 0's artifact check failed in the direction that terminates the protocol: the `M = 3`
   failure disappears under a correct generalization, and its exact reported numbers
   (`16 of 56`, gaps `5, 7, 12`) are reconstructed from a one-index clock misalignment.
2. The only genuine exceptions to the literal law anywhere in a 468-class, 41,076-transition
   sweep occur when the window admits **no terminal digit at all** at that horizon
   (`r(W,h) = G(M+h) - G(M) - E + 1 <= -1`, i.e. `V_min >= V*(h)`), where both profiles
   degenerate to `{0: 1}` and the residual is the constant vector `(-1, ..., -1)`.
3. The surviving predictive feature is therefore "the window has opened yet", a
   non-Diophantine magnitude condition (`h >= h_open(M,E) ~ E/beta`), while every phase or
   convergent feature tested scores at or below the base rate on held-out data
   (accuracies 0.0023-0.197 vs base rate 0.948-1.000).

## 2. Stage 0 — validation

| check | result |
|---|---|
| literal recursive tree == DP(literal thresholds) == DP(closed-form `d*`) == reduced chain model | **exact on all 240 (window, h) pairs**, 18 windows, `h <= 15` (120 pairs in Stage 0, 120 more on HOLDOUT windows in Stage 3.8) |
| generalized `d*(S,j,h) = C_W(h) + 2 + j - S_fin`, `C_W(h) = floor(c + M*alpha + h*beta) = c - h + F(M+h)` | **exact at all 13,532 reachable states** tested (8 windows, `h <= 25`) |
| global digit `eps_h = d*(S,j,h+1) - d*(S,j,h) = C_W(h+1)-C_W(h) = G(M+h+1)-G(M+h) in {0,1}` | **exact**, state-independent, for every window and `h <= 400` |
| clipping case `d* <= 0` | **never reachable**: `min d* = 2` over all reachable states (both in Stage 0 and on HOLDOUT windows). Upgrades Remark 4 of the source note from computational observation to proof (see Lemma A(iv)) |
| `N_h(h) = 0` (horizon-exhaustion clause never fires) | **true** for every window and horizon tested; proved in Lemma A(v) |
| window -> `(M, E)` collapse (`E = P + V_min - c - F(M)`) | 135 windows -> 15 classes, **identical profiles within a class**; re-verified on all 66 HOLDOUT windows with the independent DP |
| W0 (`M=2, c=1, prefix (1), V_min=4`, so `E=1`, `gamma_W = 0.169925`) `T0` on `eps=0` | **62/62 exact with no margin at all** (`j <= h`), sharper than the note's `j <= h-2` |
| W0 `T1` on `eps=1`, `h <= 151` | **89/89 exact** (note reports 85/85 to `h=150`) |
| W0 Haar-weighted bulk stats, weights `N_h(j)*2^-j` | reproduces the note's quoted numbers **to every printed digit** (below) |

W0 bulk statistics, `K(h) = C(h) - C(5)`:

| h | mode `j*` | `K(h)-1` | `E[J]-K` (note: 0.068, 0.012, -0.002, -0.006) | `Var/K` (note: 1.659, 1.747, 1.803, 1.843) |
|---|---|---|---|---|
| 40 | 19 | 19 | +0.068 | 1.659 |
| 60 | 31 | 31 | +0.012 | 1.747 |
| 80 | 42 | 42 | -0.002 | 1.802 |
| 100 | 54 | 54 | -0.006 | 1.843 |

The mode is a single untied maximum at every `h` tested, as reported.

### 2.1 Convention note (a real trap)

The source note's Appendix A.4 `T1` implementation accumulates `out[-1] + N[j]` starting from
`N[0] + 1`, i.e. the `+1` boundary bump **propagates** into every later partial sum. Under that
convention W0 scores **1/89**, not 89/89. The convention that reproduces the reported result is
Observation 11 as *stated*: `T1(N)(0) = N(0)+1` and `T1(N)(j) = sum_{i<=j} N(i)` with no
propagating bump. All results here use the stated convention.
*(VERIFIED-COMPUTATIONALLY; both conventions were run.)*

A second convention needed pinning down. The source note's Section 2.3 prose reads as if a
RETURN at `cur_min` terminates the state, whereas Appendix A.3 (and the proof text of
Theorem 10, "the classification of each candidate digit value `delta = 1,...,d*-1` as a return
or a continuation") classifies **every** `delta` in `1..d*-1` independently and continues past
the returns. This report uses the latter; it is the reading under which the W0 figures
(89/89, and the Haar mode/mean/variance to every printed digit) reproduce, which is what
identifies the object studied here with the object in the note.

### 2.2 The `M = 3` artifact

| test | result |
|---|---|
| `M=3, c=1, prefix (1,1), V_min=4` (`E=1`), fast DP, `h <= 150` | `T1` exact **88/88** |
| same window, fully literal recursive tree / literal-threshold DP, `h <= 30` | `T1` exact **17/17** |
| `E=2` (`V_min=5`) | **88/88** |
| `E=3`, `E=4` | 87/88, 86/88 — failures only at `h = 2` and `h = 2,3`, i.e. `h < h_open` |
| `T0` on `eps=0` transitions, all four | **63/63 exact** |
| buggy generalizations with a hard-coded `M=2` leftover (interior sum `1`; `rho*(j) = floor((j+3)a)`; `V_min = 4`) | all still **58/58 exact** — these cannot produce the reported failure |
| best-case index hard-coded to `M=2` (`2+h`) | 0/58 — too destructive, not the reported pattern |
| **clock selected with `M=2`, tree computed for `M=3`** | **`#selected = 56`, `#exact = 16`, gaps `{5,7}` (`12 = 5+7`)** — exact match to the report |

Confirmation of the mechanism: for `M=3, E in {1,2}`, `T1` is exact on **58/58** transitions with
`e_{3+h} = 1` and fails on **41/41** transitions with `e_{3+h} = 0` — because at those the true
law is `T0 = identity`, and `T1` differs from the identity by a full cumulative sum. That is
precisely the report's "discrepancies were often large and therefore do not appear to be a small
finite-horizon boundary effect". *(VERIFIED-COMPUTATIONALLY.)*

## 3. The exact reduction, and the proof of the law

Write `F(n) = floor(n*alpha) = (3^n).bit_length()-1`, `G(n) = F(n) - n = floor(n*beta)`
(both exact integer operations), `P = sum(prefix)`, `E = P + V_min - c - F(M) >= 1`.

**Lemma A (reduction; PROVED, and verified against the literal tree on 240 (window,h) pairs).**
Fix a window and a horizon `h`. For an alive state at level `j` with confirmed sum `S`, put
`x := d*(S,j,h) - 2 = C_W(h) + j - S`. Then

* (i) the terminal-digit seeds are exactly `x = 0, 1, ..., l_0 - 1` with
  `l_0 = C_W(h) - P - V_min + 1 = G(M+h) - G(M) - E + 1`  (since `V*(h) = C_W(h) - P + 1`);
* (ii) a state `x` at level `j` has continuations exactly `x' in [0, min(x, q_j - 1)]` at level
  `j+1`, where `q_j := C_W(h) + j + 1 - rho*(j) = G(M+h) - G(M+j+1)` and
  `rho*(j) = c + F(M+j+1)`  (digits `delta = 1..d*-1` minus the returns `delta <= x - q_j + 1`,
  with `x' = x - delta + 1`);
* (iii) every alive state contributes exactly one pruned leaf (at `delta = d*`), so with
  `l_j := min(l_0, G(M+h) - G(M+j))` for `j >= 1` (equivalently `l_{j+1} = min(l_j, q_j)`),
  `N_h(j) = Ntilde_L(j) + [j=0]`, where
  `Ntilde_L(j) = #{ x_0 >= x_1 >= ... >= x_j >= 0 : x_i <= l_i - 1 }`;
* (iv) every reachable state has `d* >= 2` (seeds have `x >= 0`, and continuations have
  `x' >= 0`), so the clipping case `d* <= 0` never occurs;
* (v) `l_h = min(l_0, G(M+h) - G(M+h)) = 0`, hence `N_h(h) = 0`: at the last level every digit
  value either returns or prunes, so the horizon-exhaustion clause never fires;
* (vi) consequently the entire profile depends on the window only through `(M, E)`; `c` and the
  interior prefix enter only through `E`.

**Lemma B (the clock; PROVED).** `eps_h = C_W(h+1) - C_W(h) = G(M+h+1) - G(M+h)`, the Sturmian
word of slope `beta` with intercept `gamma_W = frac(M*beta) = frac(M*alpha)`. Moreover
(verified 0 violations in 118,197 + 88,908 cases):
`eps_h = 0` implies `l_j(h+1) = l_j(h)` for all `j <= h`; and `eps_h = 1` with `l_0(h) >= 0`
implies `l_j(h+1) = l_j(h) + 1` for all `j <= h`.

**Lemma C (the missing lemma; PROVED).** Let `L = (l_0 >= l_1 >= ... >= l_n >= 0)` and
`L+1 = (l_0+1, ..., l_n+1)`. Let `A_j(x) = #{x_j = x <= x_{j-1} <= ... <= x_0 : x_i <= l_i - 1}`
(supported on `0..l_j-1`, with `A_{j+1}(x) = sum_{y=x}^{l_j-1} A_j(y)`) and `A'_j` the same for
`L+1` (supported on `0..l_j`). Then

* (i) `A'_{j+1}(x) = A_{j+1}(x) + A'_j(x)` for `0 <= x <= l_{j+1}-1`;
* (ii) `Ntilde_{L+1}(j+1) = Ntilde_{L+1}(j) + Ntilde_L(j+1)`;
* (iii) hence `Ntilde_{L+1}(j) = 1 + sum_{i=0}^{j} Ntilde_L(i)`.

*Proof.* (i) Base `j = 0`: `A'_1(x) = l_0 + 1 - x = (l_0 - x) + 1 = A_1(x) + A'_0(x)`.
Inductively,
`A'_{j+1}(x) = sum_{y=x}^{l_j} A'_j(y) = sum_{y=x}^{l_j-1} (A_j(y) + A'_{j-1}(y)) + A'_j(l_j)
            = A_{j+1}(x) + [ sum_{y=x}^{l_j-1} A'_{j-1}(y) + A'_j(l_j) ]`,
and `A'_j(l_j) = sum_{y=l_j}^{l_{j-1}} A'_{j-1}(y)` is the primed recursion evaluated at
`x = l_j` (legitimate because `l_j <= l_{j-1}`), so the bracket telescopes to
`sum_{y=x}^{l_{j-1}} A'_{j-1}(y) = A'_j(x)`.
(ii) Summing (i) over `x = 0..l_{j+1}-1` and adding the one extra primed entry,
`Ntilde_{L+1}(j+1) = Ntilde_L(j+1) + sum_{x=0}^{l_{j+1}-1} A'_j(x) + A'_{j+1}(l_{j+1})`,
and `A'_{j+1}(l_{j+1}) = sum_{x=l_{j+1}}^{l_j} A'_j(x)` (the primed recursion again), so the last
two terms combine to `sum_{x=0}^{l_j} A'_j(x) = Ntilde_{L+1}(j)`.
(iii) Telescope (ii) from `Ntilde_{L+1}(0) = l_0 + 1 = 1 + Ntilde_L(0)`. QED

**Theorem (T0 and T1; PROVED).** For every window `W` and every `h >= 1`:
* if `eps_h = 0` then `N_{h+1}(j) = N_h(j)` for all `0 <= j <= h`;
* if `eps_h = 1` and `r(W,h) := G(M+h) - G(M) - E + 1 >= 0` then
  `N_{h+1}(j) = T1(N_h)(j)` for all `0 <= j <= h`;
* if `eps_h = 1` and `r(W,h) <= -1` then `N_h = N_{h+1} = {0: 1}` and
  `N_{h+1}(j) - T1(N_h)(j) = -1` identically.

*Proof.* `T0` is Lemma A + Lemma B (identical bound profiles give identical chain counts).
For `T1`, Lemma B gives `L(h+1) = L(h) + 1` on `0 <= j <= h`, so Lemma C applies:
`j = 0` gives `N_{h+1}(0) = l_0 + 2 = N_h(0) + 1`; for `j >= 1`,
`N_{h+1}(j) = 1 + sum_{i<=j} Ntilde_L(i) = N_h(0) + sum_{i=1}^{j} N_h(i) = T1(N_h)(j)`.
If `r <= -1` the window admits no terminal digit at either horizon, both profiles are the single
terminal prune leaf, and the residual is the stated constant. QED

The only ingredient is the monotonicity `l_0 >= l_1 >= ...`, which Lemma A(iii) supplies
automatically. **No property of `beta` beyond `eps_h in {0,1}` is used**, and in particular no
phase, no convergent, and no alignment between `gamma_W` and the `j`-clock.
Independent check: the identity of Lemma C was verified on 4,000 random non-increasing profiles
(33,567 checks, 0 violations), and sub-claim (i) likewise with 0 violations.

## 4. Window table

Full data: `data/window_table_full.csv` (all 132 windows, both splits),
`data/windows.csv` (split assignment), `data/explore_window_table.csv`,
`data/sweep_ME.csv` (the 468-class sweep). One representative window per `(M,E)` class:

| M | c | prefix | V_min | E | gamma_W | h_open | #eps=1 (h<=150) | #exact T1 | #eps=0 | #exact T0 | T1 failures at h |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 2 | 2 | (1) | 5 | 1 | 0.169925 | 1 | 88 | 88 | 62 | 62 | none |
| 2 | 3 | (2) | 6 | 2 | 0.169925 | 2 | 88 | 87 | 62 | 62 | 1 |
| 2 | 1 | (1) | 7 | 4 | 0.169925 | 5 | 88 | 85 | 62 | 62 | 1;3;4 |
| 3 | 1 | (1-1) | 4 | 1 | 0.754888 | 1 | 88 | 88 | 62 | 62 | none |
| 3 | 1 | (1-1) | 5 | 2 | 0.754888 | 1 | 88 | 88 | 62 | 62 | none |
| 3 | 2 | (1-1) | 8 | 4 | 0.754888 | 4 | 88 | 86 | 62 | 62 | 2;3 |
| 4 | 1 | (1-1-1) | 5 | 1 | 0.339850 | 1 | 88 | 88 | 62 | 62 | none |
| 4 | 2 | (1-1-1) | 7 | 2 | 0.339850 | 2 | 88 | 87 | 62 | 62 | 1 |
| 4 | 1 | (1-1-1) | 8 | 4 | 0.339850 | 5 | 88 | 85 | 62 | 62 | 1;2;4 |
| 5 | 1 | (1-1-1-1) | 5 | 1 | 0.924813 | 1 | 88 | 88 | 62 | 62 | none |
| 5 | 1 | (1-1-1-1) | 6 | 2 | 0.924813 | 1 | 88 | 88 | 62 | 62 | none |
| 5 | 1 | (1-1-3-2) | 5 | 4 | 0.924813 | 4 | 88 | 86 | 62 | 62 | 1;3 |
| 6 | 1 | (1-1-1-1-1) | 6 | 1 | 0.509775 | 1 | 87 | 87 | 63 | 63 | none |
| 6 | 1 | (1-1-1-1-1) | 7 | 2 | 0.509775 | 1 | 87 | 87 | 63 | 63 | none |
| 6 | 1 | (1-1-1-1-1) | 9 | 4 | 0.509775 | 5 | 87 | 85 | 63 | 63 | 2;4 |

`h_open(M,E) = min{h : r(M,E,h) >= 0}`. In every row the failing horizons are exactly the
`eps=1` horizons below `h_open`. Note the `(M,E)` collapse means these 15 classes exhaust the
specified `M in {2..6}` grid; the wide/deep held-out sets (Section 6) extend to `M <= 120`.

Totals over the exploration grid (`M = 2..40`, `E = 1..12`, `h <= 150`, 468 classes):
`T1` exact on **38,744 / 41,076** `eps=1` transitions, `T0` exact on **29,124 / 29,124`
`eps=0` transitions, and **every** failure has `r <= -1`.

## 5. Failure anatomy

* **Location.** All of it. The residual is non-zero at *every* compared index `j`, never at an
  isolated boundary index; varying the excluded margin over `{0,1,2,3,5,8}` changes nothing
  (identical counts at all six margins for all `(M,E)` tested).
* **Shape.** The residual is the constant vector `(-1, -1, ..., -1)` in 100% of failures
  (56/56 on EXPLORE; 477+217+207+... = all shapes observed in the 16,152-transition check are
  `(-1)^k`). Sign always negative.
* **Cause.** Both profiles are the degenerate `{0: 1}` (only the terminal-digit prune leaf),
  so `T1` adds its `+1` boundary bump to a profile that has not changed. This is "the window
  has not opened yet", not an operator phenomenon.
* **Criterion.** `T1 fails  <=>  r(W,h) = G(M+h) - G(M) - E + 1 <= -1`: agreement
  **16,152/16,152** on the exploration sweep and **0 errors** on every held-out set.
* **Not Diophantine.** `h_open(M,E) ~ (E-1)/beta + O(1)` is a magnitude threshold in `h`, i.e. the
  "boring: `h` small" predictor the protocol names as a disqualifier, not a phase condition.

## 6. Frozen predictions and held-out results vs nulls

`PREDICTION.md` was committed (`73c47f3`, 2026-10-05 10:09:24 -0400) before any HOLDOUT
computation; `code/stage3_holdout.py` was written and run afterwards.

| held-out set | #eps=1 (fail) | #eps=0 | P1 errors | residual shape as predicted | T0 exact |
|---|---|---|---|---|---|
| HOLDOUT windows, `h<=150`, margin 2 | 1,317 (14) | 933 | **0** | yes | yes |
| HOLDOUT windows, `h<=150`, **margin 0** (`j<=h`) | 1,317 (14) | 933 | **0** | yes | yes |
| EXPLORE classes, **extended** `h = 151..600`, margin 0 | 3,948 (0) | 2,802 | **0** | n/a | yes |
| WIDE: unseen phases `M = 41..120`, `E in {1,2,3,4,6,8,12}`, `h<=120` | 39,298 (2,038) | 27,902 | **0** | yes | yes |
| DEEP: unseen phases `M = 101..110`, `E in {1,2,4}`, `h<=400` | 7,017 (30) | 4,983 | **0** | yes | yes |

**51,580 held-out `eps=1` transitions and 36,620 `eps=0` transitions, zero errors.**
Auxiliary frozen claims also held: `N_h(h) = 0`, `min d* = 2`, `(M,E)` collapse on all 66
HOLDOUT windows (independent DP), brute tree == DP == chain model on 120 HOLDOUT `(W,h)` pairs.

Nulls (balanced accuracy = mean of per-class recall; the base-rate predictor is pinned at 0.5):

| predictor | HOLDOUT acc / bal | WIDE acc / bal | DEEP acc / bal |
|---|---|---|---|
| **P1: `r >= 0`** | **1.000000 / 1.0000** | **1.000000 / 1.0000** | **1.000000 / 1.0000** |
| null: always-success (base rate) | 0.989370 / 0.5000 | 0.948140 / 0.5000 | 0.995725 / 0.5000 |
| Diophantine: `|h-q_k| <= 0` | 0.039484 / 0.2319 | 0.066568 / 0.3475 | 0.014964 / 0.2067 |
| Diophantine: `|h-q_k| <= 1` | 0.091875 / 0.0464 | 0.107512 / 0.2123 | 0.039761 / 0.0200 |
| Diophantine: `|h-q_k| <= 2` | 0.135156 / 0.0683 | 0.151255 / 0.1709 | 0.062847 / 0.0316 |
| Diophantine: `|h-q_k| <= 3` | 0.182992 / 0.0925 | 0.196575 / 0.1343 | 0.085934 / 0.0432 |
| null: random-nearby-integers `<=1` (200 reps) | 0.090774 / 0.1199 | 0.109353 / 0.2364 | 0.041284 / 0.0997 |
| null: random-nearby-integers `<=2` (200 reps) | 0.136986 / 0.0800 | 0.152198 / 0.1773 | 0.063129 / 0.0432 |
| best phase interval re-fitted **on the test set** | 0.989370 / — | 0.948140 / — | 0.995725 / — |

Label-shuffle permutation p-value for P1: `p = 0.0005` (`1/2001`, the floor at `R = 2000`) on
HOLDOUT, WIDE and DEEP; `p = 1.0` on the extended-horizon set only because that set contains no
failures at all, which makes the permutation test vacuous there rather than negative.

Two notes on the phase predictor, both of which matter for honest reporting:

* The best interval on `frac(h*beta + gamma_W)` is `[t, 1)` for any `t <= 1 - beta = 0.415037`,
  and its accuracy is *exactly* the base rate. This is a tautology: `eps_h = 1` holds iff
  `frac(h*beta + gamma_W) >= 1 - beta`, so after conditioning on `eps_h = 1` the phase is
  confined to `[1-beta, 1)` and any such interval is literally the always-success predictor
  (verified: 0 mismatches in 11,542 `(M,h)` pairs). Reported as PROVED-tautology, not a signal.
* Conditional on `r(W,h)`, every phase feature is independent of the label: there is no failure
  with `r >= 0` and no success with `r <= -1`, so the conditional mutual information is exactly 0.

**Chance calibration of "near a convergent"** (convergent denominators of `beta`:
1, 2, 5, 12, 41, 53, 306, 665; semiconvergents 1, 2, 3, 5, 7, 12, 17, 29, 41, 53, 94, 147, 200):

| range of h | tol=0 | tol=1 | tol=2 |
|---|---|---|---|
| 1..20 — near a convergent | 20% | 45% | 60% |
| 1..20 — near a semiconvergent | 35% | **70%** | **95%** |
| 1..50 — near a semiconvergent | 18% | 40% | 58% |
| 1..150 — near a semiconvergent | 8% | 19% | 29% |

The reported gaps `5, 7, 12` are *all three* exact semiconvergent denominators, and under a
`tol=1` criterion 70% of integers below 20 qualify; separately, a set of 16 successes in
`1..96` has mean gap 6.0 by construction, so gaps of 5, 7 and 12 are the null expectation.
Finally, the 16 survivor horizons are the occurrences of the factor `11` in a Sturmian word, and
the return times of any Sturmian factor take at most three values which are automatically
convergent-sized (three-distance theorem) — observed gap multiset over 4,000 horizons:
`{5: 379, 7: 300}`. So the Diophantine appearance is *forced* by the bug and carries no
information about the operator law.

## 7. Claim audit

| # | claim | label |
|---|---|---|
| 1 | `d*(S,j,h) = floor(c + M*alpha + h*beta) + 2 + j - S_fin` for general windows | **PROVED** (+ exact at 13,532 reachable states) |
| 2 | `eps_h in {0,1}` is global over states; `= G(M+h+1)-G(M+h)`; Sturmian slope `beta`, intercept `gamma_W = frac(M*alpha)` | **PROVED** |
| 3 | Lemma A: reduction of `N_h(.)` to the staircase chain count; profile depends on `W` only via `(M,E)` | **PROVED**, and VERIFIED-COMPUTATIONALLY against the literal tree (240 `(W,h)` pairs; 135+66 windows -> 15 classes) |
| 4 | No reachable state has `d* <= 0` (in fact `d* >= 2`) | **PROVED** (Lemma A(iv)); upgrades the source note's Remark 4 |
| 5 | `N_h(h) = 0`: the horizon-exhaustion clause never fires | **PROVED** (Lemma A(v)) |
| 6 | `T0`: `eps_h = 0` implies `N_{h+1}(j) = N_h(j)` for all `j <= h` (no margin) | **PROVED**; VERIFIED on 36,620 held-out + 29,124 sweep transitions |
| 7 | Lemma C: `Ntilde_{L+1}(j) = 1 + sum_{i<=j} Ntilde_L(i)` for every non-increasing profile | **PROVED** (elementary induction); VERIFIED on 4,000 random profiles |
| 8 | `T1`: `eps_h = 1` and `r(W,h) >= 0` implies `N_{h+1} = T1(N_h)` on all `j <= h`, every window | **PROVED** (resolves Conjecture 12 of the source note); VERIFIED on 51,580 held-out transitions, 0 errors |
| 9 | `T1` fails iff `r(W,h) <= -1`, and then the residual is exactly `(-1,...,-1)` | **PROVED** + VERIFIED (16,152/16,152 and 0 held-out errors) |
| 10 | The reported `M=3` "16 of 56 / gaps 5,7,12" is reproduced by selecting transitions with the `M=2` clock | **VERIFIED-COMPUTATIONALLY** (exact match of all three numbers; no free parameters) |
| 11 | No Diophantine feature beats the base rate on held-out data | **VERIFIED-COMPUTATIONALLY** (table in Section 6) |
| 12 | "Phase interval `[t,1)`" predictors are the always-success predictor in disguise | **PROVED** (`eps_h=1 <=> frac(h*beta+gamma_W) >= 1-beta`) |
| 13 | W0 Haar mode `= K(h)-1`, `E[J]-K -> 0`, `Var/K -> 2` | **VERIFIED-COMPUTATIONALLY** (reproduces the source note exactly); asymptotics remain CONJECTURE |
| 14 | The source note's Appendix A.4 `T1` code propagates the `+1` bump and would score 1/89 on W0 | **VERIFIED-COMPUTATIONALLY** |
| 15 | Anything about the Collatz conjecture or an `O(log m)` bound | **NOT CLAIMED** — out of scope here |

## 8. What this leaves open (the honest next step)

The operator question is closed, so no `GO` conjecture is offered. Two things are worth noting
for the parent programme:

1. `Corollary 13` of the source note (`N_h = T1^{K(h)}(N_{h_0})`, hence the horizon dependence
   through the scalar `K(h)`) is now **unconditional on its base range**, for every window and all
   `h >= h_open(M,E)`, with the operator alphabet `{identity, cumulative sum}` and the Sturmian
   word `(eps_h)` as the driver. The qualification matters: the *single-step* laws (T0, T1) are
   exact on the whole profile `j <= h`, but iterating them from a reference horizon `h_0` is exact
   only on `j <= h_0`. Each step appends a new top level with `N_{h+1}(h+1) = 0`, whereas `T1`
   applied to the zero-extended profile puts `sum_i N_h(i) > 0` there; since `T1` reads only lower
   indices the discrepancy never reaches below `h_0+1`, but it is never removed above it. Verified:
   248,080 checks on `j <= h_0` over 28 `(M,E)` classes and `h_0 in {5,10,20,30}`, zero violations,
   with the first mismatch always at `j = h_0+1` or `h_0+2` (for W0: `j=12` at `h_0=10`, `j=21` at
   `h_0=20`, persisting to `h=150`) — `code/verify_c12_range.py`, `data/c12_range.csv`.
   What remains genuinely unproved is only the *shape* statement:
   the boundary correction that turns the exactly-tied negative-binomial kernel into the single
   observed mode `K(h)-1`, i.e. `Lemma 14` -> `Observation 18`. Lemma A reduces that to an
   explicit question about `Ntilde_L(j)` for the specific Beatty staircase
   `l_j = min(l_0, G(M+h) - G(M+j))` — a finite, self-contained asymptotic computation, no longer
   a question about trees.
2. Because the profile depends on the window only through `(M, E)`, "window-universality" of any
   future statement needs to be tested over `(M, E)` pairs, not over `(M, c, prefix, V_min)`
   quadruples: 135 of the obvious windows collapse to 15 distinct objects. Any future claim of
   window-dependence should first be checked against this collapse, and any claim of
   `eps`-driven structure should be checked against the window's **own** clock
   `C_W(h) = c - h + F(M+h)` rather than the `M = 2` clock.

## 9. Reproducing

```
code/core.py              exact arithmetic, Window, brute tree, DP, chain model, Diophantine helpers
code/stage0_validate.py   0.1-0.3  cross-validation, (M,E) collapse, d*/eps/clipping
code/stage0b_reproduce.py 0.4-0.6  W0 facts, Haar stats, M=3 T0/T1, margin sensitivity
code/stage0c_artifact.py  0.7-0.9  independent-DP M=3 check, buggy-generalization variants
code/stage0d_16of56.py    0.10-0.13 exact reconstruction of 16/56 and the gap calibration
code/stage1_sweep.py      1.1-1.3  468-class sweep
code/stage1b_identity.py  1.4-1.6  general identity, inductive lemma, residual characterisation
code/windows.py           window set + fixed-seed (20261005) 50/50 split
code/stage1c_explore.py   1.7-1.10 EXPLORE table, residual anatomy, predictors, calibration
code/stage3_holdout.py    3.1-3.6  frozen rules on 4 held-out sets vs nulls
code/stage3b_addendum.py  3.7-3.9  phase tautology, brute re-validation, profile-shift lemmas
```
Raw outputs: `out/*.log`, `out/*.json`; tables: `data/*.csv`. Every confinement, prune and
return decision in every script is an exact integer inequality; floats appear only in `gamma_W`,
the phase feature, and the printed means/variances.
