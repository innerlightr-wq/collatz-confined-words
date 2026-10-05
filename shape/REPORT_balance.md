# VERDICT: **REJECTED** — balance is not what selects a tie; ties are small-integer coincidences

The mechanism predicts that ties appear where the staircase is locally uniform enough for the
profile to approach the perfectly balanced ideal kernel. The data run the other way. As `r` grows
the profile becomes **steadily more balanced** — the median balance margin falls from `0.47` at
`r < 5` to `0.0007` at `r >= 100` — and yet exact ties **stop happening entirely** above `r = 17`.
What actually decides a tie is the exact integer gap `Delta = Pi(j*+1) - 2 Pi(j*)`, whose median
magnitude explodes from `3` to `10^122` over the same range: a tie needs `Delta = 0`, and that is
an arithmetic coincidence available only while the profile entries are small integers. On held-out
data the largest `Pi(j*)` carrying a tie was **1,001**, against untied profiles reaching 101
digits; and conditioning on the **10% most balanced rows** of each `r` band — margins down to
`4e-4`, i.e. maximal balance — the tie rate is **0 for every band with `r >= 10`**.
Applying the frozen verdict rule (`REJECTED if P2 fails and P5 holds`): P2's required positive
relation between margin and top-irregularity is absent (held-out median within-cell Spearman
`-0.118`, 25% of cells positive), and P5 holds in every clause.

Two parts of the mechanism's picture do survive, and are reported as separate findings: all ties
are adjacent runs (P1), and the distance to the ideal kernel really does grow with `theta` at
matched `r` (P3's direction, `Spearman(theta, m) = +0.20`, `p = 0.0005`) — but neither is what
produces a tie.

No claim here concerns the Collatz conjecture or any `O(log m)` bound.

---

## 1. What was computed

EXPLORE: the 47 slopes already in the study x 8 windows `{(2,1),(3,1),(5,2),(8,1),(13,3),(4,2),(7,1),(11,2)}`
x `h <= 200` = **74,271 rows**. HOLDOUT (frozen before use): **8 new random-CF slopes**, seed
20261007, `theta` spread over `[0.15, 0.9]` one per band, windows `(2,1),(4,2),(7,1)`, `h <= 400`
= **9,567 rows**. Profiles advanced with the proved operators and spot-checked against direct
chain-count recomputation (0 failures). All ratios exact `Fraction`s; argmax, ties and drop counts
exact integer operations; logs are a descriptive transform of exact rationals.

Definitions fixed before computing: `rho(j) = Pi(j+1)/(2 Pi(j))`; peak `j*` = argmax (smallest
index if tied); `m = min(|log rho(j*-1)|, |log rho(j*)|)`; `D` = staircase drops in
`j in [j*-5, j*+5]`; `V` = variance of the gaps between those drops; `k_fit` = argmin TV to
`w*_k(j) = C(j+k-1,k-1)2^{-j}` over `k` within 20 of `round(mean)`; `S4` = that TV;
`Delta = Pi(j*+1) - 2 Pi(j*)`.

## 2. Balance margin and the integer gap, by `r` (EXPLORE)

| `r` band | n | tie rate | median `m` | median `\|Delta\|` | median digits of `Pi(j*)` |
|---|---|---|---|---|---|
| 1–5 | 3,953 | 0.2100 | 0.470004 | 3 | 1 |
| 5–10 | 4,962 | 0.3116 | 0.010471 | 5 | 3 |
| 10–15 | 4,966 | 0.0818 | 0.007868 | 2,682 | 6 |
| 15–20 | 4,952 | 0.0153 | 0.006057 | 1,875,119 | 9 |
| 20–30 | 9,731 | **0.0000** | 0.003978 | 2.6e10 | 13 |
| 30–50 | 14,764 | **0.0000** | 0.003420 | 4.5e18 | 21 |
| 50–100 | 22,712 | **0.0000** | 0.001895 | 3.4e46 | 40 |
| 100–400 | 8,231 | **0.0000** | 0.001130 | 3.0e82 | 67 |

The two columns move in opposite directions. This is the whole result; everything else is detail.
Figure: `figs/balance_summary.png` panels (a) and (b).

**Conditioning on balance (EXPLORE).** Tie rate among the 10% most balanced rows of each band:

| `r` band | 1–5 | 5–10 | 10–15 | 15–20 | 20–30 | 30–50 | 50–100 | 100–400 |
|---|---|---|---|---|---|---|---|---|
| most-balanced tie rate | 1.00 | 1.00 | 0.82 | 0.15 | **0.00** | **0.00** | **0.00** | **0.00** |
| margin cutoff `m <=` | 0 | 0 | 1e-6 | 1.8e-5 | 3.6e-5 | 3.2e-4 | 3.9e-4 | 2.6e-4 |

Perfect balance, zero ties, as soon as the integers are large.

## 3. Margin `m` and `S4` by `theta` and `r` (EXPLORE)

Median `m`:

| `theta` \ `r` | 1–10 | 10–20 | 20–40 | 40–80 | 80–400 |
|---|---|---|---|---|---|
| 0.13–0.25 | 0.0059 | 0.0009 | 0.0009 | 0.0036 | – |
| 0.25–0.40 | 0.0296 | 0.0051 | 0.0025 | 0.0016 | – |
| 0.40–0.55 | 0.0335 | 0.0089 | 0.0050 | 0.0025 | 0.0017 |
| 0.55–0.70 | 0.0338 | 0.0123 | 0.0052 | 0.0024 | 0.0013 |
| 0.70–0.90 | 0.0513 | 0.0152 | 0.0060 | 0.0028 | 0.0012 |

Median `S4`:

| `theta` \ `r` | 1–10 | 10–20 | 20–40 | 40–80 | 80–400 |
|---|---|---|---|---|---|
| 0.13–0.25 | 0.0479 | 0.0226 | 0.0106 | 0.0137 | – |
| 0.25–0.40 | 0.0496 | 0.0196 | 0.0101 | 0.0073 | – |
| 0.40–0.55 | 0.0519 | 0.0279 | 0.0184 | 0.0120 | 0.0095 |
| 0.55–0.70 | 0.0672 | 0.0467 | 0.0339 | 0.0210 | 0.0141 |
| 0.70–0.90 | 0.0860 | 0.0786 | 0.0719 | 0.0646 | 0.0494 |

Both rise with `theta` at matched `r` — the mechanism's one correct qualitative prediction
(figure panel (c)) — and both fall with `r`, which is the direction that sinks it.

Top-irregularity `D` is `~11*theta` by construction, so any pooled `m`-vs-`D` relation merely
restates `theta`; the test must hold `theta` and `r` fixed (figure panel (d)).

## 4. Frozen predictions vs HOLDOUT

| | prediction (frozen) | observed on HOLDOUT | verdict |
|---|---|---|---|
| **P1** | 0 non-adjacent ties | 71 ties, **0** non-adjacent; runs of length 2 (46) and 3 (25), all gaps = 1 | **HOLDS** |
| **P2** | mechanism needs Spearman(`m`,`D`) `> 0` within (slope, `r`-band); I predicted it would *fail*, median in `[-0.10,+0.10]`, fraction positive in `[0.35,0.65]` | median **−0.1176**, fraction positive **0.250** (32 cells) | mechanism's claim **REJECTED**; my own numeric bracket slightly missed — the relation is weakly *negative*, not null |
| **P3** | at `r` 10–20, `m` ratio `>= 3`, `S4` ratio `>= 1.5`, Spearman(`theta`,`m`) `> 0` at `p < 0.05` | `m` ratio **1.39** (predicted ≥3), `S4` ratio **3.34** ✓, Spearman **+0.2006**, `p = 0.0005` ✓ | **FAILS as frozen** (direction right, magnitude for `m` overstated by the EXPLORE fit) |
| **P4** | as stated, vacuous: `\|j* − (k_fit−2)\| <= 1` for `>= 99%` of both classes | tied **1.0000**, untied **0.9937** | **CONFIRMED vacuous** |
| **P4'** | `P(j* = k_fit−2 \| tied) − P(... \| untied) >= +0.15` | 0.437 − 0.357 = **+0.079** | **FAILS** (sign right, half the predicted size) |
| **P5** | no tie with `Pi(j*) >= 10^10`; tie rate 0 for all `r >= 20`; most-balanced tie rate 0 for `r >= 20`; median `\|Delta\|` monotone in `r` | largest `Pi(j*)` with a tie **1,001**; **0** ties in 8,534 rows with `r >= 20`; most-balanced tie rate **0** from `r >= 10`; `\|Delta\|` monotone, 3 → 1.6e122 | **HOLDS, every clause** |

Frozen verdict rule: *REJECTED if P2 fails and P5 holds.* Both conditions are met.

## 5. Claim audit

| # | claim | label |
|---|---|---|
| 1 | A tie requires `Pi(j*+1) = 2 Pi(j*)` exactly | **PROVED** (definition of `w(j) = Pi(j)2^{-j}`) |
| 2 | All tied sets are adjacent runs (length 2 or 3); no non-adjacent tie | **VERIFIED** (2,858 EXPLORE + 71 HOLDOUT ties, 0 non-adjacent) |
| 3 | The balance margin `m` decreases with `r` while `\|Delta\|` increases | **VERIFIED** over 83,838 rows |
| 4 | Ties vanish above `r = 17` despite margins continuing to shrink | **VERIFIED**; consistent with the separately frozen `r >= 20` threshold |
| 5 | Conditioning on maximal balance does not produce ties once `r >= 20` | **VERIFIED** (0 of the most-balanced decile in every band `r >= 20`, EXPLORE and HOLDOUT) |
| 6 | `m` and `S4` increase with `theta` at matched `r` | **VERIFIED** (EXPLORE tables; HOLDOUT Spearman `+0.20`, `p = 0.0005`) — but this is **not** the tie mechanism |
| 7 | `m` increases with top-irregularity `D` at fixed `theta`, `r` | **REJECTED** (held-out median Spearman `−0.118`, 25% of cells positive) |
| 8 | The peak sits within one index of the ideal balance point `k_fit − 2` | **VERIFIED but vacuous** — true of 99.4–100% of rows, tied and untied |
| 9 | Ties are an arithmetic small-integer effect, not a balance effect | **OBSERVED** — the reading the data support; not proved |
| 10 | Anything about Collatz or `O(log m)` | **NOT CLAIMED** |

## 6. What this means for the uniqueness conjecture

The conjecture of `REPORT.md` §7 is unaffected in substance but its *proof strategy* changes. Since
ties are exact integer coincidences rather than a balance phenomenon, a proof of uniqueness for
`r >= 20` should not be sought by showing the profile is far from balanced — it is not; it is
extremely close to balanced, and increasingly so. The statement to prove is arithmetic:

> **CONJECTURE (restated after this study).** For the staircase chain counts `Q_L` of
> `REPORT.md` §7 and `r >= 20`, the exact integer `Delta(j) = Pi(j+1) - 2 Pi(j)` is non-zero at
> the peak — equivalently, `Pi(j*+1) != 2 Pi(j*)`: a `2`-adic rather than an analytic statement.

This also sharpens the log-concavity route already recorded there: log-concavity would give
unimodality, and the separate argument needed to exclude a two-way tie is exactly a statement that
`2 Pi(j*)` is never hit on the nose. The observation that every tied set is an adjacent run of
length 2 or 3 says the coincidence, when it happens, happens on consecutive indices.

**Update (`REPORT_2adic.md`): the 2-adic route suggested here is closed.** This section proposed a
`v_2` argument as the plausible way to exclude the tie. Testing it on all 63,972 rows with
`r >= 20` found `v2(Delta)` distributed as for a random integer — `P(v2 >= k) ~ 1.1 * 2^-k`, max
**17** against `log2(n) = 16.0` expected, no growth or ceiling in `r` — so there is no bounded
`k0`, no modulus at which `Delta != 0` can be certified, and no residue class excluded. Tie
exclusion is not a fixed-modulus 2-adic statement. What the same data do give is a heuristic
account of the tie-free regime: `Delta` has median bit-length 35 to 467 across the `r` bands, so
`P(Delta = 0)` under a random-integer model is around `2^-35` summed over everything computed.

## 7. Files

```
shape/code/balance_core.py       exact rho, peak, margin, irregularity, ideal-kernel fit
shape/code/stage_b1_explore.py   B1.1-B1.6  descriptive (EXPLORE)
shape/code/stage_b1b_gap.py      B1.7-B1.9  the integer-gap diagnostic
shape/code/stage_b2_holdout.py   P1-P5 on 8 new slopes, h <= 400
shape/code/stage_b3_figs.py      figures
shape/figs/balance_summary.png   (a) margin vs r with ties  (b) m vs |Delta|
                                 (c) S4 vs theta at matched r  (d) within-slope Spearman(m,D)
shape/figs/balance_holdout.png   holdout: ties only at small r; most-balanced tie rate
shape/data/balance_explore.csv, balance_holdout.csv ; shape/out/balance_*.log
shape/PREDICTION_balance.md      pre-registration (committed before the holdout run)
```
