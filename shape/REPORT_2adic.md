# VERDICT: **ROUTE CLOSED** — v2(Delta) is unbounded and geometric; no fixed-modulus certificate exists

Over **63,972 rows with `r >= 20`** (47 study slopes + 10 re-test slopes x 8 windows, `h <= 200`;
8 holdout slopes x 3 windows, `h <= 400`; 0 ties among them), the 2-adic valuation of the peak gap
follows the law of a random integer almost exactly, reaching **max v2(Delta) = 17** against the
`log2(63972) = 16.0` expected from geometric draws. There is therefore no small `k0` with
`v2(Delta) <= k0`, and no modulus `2^(k0+1)` at which `Delta != 0` could be certified. The
tie-exclusion step cannot be made 2-adic at a fixed modulus.

## 1. Distribution of v2

`Delta = Pi(j*+1) - 2 Pi(j*)` (`< 0` off a tie), `Delta' = Pi(j*) - 2 Pi(j*-1)` (`> 0`), both exact.

| k | count | obs `P(v2 = k)` | geometric `2^-(k+1)` | obs `P(v2 >= k)` | `2^-k` | ratio |
|---|---|---|---|---|---|---|
| 0 | 29,391 | 0.45944 | 0.50000 | 1.00000 | 1.00000 | 1.00 |
| 1 | 16,706 | 0.26115 | 0.25000 | 0.54057 | 0.50000 | 1.08 |
| 2 | 8,766 | 0.13703 | 0.12500 | 0.27942 | 0.25000 | 1.12 |
| 3 | 4,454 | 0.06962 | 0.06250 | 0.14239 | 0.12500 | 1.14 |
| 4 | 2,409 | 0.03766 | 0.03125 | 0.07277 | 0.06250 | 1.16 |
| 5 | 1,138 | 0.01779 | 0.01563 | 0.03511 | 0.03125 | 1.12 |
| 6 | 540 | 0.00844 | 0.00781 | 0.01732 | 0.01563 | 1.11 |
| 7 | 325 | 0.00508 | 0.00391 | 0.00888 | 0.00781 | 1.14 |
| 8 | 118 | 0.00185 | 0.00195 | 0.00380 | 0.00391 | 0.97 |
| 9 | 61 | 0.00095 | 0.00098 | 0.00195 | 0.00195 | 1.00 |
| 10–12 | 60 | — | — | 0.00100 | 0.00098 | 1.02 |
| 14 | 2 | — | — | 0.000063 | 0.000061 | 1.03 |
| 17 | 2 | — | — | 0.000031 | 0.000008 | 4.1 |

Mean `v2(Delta) = 1.1043` (geometric mean is 1.0). `v2(Delta')` is the same picture: max **13**,
mean 1.0866. The only systematic departure is a mild constant inflation of the tail,
`P(v2 >= k) ~ 1.1 * 2^-k` for `1 <= k <= 7` — `Delta` is even slightly more often than a random
integer (54.1% vs 50%) — but the decay rate is `2^-k`, which is what matters.

## 2. Boundedness and growth in `r`

| `r` band | n | max `v2(Delta)` | max `v2(Delta')` | mean `v2` | median bits of `Delta` |
|---|---|---|---|---|---|
| 20–30 | 10,277 | 14 | 11 | 1.2214 | 35 |
| 30–50 | 15,853 | 11 | 11 | 1.0976 | 63 |
| 50–80 | 17,177 | 11 | 13 | 1.0810 | 111 |
| 80–120 | 14,279 | **17** | 12 | 1.0950 | 176 |
| 120–200 | 5,092 | 12 | 11 | 1.0147 | 257 |
| 200–400 | 1,294 | 10 | 9 | 1.0216 | 467 |

`max v2` does not grow with `r` — it fluctuates between 10 and 17 because it is the maximum of
roughly geometric draws, and the bands differ in size. Equally, it shows no sign of a ceiling:
the largest value appears in the largest band, and the mean sits at the geometric value in every
band. **There is no `k0`.**

## 3. Residues, and the mod-`2^m` dynamics

`Delta mod 2^m` is close to uniform, with the mild even-excess implied by the `k = 0` deficit:

```
Delta  mod 2 : {0: 34581, 1: 29391}      Delta  mod 4 : {0:17875, 1:14560, 2:16706, 3:14831}
Delta' mod 2 : {0: 34574, 1: 29398}      Delta' mod 4 : {0:17678, 1:14926, 2:16896, 3:14472}
```
No residue class is excluded, including `0` at every modulus tested up to `2^16` — which is the
direct statement that a fixed-modulus certificate does not exist.

The reduction itself is exact, as it must be: `T0` (identity) and `T1` (cumulative sum with `+1`
at `j = 0`) are integer affine maps, so reduction mod `2^m` commutes with them. Verified directly:
the profile evolved entirely in `Z/2^m` from its first open horizon agrees with the exact profile
reduced mod `2^m` at **9,452 of 9,452** checked horizons (`m = 8` and `m = 16`, 40 slope-window
pairs), 0 mismatches. This confirms the implementation but yields no finite-state description: the
state is the whole profile vector mod `2^m`, whose length grows with `h`, and §2 shows no fixed `m`
suffices anyway.

## 4. What this does and does not say

The negative result is specific: it closes the *fixed-modulus* route to tie exclusion. It does not
contradict the empirical tie-free behaviour at `r >= 20` — it explains it heuristically. If
`Delta` behaves like a random integer of the observed size, then `P(Delta = 0)` is about
`2^-bits`, and the median bit-length of `Delta` runs 35, 63, 111, 176, 257, 467 across the `r`
bands. Summed over every row ever computed that is of order `2^-35`. So the absence of ties above
`r = 17` is exactly what a random-integer model predicts — which is reassuring for the
*conjecture* and discouraging for any short *proof*, since beating a random-integer heuristic for a
specific family of integers is the hard kind of problem.

## 5. Claim audit

| # | claim | label |
|---|---|---|
| 1 | `Delta = 0` iff the peak ties with its upper neighbour | **PROVED** (definition) |
| 2 | Reduction mod `2^m` commutes with `T0`/`T1` | **PROVED** (integer affine maps) + VERIFIED (9,452/9,452) |
| 3 | `v2(Delta)` follows `P(v2 >= k) ~ 1.1 * 2^-k`; max 17 over 63,972 rows | **VERIFIED** |
| 4 | `v2(Delta)` is unbounded; no `k0` exists | **OBSERVED** — consistent with geometric at every band; not proved |
| 5 | No fixed-modulus 2-adic certificate for tie exclusion | **OBSERVED**, following from 4 |
| 6 | Tie-freeness at `r >= 20` is consistent with a random-integer model | **OBSERVED** (heuristic, not a proof) |
| 7 | Anything about Collatz or `O(log m)` | **NOT CLAIMED** |

Data: `shape/data/twoadic.csv`; log `shape/out/twoadic.log`; code `shape/code/stage_b4_2adic.py`.
