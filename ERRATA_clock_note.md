# Errata and status changes — *A Sturmian Complexity Clock in Accelerated Collatz Words*
### (the 26-page revision, `SturmianClock (3) (1).pdf`, built 2026-08-25)

Keyed to that note's own numbering. Every entry is supported by
`../REPORT.md` (results), `../code/` (exact-arithmetic scripts) and
`appendix_C_reduction.tex` (proofs, as Appendix C of the September capacity note).

Notation here: `G(n) = floor(n*beta)`, `E = sigma + V_min - c - floor(M*alpha)` (exit margin,
`sigma` = interior prefix sum), `r(W,h) = G(M+h) - G(M) - E + 1` (number of admissible terminal
digits at horizon `h`).

| # | item in the note | status | change |
|---|---|---|---|
| 1 | **Observation 11** (`T1` = cumulative summation on the one-digit) | Computational -> **PROVED** | Now Theorem C.11(2). Holds for **every** window, not only the `M=2` checkpoint, and on the **whole** profile `0 <= j <= h` rather than "interior `j`" — subject only to `r(W,h) >= 0`. |
| 2 | **Conjecture 12** (the missing combinatorial lemma) | Open -> **RESOLVED** | The exact argument asked for is Lemma C.7 (uniform unit shift of a descending Beatty staircase) applied to the reduction Theorem C.6. It is bijective/algebraic, uses only the closed form of Prop. 3 and the definitions, and appeals to no computation. |
| 3 | **Theorem 10** (`T0` = identity), stated for `j <= h-2` | Proved -> **range widened** | Holds for every `j <= h`. The caution about the exhaustion clause is unnecessary: the clause never fires (entry 5), so nothing has to be reconciled at `j = h-1`. |
| 4 | **Remark 4** (clipping at `delta = 1` "downgraded to a computational statement") | Computational -> **PROVED** | Theorem C.6(i): every reachable state has `d* >= 2`, seeds included, so the unclipped formula is always the true one. The global reachability invariant the remark says was missing is exactly the change of variables `x = d* - 2 >= 0`. |
| 5 | **Theorem 10, "where the two horizons first genuinely diverge"** | — | `N_h(h) = 0` identically: at level `h-1` a continuing digit would need `S + delta > c + floor((M+h)*alpha)` and `S + delta < c + floor((M+h)*alpha) + 1`. The horizon-exhaustion clause is never reached, so `N_h(h-1)` does **not** mix ordinary and exhaustion pruning. |
| 6 | **Corollary 13** (`N_h = T1^{K(h)}(N_{h_0})`) | Conditional -> **UNCONDITIONAL on its base range** | Exact for every window and every `h >= h_open(M,E) = min{h : r(W,h) >= 0}`, but only on `0 <= j <= h_0` (Corollary C.12). The iterated statement does **not** extend to `j > h_0`: each step appends a top level with `N_{h+1}(h+1) = 0` while `T1` of the zero-extended profile puts `sum_i N_h(i) > 0` there, and since `T1` reads only lower indices the discrepancy never reaches below `h_0+1` and is never removed above it (Remark C.12a). For W0 the first mismatch is at `j = 12` when `h_0 = 10` and at `j = 21` when `h_0 = 20`. The single-step laws (entry 1) are unaffected: those are exact on all `j <= h`. |
| 7 | **Appendix A.4**, the `T1` code | — | **Bug.** `out = [N[0] + boundary_bump]` followed by `out.append(out[-1] + N[j])` propagates the `+1` into every later partial sum. Under that convention W0 scores **1/89**, not 89/89. Observation 11 as *stated* (no propagating bump) is the correct convention and the one that reproduces the reported verifications. Fix: `out.append(sum(N[:j+1]))`, or subtract the bump before accumulating. |
| 8 | **Appendix C.4** (*Non-universality across windows*): "for a deeper checkpoint with `M=3, c=1` the literal identity matched only 16 of 56 tested one-digit transitions and failed in the remaining 40", "discrepancies were often large", "exact recoveries ... at a highly structured set of horizons, with gaps involving values such as 5, 7, and 12 ... close to convergent or semiconvergent scales", "the cumulative operator may depend on a Diophantine phase alignment" | **WITHDRAWN (artifact)** | The `M=3, c=1` windows with `E = 1` and `E = 2` satisfy the literal identity on **88/88** one-digit transitions to `h = 150` (and 17/17 under the fully literal recursive tree). The reported figures are reproduced exactly, with no free parameters, by selecting the `eps = 1` horizons with the hard-coded `M = 2` clock `C(h) = floor(2*alpha + h*beta)` (Definition 5, Appendix A.1) while enumerating an `M = 3` tree, whose own clock is `floor(c + 3*alpha + h*beta)` — the same Sturmian word shifted one index: `#{h <= 95 : e_{2+h} = 1} = 56`, of which `#{e_{3+h} = 1} = 16`, gaps exactly `{5, 7}` (so `12 = 5+7`). The 40 "failures" are precisely the `e_{3+h} = 0` transitions, where the true law is `T0` = identity; applying `T1` there differs from the identity by a full cumulative sum, which is why the discrepancies were large rather than boundary-sized. |
| 9 | the Diophantine reading of "gaps 5, 7, 12" | **WITHDRAWN** | The 16 surviving horizons are the occurrences of the Sturmian factor `11`; the return times of any Sturmian factor take at most three values and are automatically convergent-sized (three-distance theorem). Observed gap multiset over 4,000 horizons: `{5: 379, 7: 300}`. Calibration: at tolerance 1, 70% of integers below 20 are near a semiconvergent denominator of `beta` (1, 2, 3, 5, 7, 12, 17, 29, 41, 53, ...); and 16 successes spread over `1..96` have mean gap 6.0 by construction. Held out, no phase or convergent feature beat the base-rate predictor (accuracies 0.015-0.197 vs 0.948-1.000). |
| 10 | **Appendix C.8**, conclusion 5 / §C.4's closing sentence ("the `M = 2` checkpoint ... exhibits an unusually rigid exact operator law whose scope and arithmetic alignment deserve separate investigation") | **WITHDRAWN** | The law is not special to `M = 2`; it is universal over windows and, by Lemma C.7 + Remark C.8, uses no property of `beta` whatsoever. There is no arithmetic alignment to investigate. |
| 11 | **Section 2.3**, the prose branching rule | **AMBIGUOUS — clarify** | "Otherwise, if finalizing the current digit at `cur_min` already satisfies `R_{M+j+1} <= c`, the branch is a return" reads as though a return terminates the state. Appendix A.3, and the proof of Theorem 10 ("the classification of each candidate digit value `delta = 1,...,d*-1` as a return ... or a continuation"), classify **every** `delta` independently and continue past returns. The latter is the object analyzed, and the only one under which the reported statistics reproduce; the prose should say so. |
| 12 | **Observations 18-19** (Haar mode `K(h)-1`; `E[J]-K -> 0`; `Var/K -> 2`) | Computational — **unchanged, and reproduced** | Reproduced to every printed digit: mode `= K(h)-1` untied at `h = 40,60,80,100` with `K(h) = C(h)-C(5)`; `E[J]-K = +0.068, +0.012, -0.002, -0.006`; `Var/K = 1.659, 1.747, 1.802, 1.843`. These remain computational. By Theorem C.6 the shape question is now a self-contained asymptotic question about chain counts under the explicit staircase `l_j = min(l_0, G(M+h) - G(M+j))`, with no tree in it (Remark C.16). |
| 13 | **Section 13 claim audit** | — | Rows to update: "T1 = cumulative sum on `eps_h = 1` (Obs. 11)" Computational -> Proved; "Bulk mode/mean/variance rates (Conj. 20)" keeps its conjectural status but its *dependency* changes from "Obs. 11 (unproved half)" to "Thm. 8 + Thm. C.11 (both proved) + the unquantified boundary correction of Lem. 14", i.e. the only unproved input is now the boundary correction. |
| 14 | new material worth adding | — | **Window space collapses.** The profile depends on the window only through `(M, E)`: 135 of the obvious windows (`M in 2..6`, `c in 1..3`, several prefixes and `V_min`) are 15 distinct objects (Corollary C.13). Any future claim of window-dependence should be tested across `(M, E)` pairs, and any claim of `eps`-driven structure should use the window's own clock `C_W(h) = c - h + floor((M+h)*alpha)`, never the `M = 2` clock. |

## Degenerate branch (the only genuine exception, and it is not window-dependence)

`T1` fails exactly when `r(W,h) <= -1`, i.e. `V_min >= V*(h)`: the window admits no terminal digit
at horizon `h` and none at `h+1`, both profiles are the single terminal prune leaf `{0: 1}`, and the
residual `N_{h+1} - T1(N_h)` is identically `-1`. Equivalently, failure occurs only for
`h < h_open(M,E) ~ (E-1)/beta + O(1)`. This is a "the window has not opened yet" degeneracy — a
magnitude threshold in `h`, not a phase condition.

## Scale of verification

Exact integer arithmetic throughout, three independent implementations (literal recursive tree,
DP over `{S_fin: multiplicity}`, chain count), agreeing on all 240 tested (window, horizon) pairs.
`T0`/`T1` confirmed on 41,076 one-digit and 29,124 zero-digit transitions over 468 `(M,E)` classes
(`M <= 40`, `E <= 12`, `h <= 150`), plus 51,580 one-digit and 36,620 zero-digit transitions held
out behind a pre-registered prediction (unseen windows; `151 <= h <= 600`; unseen phases
`41 <= M <= 120` and `101 <= M <= 110` with `h <= 400`). Zero exceptions; every failure of the
literal rule satisfied `r(W,h) <= -1` with residual identically `-1`.

## Non-claims (unchanged)

Nothing here bears on the Collatz conjecture, on the EOC bound `O_c(m) = O(log_2 m)`, or on any
relation between these quantities and the bit-length of a genuine integer realizer.
