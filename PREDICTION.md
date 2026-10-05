# PREDICTION.md — frozen before any HOLDOUT computation

Frozen at the timestamp of the git commit that adds this file.
Nothing below was fitted on, or informed by, the HOLDOUT half of the window split
(`data/windows.csv`, column `split`), nor by any horizon h > 150, nor by any
window with M > 40.

## Notation (all exact integer arithmetic)

* `F(n) = floor(n*log2 3)` = `(3**n).bit_length()-1`;  `G(n) = F(n) - n = floor(n*beta)`, `beta = log2(3/2)`.
* Window `W = (M, c, prefix, V_min)`, `P = sum(prefix)`, **exit margin** `E = P + V_min - c - F(M) >= 1`.
* `eps_h = C_W(h+1) - C_W(h) = G(M+h+1) - G(M+h)`,  `C_W(h) = c - h + F(M+h)`.
* **Seed count** `r(W,h) = G(M+h) - G(M) - E + 1`  ( = number of admissible terminal digits
  `V in [V_min, V*(h)-1]` at horizon h; `r <= 0` means the window admits none).
* `T1(N)(0) = N(0)+1`, `T1(N)(j) = sum_{i<=j} N(i)` for `j >= 1` (the `+1` does **not** propagate).

## P1 — Generalized operator law (claimed EXACT, all windows, all j)

For every window `W` and every horizon `h >= 1`:

1. **(T0)** If `eps_h = 0` then `N_{h+1}(j) = N_h(j)` for **every** `0 <= j <= h`
   (no interior margin needed — stronger than the `j <= h-2` of the source note).
2. **(T1)** If `eps_h = 1` **and** `r(W,h) >= 0` then `N_{h+1}(j) = T1(N_h)(j)` for **every**
   `0 <= j <= h` (again no interior margin needed).
3. **(degenerate branch)** If `eps_h = 1` and `r(W,h) <= -1` then `N_h = N_{h+1} = {0: 1}`
   and the residual is exactly the constant vector:
   `N_{h+1}(j) - T1(N_h)(j) = -1` for every `0 <= j <= j_max`.

Expected accuracy of P1 as a decision rule ("does the literal T1 identity hold at this
`eps_h = 1` transition?"): **100.00%**, i.e. zero errors, on
(a) HOLDOUT windows, `h <= 150`; (b) EXPLORE windows at extended horizons `151 <= h <= 400`;
(c) a wide holdout of unseen phases `M = 41..90`, `E = 1..12`, `h <= 200`.
Any single counterexample falsifies P1.

Auxiliary exact claims frozen with P1:
* `N_h(h) = 0` for every window and every `h >= 1` (the horizon-exhaustion clause never fires).
* every reachable state has `d* >= 2` (no clipping).
* two windows with equal `(M, E)` have identical `N_h(.)` for every `h`.

## P2 — Negative prediction about Diophantine structure

No Diophantine feature improves on P1, and none beats the trivial base-rate predictor
by more than permutation noise on held-out data. Concretely, on each held-out set:

* the predictor `|h - q_k| <= tol` (`q_k` = convergent denominators of `beta`:
  1, 2, 5, 12, 41, 53, 306, 665), for `tol = 0,1,2,3`, will have accuracy **below** the
  base-rate predictor "always success";
* the best EXPLORE-fitted interval rule on `frac(h*beta + gamma_W)`, `gamma_W = frac(M*alpha)`,
  will equal the base-rate predictor (the fitted interval is the whole circle `[0,1)`)
  and so will carry **zero** information;
* conditional on `r(W,h)`, every phase feature is independent of success: predicted
  conditional mutual information exactly 0 (no failure with `r >= 0`, no success with `r <= -1`).

## P3 — Artifact reconstruction (frozen as a falsifiable claim)

The reported "M = 3, c = 1: literal identity matched only 16 of 56 tested one-digit
transitions" is reproduced exactly, and with no free parameters, by selecting the
`eps = 1` transitions with the **M = 2** clock `C(h) = floor(2*alpha + h*beta)` while
computing the tree for an `M = 3` window (whose own clock is `floor(c + 3*alpha + h*beta)`,
i.e. the same Sturmian word shifted by one index):

* `#{1 <= h <= 95 : e_{2+h} = 1} = 56`;
* of those, `#{h : e_{3+h} = 1 too} = 16`;
* the gaps between those 16 horizons are exactly `{5, 7}` (and `12 = 5 + 7` whenever one is
  skipped), which are the return times of the Sturmian factor `11` and are therefore forced
  to be convergent-denominator-sized by the three-distance theorem — carrying no information
  about the operator law.

The 40 "failures" are then exactly the transitions with `e_{3+h} = 0`, where the true law is
`T0 = identity`; `T1` applied there differs from the identity by a full cumulative sum, which
explains the reported "often large" discrepancies.

## Decision rule for the verdict (frozen)

* **GO** only if some held-out transition violates P1 **and** some Diophantine feature then
  beats every null on held-out data.
* **STOP / ARTIFACT** if P1 holds with zero exceptions held out and P2's negative predictions hold.
