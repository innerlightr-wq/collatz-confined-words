# PREDICTION.md — frozen before any HOLDOUT computation (shape study)

Frozen at the timestamp of the git commit adding this file. Nothing below was fitted on, or
informed by, the HOLDOUT half of the slope split (`shape/code/slopes.py:split`, seed 20261005),
nor by any horizon h > 150. Stage 0 used `beta` for the reproduction check specified by the task;
no model was fitted to it.

## Definitions (all staircases and counts exact)

Slope `theta` in (0,1), `G(n) = floor(n*theta)` exact. Window `(M,E)`:
`r(h) = G(M+h)-G(M)-E+1 = ell_0`, `ell_j = min(ell_0, G(M+h)-G(M+j))`,
`Pi_h(j) = Q_{L(h)}(j) + [j=0]`, Haar weights `w(j) = Pi_h(j) 2^{-j}`,
`K(h) = G(M+h) - G(M+5)` (the existing study code's convention).

Convention-free statistics (K depends on an arbitrary reference horizon; r does not):
`U_var = Var/r`, `U_loc = (mean - r)/sqrt(r)`, `U_tv = S4m*sqrt(r)` where `S4m` is the total
variation distance to the mean-matched negative binomial NB(round(mean), 1/2).

## P1 — The untied mode is UNIVERSAL (claimed exact)

For every HOLDOUT slope, every window in {(2,1),(3,1),(5,2),(8,1),(13,3)} and every horizon,
the Haar-weighted mode is **unique** whenever `K(h) >= 17`. Predicted: **zero ties** among all
HOLDOUT rows with K >= 17, including at extended horizons h = 160..400. Ties occur only at small
K, with a frequency falling from ~4% (K>=1) to 0% (K>=17), as on EXPLORE.
A single tie at K >= 17 falsifies P1.

## P2 — Shape depends on theta smoothly, not on Diophantine type

The per-slope means of `U_var` and `U_loc` over rows with K >= 20 are predicted by the cubic in
theta fitted on the 18 EXPLORE slopes (EXPLORE residual sd: U_var 0.0034, U_loc 0.0059):

    U_var(theta) = -3.22576 t^3 + 2.31929 t^2 - 0.83491 t + 2.05946
    U_loc(theta) = -2.85963 t^3 + 2.09005 t^2 - 0.63943 t - 0.17791

Frozen numeric predictions for each HOLDOUT slope are in `out/frozen_curve.json`; the key ones:

| slope | family | theta | U_var pred | U_loc pred |
|---|---|---|---|---|
| pi-3 | named | 0.141593 | +1.9786 | -0.2347 |
| liouville a_3=50 | liouville | 0.504891 | +1.8140 | -0.3360 |
| **beta = log2(3)-1** | named | 0.584963 | **+1.7190** | **-0.4092** |
| liouville a_7=500 | liouville | 0.615396 | +1.6722 | -0.4463 |
| rational 55/89 | rational | 0.617978 | +1.6679 | -0.4498 |
| rational 144/233 | rational | 0.618026 | +1.6679 | -0.4498 |
| **golden [0;1,1,1,...]** | bounded | 0.618034 | **+1.6678** | **-0.4498** |
| sqrt3-1 | bounded | 0.732051 | +1.4257 | -0.6478 |

Predicted out-of-sample accuracy on HOLDOUT, for slopes with theta inside the fitted range
[0.1921, 0.7183]: **residual sd <= 0.02 (U_var) and <= 0.03 (U_loc), max |residual| <= 0.08.**
Two HOLDOUT slopes extrapolate (`random#19` at theta=0.1384, `random#07` at theta=0.8150); they
are excluded from the tolerance and reported separately.

## P3 — No Diophantine feature adds anything

On HOLDOUT rows with K >= 20, after the baseline `B1t` (constant, 1/sqrt(r), 1/r, theta, theta^2,
theta/sqrt(r), and the local-phase block frac((M+h)theta), frac(M theta), r mod 2, r mod 3):

* adding the CF block (log max a_k, log sum a_k, CF depth, distance of M+h and of h to the
  nearest convergent denominator) reduces out-of-sample RMSE by **less than 5%**;
* no CF feature correlates with the per-slope residuals at p < 0.05 (4000-permutation test);
* the resonance test — mean |residual| for rows with `min_k |h - q_k| <= 2` versus the rest —
  gives p > 0.05 against a random-denominator null (on EXPLORE: p = 0.43 and 0.78, with the
  near-convergent rows having *smaller* residuals).

## P4 — The sharpest control (rational vs badly approximable at equal theta)

`golden` (theta = 0.6180340, continued fraction all 1s, the hardest irrational to approximate)
and `rational 144/233` (theta = 0.6180258, a ratio of integers, periodic staircase) differ in
theta by 8.2e-6. Predicted: their measured `U_var` means differ by **less than 0.01**, and
likewise for `rational 55/89` (theta = 0.6179775) and `liouville a_7=500` (theta = 0.6153959).
If Diophantine type entered the global shape at all, these four — spanning rational, Liouville-like
and badly-approximable — should separate. Predicted: they do not.

## P5 — Extended horizons

At h = 160..400 (beyond anything used in exploration) the same statements hold: no tie at
K >= 17, and the per-slope means stay within the P2 tolerance of the frozen curve.

## Verdict rule (frozen)

* **GO** only if some CF-feature model beats B1t on HOLDOUT by a clear margin for at least one
  statistic, the effect survives the calibration nulls, and it is not reducible to the local
  rotation phase or to the value of theta.
* **UNIVERSAL** if P1-P5 hold: shape statistics are a smooth function of (r, theta) plus a local
  phase term, with rational controls on the same curve.
* **STOP** for the Diophantine hypothesis if P3 holds with no CF signal.
