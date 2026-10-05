# Note for a future revision of Remark C.16 (September capacity note)

Not applied to `main.tex` or to the published v2 PDF; this is draft wording for whenever that
note is next revised.

> The shape question left open in Remark C.16 has been tested across slopes and is **not**
> arithmetic: over 46,911 exact-integer profiles spanning 47 slopes with
> `theta in [0.1384, 0.8791]` — badly approximable quadratics, Liouville-like slopes with partial
> quotients 50 and 500, rationals, `e-2`, `pi-3` and random continued fractions — every global
> shape statistic of the Haar-weighted profile is reproduced by a smooth function of the slope's
> *value* `theta` together with the local rotation phase, while continued-fraction features
> (largest and summed partial quotients, star discrepancy, proximity of the horizon to a
> convergent denominator) add 1.0–1.3% of out-of-sample RMSE and show no correlation with
> per-slope residuals (`p >= 0.17`); at matched `r` the golden slope and the rationals `55/89`
> and `144/233` agree in `Var/r` to five decimal places. What survives as slope-independent is
> exactly the uniqueness of the maximum: no tie was observed at `r >= 20` in 36,863 such rows
> (largest observed tie `r = 17`, at the two smallest slopes tested), so the conjecture worth
> recording is that for every `theta` in the tested range `[0.1384, 1)`, rational or irrational,
> and every window tested, `w(j) = Q_L(j) 2^{-j}` has a unique maximum once `r >= 20` — a
> statement about monotone Beatty staircases with no arithmetic input, for which log-concavity of
> `j -> Q_L(j)` is the natural lemma, though log-concavity gives unimodality only and a separate
> argument is needed to exclude a two-way tie at the summit. The *location* claimed in that
> remark should be narrowed at the same time: `mode = K(h) - 1` is not a general fact but is
> specific to `beta`, the checkpoint `W_0 = (M,E) = (2,1)` and the reference convention
> `h_0 = 5`; the convention-free form of the same observation is `mode = r - 4`, which held in
> all 136 tested rows with `K >= 20` up to `h = 400`, whereas across other slopes and windows
> `mode - r` ranges over `{-6, …, -3}`. Slopes below the tested range are not claimed: ties
> concentrate at small `theta`, so the threshold may grow as `theta -> 0`.

Status labels: the uniqueness statement is OBSERVED / CONJECTURE, not proved. The
slope-independence finding (that the dependence is on `theta`, not on Diophantine type) is
VERIFIED computationally on held-out data against pre-registered predictions. Supporting
material: `shape/REPORT.md` §7 (conjecture and scope), §7b (confirmatory re-test),
`shape/PREDICTION.md` (pre-registration), `shape/code/`, `shape/data/`.
