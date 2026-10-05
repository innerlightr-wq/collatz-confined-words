# Draft text for the Zenodo description (record 10.5281/zenodo.23167068)

Fact-checked against `shape/REPORT.md`, `shape/REPORT_balance.md`, `shape/REPORT_2adic.md` and
the committed pre-registrations. One number is changed from the draft: the slope threshold is
**0.1384**, not 0.175 — see the note below the text.

---

**Update — work after v2 (October 2026; to be included in the next revision)**

Follow-up computations on the open shape question of Remark C.16 established:

**Slope-independence (verified).** The global shape of the staircase profile does not depend on
the Diophantine type of the slope. Rational slopes, bounded-type irrationals, generic irrationals
and Liouville-type slopes give the same shape statistics once matched by slope value. The local
transition law (Theorem C.11) was re-verified exactly for every slope tested.

**Narrowed conjecture.** The weighted profile has a unique maximum for r ≥ 20 and slopes
θ ≥ 0.1384. This threshold was fixed in advance and confirmed on held-out slopes and windows with
zero exceptions. Slopes near 0 are not claimed.

**The peak location is specific to the Collatz case.** For θ = log₂3 − 1 and the reference window,
the peak sits at r − 4. This corrects the earlier "K − 1" phrasing, which depended on a choice of
convention.

**Two proof routes ruled out.** Ties are not caused by near-balance of the profile; they are
coincidences between small integers. The 2-adic valuation of the tie gap is distributed like a
random integer's, so no fixed-modulus argument can exclude ties. Unimodality remains approachable
through log-concavity, but proving the absence of ties is expected to be hard.

Code and reports: https://github.com/innerlightr-wq/collatz-confined-words

---

## Note on the threshold

The draft said **θ ≥ 0.175**; the text above says **θ ≥ 0.1384**, which is what the committed
conjecture states. `0.175` is the smallest slope in the *confirmatory re-test* only
(`new-random#00`); across all three runs the smallest slope tested is `random#19` at
**θ = 0.138380349690391**, with `pi-3` next at `0.141592654`, and both were tested at `r >= 20`
with zero ties. `shape/REPORT.md` §7 states the conjecture with
`theta_min = 0.138380349690391`, so the record and the Zenodo text now agree.

## Supporting figures, for anyone who checks

* slope-independence: at matched `r`, one cubic in θ fits 15 slopes with residual sd 0.0011–0.0092;
  golden vs rationals 55/89 and 144/233 agree in `Var/r` to five decimal places.
* threshold: 0 ties in 36,863 rows with `r >= 20`; largest observed tie at `r = 17`.
* ties not from balance: the margin falls from 0.47 to 0.0007 as `r` grows while ties stop
  entirely; among the 10% most balanced rows the tie rate is 0 for every band `r >= 20`.
* 2-adic: `P(v2 >= k) ~ 1.1 · 2^-k`, max 17 over 63,972 rows against log₂(n) = 16.0 expected.
