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

**Narrowed conjecture.** The weighted profile has a unique maximum for r ≥ 20 and θ ≥ 0.1384. The
threshold was confirmed with zero exceptions on held-out slopes and windows for θ ≥ 0.175, and
holds without exception in all data down to θ = 0.1384. Slopes near 0 are not claimed.

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

Both numbers appear, and they mean different things.

`0.175079` (`new-random#00`) is the smallest slope in the **pre-registered confirmatory re-test**,
so that is how far down the threshold is confirmed *out of sample*. `0.138380349690391`
(`random#19`, with `pi-3` next at `0.141592654`) is the smallest slope **tested at all**; both
satisfy `r >= 20` with zero ties, but both were in the dataset used to choose the threshold — they
are where the earlier `K >= 17` prediction broke — so at the low edge the evidence is consistency,
not out-of-sample confirmation. The bullet above states each separately. `shape/REPORT.md` §7
carries the same distinction in its Scope paragraph.

## Supporting figures, for anyone who checks

* slope-independence: at matched `r`, one cubic in θ fits 15 slopes with residual sd 0.0011–0.0092;
  golden vs rationals 55/89 and 144/233 agree in `Var/r` to five decimal places.
* threshold: 0 ties in 36,863 rows with `r >= 20`; largest observed tie at `r = 17`.
* ties not from balance: the margin falls from 0.47 to 0.0007 as `r` grows while ties stop
  entirely; among the 10% most balanced rows the tie rate is 0 for every band `r >= 20`.
* 2-adic: `P(v2 >= k) ~ 1.1 · 2^-k`, max 17 over 63,972 rows against log₂(n) = 16.0 expected.
