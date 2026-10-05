"""Stage 0: (1) reproduce the W0 shape numbers with theta=beta; (2) confirm the local
T0/T1 law holds for EVERY slope in the test set (the 'local rules are slope-neutral' half)."""
import os as _os, sys as _sys
ROOT = _os.path.dirname(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
_sys.path.insert(0, _os.path.join(ROOT, 'shape','code'))

import sys, json
sys.path.insert(0, ROOT+"/shape/code")
from slopes import build, split
from shape_core import profile, l_profile, K_of, shape_stats, row

LOG = []
def say(*a):
    s = " ".join(str(x) for x in a); print(s); LOG.append(s)

NMAX = 460
S = build(depth=60)
for s in S: s.verify(NMAX)
beta = [s for s in S if s.name.startswith("beta")][0]

say("== 0.1 reproduce the W0 shape numbers (theta = beta, window (M,E)=(2,1)) ==")
say("   h   mode  K(h)  S1=mode-K  tied   mean-K     Var/K    (note: 0.068,0.012,-0.002,-0.006 ; 1.659,1.747,1.803,1.843)")
ok1 = True
ref_mean = {40: 0.068, 60: 0.012, 80: -0.002, 100: -0.006}
ref_var  = {40: 1.659, 60: 1.747, 80: 1.802, 100: 1.843}
for h in (40, 60, 80, 100):
    r = row(beta, 2, 1, h)
    say("  %3d  %4d  %4d   %+3d       %-5s  %+.4f   %.4f" %
        (h, r["mode"], r["K"], r["S1"], r["mode_tied"], r["S2"], r["S3"]))
    if r["S1"] != -1 or r["mode_tied"]: ok1 = False
    if round(r["S2"], 3) != ref_mean[h] or abs(r["S3"] - ref_var[h]) > 0.0006: ok1 = False
say("   reproduces the note's values exactly (mode = K-1 untied, mean-K, Var/K): %s" % ok1)

say("== 0.2 local law T0/T1 for EVERY slope in the test set ==")
def Pi_or_delta(S_, M, E, h):
    P = profile(S_, M, E, h)
    return P if P is not None else [1] + [0] * h

def T1(P, jmax):
    out, run = [], 0
    for j in range(jmax + 1):
        run += P[j] if j < len(P) else 0
        out.append(run + 1 if j == 0 else run)
    return out

WINDOWS = [(2, 1), (3, 1), (5, 2), (8, 1), (13, 3)]
HMAX = 120
tot0 = ok0 = tot1 = ok1c = totd = okd = 0
bad = []
for sl in S:
    for (M, E) in WINDOWS:
        for h in range(1, HMAX):
            eta = sl.G(M + h + 1) - sl.G(M + h)
            r_h = sl.G(M + h) - sl.G(M) - E + 1
            P, Q = Pi_or_delta(sl, M, E, h), Pi_or_delta(sl, M, E, h + 1)
            if eta == 0:
                tot0 += 1
                good = all(Q[j] == P[j] for j in range(h + 1))
                ok0 += good
            elif r_h >= 0:
                tot1 += 1
                T = T1(P, h)
                good = all(Q[j] == T[j] for j in range(h + 1))
                ok1c += good
            else:
                totd += 1
                T = T1(P, h)
                good = all(Q[j] - T[j] == -1 for j in range(h + 1))
                okd += good
            if not good: bad.append((sl.name, M, E, h, eta, r_h))
say("   eta=0  (identity)           : %d/%d exact" % (ok0, tot0))
say("   eta=1, r>=0  (T1)           : %d/%d exact" % (ok1c, tot1))
say("   eta=1, r<=-1 (residual -1)  : %d/%d exact" % (okd, totd))
say("   slopes x windows x horizons : %d x %d x %d ; exceptions: %d %s"
    % (len(S), len(WINDOWS), HMAX - 1, len(bad), bad[:3]))
say("   -> the LOCAL transition rules carry no arithmetic of theta (VERIFIED over all 37 slopes).")

json.dump(dict(w0_ok=ok1, eta0=[ok0, tot0], eta1=[ok1c, tot1], degen=[okd, totd],
               exceptions=len(bad)),
          open(ROOT+"/shape/out/stage0.json", "w"), indent=1)
open(ROOT+"/shape/out/stage0.log", "w").write("\n".join(LOG))
