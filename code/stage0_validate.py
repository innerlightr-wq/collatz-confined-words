"""Stage 0: validation + artifact check.  Writes out/stage0.json and prints a log."""
import os as _os, sys as _sys
ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
_sys.path.insert(0, _os.path.join(ROOT, 'code'))

import json, sys, time, math
from fractions import Fraction
sys.path.insert(0, ROOT+"/code")
from core import *

LOG = []
def say(*a):
    s = " ".join(str(x) for x in a)
    print(s); LOG.append(s)

W0 = Window(2, 1, (1,), 4, name="W0")

# a few windows for cross-checking (valid: all-1 prefixes are always c-confined)
def mkwin(M, c, prefix, extra=0):
    """V_min = minimal exiting digit + extra  (so E = 1 + extra)."""
    P = sum(prefix)
    Vmin = c + F(M) + 1 - P + extra
    if Vmin < 1:
        return None
    try:
        return Window(M, c, prefix, Vmin)
    except AssertionError:
        return None

results = {}

# ---------------------------------------------------------------- 0.1 cross-check
say("== 0.1 brute tree vs DP(literal) vs DP(closed form) vs reduced chain model ==")
xcheck = []
wins = [W0, mkwin(3,1,(1,1)), mkwin(3,1,(1,2)), mkwin(4,2,(1,1,1)), mkwin(2,1,(1,),1),
        mkwin(5,1,(1,1,1,1)), mkwin(4,1,(2,1,1)), mkwin(6,3,(1,1,1,1,1),2)]
wins = [w for w in wins if w]
HMAX_BRUTE = 15
for W in wins:
    for h in range(1, HMAX_BRUTE+1):
        b  = brute_N(W, h)
        d1 = dp_N(W, h, literal_thresholds=True)
        d2 = dp_N(W, h, literal_thresholds=False)
        cm = chain_N(W.M, W.E, h)
        b  = {k:v for k,v in b.items()  if v}
        d1 = {k:v for k,v in d1.items() if v}
        d2 = {k:v for k,v in d2.items() if v}
        cm = {k:v for k,v in cm.items() if v}
        ok = (b==d1==d2==cm)
        xcheck.append(dict(window=W.tag(), M=W.M, c=W.c, E=W.E, h=h, ok=bool(ok),
                           brute_eq_dp=bool(b==d1), dp_eq_closed=bool(d1==d2),
                           dp_eq_chain=bool(d1==cm)))
        if not ok:
            say("  MISMATCH", W.tag(), "h=",h); say("   brute",b); say("   dp   ",d1)
            say("   closed",d2); say("   chain",cm)
say("  windows tested: %d ; (window,h) pairs: %d ; all four agree: %s"
    % (len(wins), len(xcheck), all(x["ok"] for x in xcheck)))
results["crosscheck"] = dict(n_windows=len(wins), n_pairs=len(xcheck),
                             all_agree=all(x["ok"] for x in xcheck),
                             detail=xcheck)

# ---------------------------------------------------------------- 0.2 (M,E) collapse
say("== 0.2 window -> (M,E) collapse ==")
groups = {}
cand = []
for M in range(2,7):
    for c in (1,2,3):
        for pf in ([1]*(M-1), [1]*(M-2)+[2] if M>=2 else None, [2]+[1]*(M-2)):
            if pf is None or len(pf)!=M-1: continue
            for extra in (0,1,2):
                W = mkwin(M,c,tuple(pf),extra)
                if W: cand.append(W)
seen = {}
collapse_ok = True
for W in cand:
    key = (W.M, W.E)
    prof = tuple(sorted(chain_N(W.M,W.E,24).items()))
    prof_dp = tuple(sorted({k:v for k,v in dp_N(W,24).items() if v}.items()))
    if prof != prof_dp: collapse_ok = False; say("  chain!=dp for", W.tag())
    if key in seen and seen[key][1] != prof_dp:
        collapse_ok = False
        say("  COLLAPSE FAILS", W.tag(), "vs", seen[key][0].tag())
    seen.setdefault(key, (W, prof_dp))
    groups.setdefault(str(key), []).append(W.tag())
say("  %d candidate windows -> %d distinct (M,E) classes; collapse exact: %s"
    % (len(cand), len(seen), collapse_ok))
results["ME_collapse"] = dict(n_windows=len(cand), n_classes=len(seen), ok=collapse_ok,
                              classes={k:v for k,v in groups.items()})

# ------------------------------------------- 0.3 d* closed form / eps / no clipping
say("== 0.3 generalized d* closed form, global eps_h, clipping ==")
dstar_ok = True; eps_ok = True; min_dstar = 10**9; n_states = 0
for W in wins:
    for h in range(1, 26):
        CW = W.CW(h)
        # walk the DP collecting reachable states
        alive = {}
        V = W.Vmin
        while confined(W.P+V+h, W.M+h, W.c):
            alive[W.P+V] = alive.get(W.P+V,0)+1; V += 1
        j = 0
        while alive and j < h:
            nxt = {}
            for S in alive:
                d = 1
                while confined(S+d+(h-j-1), W.M+h, W.c): d += 1
                n_states += 1
                if d != CW+2+j-S: dstar_ok = False
                min_dstar = min(min_dstar, d)
                # global eps: d*(S,j,h+1)-d*(S,j,h) == eps_h
                d2 = 1
                while confined(S+d2+(h+1-j-1), W.M+h+1, W.c): d2 += 1
                if d2-d != W.eps(h): eps_ok = False
                for dd in range(1,d):
                    if not confined(S+dd, W.M+j+1, W.c):
                        if j+1 < h: nxt[S+dd] = 1
            alive = nxt; j += 1
    for h in range(0, 400):
        if W.eps(h) not in (0,1): eps_ok = False
        if W.eps(h) != G(W.M+h+1)-G(W.M+h): eps_ok = False
        if W.CW(h) != W.c - h + F(W.M+h): eps_ok = False
say("  reachable states checked: %d ; d* closed form exact: %s ; eps_h global/Sturmian: %s"
    % (n_states, dstar_ok, eps_ok))
say("  min d* over all reachable states = %d  (clipping d*<=0 never occurs; d*>=2)" % min_dstar)
results["dstar"] = dict(n_states=n_states, closed_form_exact=dstar_ok,
                        eps_global=eps_ok, min_dstar=min_dstar)
json.dump(results, open(ROOT+"/out/stage0a.json","w"), indent=1)
open(ROOT+"/out/stage0a.log","w").write("\n".join(LOG))
