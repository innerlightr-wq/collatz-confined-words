"""Stage 8: confirmatory re-test of the uniqueness threshold on data that played no part
in setting it -- 10 NEW random-CF slopes (new seed), 3 NEW windows, horizons to 600.

Thresholds under test, fixed before this run (shape/REPORT.md, committed earlier):
    conjecture, in r : no tie for r >= 20
    observed, in K   : no tie for K >= 19
No threshold is adjusted below; the run only measures.

Speed: profiles are advanced with the PROVED operators (eta=0 -> identity, eta=1 -> T1),
which is O(h) per horizon instead of O(h^2), and spot-checked against a direct
recomputation of the chain count at sampled horizons.
"""
import os as _os, sys as _sys
ROOT = _os.path.dirname(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
_sys.path.insert(0, _os.path.join(ROOT, 'shape','code'))

import sys, csv, json, random, collections
sys.path.insert(0,ROOT+"/shape/code")
from slopes import Slope, build
from shape_core import profile, l_profile, K_of
LOG=[]
def say(*a):
    s=" ".join(str(x) for x in a); print(s); LOG.append(s)

HMAX=600; NMAX=HMAX+20
WINDOWS=[(4,2),(7,1),(11,2)]
SEED=20261006                      # NEW seed: none of these slopes appears in the earlier study

old=build(depth=60)
old_cfs={tuple(s.cf) for s in old if s.cf}
rng=random.Random(SEED)
NEW=[]
while len(NEW)<10:
    cf=[rng.choice([1,1,1,2,2,3,4,5,7]) for _ in range(60)]
    if tuple(cf) in old_cfs: continue
    s=Slope("new-random#%02d"%len(NEW),"random",cf)
    try: s.verify(NMAX)
    except AssertionError: continue
    NEW.append(s)
say("10 new random-CF slopes (seed %d), all exact for n <= %d, none shared with the earlier set:"%(SEED,NMAX))
for s in NEW: say("   %-16s theta=%.9f  CF head %s"%(s.name,s.float_value(),s.cf[:8]))
say("windows: %s ; horizons 1..%d"%(WINDOWS,HMAX))

def T1(P):
    out=[];run=0
    for j,v in enumerate(P):
        run+=v; out.append(run+1 if j==0 else run)
    return out

def mode_and_tie(P,h):
    W=[P[j]<<(h-j) for j in range(len(P))]
    mx=max(W); idx=[j for j,w in enumerate(W) if w==mx]
    return idx[0],len(idx)>1,len(idx)

rows=[];checks=0;check_fail=0
for sl in NEW:
    for (M,E) in WINDOWS:
        P=None
        for h in range(1,HMAX+1):
            r=sl.G(M+h)-sl.G(M)-E+1
            if r<1:
                P=None; continue
            if P is None:
                P=profile(sl,M,E,h)              # direct chain count to start
            else:
                eta_prev=sl.G(M+h)-sl.G(M+h-1)
                P=(P if eta_prev==0 else T1(P))+[0]
            if h%97==0 or h in (HMAX,):          # spot-check the incremental advance
                D=profile(sl,M,E,h); checks+=1
                if D[:len(P)]!=P[:len(D)] or sum(D)!=sum(P[:len(D)]):
                    check_fail+=1
                    say("   SPOT-CHECK FAILURE %s (%d,%d) h=%d"%(sl.name,M,E,h))
            m,tied,nm=mode_and_tie(P,h)
            rows.append(dict(slope=sl.name,theta=sl.float_value(),M=M,E=E,h=h,
                             r=r,K=K_of(sl,M,h),mode=m,tied=tied,n_modes=nm))
say("rows: %d ; incremental-advance spot-checks against the direct chain count: %d, failures: %d"
    %(len(rows),checks,check_fail))

with open(ROOT+"/shape/data/retest_shape.csv","w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

say("")
say("== 8.1 ties by window, reported in r ==")
say("   window   rows    r range     #ties   largest r with a tie   smallest tie-free r")
summary={}
for (M,E) in WINDOWS:
    sub=[d for d in rows if (d["M"],d["E"])==(M,E)]
    t=[d for d in sub if d["tied"]]
    rmax=max(d["r"] for d in t) if t else None
    free=(rmax+1) if t else 1
    summary["(%d,%d)"%(M,E)]=dict(rows=len(sub),ties=len(t),largest_tie_r=rmax,tie_free_from_r=free,
                                  r_min=min(d["r"] for d in sub),r_max=max(d["r"] for d in sub))
    say("   (%2d,%d)  %5d   %3d..%3d     %4d    %s                    r >= %s"
        %(M,E,len(sub),min(d["r"] for d in sub),max(d["r"] for d in sub),len(t),
          ("%3d"%rmax) if t else " none", free))
allt=[d for d in rows if d["tied"]]
say("   ALL      %5d   %3d..%3d     %4d    %s                    r >= %s"
    %(len(rows),min(d["r"] for d in rows),max(d["r"] for d in rows),len(allt),
      ("%3d"%max(d["r"] for d in allt)) if allt else " none",
      (max(d["r"] for d in allt)+1) if allt else 1))

say("")
say("== 8.2 the same ties expressed in K, for comparison with the earlier runs ==")
for (M,E) in WINDOWS:
    t=[d for d in rows if (d["M"],d["E"])==(M,E) and d["tied"]]
    say("   (%2d,%d) largest K with a tie: %s"%(M,E,max((d["K"] for d in t),default=None)))
say("   ALL    largest K with a tie: %s"%max((d["K"] for d in allt),default=None))

say("")
say("== 8.3 the thresholds under test (fixed before this run; NOT adjusted) ==")
for lab,key,thr in (("r >= 20 (conjecture, REPORT.md section 7)","r",20),
                    ("K >= 19 (observed in the earlier study)","K",19)):
    v=[d for d in rows if d[key]>=thr]
    bad=[d for d in v if d["tied"]]
    say("   %-42s : %5d rows, %d ties -> %s"%(lab,len(v),len(bad),"HOLDS" if not bad else "FALSIFIED"))
    if bad:
        for d in bad[:10]:
            say("        tie at %s (M,E)=(%d,%d) h=%d r=%d K=%d n_modes=%d"
                %(d["slope"],d["M"],d["E"],d["h"],d["r"],d["K"],d["n_modes"]))

say("")
say("== 8.4 tie frequency by r, pooled over the three windows ==")
for lo,hi in ((1,5),(5,10),(10,15),(15,20),(20,25),(25,40),(40,80),(80,400)):
    sub=[d for d in rows if lo<=d["r"]<hi]
    if not sub: continue
    t=sum(d["tied"] for d in sub)
    say("   r %3d-%3d : %5d rows, %4d tied (%.3f%%)"%(lo,hi,len(sub),t,100*t/len(sub)))

say("")
say("== 8.5 per-slope largest tie r (are ties concentrated in particular slopes?) ==")
for sl in NEW:
    t=[d for d in rows if d["slope"]==sl.name and d["tied"]]
    say("   %-16s theta=%.6f  ties=%4d  largest r with a tie: %s"
        %(sl.name,sl.float_value(),len(t),max((d["r"] for d in t),default=None)))
json.dump(dict(summary=summary,n_rows=len(rows),spot_checks=checks,spot_fail=check_fail,
               largest_tie_r_overall=(max(d["r"] for d in allt) if allt else None),
               largest_tie_K_overall=(max(d["K"] for d in allt) if allt else None)),
          open(ROOT+"/shape/out/stage8.json","w"),indent=1)
open(ROOT+"/shape/out/stage8.log","w").write("\n".join(LOG))
