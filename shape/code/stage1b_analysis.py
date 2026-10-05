"""Stage 1b (EXPLORE only): convention-free statistics, the K-collapse, and tie anatomy."""
import os as _os, sys as _sys
ROOT = _os.path.dirname(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
_sys.path.insert(0, _os.path.join(ROOT, 'shape','code'))

import sys, csv, json, collections, statistics
sys.path.insert(0,ROOT+"/shape/code")
from slopes import build, split
from shape_core import row2
LOG=[]
def say(*a):
    s=" ".join(str(x) for x in a); print(s); LOG.append(s)
S=build(depth=60)
for s in S: s.verify(460)
EX,HO=split(S)
WINDOWS=[(2,1),(3,1),(5,2),(8,1),(13,3)]
HMAX=150
rows=[]
for sl in EX:
    for (M,E) in WINDOWS:
        for h in range(1,HMAX+1):
            r=row2(sl,M,E,h,want_tv=(h%25==0))
            if r: rows.append(r)
say("EXPLORE rows: %d (slopes %d, windows %d, h<=%d)"%(len(rows),len(EX),len(WINDOWS),HMAX))
fields=["slope","family","M","E","h","K","r","eta","mode","mode_tied","n_modes","S1","S1r",
        "S2","S2r","S3","S3r","S4","S4m","kappa","mean","var","frac_h","frac_M","plateau",
        "mode_minus_roundmean","w_top","w_mid"]
with open(ROOT+"/shape/data/explore_shape2.csv","w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=fields,extrasaction="ignore"); w.writeheader()
    for r in rows: w.writerow(r)

big=[r for r in rows if r["K"]>=20]
say("== 1.6 does the shape collapse as a function of K (baseline B3: one universal curve)? ==")
say("    K-bin    n   S3=Var/K  mean+-sd        S2r=mean-r  mean+-sd      spread across slopes")
for lo,hi in ((20,30),(30,45),(45,60),(60,80),(80,120)):
    sub=[r for r in big if lo<=r["K"]<hi]
    if len(sub)<5: continue
    a=[r["S3"] for r in sub]; b=[r["S2r"] for r in sub]
    byslope=collections.defaultdict(list)
    for r in sub: byslope[r["slope"]].append(r["S3"])
    spread=max(statistics.mean(v) for v in byslope.values())-min(statistics.mean(v) for v in byslope.values())
    say("   %3d-%3d %5d   %.4f+-%.4f   %+.4f+-%.4f   per-slope mean S3 spread %.4f"
        %(lo,hi,len(sub),statistics.mean(a),statistics.pstdev(a),
          statistics.mean(b),statistics.pstdev(b),spread))
say("== 1.7 ties: where does the Haar mode fail to be unique? ==")
tied=[r for r in rows if r["mode_tied"]]
say("   tied rows: %d / %d (%.3f%%)"%(len(tied),len(rows),100*len(tied)/len(rows)))
say("   by family: %s"%dict(collections.Counter(r["family"] for r in tied)))
say("   by K: %s"%dict(sorted(collections.Counter(min(r["K"],30) for r in tied).items())))
say("   by slope (top 8): %s"%collections.Counter(r["slope"] for r in tied).most_common(8))
tbig=[r for r in tied if r["K"]>=20]
say("   tied rows with K>=20: %d (of %d rows with K>=20)"%(len(tbig),len(big)))
for r in tbig[:10]:
    say("      %-22s (M,E)=(%2d,%d) h=%3d K=%2d r=%2d mode=%2d w_top=%s"
        %(r["slope"],r["M"],r["E"],r["h"],r["K"],r["r"],r["mode"],r["w_top"]))
say("== 1.8 mode relative to its own mean (convention-free) ==")
c=collections.Counter(r["mode_minus_roundmean"] for r in big)
say("   mode - round(mean), K>=20: %s"%dict(sorted(c.items())))
c2=collections.Counter(r["S1r"] for r in big)
say("   S1r = mode - r,      K>=20: %s"%dict(sorted(c2.items())))
say("== 1.9 S4 / S4m (TV to the negative-binomial kernel) vs K ==")
tv=[r for r in rows if "S4m" in r and r["K"]>=10]
for lo,hi in ((10,20),(20,40),(40,70),(70,120)):
    sub=[r for r in tv if lo<=r["K"]<hi]
    if len(sub)<3: continue
    say("   K %3d-%3d  n=%3d  S4(study K) %.4f+-%.4f   S4m(mean-matched) %.4f+-%.4f"
        %(lo,hi,len(sub),statistics.mean([r["S4"] for r in sub]),statistics.pstdev([r["S4"] for r in sub]),
          statistics.mean([r["S4m"] for r in sub]),statistics.pstdev([r["S4m"] for r in sub])))
json.dump(dict(rows=len(rows),tied=len(tied),tied_bigK=len(tbig)),
          open(ROOT+"/shape/out/stage1b.json","w"),indent=1)
open(ROOT+"/shape/out/stage1b.log","w").write("\n".join(LOG))
