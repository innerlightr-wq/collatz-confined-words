"""Stage 6: is the extended-horizon drift a UNIVERSAL finite-size law in r,
or does it separate the slopes?  Fit U_var = cubic(theta) + a/sqrt(r) + b/r on EXPLORE
(h<=150, r<=~100) and predict the HOLDOUT extended rows (h=160..400, r up to ~250)."""
import os as _os, sys as _sys
ROOT = _os.path.dirname(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
_sys.path.insert(0, _os.path.join(ROOT, 'shape','code'))

import sys, csv, math, collections, statistics
import numpy as np
sys.path.insert(0,ROOT+"/shape/code")
from slopes import build, split
LOG=[]
def say(*a):
    s=" ".join(str(x) for x in a); print(s); LOG.append(s)
S=build(depth=60); 
for s in S: s.verify(460)
EX,HO=split(S); byname={s.name:s for s in S}
def load(path, extflag=False):
    out=[]
    with open(path) as f:
        for d in csv.DictReader(f):
            for k in ("K","r","h","M"): d[k]=int(d[k])
            for k in ("S2r","S3r"): d[k]=float(d[k])
            d["theta"]=byname[d["slope"]].float_value()
            d["ext"]= (d.get("ext")=="True")
            if d["K"]>=20 and d["r"]>0:
                d["U_var"]=d["S3r"]; d["U_loc"]=d["S2r"]/math.sqrt(d["r"]); out.append(d)
    return out
E=load(ROOT+"/shape/data/explore_shape2.csv")
H=load(ROOT+"/shape/data/holdout_shape.csv")
LO,HI=0.1921,0.7183
Hin=[d for d in H if LO<=d["theta"]<=HI]
def X(rows):
    return np.array([[1,d["theta"],d["theta"]**2,d["theta"]**3,
                      1/math.sqrt(d["r"]),1/d["r"]] for d in rows])
for key in ("U_var","U_loc"):
    b,*_=np.linalg.lstsq(X(E),np.array([d[key] for d in E]),rcond=None)
    say("== %s : single law  f(theta) + a/sqrt(r) + b/r  fitted on EXPLORE (h<=150) =="%key)
    say("   a = %+.4f, b = %+.4f"%(b[4],b[5]))
    for lab,sel in (("HOLDOUT h<=150",lambda d: not d["ext"]),
                    ("HOLDOUT h=160..400",lambda d: d["ext"])):
        rows=[d for d in Hin if sel(d)]
        if not rows: continue
        pred=X(rows)@b; res=np.array([d[key] for d in rows])-pred
        per=collections.defaultdict(list)
        for d,e in zip(rows,res): per[d["slope"]].append(e)
        mu=[np.mean(v) for v in per.values()]
        say("   %-20s n=%5d  rmse %.4f | per-slope mean resid: %+.4f +- %.4f (max|.| %.4f)"
            %(lab,len(rows),res.std(),np.mean(mu),np.std(mu),max(abs(m) for m in mu)))
    say("   -> one r-correction, shared by every slope, carries the extended horizons:")
say("== 6.2 r-range actually probed ==")
say("   EXPLORE r: %d..%d ; HOLDOUT h<=150 r: %d..%d ; HOLDOUT extended r: %d..%d"
    %(min(d["r"] for d in E),max(d["r"] for d in E),
      min(d["r"] for d in Hin if not d["ext"]),max(d["r"] for d in Hin if not d["ext"]),
      min(d["r"] for d in Hin if d["ext"]),max(d["r"] for d in Hin if d["ext"])))
say("== 6.3 limiting values implied by the fit (r -> infinity) ==")
b,*_=np.linalg.lstsq(X(E),np.array([d["U_var"] for d in E]),rcond=None)
for t in (0.2,0.4,0.584963,0.618034,0.7):
    say("   theta=%.6f : Var/r -> %.4f"%(t,b[0]+b[1]*t+b[2]*t*t+b[3]*t**3))
say("   (at theta -> 0 the fit gives %.4f; the idealized NB(K,1/2) value is 2)"
    %(b[0]))
open(ROOT+"/shape/out/stage6.log","w").write("\n".join(LOG))
