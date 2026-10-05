"""Fit the EXPLORE theta-curve and emit the frozen numeric predictions for HOLDOUT."""
import os as _os, sys as _sys
ROOT = _os.path.dirname(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
_sys.path.insert(0, _os.path.join(ROOT, 'shape','code'))

import sys, csv, math, collections, statistics, json
import numpy as np
sys.path.insert(0,ROOT+"/shape/code")
from slopes import build, split
S=build(depth=60)
for s in S: s.verify(460)
EX,HO=split(S); byname={s.name:s for s in S}
rows=[]
with open(ROOT+"/shape/data/explore_shape2.csv") as f:
    for d in csv.DictReader(f):
        d["K"]=int(d["K"]); d["r"]=int(d["r"]); d["S2r"]=float(d["S2r"]); d["S3r"]=float(d["S3r"])
        rows.append(d)
BIG=[d for d in rows if d["K"]>=20]
per=collections.defaultdict(list)
for d in BIG: per[d["slope"]].append(d)
x,yv,yl=[],[],[]
for n,v in per.items():
    x.append(byname[n].float_value())
    yv.append(statistics.mean([d["S3r"] for d in v]))
    yl.append(statistics.mean([d["S2r"]/math.sqrt(d["r"]) for d in v]))
x=np.array(x); cv=np.polyfit(x,np.array(yv),3); cl=np.polyfit(x,np.array(yl),3)
out=dict(cubic_U_var=list(cv), cubic_U_loc=list(cl),
         explore_resid_sd_U_var=float((np.array(yv)-np.polyval(cv,x)).std()),
         explore_resid_sd_U_loc=float((np.array(yl)-np.polyval(cl,x)).std()),
         theta_range=[float(x.min()),float(x.max())], predictions={})
for s in HO:
    t=s.float_value()
    out["predictions"][s.name]=dict(theta=t, family=s.family,
        U_var_pred=float(np.polyval(cv,t)), U_loc_pred=float(np.polyval(cl,t)))
json.dump(out,open(ROOT+"/shape/out/frozen_curve.json","w"),indent=1)
print("cubic U_var:", " ".join("%.5f"%c for c in cv))
print("cubic U_loc:", " ".join("%.5f"%c for c in cl))
print("explore residual sd: U_var %.4f  U_loc %.4f"%(out["explore_resid_sd_U_var"],out["explore_resid_sd_U_loc"]))
print("theta range fitted: %.4f .. %.4f"%tuple(out["theta_range"]))
for n,p in sorted(out["predictions"].items(), key=lambda kv: kv[1]["theta"]):
    print("  %-26s theta=%.6f  U_var_pred=%+.4f  U_loc_pred=%+.4f  [%s]"
          %(n,p["theta"],p["U_var_pred"],p["U_loc_pred"],p["family"]))
