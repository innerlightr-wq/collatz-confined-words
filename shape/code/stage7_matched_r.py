"""Stage 7: at MATCHED r, is a single cubic in theta enough for every slope?
(This removes the finite-size term without assuming its functional form.)"""
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
S=build(depth=60)
for s in S: s.verify(460)
EX,HO=split(S); byname={s.name:s for s in S}
def load(p):
    out=[]
    with open(p) as f:
        for d in csv.DictReader(f):
            for k in ("K","r","h","M"): d[k]=int(d[k])
            for k in ("S2r","S3r"): d[k]=float(d[k])
            d["theta"]=byname[d["slope"]].float_value()
            if d["K"]>=20 and d["r"]>0:
                d["U_var"]=d["S3r"]; d["U_loc"]=d["S2r"]/math.sqrt(d["r"]); out.append(d)
    return out
E=load(ROOT+"/shape/data/explore_shape2.csv")
H=[d for d in load(ROOT+"/shape/data/holdout_shape.csv")
   if 0.1921<=d["theta"]<=0.7183]
BINS=[(20,30),(30,40),(40,55),(55,75),(75,105),(105,150),(150,260)]
for key in ("U_var","U_loc"):
    say("== %s : one cubic in theta per r-bin =="%key)
    say("   r-bin     EXPLORE->HOLDOUT (out-of-sample)        HOLDOUT-only (slope agreement)")
    for lo,hi in BINS:
        e=[d for d in E if lo<=d["r"]<hi]; h=[d for d in H if lo<=d["r"]<hi]
        line="   %3d-%3d  "%(lo,hi)
        if len(e)>=60 and len(h)>=60:
            pe=collections.defaultdict(list); 
            for d in e: pe[d["slope"]].append(d[key])
            x=np.array([byname[n].float_value() for n in pe]); y=np.array([np.mean(v) for v in pe.values()])
            c=np.polyfit(x,y,min(3,len(x)-1))
            ph=collections.defaultdict(list)
            for d in h: ph[d["slope"]].append(d[key])
            xh=np.array([byname[n].float_value() for n in ph]); yh=np.array([np.mean(v) for v in ph.values()])
            res=yh-np.polyval(c,xh)
            line+="n_sl=%2d rmse %.4f max %.4f   "%(len(xh),res.std(),max(abs(res)))
        else:
            line+="  (insufficient overlap)            "
        if len(h)>=60:
            ph=collections.defaultdict(list)
            for d in h: ph[d["slope"]].append(d[key])
            if len(ph)>=6:
                xh=np.array([byname[n].float_value() for n in ph]); yh=np.array([np.mean(v) for v in ph.values()])
                c2=np.polyfit(xh,yh,3); r2=yh-np.polyval(c2,xh)
                line+="n_sl=%2d  resid sd %.4f  (df=%d)"%(len(xh),r2.std(),len(xh)-4)
        say(line)
say("== 7.2 the four near-equal-theta slopes, at matched r ==")
quad=["golden [0;1,1,1,...]","rational 144/233","rational 55/89","liouville a_7=500"]
for lo,hi in ((20,40),(40,70),(70,110),(110,260)):
    vals={}
    for n in quad:
        v=[d["U_var"] for d in H if d["slope"]==n and lo<=d["r"]<hi]
        if v: vals[n]=(statistics.mean(v),len(v))
    if len(vals)>=3:
        say("   r %3d-%3d : %s"%(lo,hi," | ".join("%s %.5f(n=%d)"%(n.split()[0][:9],a,k)
                                                   for n,(a,k) in vals.items())))
        say("              spread = %.5f"%(max(a for a,_ in vals.values())-min(a for a,_ in vals.values())))
open(ROOT+"/shape/out/stage7.log","w").write("\n".join(LOG))
