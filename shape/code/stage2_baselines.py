"""Stage 2 (EXPLORE): do the shape statistics collapse to a universal curve in the
intrinsic parameter r = ell_0, and do CF features add anything beyond local phase /
discrepancy?  Nested OLS + between-slope variance decomposition."""
import os as _os, sys as _sys
ROOT = _os.path.dirname(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
_sys.path.insert(0, _os.path.join(ROOT, 'shape','code'))

import sys, csv, json, math, collections, statistics, random
import numpy as np
sys.path.insert(0,ROOT+"/shape/code")
from slopes import build, split
LOG=[]
def say(*a):
    s=" ".join(str(x) for x in a); print(s); LOG.append(s)

S=build(depth=60)
for s in S: s.verify(460)
EX,HO=split(S)
byname={s.name:s for s in S}
R=[]
with open(ROOT+"/shape/data/explore_shape2.csv") as f:
    for d in csv.DictReader(f):
        d["K"]=int(d["K"]); d["r"]=int(d["r"]); d["h"]=int(d["h"]); d["M"]=int(d["M"])
        for k in ("S2r","S3r","frac_h","frac_M","mean"): d[k]=float(d[k])
        d["S4m"]=float(d["S4m"]) if d["S4m"] else None
        d["mode_tied"]= d["mode_tied"]=="True"
        R.append(d)
BIG=[d for d in R if d["K"]>=20]
say("rows with K>=20: %d over %d slopes"%(len(BIG),len({d['slope'] for d in BIG})))

say("== 2.1 universality of the single untied mode ==")
for lo in (1,5,10,15,17,20):
    sub=[d for d in R if d["K"]>=lo]
    t=sum(d["mode_tied"] for d in sub)
    say("   K >= %2d : %5d rows, %4d tied (%.3f%%)"%(lo,len(sub),t,100*t/len(sub)))
say("   -> no tie anywhere with K >= 17 on EXPLORE; ties are a small-K effect, not a slope effect.")

say("== 2.2 normalized statistics (intrinsic parameter r = ell_0) ==")
for d in BIG:
    d["U_loc"]=d["S2r"]/math.sqrt(d["r"])      # location offset / sqrt(r)
    d["U_var"]=d["S3r"]                        # Var / r
    if d["S4m"] is not None: d["U_tv"]=d["S4m"]*math.sqrt(d["r"])
for nm,key in (("U_loc = (mean-r)/sqrt(r)","U_loc"),("U_var = Var/r","U_var")):
    say("   %s"%nm)
    for lo,hi in ((20,40),(40,60),(60,90),(90,200)):
        sub=[d for d in BIG if lo<=d["r"]<hi]
        if len(sub)<20: continue
        per=collections.defaultdict(list)
        for d in sub: per[d["slope"]].append(d[key])
        mu=[statistics.mean(v) for v in per.values()]
        within=statistics.mean([statistics.pstdev(v) for v in per.values() if len(v)>2])
        say("      r %3d-%3d n=%4d  pooled %.4f+-%.4f | between-slope sd %.4f, within-slope sd %.4f"
            %(lo,hi,len(sub),statistics.mean([d[key] for d in sub]),
              statistics.pstdev([d[key] for d in sub]),statistics.pstdev(mu),within))

say("== 2.3 nested OLS on EXPLORE (target: U_loc, U_var, U_tv) ==")
def feats(d, kind):
    s=byname[d["slope"]]; N=d["M"]+d["h"]
    f=[1.0, 1.0/math.sqrt(d["r"]), 1.0/d["r"]]                       # B3 universal
    if kind=="B3": return f
    ph=d["frac_h"]; pm=d["frac_M"]
    f+= [math.cos(2*math.pi*ph),math.sin(2*math.pi*ph),
         math.cos(4*math.pi*ph),math.sin(4*math.pi*ph),
         math.cos(2*math.pi*pm),math.sin(2*math.pi*pm),
         float(d["r"]%2), float(d["r"]%3==0)]                        # B1 local phase
    if kind=="B1": return f
    f+=[s.star_discrepancy(min(N,400))*math.sqrt(N)]                 # B2 discrepancy
    if kind=="B2": return f
    pq,dep=s.pq_upto(N)                                              # CF features
    mx=max(pq) if pq else 1; sm=sum(pq) if pq else 1
    qs=s.conv_denoms(N); dist=min(abs(N-q) for q in qs)
    f+=[math.log(mx),math.log(sm),float(dep),math.log1p(dist),
        math.log1p(min(abs(d["h"]-q) for q in qs))]
    return f
def ols(rows,key,kind):
    X=np.array([feats(d,kind) for d in rows]); y=np.array([d[key] for d in rows])
    b,*_=np.linalg.lstsq(X,y,rcond=None); res=y-X@b
    return res, 1-res.var()/y.var()
targets=[("U_loc",BIG),("U_var",BIG),("U_tv",[d for d in BIG if "U_tv" in d])]
res_store={}
for key,rows in targets:
    line=[]
    for kind in ("B3","B1","B2","CF"):
        res,r2=ols(rows,key,kind); line.append("%s R2=%.4f rmse=%.4f"%(kind,r2,res.std()))
        res_store[(key,kind)]=(rows,res)
    say("   %-6s n=%4d | %s"%(key,len(rows)," | ".join(line)))

say("== 2.4 does slope IDENTITY still explain residual variance after each model? ==")
say("   (between-slope variance of residuals / total residual variance; 0 = fully explained)")
for key,_ in targets:
    line=[]
    for kind in ("B3","B1","B2","CF"):
        rows,res=res_store[(key,kind)]
        per=collections.defaultdict(list)
        for d,e in zip(rows,res): per[d["slope"]].append(e)
        mu=[np.mean(v) for v in per.values()]
        line.append("%s %.3f"%(kind,np.var(mu)/np.var(res)))
    say("   %-6s %s"%(key," | ".join(line)))

say("== 2.5 do CF features correlate with the per-slope residual of the LOCAL-PHASE model? ==")
for key,_ in targets:
    rows,res=res_store[(key,"B1")]
    per=collections.defaultdict(list)
    for d,e in zip(rows,res): per[d["slope"]].append(e)
    names=sorted(per); mu=np.array([np.mean(per[n]) for n in names])
    out=[]
    for fname,fn in (("log max a_k",lambda s: math.log(max(s.pq_upto(300)[0] or [1]))),
                     ("log sum a_k",lambda s: math.log(sum(s.pq_upto(300)[0] or [1]))),
                     ("discrepancy",lambda s: s.star_discrepancy(300)),
                     ("theta",lambda s: s.float_value())):
        x=np.array([fn(byname[n]) for n in names])
        if x.std()==0: continue
        c=np.corrcoef(x,mu)[0,1]
        # permutation p-value
        rng=random.Random(99); ge=0
        for _ in range(2000):
            p=list(mu); rng.shuffle(p)
            if abs(np.corrcoef(x,p)[0,1])>=abs(c): ge+=1
        out.append("%s r=%+.3f p=%.3f"%(fname,c,(ge+1)/2001))
    say("   %-6s %s"%(key," | ".join(out)))
json.dump({"n_big":len(BIG)},open(ROOT+"/shape/out/stage2.json","w"))
open(ROOT+"/shape/out/stage2.log","w").write("\n".join(LOG))
