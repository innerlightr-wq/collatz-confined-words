"""Stage 2b (EXPLORE): with theta itself in the baseline, do CF features add anything?
Plus the sharpest control: rational p/q vs an irrational of nearly the same theta."""
import os as _os, sys as _sys
ROOT = _os.path.dirname(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
_sys.path.insert(0, _os.path.join(ROOT, 'shape','code'))

import sys, csv, math, collections, statistics, random, json
import numpy as np
sys.path.insert(0,ROOT+"/shape/code")
from slopes import build, split
LOG=[]
def say(*a):
    s=" ".join(str(x) for x in a); print(s); LOG.append(s)
S=build(depth=60)
for s in S: s.verify(460)
EX,HO=split(S); byname={s.name:s for s in S}
R=[]
with open(ROOT+"/shape/data/explore_shape2.csv") as f:
    for d in csv.DictReader(f):
        d["K"]=int(d["K"]); d["r"]=int(d["r"]); d["h"]=int(d["h"]); d["M"]=int(d["M"])
        for k in ("S2r","S3r","frac_h","frac_M"): d[k]=float(d[k])
        d["S4m"]=float(d["S4m"]) if d["S4m"] else None
        R.append(d)
BIG=[d for d in R if d["K"]>=20]
for d in BIG:
    d["theta"]=byname[d["slope"]].float_value()
    d["U_loc"]=d["S2r"]/math.sqrt(d["r"]); d["U_var"]=d["S3r"]
    if d["S4m"] is not None: d["U_tv"]=d["S4m"]*math.sqrt(d["r"])

def feats(d,kind):
    s=byname[d["slope"]]; N=d["M"]+d["h"]; t=d["theta"]; ir=1.0/math.sqrt(d["r"])
    f=[1.0, ir, 1.0/d["r"], t, t*t, t*ir]                      # B3t: smooth in (r,theta)
    if kind=="B3t": return f
    ph=d["frac_h"]; pm=d["frac_M"]
    f+=[math.cos(2*math.pi*ph),math.sin(2*math.pi*ph),math.cos(4*math.pi*ph),
        math.sin(4*math.pi*ph),math.cos(2*math.pi*pm),math.sin(2*math.pi*pm),
        float(d["r"]%2),float(d["r"]%3==0)]
    if kind=="B1t": return f
    f+=[s.star_discrepancy(min(N,400))*math.sqrt(N)]
    if kind=="B2t": return f
    pq,dep=s.pq_upto(N); mx=max(pq) if pq else 1; sm=sum(pq) if pq else 1
    qs=s.conv_denoms(N)
    f+=[math.log(mx),math.log(sm),float(dep),
        math.log1p(min(abs(N-q) for q in qs)),math.log1p(min(abs(d["h"]-q) for q in qs))]
    return f
def ols(rows,key,kind):
    X=np.array([feats(d,kind) for d in rows]); y=np.array([d[key] for d in rows])
    b,*_=np.linalg.lstsq(X,y,rcond=None); res=y-X@b
    return res,1-res.var()/y.var()
say("== 2.6 nested models WITH theta in the baseline (EXPLORE) ==")
targets=[("U_loc",BIG),("U_var",BIG),("U_tv",[d for d in BIG if "U_tv" in d])]
store={}
for key,rows in targets:
    line=[]
    for kind in ("B3t","B1t","B2t","CFt"):
        res,r2=ols(rows,key,kind); store[(key,kind)]=(rows,res)
        line.append("%s R2=%.4f rmse=%.4f"%(kind,r2,res.std()))
    say("   %-6s n=%4d | %s"%(key,len(rows)," | ".join(line)))
say("   (compare with Stage 2.3, where theta was ABSENT and the CF block reached R2 0.43-0.55:")
say("    that gain was the CF block proxying for theta, not Diophantine structure.)")

say("== 2.7 CF features vs per-slope residuals of B1t (theta + local phase) ==")
for key,_ in targets:
    rows,res=store[(key,"B1t")]
    per=collections.defaultdict(list)
    for d,e in zip(rows,res): per[d["slope"]].append(e)
    names=sorted(per); mu=np.array([np.mean(per[n]) for n in names])
    out=[]
    for fn_name,fn in (("log max a_k",lambda s: math.log(max(s.pq_upto(300)[0] or [1]))),
                       ("log sum a_k",lambda s: math.log(sum(s.pq_upto(300)[0] or [1]))),
                       ("discrepancy",lambda s: s.star_discrepancy(300)),
                       ("theta",lambda s: s.float_value())):
        x=np.array([fn(byname[n]) for n in names])
        if x.std()==0: continue
        c=np.corrcoef(x,mu)[0,1]
        rng=random.Random(7); ge=0
        for _ in range(4000):
            p=list(mu); rng.shuffle(p)
            if abs(np.corrcoef(x,p)[0,1])>=abs(c): ge+=1
        out.append("%s r=%+.3f p=%.4f"%(fn_name,c,(ge+1)/4001))
    say("   %-6s %s"%(key," | ".join(out)))

say("== 2.8 sharpest control: RATIONAL vs IRRATIONAL at nearly the same theta ==")
say("   (rational p/q is maximally non-Diophantine; golden is the hardest to approximate)")
pairs=[("rational 8/13","golden [0;1,1,1,...]"),("rational 55/89","golden [0;1,1,1,...]")]
for a,b in pairs:
    if a not in {d['slope'] for d in BIG} or b not in {d['slope'] for d in BIG}:
        say("   (%s or %s is in HOLDOUT -- skipped here)"%(a,b)); continue
    for key in ("U_loc","U_var"):
        va=[d[key] for d in BIG if d["slope"]==a]; vb=[d[key] for d in BIG if d["slope"]==b]
        say("   %-16s theta=%.6f %s mean %.4f (n=%d) | %-22s theta=%.6f mean %.4f (n=%d)"
            %(a,byname[a].float_value(),key,statistics.mean(va),len(va),
              b,byname[b].float_value(),statistics.mean(vb),len(vb)))
say("== 2.9 per-slope means vs theta: is the dependence a smooth curve? ==")
for key in ("U_loc","U_var"):
    per=collections.defaultdict(list)
    for d in BIG: per[d["slope"]].append(d[key])
    pts=sorted((byname[n].float_value(),statistics.mean(v),n,byname[n].family) for n,v in per.items())
    say("   %s:"%key)
    for t,m,n,fam in pts: say("      theta=%.6f  %-24s %-9s mean=%+.4f"%(t,n,fam,m))
    x=np.array([p[0] for p in pts]); y=np.array([p[1] for p in pts])
    for deg in (1,2,3):
        c=np.polyfit(x,y,deg); r=y-np.polyval(c,x)
        say("      poly deg %d in theta: R2=%.4f  residual sd=%.4f"%(deg,1-r.var()/y.var(),r.std()))
open(ROOT+"/shape/out/stage2b.log","w").write("\n".join(LOG))
