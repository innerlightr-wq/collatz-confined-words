"""Stage 4: frozen rules P1-P5 on HOLDOUT slopes and at extended horizons."""
import os as _os, sys as _sys
ROOT = _os.path.dirname(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
_sys.path.insert(0, _os.path.join(ROOT, 'shape','code'))

import sys, csv, json, math, collections, statistics, random
import numpy as np
sys.path.insert(0,ROOT+"/shape/code")
from slopes import build, split
from shape_core import row2
LOG=[]
def say(*a):
    s=" ".join(str(x) for x in a); print(s); LOG.append(s)
FR=json.load(open(ROOT+"/shape/out/frozen_curve.json"))
cv=np.array(FR["cubic_U_var"]); cl=np.array(FR["cubic_U_loc"]); LO,HI=FR["theta_range"]
S=build(depth=60)
for s in S: s.verify(460)
EX,HO=split(S); byname={s.name:s for s in S}
say("HOLDOUT slopes (%d): %s"%(len(HO),", ".join(s.name for s in HO)))
WINDOWS=[(2,1),(3,1),(5,2),(8,1),(13,3)]
rows=[]
for sl in HO:
    for (M,E) in WINDOWS:
        for h in list(range(1,151))+list(range(160,401,10)):
            r=row2(sl,M,E,h,want_tv=(h%50==0))
            if r:
                r["theta"]=sl.float_value(); r["ext"]= h>150
                rows.append(r)
say("HOLDOUT rows: %d (%d at extended horizons h=160..400)"%(len(rows),sum(r["ext"] for r in rows)))
with open(ROOT+"/shape/data/holdout_shape.csv","w",newline="") as f:
    fl=["slope","family","theta","M","E","h","ext","K","r","mode","mode_tied","n_modes",
        "S1r","S2r","S3r","S4m","mean","var","frac_h","frac_M"]
    w=csv.DictWriter(f,fieldnames=fl,extrasaction="ignore"); w.writeheader()
    for r in rows: w.writerow(r)

say("== 4.1 P1: untied mode for K >= 17 ==")
for lo in (1,5,10,15,17,20):
    sub=[r for r in rows if r["K"]>=lo]; t=sum(r["mode_tied"] for r in sub)
    say("   K >= %2d : %5d rows, %4d tied (%.3f%%)"%(lo,len(sub),t,100*t/len(sub)))
ties=[r for r in rows if r["K"]>=17 and r["mode_tied"]]
ext=[r for r in rows if r["ext"] and r["K"]>=17]
say("   P1 verdict: ties at K>=17 : %d  -> %s"%(len(ties),"HOLDS" if not ties else "FALSIFIED"))
say("   (extended horizons h=160..400 with K>=17: %d rows, %d ties)"
    %(len(ext),sum(r["mode_tied"] for r in ext)))

BIG=[r for r in rows if r["K"]>=20]
for r in BIG:
    r["U_var"]=r["S3r"]; r["U_loc"]=r["S2r"]/math.sqrt(r["r"])
say("== 4.2 P2: frozen theta-curve vs measured, per HOLDOUT slope (K>=20) ==")
say("   slope                      theta     U_var pred  meas   resid | U_loc pred  meas   resid")
per=collections.defaultdict(list)
for r in BIG: per[r["slope"]].append(r)
inr,outr=[],[]
for n in sorted(per,key=lambda n: byname[n].float_value()):
    v=per[n]; t=byname[n].float_value()
    mv=statistics.mean([d["U_var"] for d in v]); ml=statistics.mean([d["U_loc"] for d in v])
    pv,pl=float(np.polyval(cv,t)),float(np.polyval(cl,t))
    tag="" if LO<=t<=HI else "  <-- EXTRAPOLATION"
    (inr if LO<=t<=HI else outr).append((n,mv-pv,ml-pl))
    say("   %-26s %.6f  %+.4f %+.4f %+.4f | %+.4f %+.4f %+.4f%s"
        %(n,t,pv,mv,mv-pv,pl,ml,ml-pl,tag))
rv=[d[1] for d in inr]; rl=[d[2] for d in inr]
say("   in-range slopes (%d): U_var resid sd %.4f max|.| %.4f  (predicted <=0.02 / <=0.08)"
    %(len(inr),statistics.pstdev(rv),max(abs(x) for x in rv)))
say("                        U_loc resid sd %.4f max|.| %.4f  (predicted <=0.03 / <=0.08)"
    %(statistics.pstdev(rl),max(abs(x) for x in rl)))
say("   P2 verdict: %s"%("HOLDS" if statistics.pstdev(rv)<=0.02 and max(abs(x) for x in rv)<=0.08
                         and statistics.pstdev(rl)<=0.03 and max(abs(x) for x in rl)<=0.08
                         else "FALSIFIED"))
for n,a,b in outr: say("   extrapolation slope %-22s U_var resid %+.4f  U_loc resid %+.4f"%(n,a,b))

say("== 4.3 P4: rational vs badly-approximable vs Liouville at essentially equal theta ==")
quad=["golden [0;1,1,1,...]","rational 144/233","rational 55/89","liouville a_7=500"]
vals={}
for n in quad:
    if n in per:
        vals[n]=(byname[n].float_value(),
                 statistics.mean([d["U_var"] for d in per[n]]),
                 statistics.mean([d["U_loc"] for d in per[n]]))
for n,(t,a,b) in sorted(vals.items(),key=lambda kv: -kv[1][0]):
    say("   %-24s theta=%.7f  U_var=%+.5f  U_loc=%+.5f  [%s]"%(n,t,a,b,byname[n].family))
if len(vals)>=2:
    sv=[v[1] for v in vals.values()]; sl_=[v[2] for v in vals.values()]
    say("   spread across these %d slopes: U_var %.5f, U_loc %.5f  (predicted U_var < 0.01)"
        %(len(vals),max(sv)-min(sv),max(sl_)-min(sl_)))
    say("   P4 verdict: %s"%("HOLDS" if max(sv)-min(sv)<0.01 else "FALSIFIED"))

say("== 4.4 P3: does the CF block add anything out of sample? ==")
EXR=[]
with open(ROOT+"/shape/data/explore_shape2.csv") as f:
    for d in csv.DictReader(f):
        d["K"]=int(d["K"]); d["r"]=int(d["r"]); d["h"]=int(d["h"]); d["M"]=int(d["M"])
        d["S2r"]=float(d["S2r"]); d["S3r"]=float(d["S3r"])
        d["frac_h"]=float(d["frac_h"]); d["frac_M"]=float(d["frac_M"])
        if d["K"]>=20:
            d["theta"]=byname[d["slope"]].float_value()
            d["U_var"]=d["S3r"]; d["U_loc"]=d["S2r"]/math.sqrt(d["r"]); EXR.append(d)
def feats(d,kind):
    s=byname[d["slope"]]; N=d["M"]+d["h"]; t=d["theta"]; ir=1.0/math.sqrt(d["r"])
    f=[1.0,ir,1.0/d["r"],t,t*t,t*ir,
       math.cos(2*math.pi*d["frac_h"]),math.sin(2*math.pi*d["frac_h"]),
       math.cos(4*math.pi*d["frac_h"]),math.sin(4*math.pi*d["frac_h"]),
       math.cos(2*math.pi*d["frac_M"]),math.sin(2*math.pi*d["frac_M"]),
       float(d["r"]%2),float(d["r"]%3==0)]
    if kind=="B1t": return f
    pq,dep=s.pq_upto(N); mx=max(pq) if pq else 1; sm=sum(pq) if pq else 1
    qs=s.conv_denoms(400)
    f+=[math.log(mx),math.log(sm),float(dep),
        math.log1p(min(abs(N-q) for q in qs)),math.log1p(min(abs(d["h"]-q) for q in qs))]
    return f
for key in ("U_var","U_loc"):
    line=[]
    for kind in ("B1t","CFt"):
        X=np.array([feats(d,kind) for d in EXR]); y=np.array([d[key] for d in EXR])
        b,*_=np.linalg.lstsq(X,y,rcond=None)
        Xh=np.array([feats(d,kind) for d in BIG]); yh=np.array([d[key] for d in BIG])
        rmse=float(np.sqrt(((yh-Xh@b)**2).mean())); line.append((kind,rmse))
    imp=100*(line[0][1]-line[1][1])/line[0][1]
    say("   %-6s out-of-sample RMSE  B1t %.4f -> CFt %.4f  (%+.1f%% ; predicted <5%% gain)"
        %(key,line[0][1],line[1][1],imp))
say("== 4.5 P3: CF features vs per-slope HOLDOUT residuals of B1t ==")
for key in ("U_var","U_loc"):
    X=np.array([feats(d,"B1t") for d in EXR]); y=np.array([d[key] for d in EXR])
    b,*_=np.linalg.lstsq(X,y,rcond=None)
    resid=collections.defaultdict(list)
    for d in BIG: resid[d["slope"]].append(d[key]-float(np.array(feats(d,"B1t"))@b))
    names=sorted(resid); mu=np.array([np.mean(resid[n]) for n in names])
    out=[]
    for fn_name,fn in (("log max a_k",lambda s: math.log(max(s.pq_upto(300)[0] or [1]))),
                       ("log sum a_k",lambda s: math.log(sum(s.pq_upto(300)[0] or [1]))),
                       ("discrepancy",lambda s: s.star_discrepancy(300))):
        x=np.array([fn(byname[n]) for n in names])
        if x.std()==0: continue
        c=float(np.corrcoef(x,mu)[0,1]); rng=random.Random(5); ge=0
        for _ in range(4000):
            p=list(mu); rng.shuffle(p)
            if abs(np.corrcoef(x,p)[0,1])>=abs(c): ge+=1
        out.append("%s r=%+.3f p=%.4f"%(fn_name,c,(ge+1)/4001))
    say("   %-6s %s"%(key," | ".join(out)))
say("== 4.6 P3: resonance test on HOLDOUT ==")
for key in ("U_var","U_loc"):
    X=np.array([feats(d,"B1t") for d in EXR]); y=np.array([d[key] for d in EXR])
    b,*_=np.linalg.lstsq(X,y,rcond=None)
    res=np.array([abs(d[key]-float(np.array(feats(d,"B1t"))@b)) for d in BIG])
    near,far=[],[]
    for d,e in zip(BIG,res):
        qs=byname[d["slope"]].conv_denoms(400)
        (near if min(abs(d["h"]-q) for q in qs)<=2 else far).append(e)
    obs=np.mean(near)-np.mean(far); rngl=random.Random(3); diffs=[]
    for _ in range(400):
        fake={n:[rngl.randint(1,400) for _ in byname[n].conv_denoms(400)] for n in {d["slope"] for d in BIG}}
        nn,ff=[],[]
        for d,e in zip(BIG,res):
            (nn if min(abs(d["h"]-q) for q in fake[d["slope"]])<=2 else ff).append(e)
        if nn and ff: diffs.append(np.mean(nn)-np.mean(ff))
    p=(sum(1 for v in diffs if abs(v)>=abs(obs))+1)/(len(diffs)+1)
    say("   %-6s near q_k (n=%4d) %.4f vs far (n=%4d) %.4f ; diff %+.4f ; null p=%.3f"
        %(key,len(near),np.mean(near),len(far),np.mean(far),obs,p))
open(ROOT+"/shape/out/stage4.log","w").write("\n".join(LOG))
