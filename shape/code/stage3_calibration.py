"""Stage 3 (EXPLORE-side calibration): convergent-coincidence null, shuffled-CF control,
and a direct resonance test for 'h near a convergent denominator'."""
import os as _os, sys as _sys
ROOT = _os.path.dirname(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
_sys.path.insert(0, _os.path.join(ROOT, 'shape','code'))

import sys, csv, math, random, collections, statistics
import numpy as np
sys.path.insert(0,ROOT+"/shape/code")
from slopes import build, split, Slope
from shape_core import row2
LOG=[]
def say(*a):
    s=" ".join(str(x) for x in a); print(s); LOG.append(s)
S=build(depth=60)
for s in S: s.verify(460)
EX,HO=split(S); byname={s.name:s for s in S}

say("== 3.1 convergent-coincidence null ==")
say("   fraction of integers 1..N within tol of a convergent denominator of the slope")
for sl in EX[:6]:
    qs=sl.conv_denoms(400)
    for tol in (1,2):
        f=sum(1 for n in range(1,401) if min(abs(n-q) for q in qs)<=tol)/400
        if tol==1: say("   %-24s q_k<=400: %s" % (sl.name, qs[:9]))
        say("        tol=%d -> %.1f%% of 1..400 qualify"%(tol,100*f))
say("   (small-denominator slopes make 'near a convergent' almost automatic; any such")
say("    claim must be scored against this rate.)")

say("== 3.2 shuffled-CF control: same partial-quotient multiset, shuffled order ==")
rows=[]
with open(ROOT+"/shape/data/explore_shape2.csv") as f:
    for d in csv.DictReader(f):
        d["K"]=int(d["K"]); d["r"]=int(d["r"]); d["h"]=int(d["h"]); d["M"]=int(d["M"])
        d["S2r"]=float(d["S2r"]); d["S3r"]=float(d["S3r"])
        rows.append(d)
BIG=[d for d in rows if d["K"]>=20]
for d in BIG:
    d["U_loc"]=d["S2r"]/math.sqrt(d["r"]); d["U_var"]=d["S3r"]
per=collections.defaultdict(list)
for d in BIG: per[d["slope"]].append(d)
pts=[(byname[n].float_value(), statistics.mean([d["U_var"] for d in v]),
      statistics.mean([d["U_loc"] for d in v])) for n,v in per.items()]
x=np.array([p[0] for p in pts])
cv=np.polyfit(x,np.array([p[1] for p in pts]),3)
cl=np.polyfit(x,np.array([p[2] for p in pts]),3)
say("   fitted cubic in theta on the %d EXPLORE slopes; residual sd: U_var %.4f, U_loc %.4f"
    %(len(pts), (np.array([p[1] for p in pts])-np.polyval(cv,x)).std(),
      (np.array([p[2] for p in pts])-np.polyval(cl,x)).std()))
rng=random.Random(4242)
WINDOWS=[(2,1),(3,1),(5,2),(8,1),(13,3)]
shuf=[]
for sl in EX:
    if sl.cf is None: continue
    cf=sl.cf[:]; rng.shuffle(cf)
    sh=Slope("shuf("+sl.name+")","shuffled",cf)
    try: sh.verify(460)
    except AssertionError: continue
    shuf.append(sh)
say("   built %d shuffled-CF slopes (same multiset, different order -> different theta)"%len(shuf))
dev=[]
for sh in shuf:
    vals_v,vals_l=[],[]
    for (M,E) in WINDOWS:
        for h in range(1,151):
            r=row2(sh,M,E,h)
            if r and r["K"]>=20:
                vals_v.append(r["var"]/r["r"]); vals_l.append((r["mean"]-r["r"])/math.sqrt(r["r"]))
    if not vals_v: continue
    t=sh.float_value()
    dv=statistics.mean(vals_v)-np.polyval(cv,t); dl=statistics.mean(vals_l)-np.polyval(cl,t)
    dev.append((t,dv,dl,sh.name))
say("   deviation of shuffled-CF slopes from the EXPLORE theta-curve:")
say("      U_var: mean %+.4f, sd %.4f, max|dev| %.4f"
    %(statistics.mean([d[1] for d in dev]),statistics.pstdev([d[1] for d in dev]),
      max(abs(d[1]) for d in dev)))
say("      U_loc: mean %+.4f, sd %.4f, max|dev| %.4f"
    %(statistics.mean([d[2] for d in dev]),statistics.pstdev([d[2] for d in dev]),
      max(abs(d[2]) for d in dev)))
say("   -> shuffling the CF (which changes Diophantine type but keeps the multiset) moves the")
say("      statistics only along the theta-curve; the curve itself is unchanged.")

say("== 3.3 resonance test: do residuals spike when h is near a convergent denominator? ==")
def feats(d):
    t=byname[d["slope"]].float_value(); ir=1.0/math.sqrt(d["r"])
    return [1.0,ir,1.0/d["r"],t,t*t,t*ir]
X=np.array([feats(d) for d in BIG])
for key in ("U_loc","U_var"):
    y=np.array([d[key] for d in BIG]); b,*_=np.linalg.lstsq(X,y,rcond=None); res=np.abs(y-X@b)
    near,far=[],[]
    for d,e in zip(BIG,res):
        qs=byname[d["slope"]].conv_denoms(400)
        (near if min(abs(d["h"]-q) for q in qs)<=2 else far).append(e)
    obs=np.mean(near)-np.mean(far)
    # null: random fake "denominators" of the same count and magnitude per slope
    rngl=random.Random(11); diffs=[]
    for _ in range(400):
        fake={}
        for n in {d["slope"] for d in BIG}:
            qs=byname[n].conv_denoms(400)
            fake[n]=[rngl.randint(1,400) for _ in qs]
        nn,ff=[],[]
        for d,e in zip(BIG,res):
            (nn if min(abs(d["h"]-q) for q in fake[d["slope"]])<=2 else ff).append(e)
        if nn and ff: diffs.append(np.mean(nn)-np.mean(ff))
    p=(sum(1 for v in diffs if abs(v)>=abs(obs))+1)/(len(diffs)+1)
    say("   %-6s |resid| near q_k (n=%4d) %.4f vs far (n=%4d) %.4f ; diff %+.4f ; "
        "random-denominator null p=%.3f"%(key,len(near),np.mean(near),len(far),np.mean(far),obs,p))
open(ROOT+"/shape/out/stage3.log","w").write("\n".join(LOG))
