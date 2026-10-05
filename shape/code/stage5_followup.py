"""Stage 5: the two things the holdout flagged -- (a) the tie threshold, (b) the
significant U_loc resonance at p=0.007: is it a Diophantine effect or an h-confound?"""
import os as _os, sys as _sys
ROOT = _os.path.dirname(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
_sys.path.insert(0, _os.path.join(ROOT, 'shape','code'))

import sys, csv, math, collections, statistics, random
import numpy as np
sys.path.insert(0,ROOT+"/shape/code")
from slopes import build, split
LOG=[]
def say(*a):
    s=" ".join(str(x) for x in a); print(s); LOG.append(s)
S=build(depth=60)
for s in S: s.verify(460)
EX,HO=split(S); byname={s.name:s for s in S}
H=[]
with open(ROOT+"/shape/data/holdout_shape.csv") as f:
    for d in csv.DictReader(f):
        for k in ("K","r","h","M"): d[k]=int(d[k])
        for k in ("S2r","S3r","theta","frac_h","frac_M"): d[k]=float(d[k])
        d["ext"]=d["ext"]=="True"; d["mode_tied"]=d["mode_tied"]=="True"
        H.append(d)
E=[]
with open(ROOT+"/shape/data/explore_shape2.csv") as f:
    for d in csv.DictReader(f):
        for k in ("K","r","h","M"): d[k]=int(d[k])
        for k in ("S2r","S3r","frac_h","frac_M"): d[k]=float(d[k])
        d["mode_tied"]=d["mode_tied"]=="True"
        d["theta"]=byname[d["slope"]].float_value(); E.append(d)
for d in H+E:
    if d["r"]>0: d["U_var"]=d["S3r"]; d["U_loc"]=d["S2r"]/math.sqrt(d["r"])

say("== 5.1 the tie threshold: P1 used K>=17, fitted on EXPLORE; what holds on both halves? ==")
for lo in (15,16,17,18,19,20,21,25):
    e=[d for d in E if d["K"]>=lo]; h_=[d for d in H if d["K"]>=lo]
    say("   K>=%2d : EXPLORE %4d ties / %5d ; HOLDOUT %4d ties / %5d"
        %(lo,sum(d["mode_tied"] for d in e),len(e),sum(d["mode_tied"] for d in h_),len(h_)))
t=[d for d in H if d["K"]>=17 and d["mode_tied"]]
say("   the 28 HOLDOUT ties with K in 17..19:")
say("     slopes: %s"%dict(collections.Counter(d["slope"] for d in t)))
say("     K values: %s ; windows: %s"
    %(dict(sorted(collections.Counter(d["K"] for d in t).items())),
      dict(collections.Counter((d["M"],d["E"]) for d in t))))
say("   -> K>=17 was an EXPLORE artifact; K>=20 is tie-free on BOTH halves (0/8663 and 0/11189).")

say("== 5.2 P2 residuals: is the positive bias caused by the extended horizons? ==")
import json
FR=json.load(open(ROOT+"/shape/out/frozen_curve.json"))
cv=np.array(FR["cubic_U_var"]); LO,HI=FR["theta_range"]
for lab,sel in (("h<=150 only (matches the EXPLORE range)",lambda d: not d["ext"]),
                ("h=160..400 only (extended)",lambda d: d["ext"]),
                ("all", lambda d: True)):
    per=collections.defaultdict(list)
    for d in H:
        if d["K"]>=20 and sel(d) and LO<=d["theta"]<=HI: per[d["slope"]].append(d["U_var"])
    res=[statistics.mean(v)-float(np.polyval(cv,byname[n].float_value())) for n,v in per.items()]
    if res: say("   %-38s n_slopes=%2d  mean resid %+.4f  sd %.4f  max|.| %.4f"
                %(lab,len(res),statistics.mean(res),statistics.pstdev(res),max(abs(x) for x in res)))
say("   -> the bias is a finite-size (r) effect: the curve was fitted at h<=150 and the")
say("      statistics still drift with r, so extended horizons sit slightly off it.")

say("== 5.3 the U_loc resonance (p=0.007): confounded with h? ==")
def feats(d):
    t=d["theta"]; ir=1.0/math.sqrt(d["r"])
    return [1.0,ir,1.0/d["r"],t,t*t,t*ir,
            math.cos(2*math.pi*d["frac_h"]),math.sin(2*math.pi*d["frac_h"]),
            math.cos(4*math.pi*d["frac_h"]),math.sin(4*math.pi*d["frac_h"]),
            math.cos(2*math.pi*d["frac_M"]),math.sin(2*math.pi*d["frac_M"]),
            float(d["r"]%2),float(d["r"]%3==0)]
EB=[d for d in E if d["K"]>=20]; HB=[d for d in H if d["K"]>=20]
X=np.array([feats(d) for d in EB]); y=np.array([d["U_loc"] for d in EB])
b,*_=np.linalg.lstsq(X,y,rcond=None)
for d in HB: d["res"]=abs(d["U_loc"]-float(np.array(feats(d))@b))
def nearq(d,tol=2):
    return min(abs(d["h"]-q) for q in byname[d["slope"]].conv_denoms(400))<=tol
near=[d for d in HB if nearq(d)]; far=[d for d in HB if not nearq(d)]
say("   raw: near %.4f (n=%d) vs far %.4f (n=%d)"
    %(np.mean([d["res"] for d in near]),len(near),np.mean([d["res"] for d in far]),len(far)))
say("   h distribution: near mean h=%.1f, far mean h=%.1f"
    %(np.mean([d["h"] for d in near]),np.mean([d["h"] for d in far])))
say("   residual vs h (all HOLDOUT rows): corr = %+.3f"
    %np.corrcoef([d["h"] for d in HB],[d["res"] for d in HB])[0,1])
say("   within h-bins (removes the h confound):")
tot_n=tot_f=0; diffs=[]
for lo,hi in ((20,60),(60,100),(100,150),(150,250),(250,401)):
    n=[d["res"] for d in near if lo<=d["h"]<hi]; f=[d["res"] for d in far if lo<=d["h"]<hi]
    if len(n)<5: continue
    diffs.append(np.mean(n)-np.mean(f)); tot_n+=len(n); tot_f+=len(f)
    say("      h %3d-%3d : near %.4f (n=%3d) vs far %.4f (n=%4d)  diff %+.4f"
        %(lo,hi,np.mean(n),len(n),np.mean(f),len(f),np.mean(n)-np.mean(f)))
say("   mean within-bin difference: %+.4f (raw difference was %+.4f)"
    %(np.mean(diffs),np.mean([d["res"] for d in near])-np.mean([d["res"] for d in far])))
# stratified permutation null: shuffle near/far labels WITHIN h-bins
rng=random.Random(17); obs=np.mean(diffs); null=[]
bins=[(20,60),(60,100),(100,150),(150,250),(250,401)]
for _ in range(2000):
    ds=[]
    for lo,hi in bins:
        pool=[d["res"] for d in HB if lo<=d["h"]<hi]
        k=len([d for d in near if lo<=d["h"]<hi])
        if k<5 or len(pool)-k<5: continue
        rng.shuffle(pool); ds.append(np.mean(pool[:k])-np.mean(pool[k:]))
    if ds: null.append(np.mean(ds))
p=(sum(1 for v in null if abs(v)>=abs(obs))+1)/(len(null)+1)
say("   h-stratified permutation p = %.4f  (unstratified was 0.007)"%p)
say("== 5.4 is the effect carried by one slope? ==")
rowsbys=collections.defaultdict(lambda:[[],[]])
for d in HB: rowsbys[d["slope"]][0 if nearq(d) else 1].append(d["res"])
for n in sorted(rowsbys,key=lambda n: -(np.mean(rowsbys[n][0]) - np.mean(rowsbys[n][1]) if rowsbys[n][0] else -9))[:6]:
    a,b_=rowsbys[n]
    if a: say("   %-24s near %.4f (n=%2d) far %.4f (n=%4d) diff %+.4f"
              %(n,np.mean(a),len(a),np.mean(b_),len(b_),np.mean(a)-np.mean(b_)))
say("== 5.5 multiple-comparison accounting ==")
say("   resonance tests run: 2 statistics x 2 halves = 4.  Results: EXPLORE U_loc p=0.43,")
say("   EXPLORE U_var p=0.78, HOLDOUT U_var p=0.137, HOLDOUT U_loc p=0.007 (unstratified).")
say("   The EXPLORE U_loc test had the OPPOSITE sign (near < far).")
open(ROOT+"/shape/out/stage5.log","w").write("\n".join(LOG))
