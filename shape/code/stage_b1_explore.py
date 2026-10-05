"""Stage 1 (EXPLORE = the 47 slopes already studied): descriptive balance statistics."""
import os as _os, sys as _sys
ROOT = _os.path.dirname(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
_sys.path.insert(0, _os.path.join(ROOT, 'shape','code'))

import sys, csv, json, math, collections, statistics, random
sys.path.insert(0,ROOT+"/shape/code")
from slopes import build, Slope
from balance_core import iter_profiles, peak_and_margin, top_irregularity, best_ideal_fit, mean_of
LOG=[]
def say(*a):
    s=" ".join(str(x) for x in a); print(s); LOG.append(s)

HMAX=200
WINDOWS=[(2,1),(3,1),(5,2),(8,1),(13,3),(4,2),(7,1),(11,2)]
S=build(depth=60)
rng=random.Random(20261006); old={tuple(s.cf) for s in S if s.cf}
NEW=[]
while len(NEW)<10:
    cf=[rng.choice([1,1,1,2,2,3,4,5,7]) for _ in range(60)]
    if tuple(cf) in old: continue
    sl=Slope("new-random#%02d"%len(NEW),"random",cf)
    try: sl.verify(HMAX+20)
    except AssertionError: continue
    NEW.append(sl)
ALL=S+NEW
for s in ALL: s.verify(HMAX+20)
say("EXPLORE: %d slopes x %d windows x h<=%d"%(len(ALL),len(WINDOWS),HMAX))

S4EVERY=10          # pre-specified subsample for the (expensive) ideal-kernel fit
SPAN=20             # pre-specified k-search half-width
rows=[]; checks=fails=0
for sl in ALL:
    th=sl.float_value()
    for (M,E) in WINDOWS:
        for h,r,L,Pi,c,f in iter_profiles(sl,M,E,HMAX,check_every=50):
            checks,fails=c,f
            pm=peak_and_margin(Pi,h)
            D,V,wid=top_irregularity(L,pm["peak"],h)
            rec=dict(slope=sl.name,theta=th,M=M,E=E,h=h,r=r,
                     peak=pm["peak"],tied=pm["tied"],n_modes=pm["n_modes"],
                     nonadj=pm["nonadjacent_tie"],margin=pm["margin"],
                     log_lo=pm["log_lo"],log_hi=pm["log_hi"],tie_gap=pm["tie_gap"],
                     D=D,V=V,win=wid)
            if pm["tied"] or h%S4EVERY==0:
                mu=mean_of(Pi,h)
                k,tv=best_ideal_fit(Pi,h,mu,SPAN)
                rec["k_fit"]=k; rec["S4"]=tv; rec["balance_point"]=(k-2) if k else None
            rows.append(rec)
say("rows: %d | incremental spot-checks %d, failures %d"%(len(rows),checks,fails))
with open(ROOT+"/shape/data/balance_explore.csv","w",newline="") as f:
    fl=["slope","theta","M","E","h","r","peak","tied","n_modes","nonadj","margin",
        "log_lo","log_hi","tie_gap","D","V","win","k_fit","S4","balance_point"]
    w=csv.DictWriter(f,fieldnames=fl,extrasaction="ignore"); w.writeheader(); w.writerows(rows)

T=[d for d in rows if d["tied"]]
say("")
say("== B1.1 are all ties adjacent pairs? ==")
say("   tied rows: %d of %d (%.2f%%)"%(len(T),len(rows),100*len(T)/len(rows)))
say("   non-adjacent ties: %d"%sum(d["nonadj"] for d in T))
say("   number of tied indices: %s"%dict(sorted(collections.Counter(d["n_modes"] for d in T).items())))
say("   gap between the two tied indices: %s"%dict(sorted(collections.Counter(d["tie_gap"] for d in T).items())))

say("")
say("== B1.2 distribution of the balance margin m, by theta band and r band ==")
say("   (m = 0 exactly at a tie; small m > 0 = near-tie)")
def band(x,edges):
    for i in range(len(edges)-1):
        if edges[i]<=x<edges[i+1]: return i
    return None
TH=[0.13,0.25,0.40,0.55,0.70,0.90]
RB=[1,10,20,40,80,400]
say("   median m   |" + "".join("  r %3d-%3d"%(RB[i],RB[i+1]) for i in range(len(RB)-1)))
for t in range(len(TH)-1):
    line="   th %.2f-%.2f|"%(TH[t],TH[t+1])
    for rb in range(len(RB)-1):
        sub=[d["margin"] for d in rows if TH[t]<=d["theta"]<TH[t+1] and RB[rb]<=d["r"]<RB[rb+1]
             and math.isfinite(d["margin"])]
        line+="   %7.4f"%statistics.median(sub) if sub else "       -   "
    say(line)
say("   fraction of rows with m < 0.05 (near-ties), same bands:")
for t in range(len(TH)-1):
    line="   th %.2f-%.2f|"%(TH[t],TH[t+1])
    for rb in range(len(RB)-1):
        sub=[d for d in rows if TH[t]<=d["theta"]<TH[t+1] and RB[rb]<=d["r"]<RB[rb+1]]
        line+="   %6.3f "%(sum(1 for d in sub if d["margin"]<0.05)/len(sub)) if sub else "       -   "
    say(line)

say("")
say("== B1.3 is there a near-tie continuum, or are ties isolated? ==")
say("   distribution of m in the r band where ties live (r<20), by theta band:")
for t in range(len(TH)-1):
    sub=[d["margin"] for d in rows if TH[t]<=d["theta"]<TH[t+1] and d["r"]<20 and math.isfinite(d["margin"])]
    if not sub: continue
    sub.sort()
    q=lambda p: sub[min(len(sub)-1,int(p*len(sub)))]
    say("   theta %.2f-%.2f : n=%4d  min %.4f  q10 %.4f  q25 %.4f  median %.4f  (exact zeros: %d)"
        %(TH[t],TH[t+1],len(sub),sub[0],q(.10),q(.25),q(.50),sum(1 for x in sub if x==0)))

say("")
say("== B1.4 S4 (distance to the best-fit ideal kernel) vs theta at matched r ==")
s4=[d for d in rows if d.get("S4") is not None]
say("   median S4    |" + "".join("  r %3d-%3d"%(RB[i],RB[i+1]) for i in range(len(RB)-1)))
for t in range(len(TH)-1):
    line="   th %.2f-%.2f |"%(TH[t],TH[t+1])
    for rb in range(len(RB)-1):
        sub=[d["S4"] for d in s4 if TH[t]<=d["theta"]<TH[t+1] and RB[rb]<=d["r"]<RB[rb+1]]
        line+="   %7.4f"%statistics.median(sub) if sub else "       -   "
    say(line)

say("")
say("== B1.5 does the peak sit at the ideal balance point k_fit - 2 ? ==")
for lab,sub in (("tied rows",[d for d in s4 if d["tied"]]),
                ("untied rows",[d for d in s4 if not d["tied"]])):
    if not sub: continue
    off=collections.Counter(d["peak"]-d["balance_point"] for d in sub)
    within1=sum(v for k,v in off.items() if abs(k)<=1)/len(sub)
    say("   %-12s n=%5d  peak - (k_fit-2): %s  | within 1: %.3f"
        %(lab,len(sub),dict(sorted(off.items())[:9]),within1))

say("")
say("== B1.6 margin m vs top-irregularity D (the mechanism's core claim) ==")
import numpy as np
fin=[d for d in rows if math.isfinite(d["margin"])]
def spearman(x,y):
    rx=np.argsort(np.argsort(x)); ry=np.argsort(np.argsort(y))
    return float(np.corrcoef(rx,ry)[0,1])
say("   pooled over everything: spearman(m, D) = %+.4f (n=%d)"
    %(spearman([d["D"] for d in fin],[d["margin"] for d in fin]),len(fin)))
say("   NOTE: D is ~ 11*theta by construction, so the pooled figure largely restates theta.")
say("   within-slope (theta fixed) spearman(m, D):")
per=collections.defaultdict(list)
for d in fin: per[d["slope"]].append(d)
cs=[]
for n,v in per.items():
    if len(v)<50: continue
    c=spearman([d["D"] for d in v],[d["margin"] for d in v]); cs.append(c)
say("   %d slopes: median %+.4f, mean %+.4f, sd %.4f, fraction > 0: %.3f"
    %(len(cs),statistics.median(cs),statistics.mean(cs),statistics.pstdev(cs),
      sum(1 for c in cs if c>0)/len(cs)))
say("   within (slope, r-band) spearman(m, D) -- removes r as well:")
cs2=[]
for n,v in per.items():
    for rb in range(len(RB)-1):
        s_=[d for d in v if RB[rb]<=d["r"]<RB[rb+1]]
        if len(s_)<30: continue
        if len({d["D"] for d in s_})<2: continue
        cs2.append(spearman([d["D"] for d in s_],[d["margin"] for d in s_]))
say("   %d (slope,r-band) cells: median %+.4f, fraction > 0: %.3f"
    %(len(cs2),statistics.median(cs2),sum(1 for c in cs2 if c>0)/len(cs2)))
json.dump(dict(rows=len(rows),ties=len(T),nonadj=sum(d["nonadj"] for d in T)),
          open(ROOT+"/shape/out/balance_b1.json","w"),indent=1)
open(ROOT+"/shape/out/balance_b1.log","w").write("\n".join(LOG))
