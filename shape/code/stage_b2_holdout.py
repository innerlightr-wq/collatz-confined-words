"""Stage 2: the frozen predictions P1-P5 on 8 NEW slopes (seed 20261007), theta spread over
[0.15,0.9], windows (2,1),(4,2),(7,1), horizons h <= 400."""
import os as _os, sys as _sys
ROOT = _os.path.dirname(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
_sys.path.insert(0, _os.path.join(ROOT, 'shape','code'))

import sys, csv, json, math, random, collections, statistics
import numpy as np
sys.path.insert(0,ROOT+"/shape/code")
from slopes import build, Slope
from balance_core import iter_profiles, peak_and_margin, top_irregularity, best_ideal_fit, mean_of
LOG=[]
def say(*a):
    s=" ".join(str(x) for x in a); print(s); LOG.append(s)
HMAX=400; WINDOWS=[(2,1),(4,2),(7,1)]
old={tuple(s.cf) for s in build(depth=60) if s.cf}
rng=random.Random(20261006)                       # reproduce the 10 retest CFs to exclude them
for _ in range(400):
    cf=[rng.choice([1,1,1,2,2,3,4,5,7]) for _ in range(60)]; old.add(tuple(cf))
rng=random.Random(20261007)
BANDS=[(0.15,0.24),(0.24,0.33),(0.33,0.42),(0.42,0.51),(0.51,0.60),(0.60,0.70),(0.70,0.80),(0.80,0.90)]
HO=[]
for lo,hi in BANDS:
    for _ in range(20000):
        cf=[rng.choice([1,1,1,2,2,3,4,5,7]) for _ in range(60)]
        if tuple(cf) in old: continue
        s=Slope("hold-random#%d"%len(HO),"random",cf)
        if not (lo<=s.float_value()<hi): continue
        try: s.verify(HMAX+20)
        except AssertionError: continue
        HO.append(s); old.add(tuple(cf)); break
say("HOLDOUT slopes (%d), theta spread over [0.15,0.9]:"%len(HO))
for s in HO: say("   %-16s theta=%.9f"%(s.name,s.float_value()))
say("windows %s ; h <= %d"%(WINDOWS,HMAX))

S4EVERY=10; SPAN=20
rows=[]
for sl in HO:
    th=sl.float_value()
    for (M,E) in WINDOWS:
        for h,r,L,Pi,c,f in iter_profiles(sl,M,E,HMAX,check_every=100):
            pm=peak_and_margin(Pi,h); js=pm["peak"]
            D,V,wid=top_irregularity(L,js,h)
            delta = Pi[js+1]-2*Pi[js] if js+1<len(Pi) else None
            rec=dict(slope=sl.name,theta=th,M=M,E=E,h=h,r=r,peak=js,tied=pm["tied"],
                     n_modes=pm["n_modes"],nonadj=pm["nonadjacent_tie"],margin=pm["margin"],
                     tie_gap=pm["tie_gap"],D=D,V=V,delta=(abs(delta) if delta is not None else None),
                     pi_peak=Pi[js],digits=len(str(Pi[js])))
            if pm["tied"] or h%S4EVERY==0:
                k,tv=best_ideal_fit(Pi,h,mean_of(Pi,h),SPAN)
                rec["k_fit"]=k; rec["S4"]=tv; rec["bp"]=(k-2) if k else None
            rows.append(rec)
say("HOLDOUT rows: %d"%len(rows))
with open(ROOT+"/shape/data/balance_holdout.csv","w",newline="") as f:
    fl=["slope","theta","M","E","h","r","peak","tied","n_modes","nonadj","margin","tie_gap",
        "D","V","delta","digits","k_fit","S4","bp"]
    w=csv.DictWriter(f,fieldnames=fl,extrasaction="ignore"); w.writeheader(); w.writerows(rows)
T=[d for d in rows if d["tied"]]
def spearman(x,y):
    rx=np.argsort(np.argsort(x)); ry=np.argsort(np.argsort(y)); return float(np.corrcoef(rx,ry)[0,1])

say(""); say("== P1: all ties adjacent ==")
say("   ties %d / %d ; non-adjacent: %d ; run lengths: %s ; gaps: %s"
    %(len(T),len(rows),sum(d["nonadj"] for d in T),
      dict(sorted(collections.Counter(d["n_modes"] for d in T).items())),
      dict(sorted(collections.Counter(d["tie_gap"] for d in T).items()))))
say("   P1 verdict: %s"%("HOLDS" if sum(d["nonadj"] for d in T)==0 else "FALSIFIED"))

say(""); say("== P2: margin m vs top-irregularity D, within (slope, r-band) ==")
RB=[1,10,20,40,80,400]
fin=[d for d in rows if math.isfinite(d["margin"])]
cells=[]
per=collections.defaultdict(list)
for d in fin: per[d["slope"]].append(d)
for n,v in per.items():
    for i in range(len(RB)-1):
        s_=[d for d in v if RB[i]<=d["r"]<RB[i+1]]
        if len(s_)<30 or len({d["D"] for d in s_})<2: continue
        cells.append(spearman([d["D"] for d in s_],[d["margin"] for d in s_]))
med=statistics.median(cells); frac=sum(1 for c in cells if c>0)/len(cells)
say("   %d cells: median Spearman %+.4f, fraction positive %.3f"%(len(cells),med,frac))
say("   frozen prediction: FAILS, median in [-0.10,+0.10] and fraction in [0.35,0.65]")
say("   P2 verdict: prediction %s (mechanism's own claim of a positive relation: %s)"
    %("CONFIRMED" if (-0.10<=med<=0.10 and 0.35<=frac<=0.65) else "NOT CONFIRMED",
      "SUPPORTED" if med>0.10 else "REJECTED"))

say(""); say("== P3: at matched r (band 10-20), m and S4 vs theta ==")
TH=[0.15,0.30,0.45,0.60,0.75,0.90]
band=[d for d in rows if 10<=d["r"]<20]
ms,s4s,ths=[],[],[]
for i in range(len(TH)-1):
    sub=[d for d in band if TH[i]<=d["theta"]<TH[i+1]]
    s4sub=[d["S4"] for d in sub if d.get("S4") is not None]
    if not sub: continue
    mm=statistics.median([d["margin"] for d in sub if math.isfinite(d["margin"])])
    ss=statistics.median(s4sub) if s4sub else float("nan")
    ms.append(mm); s4s.append(ss); ths.append((TH[i],TH[i+1]))
    say("   theta %.2f-%.2f : n=%4d  median m %.6f   median S4 %.4f"%(TH[i],TH[i+1],len(sub),mm,ss))
rm=ms[-1]/ms[0] if ms[0]>0 else float("inf"); rs=s4s[-1]/s4s[0] if s4s[0]>0 else float("inf")
sp=spearman([d["theta"] for d in band if math.isfinite(d["margin"])],
            [d["margin"] for d in band if math.isfinite(d["margin"])])
rngp=random.Random(1); x=[d["theta"] for d in band if math.isfinite(d["margin"])]
y=[d["margin"] for d in band if math.isfinite(d["margin"])]; ge=0
for _ in range(2000):
    p=list(y); rngp.shuffle(p)
    if abs(spearman(x,p))>=abs(sp): ge+=1
say("   m ratio (highest/lowest theta band) = %.2f (predicted >= 3)"%rm)
say("   S4 ratio                            = %.2f (predicted >= 1.5)"%rs)
say("   Spearman(theta, m) = %+.4f, permutation p = %.4f (predicted > 0, p < 0.05)"%(sp,(ge+1)/2001))
say("   P3 verdict: %s"%("HOLDS" if rm>=3 and rs>=1.5 and sp>0 and (ge+1)/2001<0.05 else "FAILS"))

say(""); say("== P4 / P4': the ideal balance point ==")
s4rows=[d for d in rows if d.get("bp") is not None]
for lab,sub in (("tied",[d for d in s4rows if d["tied"]]),("untied",[d for d in s4rows if not d["tied"]])):
    if not sub: continue
    w1=sum(1 for d in sub if abs(d["peak"]-d["bp"])<=1)/len(sub)
    ex=sum(1 for d in sub if d["peak"]==d["bp"])/len(sub)
    say("   %-6s n=%5d : within 1 of k_fit-2 = %.4f ; exactly k_fit-2 = %.4f"%(lab,len(sub),w1,ex))
tt=[d for d in s4rows if d["tied"]]; uu=[d for d in s4rows if not d["tied"]]
if tt and uu:
    a=sum(1 for d in tt if d["peak"]==d["bp"])/len(tt); b=sum(1 for d in uu if d["peak"]==d["bp"])/len(uu)
    say("   P4 (as stated) is vacuous if both 'within 1' rates are >= 0.99: %s"
        %("CONFIRMED vacuous" if min(sum(1 for d in tt if abs(d["peak"]-d["bp"])<=1)/len(tt),
                                     sum(1 for d in uu if abs(d["peak"]-d["bp"])<=1)/len(uu))>=0.99 else "not vacuous"))
    say("   P4' : tied %.3f - untied %.3f = %+.3f (predicted >= +0.15) -> %s"
        %(a,b,a-b,"HOLDS" if a-b>=0.15 else "FAILS"))

say(""); say("== P5: ties as a small-integer coincidence ==")
mx=max((d["pi_peak"] for d in T),default=0)
say("   largest Pi(peak) carrying a tie: %s (%d digits) ; predicted < 10^10 -> %s"
    %("{:,}".format(mx),len(str(mx)) if mx else 0,"HOLDS" if mx<10**10 else "FALSIFIED"))
big=[d for d in rows if d["r"]>=20]
say("   rows with r >= 20: %d ; ties there: %d -> %s"
    %(len(big),sum(d["tied"] for d in big),"HOLDS" if not sum(d["tied"] for d in big) else "FALSIFIED"))
RB2=[1,5,10,15,20,30,50,100,400]
say("   r band    n     tie rate   median m     median |Delta|            10%-most-balanced tie rate")
for i in range(len(RB2)-1):
    sub=[d for d in rows if RB2[i]<=d["r"]<RB2[i+1] and math.isfinite(d["margin"]) and d["delta"] is not None]
    if len(sub)<20: continue
    sub2=sorted(sub,key=lambda d:d["margin"]); top=sub2[:max(1,len(sub2)//10)]
    say("   %3d-%3d %6d   %.4f    %.6f    %-24s %.4f"
        %(RB2[i],RB2[i+1],len(sub),sum(d["tied"] for d in sub)/len(sub),
          statistics.median([d["margin"] for d in sub]),
          "{:,}".format(int(statistics.median([d["delta"] for d in sub]))),
          sum(d["tied"] for d in top)/len(top)))
json.dump(dict(rows=len(rows),ties=len(T),nonadj=sum(d["nonadj"] for d in T),
               max_pi_peak_tied=str(mx)),
          open(ROOT+"/shape/out/balance_b2.json","w"),indent=1)
open(ROOT+"/shape/out/balance_b2.log","w").write("\n".join(LOG))
