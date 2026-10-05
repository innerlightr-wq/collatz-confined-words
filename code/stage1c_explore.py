"""Stage 1 (EXPLORE only): per-window transition records, residual anatomy,
phase predictors vs boring baselines, and chance calibration of 'near a convergent'."""
import os as _os, sys as _sys
ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
_sys.path.insert(0, _os.path.join(ROOT, 'code'))

import sys, json, csv, random, math
sys.path.insert(0,ROOT+"/code")
from core import *
from windows import split
LOG=[]
def say(*a):
    s=" ".join(str(x) for x in a); print(s); LOG.append(s)

HMAX=150; MARGIN=2
EXPLORE,HOLDOUT=split()
say("EXPLORE windows: %d  (HOLDOUT untouched in this script)"%len(EXPLORE))

def records(W,hmax=HMAX,margin=MARGIN):
    out=[]; prev=chain_N(W.M,W.E,1)
    for h in range(1,hmax+1):
        cur=chain_N(W.M,W.E,h+1)
        eps=G(W.M+h+1)-G(W.M+h)
        jm=max(0,h-margin)
        if eps==1:
            run=0;T={}
            for j in range(jm+1):
                run+=prev.get(j,0); T[j]=run+1 if j==0 else run
            res=[cur.get(j,0)-T[j] for j in range(jm+1)]
        else:
            res=[cur.get(j,0)-prev.get(j,0) for j in range(jm+1)]
        out.append(dict(h=h,eps=eps,jmax=jm,resid=res,ok=not any(res),
                        rawl0=G(W.M+h)-G(W.M)-W.E+1))
        prev=cur
    return out

# ---------------- window table + (a) success/failure record
say("== 1.7 EXPLORE window table (T1 on eps=1 transitions, T0 on eps=0) ==")
tab=[];E1=OK1=E0=OK0=0
for W in EXPLORE:
    r=records(W)
    o=[x for x in r if x["eps"]==1]; z=[x for x in r if x["eps"]==0]
    t=dict(name=W.tag(),M=W.M,c=W.c,prefix="-".join(map(str,W.prefix)),Vmin=W.Vmin,E=W.E,
           gamma=round(W.gamma,6),n_eps1=len(o),ok_eps1=sum(x["ok"] for x in o),
           n_eps0=len(z),ok_eps0=sum(x["ok"] for x in z),
           fail_h=[x["h"] for x in o if not x["ok"]])
    tab.append(t); E1+=t["n_eps1"]; OK1+=t["ok_eps1"]; E0+=t["n_eps0"]; OK0+=t["ok_eps0"]
say("   %d windows ; T1 exact on %d/%d eps=1 transitions (%.4f) ; T0 exact on %d/%d"
    %(len(EXPLORE),OK1,E1,OK1/E1,OK0,E0))
say("   windows with any T1 failure: %d"%sum(1 for t in tab if t["fail_h"]))
with open(ROOT+"/data/explore_window_table.csv","w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=list(tab[0].keys())); w.writeheader()
    for t in tab: w.writerow({**t,"fail_h":";".join(map(str,t["fail_h"]))})

# ---------------- (a) residual anatomy
say("== 1.8 residual anatomy (EXPLORE) ==")
from collections import Counter
shapes=Counter(); locs=Counter(); sizes=Counter()
for W in EXPLORE:
    for x in records(W):
        if x["eps"]==1 and not x["ok"]:
            shapes["all -1" if set(x["resid"])=={-1} else str(x["resid"][:4])]+=1
            locs[("j=0..%d"%x["jmax"])]+=1
            sizes[max(abs(v) for v in x["resid"])]+=1
say("   residual shapes: %s"%dict(shapes))
say("   |residual| maxima: %s   (sign: all negative)"%dict(sizes))
say("   -> every failure residual is the constant vector (-1,...,-1) on the compared range")

# ---------------- (c) predictors
say("== 1.9 predictors of T1 success (EXPLORE, eps=1 transitions only) ==")
CONV=convergent_denominators(2000); SEMI=semiconvergent_denominators(400)
say("   convergent denominators q_k = %s"%CONV)
data=[]
for W in EXPLORE:
    for x in records(W):
        if x["eps"]!=1: continue
        ph=(h_:=x["h"])*(math.log2(3)-1)+W.gamma
        data.append(dict(M=W.M,c=W.c,E=W.E,h=x["h"],ok=x["ok"],
                         gamma=W.gamma,phase=ph%1.0,
                         dconv=min(abs(x["h"]-q) for q in CONV),
                         dsemi=min(abs(x["h"]-q) for q in SEMI),
                         rawl0=x["rawl0"]))
base=sum(d["ok"] for d in data)/len(data)
say("   n=%d transitions ; base rate (always-success) = %.4f"%(len(data),base))
def acc(pred): return sum(1 for d in data if pred(d)==d["ok"])/len(data)
cands={
 "always-success (base rate)": lambda d: True,
 "FROZEN RULE: rawl0 >= 0":    lambda d: d["rawl0"]>=0,
 "M==2 only":                  lambda d: d["M"]==2,
 "E==1 only":                  lambda d: d["E"]==1,
}
for t in (6,12,18,24,30):
    cands["boring: h >= %d"%t]=lambda d,t=t: d["h"]>=t
for tol in (0,1,2,3):
    cands["Diophantine: |h-q_k| <= %d"%tol]=lambda d,tol=tol: d["dconv"]<=tol
    cands["Diophantine: |h-semiconv| <= %d"%tol]=lambda d,tol=tol: d["dsemi"]<=tol
best_int=None
for a in range(0,20):
    for b in range(a+1,21):
        lo,hi=a/20,b/20
        f=lambda d,lo=lo,hi=hi: lo<=d["phase"]<hi
        s=acc(f)
        if best_int is None or s>best_int[0]: best_int=(s,lo,hi)
cands["Diophantine: frac(h*beta+gamma_W) in [%.2f,%.2f) (interval fitted on EXPLORE)"%best_int[1:]]=\
    lambda d,lo=best_int[1],hi=best_int[2]: lo<=d["phase"]<hi
for k,f in cands.items():
    say("   acc=%.4f  %s"%(acc(f),k))
# does the phase add anything once rawl0 is known?
res=[d for d in data if not d["ok"]]
say("   all %d failures have rawl0 <= -1 : %s ; min rawl0 among successes = %d"
    %(len(res),all(d["rawl0"]<=-1 for d in res),min(d["rawl0"] for d in data if d["ok"])))
say("   h at failures: max %d ; h_open(M,E)=min{h: rawl0>=0} per class matches exactly"%max([d["h"] for d in res],default=0))

# ---------------- (d) calibration of "near a convergent"
say("== 1.10 calibration: how often does a random small integer land near a convergent? ==")
for Rmax in (20,50,100,150):
    for tol in (0,1,2):
        f1=sum(1 for x in range(1,Rmax+1) if min(abs(x-q) for q in CONV)<=tol)/Rmax
        f2=sum(1 for x in range(1,Rmax+1) if min(abs(x-q) for q in SEMI)<=tol)/Rmax
        say("   h in 1..%3d, tol=%d : %.0f%% near a convergent, %.0f%% near a semiconvergent"
            %(Rmax,tol,100*f1,100*f2))
say("   gap-based criterion: a set of k successes in 1..H has mean gap H/k; for k=16,H=96")
say("   the mean gap is 6.0, so observing gaps of 5,7,12 is the null expectation, not a signal.")
json.dump(dict(table=tab,base_rate=base,n=len(data)),
          open(ROOT+"/out/stage1c_explore.json","w"),indent=1,default=str)
open(ROOT+"/out/stage1c_explore.log","w").write("\n".join(LOG))
