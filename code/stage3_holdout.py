"""Stage 3: frozen rules P1/P2 on HOLDOUT windows, on EXPLORE windows at extended
horizons, and on a wide holdout of unseen phases.  Compared against nulls."""
import os as _os, sys as _sys
ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
_sys.path.insert(0, _os.path.join(ROOT, 'code'))

import sys, json, csv, math, random
sys.path.insert(0,ROOT+"/code")
from core import *
from windows import split
LOG=[]
def say(*a):
    s=" ".join(str(x) for x in a); print(s); LOG.append(s)

CONV=convergent_denominators(5000)
rawl0=lambda M,E,h: G(M+h)-G(M)-E+1

def transitions(M,E,hlo,hhi,margin):
    """yield dicts for every transition h in [hlo,hhi]."""
    out=[]
    prev=chain_N(M,E,hlo)
    for h in range(hlo,hhi+1):
        cur=chain_N(M,E,h+1)
        eps=G(M+h+1)-G(M+h); jm=max(0,h-margin)
        if eps==1:
            run=0;T={}
            for j in range(jm+1):
                run+=prev.get(j,0); T[j]=run+1 if j==0 else run
            res=[cur.get(j,0)-T[j] for j in range(jm+1)]
        else:
            res=[cur.get(j,0)-prev.get(j,0) for j in range(jm+1)]
        out.append(dict(M=M,E=E,h=h,eps=eps,jmax=jm,ok=not any(res),
                        resid_const_minus1=(set(res)=={-1}) if any(res) else None,
                        r=rawl0(M,E,h)))
        prev=cur
    return out

def evaluate(name,recs):
    """P1 as a decision rule on the eps=1 transitions of recs, + nulls."""
    d=[x for x in recs if x["eps"]==1]
    z=[x for x in recs if x["eps"]==0]
    n=len(d); pos=sum(x["ok"] for x in d); neg=n-pos
    P1=lambda x: x["r"]>=0
    err=[x for x in d if P1(x)!=x["ok"]]
    # degenerate-branch shape check
    shape_ok=all(x["resid_const_minus1"] for x in d if not x["ok"])
    T0=all(x["ok"] for x in z)
    say("  [%s]  eps=1 transitions n=%d (success %d, failure %d) ; eps=0 transitions %d"
        %(name,n,pos,neg,len(z)))
    say("    P1 (T1 law + r>=0 branch) : errors = %d   accuracy = %.6f"%(len(err),1-len(err)/n))
    say("    P1 degenerate residual == (-1,...,-1) at every failure : %s"%shape_ok)
    say("    P2/T0 identity exact on every eps=0 transition (margin as tested) : %s"%T0)
    def acc(f): return sum(1 for x in d if f(x)==x["ok"])/n
    def bal(f):
        tp=sum(1 for x in d if x["ok"] and f(x)); fn=pos-tp
        tn=sum(1 for x in d if (not x["ok"]) and (not f(x))); fp=neg-tn
        s1=tp/pos if pos else 0; s2=tn/neg if neg else 0
        return 0.5*(s1+s2)
    rows=[("P1: r>=0",acc(P1),bal(P1)),
          ("null: always-success (base rate)",acc(lambda x:True),bal(lambda x:True))]
    for tol in (0,1,2,3):
        f=lambda x,tol=tol: min(abs(x["h"]-q) for q in CONV)<=tol
        rows.append(("null-Diophantine: |h-q_k|<=%d"%tol,acc(f),bal(f)))
    rng=random.Random(777)
    # "random nearby integers" version of the convergent rule
    for tol in (1,2):
        accs=[];bals=[]
        for rep in range(200):
            fake=[rng.randint(max(1,q-3),q+3) for q in CONV if q<=5000]
            f=lambda x,fake=fake,tol=tol: min(abs(x["h"]-q) for q in fake)<=tol
            accs.append(acc(f)); bals.append(bal(f))
        rows.append(("null: random-nearby-integers |h-q~|<=%d (mean of 200)"%tol,
                     sum(accs)/len(accs),sum(bals)/len(bals)))
    # best EXPLORE-fitted phase interval was [0,1) == always-success; also test a grid here
    best=max(((acc(lambda x,lo=a/20,hi=b/20: lo<=((x["h"]*(math.log2(3)-1)+ (x["M"]*math.log2(3))%1)%1)<hi),a/20,b/20)
              for a in range(20) for b in range(a+1,21)))
    rows.append(("null-Diophantine: BEST phase interval re-fitted ON THIS SET [%.2f,%.2f)"%best[1:],
                 best[0],float('nan')))
    for nm,a,b in rows:
        say("    acc=%.6f  bal.acc=%s  %s"%(a,("%.4f"%b if b==b else " n/a "),nm))
    # permutation p-value for P1's balanced accuracy
    labels=[x["ok"] for x in d]; preds=[P1(x) for x in d]
    obs=sum(1 for p,l in zip(preds,labels) if p==l)
    rng=random.Random(12345); ge=0; R=2000
    for _ in range(R):
        sh=labels[:]; rng.shuffle(sh)
        if sum(1 for p,l in zip(preds,sh) if p==l)>=obs: ge+=1
    say("    permutation p-value for P1 (label shuffle, R=%d): p = %.5f"%(R,(ge+1)/(R+1)))
    return dict(name=name,n=n,pos=pos,neg=neg,P1_errors=len(err),
                shape_ok=bool(shape_ok),T0=bool(T0),rows=rows,
                perm_p=(ge+1)/(R+1))

out=[]
EX,HO=split()
say("=== 3.1 HOLDOUT windows (66), h<=150, margin 2 ===")
recs=[]
seen=set()
for W in HO:
    if (W.M,W.E) in seen: continue      # identical (M,E) => identical profiles (P1 aux claim)
    seen.add((W.M,W.E))
    recs+=transitions(W.M,W.E,1,150,2)
say("  distinct (M,E) classes in HOLDOUT: %d"%len(seen))
out.append(evaluate("HOLDOUT h<=150 margin2",recs))

say("=== 3.2 HOLDOUT windows, margin 0  (j <= h : the stronger P1 claim) ===")
recs0=[]
for M,E in sorted(seen): recs0+=transitions(M,E,1,150,0)
out.append(evaluate("HOLDOUT h<=150 margin0",recs0))

say("=== 3.3 EXPLORE classes at EXTENDED horizons h=151..600, margin 0 ===")
rece=[]
for M,E in sorted({(W.M,W.E) for W in EX}): rece+=transitions(M,E,151,600,0)
out.append(evaluate("EXPLORE extended h=151..600 margin0",rece))

say("=== 3.4 WIDE holdout: unseen phases M=41..120, E in {1,2,3,4,6,8,12}, h<=120, margin 0 ===")
recw=[]
for M in range(41,121):
    for E in (1,2,3,4,6,8,12):
        recw+=transitions(M,E,1,120,0)
out.append(evaluate("WIDE M=41..120 h<=120 margin0",recw))

say("=== 3.5 DEEP unseen phases M=101..110, E in {1,2,4}, h<=400, margin 0 ===")
recd=[]
for M in range(101,111):
    for E in (1,2,4):
        recd+=transitions(M,E,1,400,0)
out.append(evaluate("DEEP M=101..110 h<=400 margin0",recd))

say("=== 3.6 auxiliary frozen claims on all held-out classes ===")
nhh=all(chain_N(M,E,h).get(h,0)==0 for M in (41,67,103) for E in (1,2,4) for h in range(1,120))
say("  N_h(h)=0 on held-out classes: %s"%nhh)
# (M,E) collapse on holdout windows, verified with the independent DP
coll=True
byclass={}
for W in HO:
    p=tuple(sorted({k:v for k,v in dp_N(W,22).items() if v}.items()))
    if (W.M,W.E) in byclass and byclass[(W.M,W.E)]!=p: coll=False
    byclass[(W.M,W.E)]=p
say("  (M,E) collapse on the 66 HOLDOUT windows (independent DP, h=22): %s"%coll)
mind=10**9
for W in HO[:12]:
    for h in range(1,16):
        alive={}; V=W.Vmin
        while confined(W.P+V+h,W.M+h,W.c): alive[W.P+V]=1; V+=1
        j=0
        while alive and j<h:
            nxt={}
            for S in alive:
                d=1
                while confined(S+d+(h-j-1),W.M+h,W.c): d+=1
                mind=min(mind,d)
                for dd in range(1,d):
                    if not confined(S+dd,W.M+j+1,W.c) and j+1<h: nxt[S+dd]=1
            alive=nxt; j+=1
say("  min d* over reachable states of 12 HOLDOUT windows (no clipping iff >=2): %d"%mind)
json.dump(out,open(ROOT+"/out/stage3_holdout.json","w"),indent=1,default=str)
open(ROOT+"/out/stage3_holdout.log","w").write("\n".join(LOG))
