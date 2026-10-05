"""The exploration window set (M in 2..6, c in 1..3, several prefixes and V_min)
and the fixed-seed 50/50 EXPLORE / HOLDOUT split."""
import os as _os, sys as _sys
ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
_sys.path.insert(0, _os.path.join(ROOT, 'code'))

import sys, random, itertools, csv
sys.path.insert(0,ROOT+"/code")
from core import *

SEED = 20261005

def valid_prefixes(M, c, maxdigit=4):
    out=[]
    for pf in itertools.product(range(1,maxdigit+1), repeat=M-1):
        S=0; ok=True
        for j,d in enumerate(pf):
            S+=d
            if not confined(S,j+1,c): ok=False; break
        if ok: out.append(pf)
    return out

def build():
    rng=random.Random(SEED)
    wins=[]
    for M in range(2,7):
        for c in (1,2,3):
            pfs=valid_prefixes(M,c)
            pick=pfs[:1]+rng.sample(pfs[1:], min(2,len(pfs)-1)) if len(pfs)>1 else pfs
            for pf in pick:
                P=sum(pf)
                for extra in (0,1,3):            # E = 1, 2, 4
                    Vmin=c+F(M)+1-P+extra
                    if Vmin<1: continue
                    try: wins.append(Window(M,c,pf,Vmin))
                    except AssertionError: pass
    # de-duplicate identical (M,c,prefix,Vmin)
    seen=set(); out=[]
    for W in wins:
        k=(W.M,W.c,W.prefix,W.Vmin)
        if k in seen: continue
        seen.add(k); out.append(W)
    return out

def split():
    wins=build()
    rng=random.Random(SEED+1)
    idx=list(range(len(wins))); rng.shuffle(idx)
    half=len(idx)//2
    ex=set(idx[:half]); 
    explore=[wins[i] for i in sorted(ex)]
    holdout=[wins[i] for i in sorted(set(idx)-ex)]
    return explore, holdout

if __name__=="__main__":
    ex,ho=split()
    print("total windows %d : EXPLORE %d, HOLDOUT %d"%(len(ex)+len(ho),len(ex),len(ho)))
    with open(ROOT+"/data/windows.csv","w",newline="") as f:
        w=csv.writer(f)
        w.writerow(["split","name","M","c","prefix","P","V_min","E","gamma_W=frac(M*alpha)"])
        for tag,lst in (("EXPLORE",ex),("HOLDOUT",ho)):
            for W in lst:
                w.writerow([tag,W.tag(),W.M,W.c,"-".join(map(str,W.prefix)),W.P,W.Vmin,W.E,
                            "%.6f"%W.gamma])
    from collections import Counter
    print("EXPLORE (M,E) classes:",sorted(Counter((W.M,W.E) for W in ex).items()))
    print("HOLDOUT (M,E) classes:",sorted(Counter((W.M,W.E) for W in ho).items()))
