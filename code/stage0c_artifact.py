"""Stage 0 (cont.): is the M=3 T1 failure real?  (a) independent-implementation check,
(b) can plausible buggy generalizations reproduce '16 of 56'?"""
import os as _os, sys as _sys
ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
_sys.path.insert(0, _os.path.join(ROOT, 'code'))

import sys, json
sys.path.insert(0,ROOT+"/code")
from core import *
LOG=[]
def say(*a):
    s=" ".join(str(x) for x in a); print(s); LOG.append(s)

def T1c(N,jmax):
    out={};run=0
    for j in range(jmax+1):
        run+=N.get(j,0); out[j]=run+1 if j==0 else run
    return out

def mkwin(M,c,prefix,extra=0):
    P=sum(prefix); Vmin=c+F(M)+1-P+extra
    return Window(M,c,prefix,Vmin)

say("== 0.7 M=3 T1 via the INDEPENDENT DP (no chain-model, no closed forms) ==")
for E in (1,2,3,4):
    W=mkwin(3,1,(1,1),E-1)
    assert W.E==E
    tot=okc=0; fails=[]
    for h in range(1,31):
        if G(3+h+1)-G(3+h)!=1: continue
        Nh ={k:v for k,v in dp_N(W,h,literal_thresholds=True).items()}
        Nh1={k:v for k,v in dp_N(W,h+1,literal_thresholds=True).items()}
        jm=max(0,h-2); T=T1c(Nh,jm)
        bad=[j for j in range(jm+1) if Nh1.get(j,0)!=T[j]]
        tot+=1; okc+= (not bad)
        if bad: fails.append((h,bad))
    say("  M=3 E=%d (c=1, prefix (1,1), V_min=%d): literal-DP T1 exact on %d/%d eps=1 transitions, h<=30  fails=%s"
        %(E,W.Vmin,okc,tot,fails))

say("== 0.8 do buggy generalizations reproduce '16 of 56'? ==")
# Each variant computes N_h(j) for an M=3,c=1 window with one hard-coded M=2 leftover.
def buggy_N(h,variant,M=3,c=1,P=2,Vmin=3):
    """Mirror of dp_N with one deliberately-wrong ingredient."""
    N={}; alive={}
    V=Vmin
    while True:
        if variant=="interior1":   ok=confined(1+V+h,M+h,c)       # S_INTERIOR=1 left from M=2
        elif variant=="Mh_is_2h":  ok=confined(P+V+h,2+h,c)       # best-case index uses M=2
        else:                      ok=confined(P+V+h,M+h,c)
        if not ok: break
        base = (1+V) if variant=="interior1" else (P+V)
        alive[base]=alive.get(base,0)+1; V+=1
    N[0]=1; j=0
    while alive:
        nxt={}; pruned=0
        for S,cnt in alive.items():
            d=1
            while (confined(S+d+(h-j-1),2+h,c) if variant=="Mh_is_2h"
                   else confined(S+d+(h-j-1),M+h,c)): d+=1
            pruned+=cnt
            for dd in range(1,d):
                idx = (2+j+1) if variant=="rho_M2" else (M+j+1)
                if confined(S+dd,idx,c): continue
                if j+1==h: N[h]=N.get(h,0)+cnt
                else: nxt[S+dd]=nxt.get(S+dd,0)+cnt
        N[j]=N.get(j,0)+pruned
        if j+1==h: break
        alive=nxt; j+=1
    return {k:v for k,v in N.items() if v}

for variant in ("correct","interior1","Mh_is_2h","rho_M2","Vmin4"):
    kw={}
    if variant=="Vmin4": f=lambda h: dp_N(Window(3,1,(1,1),4),h)   # V_min=4 copied from W0 -> E=2
    else: f=lambda h,v=variant: buggy_N(h,v)
    tot=0; ok=0; okh=[]
    for h in range(1,100):
        if G(3+h+1)-G(3+h)!=1: continue
        try:
            Nh=f(h); Nh1=f(h+1)
        except Exception as e:
            say("  %s: error %s"%(variant,e)); break
        jm=max(0,h-2); T=T1c(Nh,jm)
        bad=[j for j in range(jm+1) if Nh1.get(j,0)!=T[j]]
        tot+=1
        if not bad: ok+=1; okh.append(h)
    gaps=[okh[i+1]-okh[i] for i in range(len(okh)-1)]
    say("  variant %-10s : T1 exact %2d/%2d   success horizons %s"%(variant,ok,tot,okh[:14]))
    say("                      gaps between successes: %s"%gaps[:14])
# also: how big are the discrepancies for the best-matching variant?
say("== 0.9 size of discrepancies for variant 'interior1' (M=3) ==")
for h in range(1,40):
    if G(3+h+1)-G(3+h)!=1: continue
    Nh=buggy_N(h,"interior1"); Nh1=buggy_N(h+1,"interior1")
    jm=max(0,h-2); T=T1c(Nh,jm)
    bad={j:Nh1.get(j,0)-T[j] for j in range(jm+1) if Nh1.get(j,0)!=T[j]}
    if bad and h<30:
        ks=sorted(bad)
        say("   h=%2d #bad=%d  j in [%d,%d]  max|resid|=%d  resid[:5]=%s"
            %(h,len(ks),ks[0],ks[-1],max(abs(v) for v in bad.values()),[bad[k] for k in ks[:5]]))
open(ROOT+"/out/stage0c.log","w").write("\n".join(LOG))
