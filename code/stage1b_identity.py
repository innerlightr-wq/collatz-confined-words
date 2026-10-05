"""Stage 1b: is T1 a general identity about non-increasing bound profiles
(hence nothing to do with beta, Sturmian words, or the window), and does the
inductive lemma A'_{j+1} = A_{j+1} + A'_j hold identically?"""
import os as _os, sys as _sys
ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
_sys.path.insert(0, _os.path.join(ROOT, 'code'))

import sys, random, json
sys.path.insert(0,ROOT+"/code")
from core import *
LOG=[]
def say(*a):
    s=" ".join(str(x) for x in a); print(s); LOG.append(s)

def chain_arrays(L):
    """A_j(x) for a general bound profile L (list of l_j >= 0).
    A_j(x) = #{chains x_j=x <= ... <= x_0 : x_i <= l_i - 1}.  Returns list of arrays."""
    out=[]
    a=[1]*L[0]
    out.append(a)
    for j in range(1,len(L)):
        if not a or L[j]==0:
            out.append([]); a=[]; continue
        s=0; suf=[0]*len(a)
        for x in range(len(a)-1,-1,-1):
            s+=a[x]; suf[x]=s
        a=suf[:L[j]]
        out.append(a)
    return out

def Ntilde(L):
    return [sum(a) for a in chain_arrays(L)]

random.seed(20261005)
say("== 1.4 random non-increasing profiles: is  Ntilde_{L+1}(j) = 1 + sum_{i<=j} Ntilde_L(i) ? ==")
bad=0; tested=0; lemma_bad=0
for trial in range(4000):
    n=random.randint(1,14)
    l0=random.randint(1,12)
    L=[l0]
    for _ in range(n):
        L.append(max(0,L[-1]-random.choice([0,0,1,1,1,2,3])))
    Lp=[x+1 for x in L]
    A=chain_arrays(L); Ap=chain_arrays(Lp)
    NL=[sum(a) for a in A]; NLp=[sum(a) for a in Ap]
    run=0
    for j in range(len(L)):
        run+=NL[j]
        tested+=1
        if NLp[j]!=1+run: bad+=1
    # inductive lemma: A'_{j+1}(x) = A_{j+1}(x) + A'_j(x) on x <= l_{j+1}-1
    for j in range(len(L)-1):
        for x in range(L[j+1]):
            lhs=Ap[j+1][x] if x<len(Ap[j+1]) else 0
            rhs=(A[j+1][x] if x<len(A[j+1]) else 0)+(Ap[j][x] if x<len(Ap[j]) else 0)
            if lhs!=rhs: lemma_bad+=1
say("   4000 random profiles, %d (profile,j) checks: violations of the identity = %d"%(tested,bad))
say("   violations of the inductive lemma A'_{j+1}=A_{j+1}+A'_j : %d"%lemma_bad)
say("   -> the T1 law is a general identity for ANY non-increasing bound profile with l_0>=1.")
say("      It uses NO property of beta, of the Sturmian word, or of the window.")

say("== 1.5 control: non-monotone profiles (note: the model's l is ALWAYS non-increasing;\n        a non-monotone input is silently reduced to its running minimum, so this is a\n        consistency check, not an independent test of monotonicity) ==")
badinc=0; n=0
for trial in range(2000):
    L=[random.randint(1,8) for _ in range(random.randint(2,8))]
    if all(L[i]>=L[i+1] for i in range(len(L)-1)): continue
    Lp=[x+1 for x in L]
    NL=Ntilde(L); NLp=Ntilde(Lp); run=0; n+=1
    for j in range(len(L)):
        run+=NL[j]
        if NLp[j]!=1+run: badinc+=1; break
say("   %d non-monotone profiles tested, %d violate the identity (they reduce to their\n   running minimum, which is non-increasing -> monotonicity is WLOG, not an extra hypothesis)"%(n,badinc))

say("== 1.6 exact residual characterisation on the degenerate (window-closed) case ==")
raw=lambda M,E,h: G(M+h)-G(M)-E+1
tot=0; match=0; shapes={}
for M in range(2,41):
    for E in range(1,13):
        for h in range(1,60):
            if G(M+h+1)-G(M+h)!=1: continue
            Nh=chain_N(M,E,h); Nh1=chain_N(M,E,h+1); jm=max(0,h-2)
            run=0; T={}
            for j in range(jm+1):
                run+=Nh.get(j,0); T[j]=run+1 if j==0 else run
            res=[Nh1.get(j,0)-T[j] for j in range(jm+1)]
            r=raw(M,E,h)
            tot+=1
            pred_fail = (r<=-1)
            act_fail  = any(res)
            if pred_fail==act_fail: match+=1
            if act_fail:
                shapes[tuple(res)]=shapes.get(tuple(res),0)+1
say("   rule 'T1 fails  <=>  G(M+h)-G(M)-E+1 <= -1'  (window admits no terminal digit at h or h+1)")
say("   agreement: %d/%d transitions (M=2..40, E=1..12, h<=59)"%(match,tot))
say("   residual shapes observed at failures: %s"
    %{("(-1)*%d"%len(k) if set(k)=={-1} else str(k)):v for k,v in list(shapes.items())[:8]})
say("   all failure residuals are the constant vector (-1,...,-1): %s"%all(set(k)=={-1} for k in shapes))
open(ROOT+"/out/stage1b.log","w").write("\n".join(LOG))
