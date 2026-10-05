"""Stage 3 addendum: (i) why the 'best phase interval' is a tautology,
(ii) brute-force re-validation on HOLDOUT windows,
(iii) the two profile-shift lemmas that the T1 theorem rests on."""
import os as _os, sys as _sys
ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
_sys.path.insert(0, _os.path.join(ROOT, 'code'))

import sys, math
sys.path.insert(0,ROOT+"/code")
from core import *
from windows import split
LOG=[]
def say(*a):
    s=" ".join(str(x) for x in a); print(s); LOG.append(s)
EX,HO=split()

say("== 3.7 the 'best fitted phase interval' is a tautology ==")
beta=math.log2(3)-1
say("   eps_h = 1  <=>  floor((h+1)beta+gamma) > floor(h beta+gamma)  <=>  frac(h beta+gamma) >= 1-beta")
say("   1-beta = %.6f ; observed phase range at eps=1 transitions = [0.4181, 0.9895]"%(1-beta))
say("   so ANY interval [t,1) with t <= 1-beta contains every eps=1 transition and is")
say("   literally the always-success predictor.  Zero information, as predicted (P2).")
bad=0;n=0
for M in range(2,60):
    for h in range(1,200):
        ph=(h*beta+(M*math.log2(3))%1)%1
        e=G(M+h+1)-G(M+h); n+=1
        if e!=(1 if ph>=1-beta-1e-12 else 0): bad+=1
say("   check over %d (M,h) pairs: mismatches = %d"%(n,bad))

say("== 3.8 brute-force recursive tree vs DP vs chain model on HOLDOUT windows ==")
ok=True;cnt=0
for W in HO[:10]:
    for h in range(1,13):
        b={k:v for k,v in brute_N(W,h).items() if v}
        d={k:v for k,v in dp_N(W,h,literal_thresholds=True).items() if v}
        c={k:v for k,v in chain_N(W.M,W.E,h).items() if v}
        cnt+=1
        if not (b==d==c): ok=False; say("   MISMATCH",W.tag(),h,b,d,c)
say("   %d (HOLDOUT window, h) pairs: brute == DP == chain model : %s"%(cnt,ok))

say("== 3.9 the two profile-shift lemmas underlying the theorem ==")
b1=b0=0;n1=n0=0
for M in range(2,121):
    for E in range(1,13):
        for h in range(1,151):
            e=G(M+h+1)-G(M+h)
            L=l_profile(M,E,h); L2=l_profile(M,E,h+1)
            r=G(M+h)-G(M)-E+1
            if e==1 and r>=0:
                n1+=1
                if any(L2[j]!=L[j]+1 for j in range(h+1)): b1+=1
            if e==0:
                n0+=1
                if any(L2[j]!=L[j] for j in range(h+1)): b0+=1
say("   eps_h=1 and r>=0  =>  l_j(h+1) = l_j(h)+1 for all j<=h : violations %d / %d"%(b1,n1))
say("   eps_h=0           =>  l_j(h+1) = l_j(h)   for all j<=h : violations %d / %d"%(b0,n0))
open(ROOT+"/out/stage3b.log","w").write("\n".join(LOG))
