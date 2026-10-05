"""Stage 0 (cont.): exact reconstruction of the reported '16 of 56 / gaps 5,7,12'."""
import os as _os, sys as _sys
ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
_sys.path.insert(0, _os.path.join(ROOT, 'code'))

import sys, json, math
sys.path.insert(0,ROOT+"/code")
from core import *
LOG=[]
def say(*a):
    s=" ".join(str(x) for x in a); print(s); LOG.append(s)

def letter(M,h): return G(M+h+1)-G(M+h)      # eps_h for window with that M (Sturmian letter)

say("== 0.10 hypothesis: transitions were selected with the M=2 clock, tree computed for M=3 ==")
say("   (i.e. eps was taken as eps^{M=2}_h = e_{2+h} while the window needs e_{3+h})")
hits=[]
for H in range(60,140):
    sel=[h for h in range(1,H+1) if letter(2,h)==1]                 # 'tested one-digit transitions'
    suc=[h for h in sel if letter(3,h)==1]                          # T1 can only hold if the real digit is 1
    if len(sel)==56:
        gaps=[suc[i+1]-suc[i] for i in range(len(suc)-1)]
        say("   H=%d : selected %d transitions, of which %d have the M=3 digit also =1"
            %(H,len(sel),len(suc)))
        say("        success horizons: %s"%suc)
        say("        gaps between successes: %s   (distinct: %s)"%(gaps,sorted(set(gaps))))
        hits.append((H,len(sel),len(suc),sorted(set(gaps))))
say("== 0.11 verification that T1 for the M=3 window holds exactly at e_{3+h}=1 and")
say("        necessarily fails at e_{3+h}=0 (where the true transition is T0=identity) ==")
for E in (1,2):
    nz=ok0=n1=ok1=0
    for h in range(1,100):
        Nh=chain_N(3,E,h); Nh1=chain_N(3,E,h+1); jm=max(0,h-2)
        run=0; T={}
        for j in range(jm+1):
            run+=Nh.get(j,0); T[j]=run+1 if j==0 else run
        exact=all(Nh1.get(j,0)==T[j] for j in range(jm+1))
        if letter(3,h)==1: n1+=1; ok1+=exact
        else:              nz+=1; ok0+=exact
    say("   M=3,E=%d : at e_{3+h}=1 T1 exact %d/%d ; at e_{3+h}=0 T1 exact %d/%d (T0 regime)"
        %(E,ok1,n1,ok0,nz))
say("== 0.12 calibration: are gaps {5,7,12} evidence of anything? ==")
say("   Sturmian return times of the factor '11' (slope beta=log2(3/2)):")
pos=[n for n in range(2,4000) if letter(2,n)==1 and letter(3,n)==1]
g=[pos[i+1]-pos[i] for i in range(len(pos)-1)]
from collections import Counter
say("   gap multiset over 4000 horizons: %s"%dict(Counter(g)))
say("   convergent denominators of beta=[0;1,1,2,2,3,1,5,2,23,...]: 1,2,5,12,41,53,...")
say("   -> by the three-distance theorem the return times of ANY Sturmian factor take")
say("      at most 3 values, which are automatically of the form q_k, q_k+q_{k-1}, q_{k+1}.")
say("   So 'gaps near convergents' is forced for any such set, carrying no information.")
# chance calibration for the stated criterion
say("== 0.13 chance calibration: P(random small integer is 'near a convergent') ==")
conv=[1,2,5,12,41,53,147,347]
semis=set()
for k in range(len(conv)-1):
    for t in range(1,6):
        semis.add(conv[k]*t+ (conv[k-1] if k>0 else 0))
near=lambda x,tol: any(abs(x-q)<=tol for q in list(conv)+sorted(semis))
for tol in (0,1,2):
    for Rmax in (15,30):
        frac=sum(near(x,tol) for x in range(1,Rmax+1))/Rmax
        say("   tol=%d, integers 1..%d : %.0f%% are 'near a convergent/semiconvergent'"%(tol,Rmax,100*frac))
open(ROOT+"/out/stage0d.log","w").write("\n".join(LOG))
