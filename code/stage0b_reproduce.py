"""Stage 0 (cont.): reproduce the W0 reference facts; T0 in M=3; the M=3 T1 failure."""
import os as _os, sys as _sys
ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
_sys.path.insert(0, _os.path.join(ROOT, 'code'))

import json, sys, math
from fractions import Fraction
sys.path.insert(0,ROOT+"/code")
from core import *

LOG=[]
def say(*a):
    s=" ".join(str(x) for x in a); print(s); LOG.append(s)

HMAX=152
def profiles(M,E,hmax=HMAX):
    return {h: chain_N(M,E,h) for h in range(1,hmax+1)}

def T1_variant(N,jmax,propagate=False):
    out={}; run=0
    for j in range(jmax+1):
        run+=N.get(j,0)
        out[j]= run+1 if (j==0 or propagate) else run
    return out

def test_transitions(M,E,hmax=HMAX,margin=2,propagate=False):
    """Returns per-h records for eps=0 (T0) and eps=1 (T1) transitions."""
    pr=profiles(M,E,hmax)
    rec=[]
    for h in range(1,hmax):
        eps=G(M+h+1)-G(M+h)
        Nh,Nh1=pr[h],pr[h+1]
        jmax=max(0,h-margin)          # interior: j <= h-margin
        if eps==0:
            bad=[j for j in range(jmax+1) if Nh.get(j,0)!=Nh1.get(j,0)]
            rec.append(dict(h=h,eps=0,ok=not bad,bad=bad))
        else:
            T=T1_variant(Nh,jmax,propagate)
            res={j:(Nh1.get(j,0)-T[j]) for j in range(jmax+1)}
            bad=[j for j in res if res[j]!=0]
            rec.append(dict(h=h,eps=1,ok=not bad,bad=bad,
                            resid={j:res[j] for j in bad}))
    return rec

res={}
# ---------------------------------------------------------------- W0 facts
say("== 0.4 W0 (M=2,c=1,prefix (1),V_min=4 -> E=1): T0 and T1 ==")
W0=Window(2,1,(1,),4)
say("  E(W0) =",W0.E,"  gamma_W = frac(M*alpha) =","%.6f"%W0.gamma)
for prop in (False,True):
    rec=test_transitions(2,1,HMAX,margin=2,propagate=prop)
    z=[r for r in rec if r["eps"]==0]; o=[r for r in rec if r["eps"]==1]
    say("  T1 convention propagate_bump=%s : eps=1 transitions %d, exact %d"
        %(prop,len(o),sum(r["ok"] for r in o)))
    if not prop:
        say("  T0 (eps=0) transitions %d, exact %d  [margin 2]"%(len(z),sum(r["ok"] for r in z)))
        res["W0"]=dict(n_eps1=len(o),ok_eps1=sum(r["ok"] for r in o),
                       n_eps0=len(z),ok_eps0=sum(r["ok"] for r in z))
# T0 at full range j<=h  (sharper than the note's j<=h-2)
rec=test_transitions(2,1,HMAX,margin=0)
z=[r for r in rec if r["eps"]==0]
say("  T0 with NO margin (j<=h): %d/%d exact"%(sum(r['ok'] for r in z),len(z)))
res["W0"]["T0_no_margin"]=[sum(r['ok'] for r in z),len(z)]
# N_h(h)==0 claim
nz=all(chain_N(2,1,h).get(h,0)==0 for h in range(1,HMAX))
say("  N_h(h) == 0 for all h (exhaustion clause never fires):",nz)
res["W0"]["N_h_h_zero"]=bool(nz)

# ---------------------------------------------------------------- Haar stats
say("== 0.5 W0 Haar-weighted bulk statistics (weights N_h(j)*2^-j) ==")
def bulk(M,E,h):
    N=chain_N(M,E,h); jm=max(N)
    w=[Fraction(N.get(j,0),2**j) for j in range(jm+1)]
    tot=sum(w); 
    jstar=max(range(jm+1),key=lambda j:(w[j],-j))
    ties=[j for j in range(jm+1) if w[j]==w[jstar]]
    EJ=sum(j*w[j] for j in range(jm+1))/tot
    EJ2=sum(j*j*w[j] for j in range(jm+1))/tot
    return jstar,ties,float(EJ),float(EJ2-EJ*EJ)
rows=[]
for h in (40,60,80,100):
    jstar,ties,EJ,VJ=bulk(2,1,h)
    C=W0.CW(h)                       # = floor(c+M*alpha+h*beta)
    Cp=F(2+h)-h                      # = floor(2*alpha+h*beta)  (note's C(h))
    rows.append(dict(h=h,jstar=jstar,ties=ties,EJ=EJ,VJ=VJ,CW=C,C_note=Cp))
    say("  h=%3d  mode j*=%2d (ties %s)  E[J]=%.4f  Var=%.4f  C_W=%d  C_note=%d"
        %(h,jstar,ties,EJ,VJ,C,Cp))
# find the reference h0 that makes mode = K(h)-1 with K(h)=C(h)-C(h0)
for name,Cf in (("C_note",lambda h: F(2+h)-h),("C_W",lambda h: W0.CW(h))):
    for h0 in range(0,20):
        if all(r["jstar"]==Cf(r["h"])-Cf(h0)-1 for r in rows):
            say("  mode j*(h) = K(h)-1 holds with K(h)=%s(h)-%s(%d)"%(name,name,h0))
            for r in rows: r["K"]=Cf(r["h"])-Cf(h0)
            say("  E[J]-K =", ["%+.3f"%(r["EJ"]-r["K"]) for r in rows],
                "  Var/K =", ["%.3f"%(r["VJ"]/r["K"]) for r in rows])
            break
res["haar"]=rows

# ---------------------------------------------------------------- M=3
say("== 0.6 M=3 windows: T0, and the T1 failure ==")
m3=[]
for E in (1,2,3,4):
    rec=test_transitions(3,E,HMAX,margin=2)
    z=[r for r in rec if r["eps"]==0]; o=[r for r in rec if r["eps"]==1]
    say("  M=3,E=%d : T0 %d/%d exact | T1 %d/%d exact"
        %(E,sum(r['ok'] for r in z),len(z),sum(r['ok'] for r in o),len(o)))
    m3.append(dict(E=E,T0=[sum(r['ok'] for r in z),len(z)],
                   T1=[sum(r['ok'] for r in o),len(o)],
                   ok_h=[r["h"] for r in o if r["ok"]]))
res["M3"]=m3
# margin sensitivity + residual anatomy for M=3,E=1
say("  margin sensitivity (M=3): #exact T1 out of #eps=1, by excluded margin")
for M in (2,3,4,5,6):
    for E in (1,2):
        line=[]
        for mg in (0,1,2,3,5,8):
            rec=[r for r in test_transitions(M,E,80,margin=mg) if r["eps"]==1]
            line.append("%d/%d"%(sum(r['ok'] for r in rec),len(rec)))
        say("    M=%d E=%d  margins 0,1,2,3,5,8 -> %s"%(M,E," ".join(line)))
res["margin"]="see log"
say("  residual anatomy, M=3 E=1, first failing transitions (margin 2):")
rec=[r for r in test_transitions(3,1,60,margin=2) if r["eps"]==1]
for r in rec[:8]:
    if not r["ok"]:
        ks=sorted(r["resid"]); 
        say("    h=%2d  #bad j=%d  j range [%d,%d]  residuals %s"
            %(r["h"],len(ks),ks[0],ks[-1],[r["resid"][k] for k in ks[:6]]))
    else:
        say("    h=%2d  EXACT"%r["h"])
json.dump(res,open(ROOT+"/out/stage0b.json","w"),indent=1,default=str)
open(ROOT+"/out/stage0b.log","w").write("\n".join(LOG))
