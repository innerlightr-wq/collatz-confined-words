"""Stage 1: systematic sweep of the TRUE window space (M,E) for genuine T1 failures."""
import os as _os, sys as _sys
ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
_sys.path.insert(0, _os.path.join(ROOT, 'code'))

import sys, json, csv
sys.path.insert(0,ROOT+"/code")
from core import *
LOG=[]
def say(*a):
    s=" ".join(str(x) for x in a); print(s); LOG.append(s)

HMAX=150
def transitions(M,E,hmax=HMAX,margin=2):
    out=[]
    prev=None
    for h in range(1,hmax+2):
        cur=chain_N(M,E,h)
        if prev is not None:
            hh=h-1
            eps=G(M+hh+1)-G(M+hh)
            jm=max(0,hh-margin)
            if eps==1:
                run=0; T={}
                for j in range(jm+1):
                    run+=prev.get(j,0); T[j]=run+1 if j==0 else run
                resid={j:cur.get(j,0)-T[j] for j in range(jm+1)}
            else:
                resid={j:cur.get(j,0)-prev.get(j,0) for j in range(jm+1)}
            bad={j:v for j,v in resid.items() if v}
            out.append(dict(h=hh,eps=eps,jmax=jm,ok=not bad,bad=bad,
                            l0=l_profile(M,E,hh)[0],
                            nlev=max(prev) if prev else 0))
        prev=cur
    return out

rows=[]
say("== 1.1 sweep M=2..40, E=1..12, h<=%d, margin=2 =="%HMAX)
say("   (reporting only classes with at least one failure)")
for M in range(2,41):
    for E in range(1,13):
        tr=transitions(M,E,HMAX,margin=2)
        o=[r for r in tr if r["eps"]==1]; z=[r for r in tr if r["eps"]==0]
        f1=[r for r in o if not r["ok"]]; f0=[r for r in z if not r["ok"]]
        rows.append(dict(M=M,E=E,n_eps1=len(o),ok_eps1=len(o)-len(f1),
                         n_eps0=len(z),ok_eps0=len(z)-len(f0),
                         fail_h=[r["h"] for r in f1],
                         max_fail_h=max([r["h"] for r in f1],default=0),
                         fail_l0=[r["l0"] for r in f1],
                         T0_fail=[r["h"] for r in f0]))
        if f1 or f0:
            say("   M=%2d E=%2d : T1 %3d/%3d  (fail h=%s, l0 there=%s)  T0 %d/%d %s"
                %(M,E,len(o)-len(f1),len(o),[r['h'] for r in f1],[r['l0'] for r in f1],
                  len(z)-len(f0),len(z),[r['h'] for r in f0] or ""))
tot1=sum(r["n_eps1"] for r in rows); ok1=sum(r["ok_eps1"] for r in rows)
tot0=sum(r["n_eps0"] for r in rows); ok0=sum(r["ok_eps0"] for r in rows)
say("   TOTAL: T1 exact on %d/%d eps=1 transitions (%d classes); T0 exact on %d/%d"
    %(ok1,tot1,len(rows),ok0,tot0))
allf=[h for r in rows for h in r["fail_h"]]
say("   every T1 failure occurs at h in %s ; max failing h = %d"
    %(sorted(set(allf)),max(allf,default=0)))
say("   l_0 values at failures: %s"%sorted(set(l for r in rows for l in r["fail_l0"])))
# what characterises the failures?
say("== 1.2 characterisation of the failing transitions ==")
fl=[(r["M"],r["E"],h,l) for r in rows for h,l in zip(r["fail_h"],r["fail_l0"])]
say("   #failures=%d"%len(fl))
for M,E,h,l in fl[:40]:
    L=l_profile(M,E,h); L2=l_profile(M,E,h+1)
    say("    M=%2d E=%2d h=%2d : l-profile(h)=%s  l(h+1)=%s"%(M,E,h,L[:8],L2[:8]))
say("== 1.3 same sweep with margin 0 (j<=h, no exclusion) ==")
bad0=[]
for M in range(2,41):
    for E in range(1,13):
        tr=transitions(M,E,60,margin=0)
        o=[r for r in tr if r["eps"]==1]
        f1=[r for r in o if not r["ok"]]
        if f1: bad0.append((M,E,[r["h"] for r in f1]))
say("   classes with any eps=1 failure at margin 0: %d"%len(bad0))
say("   e.g. %s"%bad0[:6])
with open(ROOT+"/data/sweep_ME.csv","w",newline="") as f:
    w=csv.writer(f); w.writerow(["M","E","n_eps1","ok_eps1","n_eps0","ok_eps0","fail_h","max_fail_h"])
    for r in rows: w.writerow([r["M"],r["E"],r["n_eps1"],r["ok_eps1"],r["n_eps0"],r["ok_eps0"],
                               ";".join(map(str,r["fail_h"])),r["max_fail_h"]])
json.dump(rows,open(ROOT+"/out/sweep_ME.json","w"),indent=0,default=str)
open(ROOT+"/out/stage1_sweep.log","w").write("\n".join(LOG))
