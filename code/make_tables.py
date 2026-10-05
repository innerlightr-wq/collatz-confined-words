import os as _os, sys as _sys
ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
_sys.path.insert(0, _os.path.join(ROOT, 'code'))
import sys,csv,math
sys.path.insert(0,ROOT+"/code")
from core import *
from windows import split
EX,HO=split()
rows=[]
cache={}
def counts(M,E,hmax=150,margin=2):
    if (M,E) in cache: return cache[(M,E)]
    n1=ok1=n0=ok0=0; fails=[]
    prev=chain_N(M,E,1)
    for h in range(1,hmax+1):
        cur=chain_N(M,E,h+1); eps=G(M+h+1)-G(M+h); jm=max(0,h-margin)
        if eps==1:
            run=0;T={}
            for j in range(jm+1):
                run+=prev.get(j,0); T[j]=run+1 if j==0 else run
            bad=any(cur.get(j,0)!=T[j] for j in range(jm+1))
            n1+=1; ok1+= not bad
            if bad: fails.append(h)
        else:
            bad=any(cur.get(j,0)!=prev.get(j,0) for j in range(jm+1))
            n0+=1; ok0+= not bad
        prev=cur
    cache[(M,E)]=(n1,ok1,n0,ok0,fails)
    return cache[(M,E)]
for tag,lst in (("EXPLORE",EX),("HOLDOUT",HO)):
    for W in lst:
        n1,ok1,n0,ok0,f=counts(W.M,W.E)
        hopen=min([h for h in range(1,200) if G(W.M+h)-G(W.M)-W.E+1>=0],default=None)
        rows.append(dict(split=tag,M=W.M,c=W.c,prefix="-".join(map(str,W.prefix)),
                         V_min=W.Vmin,E=W.E,gamma_W="%.6f"%W.gamma,
                         n_eps1=n1,T1_exact=ok1,n_eps0=n0,T0_exact=ok0,
                         h_open=hopen,fail_h=";".join(map(str,f))))
with open(ROOT+"/data/window_table_full.csv","w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print("wrote %d rows"%len(rows))
# markdown: one row per (M,E) class with a representative window
seen={}
for r in rows: seen.setdefault((r["M"],r["E"]),r)
hdr="| M | c | prefix | V_min | E | gamma_W | h_open | #eps=1 (h<=150) | #exact T1 | #eps=0 | #exact T0 | T1 failures at h |"
print(hdr); print("|"+"---|"*12)
md=[hdr,"|"+"---|"*12]
for k in sorted(seen):
    r=seen[k]
    line="| %d | %d | (%s) | %d | %d | %s | %s | %d | %d | %d | %d | %s |"%(
        r["M"],r["c"],r["prefix"],r["V_min"],r["E"],r["gamma_W"],r["h_open"],
        r["n_eps1"],r["T1_exact"],r["n_eps0"],r["T0_exact"],r["fail_h"] or "none")
    print(line); md.append(line)
open(ROOT+"/out/window_table.md","w").write("\n".join(md))
