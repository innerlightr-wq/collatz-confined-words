"""2-adic test for tie exclusion.  For every row with r >= 20 compute, exactly, at the peak j*:
      Delta  = Pi(j*+1) - 2 Pi(j*)      (<= 0; = 0 exactly at a tie)
      Delta' = Pi(j*)   - 2 Pi(j*-1)    (>  0 with the smallest-index convention)
and their 2-adic valuations.  Compare v2 with the geometric law P(v2 >= k) = 2^-k expected of a
random integer, test boundedness, and check the mod-2^m reduction of the operator dynamics."""
import os as _os, sys as _sys
ROOT = _os.path.dirname(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
_sys.path.insert(0, _os.path.join(ROOT, 'shape','code'))

import sys, csv, json, math, random, collections, statistics
sys.path.insert(0,ROOT+"/shape/code")
from slopes import build, Slope
from balance_core import iter_profiles, peak_and_margin, T1
LOG=[]
def say(*a):
    s=" ".join(str(x) for x in a); print(s); LOG.append(s)

def v2(n):
    if n==0: return None
    n=abs(n); k=0
    while n & 1 == 0: n>>=1; k+=1
    return k

# ---- the slope sets already computed in this study ----
S=build(depth=60)
rng=random.Random(20261006); old={tuple(s.cf) for s in S if s.cf}
RETEST=[]
while len(RETEST)<10:
    cf=[rng.choice([1,1,1,2,2,3,4,5,7]) for _ in range(60)]
    if tuple(cf) in old: continue
    sl=Slope("new-random#%02d"%len(RETEST),"random",cf)
    try: sl.verify(420); RETEST.append(sl); old.add(tuple(cf))
    except AssertionError: pass
rng2=random.Random(20261007)
BANDS=[(0.15,0.24),(0.24,0.33),(0.33,0.42),(0.42,0.51),(0.51,0.60),(0.60,0.70),(0.70,0.80),(0.80,0.90)]
HOLD=[]
for lo,hi in BANDS:
    for _ in range(20000):
        cf=[rng2.choice([1,1,1,2,2,3,4,5,7]) for _ in range(60)]
        if tuple(cf) in old: continue
        s2=Slope("hold-random#%d"%len(HOLD),"random",cf)
        if not (lo<=s2.float_value()<hi): continue
        try: s2.verify(420)
        except AssertionError: continue
        HOLD.append(s2); old.add(tuple(cf)); break
JOBS=([(sl,(M,E),200) for sl in S+RETEST
       for (M,E) in [(2,1),(3,1),(5,2),(8,1),(13,3),(4,2),(7,1),(11,2)]]
    + [(sl,(M,E),400) for sl in HOLD for (M,E) in [(2,1),(4,2),(7,1)]])
say("slopes: %d study + %d holdout ; (slope,window,hmax) jobs: %d"%(len(S+RETEST),len(HOLD),len(JOBS)))

rows=[]
for sl,(M,E),HM in JOBS:
    for h,r,L,Pi,c,f in iter_profiles(sl,M,E,HM):
        if r<20: continue
        pm=peak_and_margin(Pi,h); js=pm["peak"]
        if js<1 or js+1>=len(Pi): continue
        D  = Pi[js+1]-2*Pi[js]
        Dp = Pi[js]  -2*Pi[js-1]
        rows.append(dict(slope=sl.name,theta=sl.float_value(),M=M,E=E,h=h,r=r,peak=js,
                         tied=pm["tied"],v2=v2(D),v2p=v2(Dp),
                         dmod=(D % 256),dpmod=(Dp % 256),
                         dbits=D.bit_length(),dpbits=Dp.bit_length()))
say("rows with r >= 20: %d  (ties among them: %d)"%(len(rows),sum(d["tied"] for d in rows)))
with open(ROOT+"/shape/data/twoadic.csv","w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)

say("")
say("== 1. distribution of v2(Delta) and v2(Delta') ==")
for nm,key in (("v2(Delta)","v2"),("v2(Delta')","v2p")):
    c=collections.Counter(d[key] for d in rows if d[key] is not None)
    n=sum(c.values()); mx=max(c)
    say("   %s : n=%d, max = %d, mean = %.4f"%(nm,n,mx,sum(k*v for k,v in c.items())/n))
    say("      k :  count    observed P(v2 = k)   geometric 2^-(k+1)   observed P(v2 >= k)   2^-k")
    cum=n
    for k in range(0,max(12,mx+1)):
        cnt=c.get(k,0)
        say("     %2d : %7d      %.6f            %.6f           %.6f        %.6f"
            %(k,cnt,cnt/n,2.0**-(k+1),cum/n,2.0**-k))
        cum-=cnt
        if cum==0 and k>=mx: break
say("")
say("== 2. is max v2 bounded, and does it grow with r? ==")
RB=[20,30,50,80,120,200,400]
say("   r band      n     max v2(Delta)  max v2(Delta')   mean v2   median bits of Delta")
for i in range(len(RB)-1):
    sub=[d for d in rows if RB[i]<=d["r"]<RB[i+1]]
    if not sub: continue
    say("   %3d-%3d %7d   %6d        %6d         %.4f    %d"
        %(RB[i],RB[i+1],len(sub),max(d["v2"] for d in sub),max(d["v2p"] for d in sub),
          statistics.mean([d["v2"] for d in sub]),
          int(statistics.median([d["dbits"] for d in sub]))))
K0=max(d["v2"] for d in rows); K0p=max(d["v2p"] for d in rows)
say("   overall max v2(Delta) = %d ; max v2(Delta') = %d"%(K0,K0p))
# expected max under the geometric law
N=len(rows); exp_max=math.log2(N)
say("   under the geometric law the expected max over n=%d draws is about log2(n) = %.1f"%(N,exp_max))
say("   -> observed max %d vs %.1f expected: %s"
    %(K0,exp_max,"consistent with random" if K0>=exp_max-2 else "SMALLER than random"))
say("")
say("== 3. residues of Delta modulo small powers of 2 ==")
for m in (2,4,8,16):
    c=collections.Counter(d["dmod"]%m for d in rows)
    say("   Delta mod %2d : %s"%(m,dict(sorted(c.items()))))
for m in (2,4,8,16):
    c=collections.Counter(d["dpmod"]%m for d in rows)
    say("   Delta' mod %2d: %s"%(m,dict(sorted(c.items()))))
say("")
say("== 3b. does the profile mod 2^m evolve exactly under the reduced operators? ==")
say("   (T0 = identity, T1 = cumulative sum with +1 at j=0 are integer affine maps, so")
say("    reduction mod 2^m must commute with them; verified directly below.)")
ok=bad=0
for sl,(M,E),HM in JOBS[:40]:
    for m in (8,16):
        MOD=1<<m; red=None
        for h,r,L,Pi,c,f in iter_profiles(sl,M,E,min(HM,120)):
            if red is None:
                red=[x%MOD for x in Pi]
            else:
                eta=sl.G(M+h)-sl.G(M+h-1)
                if eta==0: red=red+[0]
                else:
                    out=[];run=0
                    for j,v in enumerate(red):
                        run=(run+v)%MOD; out.append((run+1)%MOD if j==0 else run)
                    red=out+[0]
            tru=[x%MOD for x in Pi]
            if red[:len(tru)]==tru: ok+=1
            else: bad+=1
say("   checks: %d matched, %d mismatched"%(ok,bad))
json.dump(dict(rows=len(rows),max_v2=K0,max_v2p=K0p),
          open(ROOT+"/shape/out/twoadic.json","w"),indent=1)
open(ROOT+"/shape/out/twoadic.log","w").write("\n".join(LOG))
