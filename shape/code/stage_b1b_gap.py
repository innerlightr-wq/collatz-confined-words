"""Stage 1b: the decisive diagnostic.  A tie needs Pi(j*+1) = 2 Pi(j*) EXACTLY.
Write Delta = Pi(j*+1) - 2 Pi(j*) (an exact integer).  Ties are Delta = 0.
 - balance mechanism  -> ties should track the RELATIVE margin m (m small => tie likely);
 - coincidence        -> ties track the ABSOLUTE integer gap |Delta| (small integers => tie likely),
   which grows with r even while m shrinks."""
import os as _os, sys as _sys
ROOT = _os.path.dirname(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
_sys.path.insert(0, _os.path.join(ROOT, 'shape','code'))

import sys, csv, math, collections, statistics
sys.path.insert(0,ROOT+"/shape/code")
from slopes import build, Slope
from balance_core import iter_profiles, peak_and_margin
import random
LOG=[]
def say(*a):
    s=" ".join(str(x) for x in a); print(s); LOG.append(s)
HMAX=200
WINDOWS=[(2,1),(3,1),(5,2),(8,1),(13,3),(4,2),(7,1),(11,2)]
S=build(depth=60)
rng=random.Random(20261006); old={tuple(s.cf) for s in S if s.cf}
NEW=[]
while len(NEW)<10:
    cf=[rng.choice([1,1,1,2,2,3,4,5,7]) for _ in range(60)]
    if tuple(cf) in old: continue
    sl=Slope("new-random#%02d"%len(NEW),"random",cf)
    try: sl.verify(HMAX+20); NEW.append(sl)
    except AssertionError: pass
ALL=S+NEW
rows=[]
for sl in ALL:
    for (M,E) in WINDOWS:
        for h,r,L,Pi,c,f in iter_profiles(sl,M,E,HMAX):
            pm=peak_and_margin(Pi,h); js=pm["peak"]
            if js+1>=len(Pi): continue
            D=Pi[js+1]-2*Pi[js]
            rows.append(dict(theta=sl.float_value(),r=r,margin=pm["margin"],tied=pm["tied"],
                             delta=D,pi_peak=Pi[js],
                             digits=len(str(Pi[js])) ))
say("rows: %d"%len(rows))
RB=[1,5,10,15,20,30,50,100,400]
say("")
say("== B1.7 relative margin m vs absolute integer gap |Delta| at the peak, by r ==")
say("   r band    n     tie rate   median m    median |Delta|   median digits of Pi(peak)")
for i in range(len(RB)-1):
    sub=[d for d in rows if RB[i]<=d["r"]<RB[i+1]]
    if not sub: continue
    say("   %3d-%3d %6d   %.4f    %.6f    %-16s %d"
        %(RB[i],RB[i+1],len(sub),sum(d["tied"] for d in sub)/len(sub),
          statistics.median([d["margin"] for d in sub if math.isfinite(d["margin"])]),
          "{:,}".format(int(statistics.median([abs(d["delta"]) for d in sub]))),
          int(statistics.median([d["digits"] for d in sub]))))
say("")
say("   -> the relative margin SHRINKS with r (profile gets more balanced) while the")
say("      absolute integer gap EXPLODES.  Ties follow |Delta|, not m.")
say("")
say("== B1.8 condition on near-balance: among the most balanced rows, how many tie? ==")
for i in range(len(RB)-1):
    sub=[d for d in rows if RB[i]<=d["r"]<RB[i+1] and math.isfinite(d["margin"])]
    if len(sub)<50: continue
    sub.sort(key=lambda d: d["margin"])
    top=sub[:max(1,len(sub)//10)]        # the 10% most balanced rows in this r band
    say("   r %3d-%3d : 10%% most balanced (m <= %.6f) -> tie rate %.4f   (band tie rate %.4f)"
        %(RB[i],RB[i+1],top[-1]["margin"],sum(d["tied"] for d in top)/len(top),
          sum(d["tied"] for d in sub)/len(sub)))
say("")
say("== B1.9 the same, expressed as: does a tie need a SMALL profile? ==")
T=[d for d in rows if d["tied"]]; U=[d for d in rows if not d["tied"]]
say("   tied rows   : median Pi(peak) digits %d, max %d"
    %(statistics.median([d["digits"] for d in T]),max(d["digits"] for d in T)))
say("   untied rows : median Pi(peak) digits %d, max %d"
    %(statistics.median([d["digits"] for d in U]),max(d["digits"] for d in U)))
say("   largest Pi(peak) at which ANY tie occurs: %s (%d digits)"
    %("{:,}".format(max(d["pi_peak"] for d in T)),max(d["digits"] for d in T)))
open(ROOT+"/shape/out/balance_b1b.log","w").write("\n".join(LOG))
