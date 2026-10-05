"""Stage 1 (EXPLORE slopes only): shape statistics across slope families and windows."""
import os as _os, sys as _sys
ROOT = _os.path.dirname(_os.path.dirname(_os.path.dirname(_os.path.abspath(__file__))))
_sys.path.insert(0, _os.path.join(ROOT, 'shape','code'))

import sys, csv, json, collections
sys.path.insert(0, ROOT+"/shape/code")
from slopes import build, split
from shape_core import row, K_of

LOG = []
def say(*a):
    s = " ".join(str(x) for x in a); print(s); LOG.append(s)

S = build(depth=60)
for s in S: s.verify(460)
EX, HO = split(S)
say("EXPLORE slopes (%d): %s" % (len(EX), ", ".join(s.name for s in EX)))
say("HOLDOUT is untouched by this script.")

WINDOWS = [(2, 1), (3, 1), (5, 2), (8, 1), (13, 3)]
HMAX = 150
TVEVERY = 25
rows = []
for sl in EX:
    for (M, E) in WINDOWS:
        for h in range(1, HMAX + 1):
            r = row(sl, M, E, h, want_tv=(h % TVEVERY == 0))
            if r is None: continue
            r["S1r"] = r["mode"] - r["r"]          # convention-free: mode relative to ell_0
            rows.append(r)
say("rows: %d" % len(rows))

fields = ["slope","family","M","E","h","K","r","eta","mode","mode_tied","n_modes",
          "S1","S1r","S2","S3","S4","S5","mean","var","frac_h","frac_M","top_incr"]
with open(ROOT+"/shape/data/explore_shape.csv","w",newline="") as f:
    w=csv.DictWriter(f,fieldnames=fields,extrasaction="ignore"); w.writeheader()
    for r in rows: w.writerow(r)

say("== 1.1 S1 = mode - K(h), and the convention-free S1r = mode - r(h) ==")
say("   (K(h) = G(M+h)-G(M+5) is the study code's definition; it ignores E, so S1 shifts with E)")
for (M,E) in WINDOWS:
    sub=[r for r in rows if (r["M"],r["E"])==(M,E) and r["h"]>=30]
    s1=collections.Counter(r["S1"] for r in sub); s1r=collections.Counter(r["S1r"] for r in sub)
    say("   window (M,E)=(%2d,%d): S1 values %s | S1r values %s"
        %(M,E,dict(s1.most_common(4)),dict(s1r.most_common(4))))
say("== 1.2 is the mode ever tied, or S1r non-constant, for h >= 30? ==")
tied=[r for r in rows if r["mode_tied"] and r["h"]>=30]
say("   tied modes: %d of %d rows with h>=30"%(len(tied),len([r for r in rows if r["h"]>=30])))
byc=collections.defaultdict(set)
for r in rows:
    if r["h"]>=30: byc[(r["slope"],r["M"],r["E"])].add(r["S1r"])
nonconst={k:sorted(v) for k,v in byc.items() if len(v)>1}
say("   (slope,window) classes with non-constant S1r on h>=30: %d of %d"%(len(nonconst),len(byc)))
for k,v in list(nonconst.items())[:8]: say("      %s -> %s"%(k,v))
say("== 1.3 S1r by family and window (h>=30) ==")
for fam in sorted(set(r["family"] for r in rows)):
    line=[]
    for (M,E) in WINDOWS:
        v=sorted({r["S1r"] for r in rows if r["family"]==fam and (r["M"],r["E"])==(M,E) and r["h"]>=30})
        line.append("(%d,%d):%s"%(M,E,v))
    say("   %-10s %s"%(fam," ".join(line)))
say("== 1.4 S2 = mean - K and S3 = Var/K, by family, at h = 50,100,150 ==")
for fam in sorted(set(r["family"] for r in rows)):
    for h in (50,100,150):
        sub=[r for r in rows if r["family"]==fam and r["h"]==h and (r["M"],r["E"])==(2,1)]
        if not sub: continue
        say("   %-10s h=%3d  S2 in [%+.4f,%+.4f]  S3 in [%.4f,%.4f]  (n=%d)"
            %(fam,h,min(r["S2"] for r in sub),max(r["S2"] for r in sub),
              min(r["S3"] for r in sub),max(r["S3"] for r in sub),len(sub)))
say("== 1.5 S4 (TV vs ideal kernel) and S5 (mass above h0=5), window (2,1) ==")
for fam in sorted(set(r["family"] for r in rows)):
    sub=[r for r in rows if r["family"]==fam and (r["M"],r["E"])==(2,1) and "S4" in r and r["h"]==150]
    if sub:
        say("   %-10s h=150  S4 in [%.4f,%.4f]  S5 in [%.6f,%.6f]"
            %(fam,min(r["S4"] for r in sub),max(r["S4"] for r in sub),
              min(r["S5"] for r in sub),max(r["S5"] for r in sub)))
json.dump(dict(n_rows=len(rows),nonconst=len(nonconst),tied=len(tied)),
          open(ROOT+"/shape/out/stage1.json","w"),indent=1)
open(ROOT+"/shape/out/stage1.log","w").write("\n".join(LOG))
