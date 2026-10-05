"""verify_c12_range.py -- scope of the iterated scalar reduction (Corollary C.12).

Claim under test:  fix a window class (M,E) and a reference horizon h0 with r(W,h0) >= 0.
Iterate the single-step operators of Theorem C.11 from Pi_{h0}:
    Pi^iter_{h+1} = Pi^iter_h           if eta_h = 0
    Pi^iter_{h+1} = T1(Pi^iter_h)       if eta_h = 1
(so Pi^iter_h = T1^{K(h)} Pi_{h0}).  Then
    (i)  Pi^iter_h(j) = Pi_h(j) for every 0 <= j <= h0 and every h,
    (ii) the identity FAILS above h0, and the first failing index is recorded.
Exact integer arithmetic throughout (chain counts are exact big-int sums).
"""
import os as _os, sys as _sys
ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
_sys.path.insert(0, _os.path.join(ROOT, 'code'))

import sys, json, csv
sys.path.insert(0, ROOT+"/code")
from core import G, chain_N, l_profile

HMAX = 150
JTOP = HMAX + 4

def eta(M, h):
    return G(M + h + 1) - G(M + h)

def r_raw(M, E, h):
    return G(M + h) - G(M) - E + 1

def h_open(M, E):
    h = 1
    while h < 500 and r_raw(M, E, h) < 0:
        h += 1
    return h

def profile_vec(M, E, h, jtop=JTOP):
    """true Pi_h as a zero-extended vector on 0..jtop."""
    N = chain_N(M, E, h)
    return [N.get(j, 0) for j in range(jtop + 1)]

def T1(v):
    out = [0] * len(v)
    run = 0
    for j in range(len(v)):
        run += v[j]
        out[j] = run + 1 if j == 0 else run
    return out

def run_class(M, E, h0, hmax=HMAX):
    """returns dict with agreement on j<=h0 and the first mismatch index above h0."""
    if r_raw(M, E, h0) < 0:
        return None
    it = profile_vec(M, E, h0)
    rec = dict(M=M, E=E, h0=h0, h_open=h_open(M, E),
               n_h=0, n_checks_le_h0=0, bad_le_h0=[], first_mismatch={}, )
    first_overall = None
    for h in range(h0, hmax):
        e = eta(M, h)
        it = it if e == 0 else T1(it)
        true = profile_vec(M, E, h + 1)
        rec["n_h"] += 1
        # (i) agreement on j <= h0
        for j in range(h0 + 1):
            rec["n_checks_le_h0"] += 1
            if it[j] != true[j]:
                rec["bad_le_h0"].append((h + 1, j))
        # (ii) first mismatch anywhere
        fm = next((j for j in range(len(true)) if it[j] != true[j]), None)
        if fm is not None:
            rec["first_mismatch"][h + 1] = fm
            if first_overall is None:
                first_overall = (h + 1, fm)
    rec["first_overall"] = first_overall
    fm = sorted(set(rec["first_mismatch"].values()))
    rec["first_mismatch_indices"] = fm
    rec["min_first_mismatch"] = min(fm) if fm else None
    rec["mismatch_persists_to_h"] = max(rec["first_mismatch"]) if rec["first_mismatch"] else None
    return rec

def main():
    rows = []
    classes = [(2, 1)] + [(M, E) for M in range(2, 9) for E in (1, 2, 3, 4)]
    seen = set(); clean = []
    for c in classes:
        if c not in seen:
            seen.add(c); clean.append(c)
    print("classes tested: %d  (W0 = (M,E) = (2,1))" % len(clean))
    for (M, E) in clean:
        for h0 in (5, 10, 20, 30):
            rec = run_class(M, E, h0)
            if rec is None:
                print("  skip M=%d E=%d h0=%d : r(W,h0) < 0 (window not open; h_open=%d)"
                      % (M, E, h0, h_open(M, E)))
                continue
            rows.append(rec)
    tot = sum(r["n_checks_le_h0"] for r in rows)
    bad = sum(len(r["bad_le_h0"]) for r in rows)
    print("\n(i)  agreement on j <= h0 : %d checks over %d (class,h0) runs, violations = %d"
          % (tot, len(rows), bad))
    nofail = [r for r in rows if r["min_first_mismatch"] is None]
    print("(ii) runs with NO mismatch anywhere up to h=%d : %d" % (HMAX, len(nofail)))
    print("     first mismatch index relative to h0 (min over h), by h0:")
    for h0 in (5, 10, 20, 30):
        sub = [r for r in rows if r["h0"] == h0 and r["min_first_mismatch"] is not None]
        if not sub: continue
        rel = sorted({r["min_first_mismatch"] - h0 for r in sub})
        print("       h0=%2d : n=%2d classes, first mismatch at j-h0 in %s"
              % (h0, len(sub), rel))
    print("\nW0 (M=2,E=1) detail:")
    for r in [x for x in rows if (x["M"], x["E"]) == (2, 1)]:
        print("   h0=%2d : first mismatch at j=%s (first seen at h=%s), persists to h=%s ;"
              " violations on j<=h0: %d"
              % (r["h0"], r["min_first_mismatch"], r["first_overall"][0] if r["first_overall"] else None,
                 r["mismatch_persists_to_h"], len(r["bad_le_h0"])))
    with open(ROOT+"/data/c12_range.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["M","E","h0","h_open","n_h","n_checks_le_h0","violations_le_h0",
                    "first_mismatch_index","first_mismatch_at_h","persists_to_h"])
        for r in rows:
            w.writerow([r["M"],r["E"],r["h0"],r["h_open"],r["n_h"],r["n_checks_le_h0"],
                        len(r["bad_le_h0"]),r["min_first_mismatch"],
                        r["first_overall"][0] if r["first_overall"] else "",
                        r["mismatch_persists_to_h"]])
    json.dump(rows, open(ROOT+"/out/c12_range.json","w"),
              indent=0, default=str)
    return rows, tot, bad

if __name__ == "__main__":
    main()
