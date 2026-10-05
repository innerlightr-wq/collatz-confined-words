"""shape_core.py -- staircase chain counts and Haar-weighted shape statistics,
for an arbitrary slope.  Exact integer arithmetic for every staircase and count.

For a slope theta (exact G(n) = floor(n*theta)) and a "window" (M,E):
    r(h)    = G(M+h) - G(M) - E + 1            (= ell_0, the number of seeds)
    ell_j   = min(ell_0, G(M+h) - G(M+j)),  1 <= j <= h
    Q_L(j)  = #{x_0 >= ... >= x_j >= 0 : x_i <= ell_i - 1}
    Pi_h(j) = Q_L(j) + [j=0]
    eta_h   = G(M+h+1) - G(M+h)
NOTE: for theta != log2(3)-1 these are the same staircase combinatorics with a different
slope; they are NOT accelerated-Collatz windows and no Collatz claim is attached to them.
"""
from fractions import Fraction
from itertools import accumulate

H0_REF = 5   # reference horizon of the existing study code: K(h) = G(M+h) - G(M+h0)

def l_profile(S, M, E, h):
    B = S.G(M + h)
    l0 = B - S.G(M) - E + 1
    if l0 < 0:
        return None
    out = [l0]
    for j in range(1, h + 1):
        out.append(max(0, min(l0, B - S.G(M + j))))
    return out

def profile(S, M, E, h):
    """Pi_h(j) for j=0..h as a list of exact integers (None if the window is closed)."""
    L = l_profile(S, M, E, h)
    if L is None:
        return None
    Pi = [L[0] + 1]
    a = [1] * L[0]
    for j in range(1, h + 1):
        if not a or L[j] == 0:
            Pi.extend([0] * (h + 1 - j))
            break
        suf = list(accumulate(reversed(a)))
        suf.reverse()
        a = suf[:L[j]]
        Pi.append(sum(a))
    return Pi[:h + 1] + [0] * max(0, h + 1 - len(Pi))

def K_of(S, M, h, h0=H0_REF):
    return S.G(M + h) - S.G(M + h0)

def shape_stats(Pi, h, K, want_tv=False):
    """Haar weights w(j) = Pi(j) 2^{-j}; all comparisons exact, summaries as floats."""
    W = [Pi[j] << (h - j) for j in range(h + 1)]          # exact integers
    tot = sum(W)
    if tot == 0:
        return None
    mx = max(W)
    modes = [j for j, w in enumerate(W) if w == mx]
    s1 = sum(j * w for j, w in enumerate(W))
    s2 = sum(j * j * w for j, w in enumerate(W))
    mean = Fraction(s1, tot)
    var = Fraction(s2 * tot - s1 * s1, tot * tot)
    out = dict(mode=modes[0], mode_tied=len(modes) > 1, n_modes=len(modes),
               S1=modes[0] - K, mean=float(mean), S2=float(mean) - K,
               var=float(var), S3=(float(var) / K if K else float('nan')),
               mass_le_h0=float(Fraction(sum(W[:H0_REF + 1]), tot)))
    out["S5"] = 1.0 - out["mass_le_h0"]
    if want_tv and K >= 1:
        # ideal kernel binom(j+K-1,K-1) 2^{-j}, total mass 2^K
        from math import comb
        jm = h
        ideal = [comb(j + K - 1, K - 1) << (jm - j) for j in range(jm + 1)]
        it = sum(ideal)
        num = sum(abs(Fraction(W[j], tot) - Fraction(ideal[j], it)) for j in range(jm + 1))
        out["S4"] = float(num) / 2.0
    return out

def row(S, M, E, h, want_tv=False):
    Pi = profile(S, M, E, h)
    if Pi is None:
        return None
    K = K_of(S, M, h)
    if K <= 0:
        return None
    st = shape_stats(Pi, h, K, want_tv)
    if st is None:
        return None
    L = l_profile(S, M, E, h)
    st.update(slope=S.name, family=S.family, M=M, E=E, h=h, K=K,
              r=L[0], eta=S.G(M + h + 1) - S.G(M + h),
              frac_h=float(Fraction((M + h) * S.p % S.q, S.q)),
              frac_M=float(Fraction(M * S.p % S.q, S.q)),
              top_incr="".join(str(L[j] - L[j + 1]) for j in range(max(0, h - 10), h)))
    return st

# ---------------- convention-free statistics -------------------------------
def shape_stats2(Pi, h, K, r, want_tv=False):
    """As shape_stats, plus statistics that do not depend on the arbitrary additive
    convention in K(h): S2r = mean - r, S3r = Var/r, and a mean-matched NB comparison."""
    from math import comb
    W = [Pi[j] << (h - j) for j in range(h + 1)]
    tot = sum(W)
    if tot == 0: return None
    mx = max(W)
    modes = [j for j, w in enumerate(W) if w == mx]
    s1 = sum(j * w for j, w in enumerate(W))
    s2 = sum(j * j * w for j, w in enumerate(W))
    mean = Fraction(s1, tot); var = Fraction(s2 * tot - s1 * s1, tot * tot)
    fm, fv = float(mean), float(var)
    out = dict(mode=modes[0], mode_tied=len(modes) > 1, n_modes=len(modes),
               mean=fm, var=fv, S1=modes[0] - K, S1r=modes[0] - r,
               S2=fm - K, S2r=fm - r, S3=(fv / K if K else float("nan")),
               S3r=(fv / r if r else float("nan")),
               mode_minus_roundmean=modes[0] - round(fm))
    if want_tv:
        if K >= 1:
            ideal = [comb(j + K - 1, K - 1) << (h - j) for j in range(h + 1)]
            it = sum(ideal)
            out["S4"] = float(sum(abs(Fraction(W[j], tot) - Fraction(ideal[j], it))
                                  for j in range(h + 1))) / 2.0
        kap = max(1, int(round(fm)))          # mean-matched negative binomial
        ideal = [comb(j + kap - 1, kap - 1) << (h - j) for j in range(h + 1)]
        it = sum(ideal)
        out["S4m"] = float(sum(abs(Fraction(W[j], tot) - Fraction(ideal[j], it))
                               for j in range(h + 1))) / 2.0
        out["kappa"] = kap
    return out

def row2(S, M, E, h, want_tv=False):
    Pi = profile(S, M, E, h)
    if Pi is None: return None
    L = l_profile(S, M, E, h); r = L[0]
    K = K_of(S, M, h)
    if K <= 0 or r <= 0: return None
    st = shape_stats2(Pi, h, K, r, want_tv)
    if st is None: return None
    incr = [L[j] - L[j + 1] for j in range(h)]          # mechanical increment word
    st.update(slope=S.name, family=S.family, M=M, E=E, h=h, K=K, r=r,
              eta=S.G(M + h + 1) - S.G(M + h),
              frac_h=float(Fraction((M + h) * S.p % S.q, S.q)),
              frac_M=float(Fraction(M * S.p % S.q, S.q)),
              plateau=sum(1 for j in range(h) if L[j] == L[0]),
              w_top="".join(map(str, incr[-8:])),
              w_mid="".join(map(str, incr[max(0, r - 8):r])))
    return st
