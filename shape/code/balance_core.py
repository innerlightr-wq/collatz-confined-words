"""balance_core.py -- balance statistics of the Haar-weighted staircase profile.

w(j) = Pi_h(j) 2^{-j},  Pi_h(j) = Q_{L(h)}(j) + [j=0].
rho(j) = w(j+1)/w(j) = Pi(j+1) / (2 Pi(j))          <- exact Fraction
peak j* = argmax w (smallest index if tied)         <- exact integer comparison
balance margin m = min(|log rho(j*-1)|, |log rho(j*)|)   (0 iff the peak ties with a neighbour)
Every decision (argmax, tie, drop counting) is exact integer arithmetic; logs are a
descriptive transform of exact rationals, used only for reporting and correlations.

Profiles are advanced with the PROVED operators (eta=0 identity, eta=1 T1), spot-checked
against a direct chain-count recomputation.
"""
import math
from fractions import Fraction
from shape_core import profile, l_profile, K_of

def T1(P):
    out=[]; run=0
    for j,v in enumerate(P):
        run+=v; out.append(run+1 if j==0 else run)
    return out

def iter_profiles(S, M, E, hmax, check_every=0):
    """yield (h, r, L, Pi) for h=1..hmax, skipping horizons where the window is closed."""
    P=None; fails=0; checks=0
    for h in range(1, hmax+1):
        r = S.G(M+h) - S.G(M) - E + 1
        if r < 1:
            P=None; continue
        if P is None:
            P = profile(S, M, E, h)
        else:
            eta_prev = S.G(M+h) - S.G(M+h-1)
            P = (P if eta_prev==0 else T1(P)) + [0]
        if check_every and h % check_every == 0:
            D = profile(S, M, E, h); checks += 1
            if D[:len(P)] != P[:len(D)]: fails += 1
        yield h, r, l_profile(S, M, E, h), P, checks, fails

def peak_and_margin(Pi, h):
    """exact peak, tie structure, and the balance margin."""
    W = [Pi[j] << (h-j) for j in range(len(Pi))]
    mx = max(W)
    idx = [j for j,w in enumerate(W) if w == mx]
    js = idx[0]
    tied = len(idx) > 1
    nonadj = tied and any(idx[t+1]-idx[t] > 1 for t in range(len(idx)-1))
    def logratio(j):           # |log rho(j)| = |log( Pi(j+1) / (2 Pi(j)) )|
        if j < 0 or j+1 >= len(Pi): return math.inf
        a, b = Pi[j+1], 2*Pi[j]
        if a == 0: return math.inf
        if b == 0: return math.inf
        return abs(math.log(float(Fraction(a, b))))
    lo = logratio(js-1) if js >= 1 else math.inf
    hi = logratio(js)
    m = min(lo, hi)
    return dict(peak=js, tied=tied, n_modes=len(idx), nonadjacent_tie=nonadj,
                margin=m, log_lo=lo, log_hi=hi, tie_gap=(idx[1]-idx[0]) if tied else 0)

def top_irregularity(L, js, h, half=5):
    """D = number of staircase drops in j in [j*-5, j*+5]; V = variance of the gaps
    between consecutive drops in that window (0 if fewer than two drops)."""
    lo, hi = max(0, js-half), min(h-1, js+half)
    pos = [j for j in range(lo, hi+1) if j+1 < len(L) and L[j]-L[j+1] > 0]
    D = len(pos)
    if D >= 2:
        g = [pos[t+1]-pos[t] for t in range(D-1)]
        mu = sum(g)/len(g)
        V = sum((x-mu)**2 for x in g)/len(g)
    else:
        V = 0.0
    return D, V, (hi-lo+1)

def best_ideal_fit(Pi, h, mean_est, span=20):
    """k_fit = argmin_k TV(normalised w, normalised w*_k), w*_k(j) = C(j+k-1,k-1) 2^{-j}.
    k searched over [max(1, round(mean)-span), round(mean)+span] (pre-specified span)."""
    W = [Pi[j] << (h-j) for j in range(len(Pi))]
    totW = sum(W)
    if totW == 0: return None, None
    n = len(W)
    k0 = max(1, int(round(mean_est)))
    best = (None, None)
    for k in range(max(1, k0-span), k0+span+1):
        I = [math.comb(j+k-1, k-1) << (h-j) for j in range(n)]
        totI = sum(I)
        num = sum(abs(W[j]*totI - I[j]*totW) for j in range(n))
        tv = float(Fraction(num, 2*totW*totI))
        if best[1] is None or tv < best[1]: best = (k, tv)
    return best

def mean_of(Pi, h):
    W = [Pi[j] << (h-j) for j in range(len(Pi))]
    tot = sum(W)
    if tot == 0: return 0.0
    return float(Fraction(sum(j*w for j,w in enumerate(W)), tot))
