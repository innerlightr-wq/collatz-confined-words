"""
core.py -- exact machinery for the accelerated-Collatz confinement-window
finite-horizon decision tree, generalized to arbitrary windows W=(M,c,prefix,V_min).

EVERY confinement / prune / return decision is an exact integer inequality.
Floats appear only in explicitly-labelled descriptive features (gamma_W etc).

Conventions follow the reference implementation of the companion note
("A Sturmian Complexity Clock in Accelerated Collatz Words", App. A.2/A.3):

  state = (j post-window digits finalized, confirmed sum S_fin, current digit
           lower bound cur)
  PRUNE  : 2^(S_fin+cur+h-j-1) >  2^c * 3^(M+h)      [best-case continuation fails]
  RETURN : 2^(S_fin+cur)       <= 2^c * 3^(M+j+1)    [re-enters corridor]
  else   : CONTINUE to level j+1 with S_fin+cur (and cur is also incremented,
           i.e. every cur in 1..d*-1 is classified independently)
  terminal seeding: V = V_min .. V*(h)-1,  V*(h) = first V failing the best-case
           check with sum P+V+h at index M+h.  N_h(0) includes the V=V*(h) prune leaf.
  N_h(j)  = number of pruned leaves settled after exactly j post-window digits
          = (number of alive states at level j, with multiplicity) [+1 if j==0]
"""
from fractions import Fraction

# ---------------------------------------------------------------- exact floors
_P3 = [1]
def pow3(n):
    while len(_P3) <= n:
        _P3.append(_P3[-1] * 3)
    return _P3[n]

def F(n):
    """floor(n*log2 3), exact: the unique k with 2^k <= 3^n < 2^(k+1)."""
    return pow3(n).bit_length() - 1

def G(n):
    """floor(n*beta), beta = log2(3/2), exact.  (floor(n*alpha) = n + floor(n*beta))"""
    return F(n) - n

def confined(S, n, c):
    """R_n <= c  <=>  2^S <= 2^c * 3^n   (exact integer test)."""
    return (1 << S) <= (1 << c) * pow3(n)

# ---------------------------------------------------------------- the window
class Window:
    """W = (M, c, interior prefix (d_0..d_{M-2}), V_min)."""
    def __init__(self, M, c, prefix, Vmin, name=None):
        assert M >= 1 and c >= 0
        prefix = tuple(prefix)
        assert len(prefix) == M - 1, "prefix must have M-1 digits (indices 0..M-2)"
        assert all(d >= 1 for d in prefix)
        self.M, self.c, self.prefix, self.Vmin = M, c, prefix, Vmin
        self.P = sum(prefix)
        self.name = name or self.tag()
        # interior prefix must be confined for j <= M-1
        S = 0
        for j, d in enumerate(prefix):
            S += d
            assert confined(S, j + 1, c), "interior prefix not c-confined at index %d" % (j + 1)
        # terminal digit must exit the corridor at index M
        assert not confined(self.P + Vmin, M, c), "V_min does not exit the corridor at index M"
        # exit margin E >= 1  (E = 1 means V_min is the *minimal* exiting digit)
        self.E = self.P + Vmin - c - F(M)
        assert self.E >= 1

    def tag(self):
        return "M%d_c%d_p%s_V%d" % (self.M, self.c, "-".join(map(str, self.prefix)), self.Vmin)

    # descriptive (float) features -- never used in a decision
    @property
    def gamma(self):
        import math
        return (self.M * math.log2(3)) % 1.0

    def CW(self, h):
        """C_W(h) = floor(c + M*alpha + h*beta), by exact integer search."""
        return self.c - h + F(self.M + h)

    def rho(self, j):
        """rho*(j) = floor(c + (M+j+1)*alpha), exact."""
        return self.c + F(self.M + j + 1)

    def eps(self, h):
        return self.CW(h + 1) - self.CW(h)

# ------------------------------------------------- 1. literal brute-force tree
def brute_N(W, h, exhaust_at_level_h=True):
    """Fully literal recursive tree.  Only raw integer confinement tests are used:
    no closed-form threshold anywhere.  O(#leaves) -- use for small h only."""
    M, c, P, Vmin = W.M, W.c, W.P, W.Vmin
    N = {}
    def bump(j):
        N[j] = N.get(j, 0) + 1

    def best_case_ok(S_fin, cur, j):
        # current digit = cur, all later digits = 1, evaluated at index M+h
        return confined(S_fin + cur + (h - j - 1), M + h, c)

    def rec(j, S_fin):
        bump_done = False
        cur = 1
        while True:
            if not best_case_ok(S_fin, cur, j):
                bump(j)                      # the single prune leaf of this state
                return
            if confined(S_fin + cur, M + j + 1, c):
                pass                         # RETURN leaf (not counted in N)
            elif j + 1 == h:
                if exhaust_at_level_h:
                    bump(h)                  # horizon exhaustion
                else:
                    bump(j)
            else:
                rec(j + 1, S_fin + cur)
            cur += 1

    # terminal digit seeding
    V = Vmin
    while confined(P + V + h, M + h, c):     # best-case check for terminal digit V
        assert not confined(P + V, M, c), "seed V does not exit the corridor"
        rec(0, P + V)
        V += 1
    bump(0)                                  # the V = V*(h) prune leaf
    return {j: v for j, v in sorted(N.items())}

# ---------------------------------------- 2. DP over {S_fin: multiplicity}
def dp_N(W, h, literal_thresholds=True):
    """O(h^2)-ish DP, mirroring App. A.3 but for a general window.
    literal_thresholds=True finds d* by the same while-loop as the brute tree
    (no closed form); False uses the derived closed form d* = C_W(h)+2+j-S_fin."""
    M, c, P, Vmin = W.M, W.c, W.P, W.Vmin
    N = {}
    alive = {}
    V = Vmin
    while confined(P + V + h, M + h, c):
        alive[P + V] = alive.get(P + V, 0) + 1
        V += 1
    N[0] = 1                                 # terminal-digit prune leaf
    CW = W.CW(h)
    j = 0
    while alive:
        nxt = {}
        pruned = 0
        for S_fin, cnt in alive.items():
            if literal_thresholds:
                d = 1
                while confined(S_fin + d + (h - j - 1), M + h, c):
                    d += 1
                dstar = d
            else:
                dstar = CW + 2 + j - S_fin
            pruned += cnt
            for d in range(1, dstar):
                if confined(S_fin + d, M + j + 1, c):
                    continue                 # RETURN
                if j + 1 == h:
                    N[h] = N.get(h, 0) + cnt # exhaustion
                else:
                    nxt[S_fin + d] = nxt.get(S_fin + d, 0) + cnt
        N[j] = N.get(j, 0) + pruned
        if j + 1 == h:
            break
        alive = nxt
        j += 1
    return {j: v for j, v in sorted(N.items()) if v}

# ------------------------------------- 3. reduced lattice-chain model (M,E)
def l_profile(M, E, h):
    """Exact bound profile of the reduced model.
         l_0 = G(M+h) - G(M) - E + 1
         l_j = min(l_0, G(M+h) - G(M+j))      (j >= 1)
    Levels j = 0..h.  (Derived in REPORT.md Sec.2; verified against brute_N.)"""
    B = G(M + h)
    l0 = B - G(M) - E + 1
    out = [max(0, l0)]
    for j in range(1, h + 1):
        out.append(max(0, min(l0, B - G(M + j))))
    return out

def chain_N(M, E, h):
    """N_h(j) via the reduced model:
         N_h(j) = #{ x_0 >= x_1 >= ... >= x_j >= 0 : x_i <= l_i - 1 }  (+1 if j==0)
    computed by truncated suffix sums.  Exact big-int arithmetic."""
    L = l_profile(M, E, h)
    N = {0: L[0] + 1}
    a = [1] * L[0]
    for j in range(1, h + 1):
        if L[j] == 0 or not a:
            break
        s = 0
        suf = [0] * len(a)
        for x in range(len(a) - 1, -1, -1):
            s += a[x]
            suf[x] = s
        a = suf[:L[j]]
        tot = sum(a)
        if tot == 0:
            break
        N[j] = tot
    return N

def chain_N_window(W, h):
    return chain_N(W.M, W.E, h)

# ---------------------------------------------------------------- operators
def T1(N, jmax):
    """T1(N)(0) = N(0)+1 ; T1(N)(j) = sum_{i<=j} N(i).  True-index dict."""
    out = {}
    run = 0
    for j in range(jmax + 1):
        run += N.get(j, 0)
        out[j] = run + 1 if j == 0 else run
    return out

def profile_list(N, jmax):
    return [N.get(j, 0) for j in range(jmax + 1)]

# ------------------------------------- exact Diophantine helpers for beta
def _cmp_dbeta_n(d, n):
    """sign of d*beta - n for integers d>0, exact:  d*beta > n <=> 3^d > 2^(n+d)."""
    a, b = pow3(d), (1 << (n + d)) if n + d >= 0 else None
    if b is None:
        return 1
    return (a > b) - (a < b)

def frac_dist(q):
    """two-sided distance ||q*beta|| as an exact Fraction-free comparable key:
    returns (which, q, p) where the distance is  q*beta-p  (which=+1, = frac)
    or p-q*beta (which=-1, = 1-frac).  Compare two such keys with dist_cmp."""
    p = G(q)
    # frac = q*beta - p in (0,1); compare frac vs 1/2  <=> 2*q*beta vs 2*p+1
    lo = _cmp_dbeta_n(2 * q, 2 * p + 1)
    return (1, q, p) if lo < 0 else (-1, q, p + 1)

def dist_cmp(k1, k2):
    """compare the two distances represented by keys k1,k2 (exact)."""
    s1, q1, p1 = k1
    s2, q2, p2 = k2
    # distance_i = s_i*(q_i*beta - p_i);  compare d1 < d2
    # d1 - d2 = (s1*q1 - s2*q2)*beta - (s1*p1 - s2*p2)
    d = s1 * q1 - s2 * q2
    n = s1 * p1 - s2 * p2
    if d == 0:
        return (0 > n) - (0 < n) if n else 0
    if d > 0:
        return _cmp_dbeta_n(d, n)
    return -_cmp_dbeta_n(-d, -n)

def convergent_denominators(qmax):
    """denominators q of the best two-sided rational approximations to beta
    (= continued-fraction convergents), by exact integer comparison only."""
    out = []
    best = None
    for q in range(1, qmax + 1):
        k = frac_dist(q)
        if best is None or dist_cmp(k, best) < 0:
            best = k
            out.append(q)
    return out

def semiconvergent_denominators(qmax):
    """one-sided record denominators (convergents + intermediate fractions)."""
    out = []
    bl = br = None
    for q in range(1, qmax + 1):
        k = frac_dist(q)
        tgt = 'l' if k[0] == 1 else 'r'
        cur = bl if tgt == 'l' else br
        if cur is None or dist_cmp(k, cur) < 0:
            if tgt == 'l':
                bl = k
            else:
                br = k
            out.append(q)
    return sorted(set(out))
