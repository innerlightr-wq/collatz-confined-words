"""slopes.py -- irrational/rational slopes in (0,1) defined by continued fractions,
with EXACT floor(n*theta).

Exactness argument (PROVED, asserted at construction):
theta lies strictly between any two successive convergents p_k/q_k and p_{k+1}/q_{k+1}.
Hence n*theta lies strictly between n*p_k/q_k and n*p_{k+1}/q_{k+1}; if those two numbers
have the same floor, every number between them has that floor, so floor(n*theta) equals it.
`Slope.verify(N)` checks this for every 1 <= n <= N and refuses to proceed otherwise.
"""
from fractions import Fraction

def convergents(cf):
    """cf = [a1,a2,...] for theta = [0;a1,a2,...]; yields (p,q) with p/q -> theta."""
    pm1, qm1, pm2, qm2 = 0, 1, 1, 0      # p_{-1}/q_{-1}=0/1 (a0=0), p_{-2}/q_{-2}=1/0
    out = []
    for a in cf:
        p, q = a * pm1 + pm2, a * qm1 + qm2
        out.append((p, q))
        pm2, qm2, pm1, qm1 = pm1, qm1, p, q
    return out

class Slope:
    def __init__(self, name, family, cf=None, rational=None, note=""):
        self.name, self.family, self.note = name, family, note
        self.rational = rational
        if rational is not None:
            self.p, self.q = rational
            self.cf = None
            self.depth = None
        else:
            self.cf = list(cf)
            self.conv = convergents(self.cf)
            self.p, self.q = self.conv[-1]
            self.p2, self.q2 = self.conv[-2]
            self.depth = len(self.cf)
        self._cache = {}

    # ---- exact floor(n*theta) ----
    def G(self, n):
        v = self._cache.get(n)
        if v is None:
            v = (n * self.p) // self.q
            self._cache[n] = v
        return v

    def verify(self, N):
        """assert floor(n*p/q) agrees for the last two convergents for all n<=N."""
        if self.rational is not None:
            return True
        for n in range(0, N + 1):
            if (n * self.p) // self.q != (n * self.p2) // self.q2:
                raise AssertionError(
                    "slope %s: convergent depth %d insufficient at n=%d (q=%d)"
                    % (self.name, self.depth, n, self.q))
        return True

    def float_value(self):      # descriptive only
        return self.p / self.q

    # ---- Diophantine features (descriptive covariates) ----
    def pq_upto(self, N):
        """partial quotients a_k with q_k <= N; returns (list, depth)."""
        if self.rational is not None:
            return [], 0
        out = []
        for k, (p, q) in enumerate(self.conv):
            if q > N:
                break
            out.append(self.cf[k])
        return out, len(out)

    def conv_denoms(self, N):
        if self.rational is not None:
            return [self.q]
        return [q for (p, q) in self.conv if q <= N] or [1]

    def star_discrepancy(self, N):
        """exact star-discrepancy of {n*theta}, 1<=n<=N, returned as a float.
        Uses frac(n*theta) = ((n*p) mod q)/q, exact for the verified range."""
        xs = sorted(((n * self.p) % self.q) for n in range(1, N + 1))
        q = self.q
        best = 0
        for k, x in enumerate(xs, start=1):
            a = abs(Fraction(x, q) - Fraction(k, N))
            b = abs(Fraction(x, q) - Fraction(k - 1, N))
            m = max(a, b)
            if m > best:
                best = m
        return float(best)

# ------------------------------------------------------------------ the test set
def periodic(a, n):     return [a] * n
def alternating(a, b, n): return [a if i % 2 == 0 else b for i in range(n)]

def build(depth=40, seed=20261005, n_random=20):
    import random
    S = []
    # bounded type
    S.append(Slope("golden [0;1,1,1,...]", "bounded", periodic(1, depth)))
    S.append(Slope("sqrt2-1 [0;2,2,2,...]", "bounded", periodic(2, depth)))
    S.append(Slope("[0;3,3,3,...]", "bounded", periodic(3, depth)))
    S.append(Slope("sqrt3-1 [0;1,2,1,2,...]", "bounded", alternating(1, 2, depth)))
    # named
    S.append(Slope("beta = log2(3)-1", "named",
                   [1, 1, 2, 2, 3, 1, 5, 2, 23, 2, 2, 1, 1, 55, 1, 4, 3, 1, 1, 4, 1, 2,
                    1, 1, 2, 1, 1, 2, 1, 1, 3, 1, 2, 2, 5, 1, 1, 1, 2, 1, 3, 1, 2, 1, 1,
                    1, 1, 1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 1, 1]))
    S.append(Slope("e-2", "named",
                   [1, 2, 1, 1, 4, 1, 1, 6, 1, 1, 8, 1, 1, 10, 1, 1, 12, 1, 1, 14, 1, 1,
                    16, 1, 1, 18, 1, 1, 20, 1, 1, 22, 1, 1, 24, 1, 1, 26, 1, 1]))
    S.append(Slope("pi-3", "named",
                   [7, 15, 1, 292, 1, 1, 1, 2, 1, 3, 1, 14, 2, 1, 1, 2, 2, 2, 2, 1, 84,
                    2, 1, 1, 15, 3, 13, 1, 4, 2, 6, 6, 99, 1, 2, 2, 6, 3, 5, 1]))
    # Liouville-like: a huge partial quotient inserted at a given depth
    for A in (50, 500):
        for d in (3, 5, 7):
            cf = periodic(1, depth)
            cf[d - 1] = A
            S.append(Slope("liouville a_%d=%d" % (d, A), "liouville", cf,
                           note="all-ones CF with a_%d replaced by %d" % (d, A)))
    # rational controls near the golden slope (Fibonacci convergents)
    for p, q in ((8, 13), (21, 34), (55, 89), (144, 233)):
        S.append(Slope("rational %d/%d" % (p, q), "rational", rational=(p, q)))
    # random CFs, fixed seed
    rng = random.Random(seed)
    for i in range(n_random):
        cf = [rng.choice([1, 1, 1, 2, 2, 3, 4, 5, 7]) for _ in range(depth)]
        S.append(Slope("random#%02d" % i, "random", cf))
    return S

def split(slopes, seed=20261005):
    """50/50 EXPLORE/HOLDOUT, stratified so every family appears in both halves."""
    import random, collections
    rng = random.Random(seed + 7)
    byfam = collections.defaultdict(list)
    for s in slopes:
        byfam[s.family].append(s)
    ex, ho = [], []
    for fam in sorted(byfam):
        g = byfam[fam][:]
        rng.shuffle(g)
        half = len(g) // 2
        ex += g[:half] if half else g[:1]
        ho += g[half:] if half else g[1:]
    return ex, ho
