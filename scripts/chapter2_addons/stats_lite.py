import math
import numpy as np

# ---- regularized incomplete gamma (Numerical Recipes) ----
def _gser(a, x, itmax=500, eps=3e-14):
    if x <= 0: return 0.0
    ap = a; s = 1.0/a; d = s
    for _ in range(itmax):
        ap += 1.0; d *= x/ap; s += d
        if abs(d) < abs(s)*eps: break
    return s*math.exp(-x + a*math.log(x) - math.lgamma(a))

def _gcf(a, x, itmax=500, eps=3e-14):
    FPMIN = 1e-300
    b = x + 1.0 - a; c = 1.0/FPMIN; d = 1.0/b; h = d
    for i in range(1, itmax+1):
        an = -i*(i-a); b += 2.0
        d = an*d + b
        if abs(d) < FPMIN: d = FPMIN
        c = b + an/c
        if abs(c) < FPMIN: c = FPMIN
        d = 1.0/d; de = d*c; h *= de
        if abs(de-1.0) < eps: break
    return math.exp(-x + a*math.log(x) - math.lgamma(a))*h

def gammaq(a, x):
    """upper regularized incomplete gamma Q(a,x)"""
    if x < 0 or a <= 0: return float('nan')
    if x == 0: return 1.0
    if x < a + 1.0: return 1.0 - _gser(a, x)
    return _gcf(a, x)

def chi2_sf(x, df):
    if x <= 0: return 1.0
    return max(min(gammaq(df/2.0, x/2.0), 1.0), 0.0)

# ---- Kruskal-Wallis (ties-corrected) ----
def kruskal(groups):
    groups = [np.asarray(g, dtype=float) for g in groups if len(g) > 0]
    k = len(groups)
    if k < 2: return float('nan'), float('nan'), 0, float('nan')
    allv = np.concatenate(groups); n = allv.size
    order = np.argsort(allv, kind='mergesort')
    ranks = np.empty(n, float); sv = allv[order]
    i = 0; tie_sum = 0.0
    while i < n:
        j = i
        while j+1 < n and sv[j+1] == sv[i]: j += 1
        r = (i + j)/2.0 + 1.0
        ranks[order[i:j+1]] = r
        t = j - i + 1
        if t > 1: tie_sum += t**3 - t
        i = j + 1
    idx = 0; H = 0.0
    for g in groups:
        m = g.size; rs = ranks[idx:idx+m].sum(); idx += m
        H += rs*rs/m
    H = 12.0/(n*(n+1))*H - 3.0*(n+1)
    corr = 1.0 - tie_sum/(n**3 - n) if n > 1 else 1.0
    if corr > 0: H = H/corr
    p = chi2_sf(H, k-1)
    eps2 = (H - k + 1)/(n - k) if n > k else float('nan')   # epsilon-squared effect size
    return H, p, k, max(eps2, 0.0)

# ---- Benjamini-Hochberg ----
def bh(pvals):
    p = np.asarray(pvals, float); n = p.size
    ok = ~np.isnan(p); out = np.full(n, np.nan)
    pp = p[ok]; m = pp.size
    if m == 0: return out
    o = np.argsort(pp); ranked = pp[o]
    adj = ranked*m/np.arange(1, m+1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    res = np.empty(m); res[o] = np.minimum(adj, 1.0)
    out[ok] = res
    return out

# ---- Fisher right-tail (hypergeometric) ----
def _lchoose(n, k):
    if k < 0 or k > n: return -math.inf
    return math.lgamma(n+1) - math.lgamma(k+1) - math.lgamma(n-k+1)

def hyper_sf(k, N, K, n):
    """P(X >= k) for X~Hypergeom(N pop, K successes, n draws)"""
    lo, hi = max(0, n+K-N), min(K, n)
    if k <= lo: return 1.0
    if k > hi: return 0.0
    d = _lchoose(N, n)
    tot = 0.0
    for i in range(int(k), int(hi)+1):
        tot += math.exp(_lchoose(K, i) + _lchoose(N-K, n-i) - d)
    return min(max(tot, 0.0), 1.0)

def overlap_stats(A, B, N):
    A = set(A); B = set(B)
    a, b = len(A), len(B); ov = len(A & B)
    exp = a*b/N if N else float('nan')
    fold = ov/exp if exp > 0 else float('nan')
    p = hyper_sf(ov, N, a, b) if a and b else float('nan')
    jac = ov/len(A | B) if (A or B) else float('nan')
    return dict(n_A=a, n_B=b, n_shared=ov, expected=exp, fold=fold, p=p, jaccard=jac)
