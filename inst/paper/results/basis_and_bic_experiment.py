"""
Two experiments.

(2) BIC-SELECTED ORDERS.  Compute the whitened cross-moment array once at
    (pmax, qmax); every P_n(p,q) is then a nested partial sum of its squares,
    so the whole (p,q) grid costs nothing extra.  Select
        (p^,q^) = argmax { P_n(p,q) - p q log n }
    and report P_n(p^,q^).  Question asked: is the null anywhere near chi^2_1,
    and does the selected test gain power over rule (12)?

(3) HAAR BASIS.  A Haar system at two levels gives exactly three functions
    per margin, matching p = q = 3 from rule (12) at n = 100, so the two
    bases are compared at EQUAL degrees of freedom -- any difference is the
    basis, not the budget.  A Haar basis is also, in effect, the binary
    expansion that BET uses, so this is the binary-versus-smooth comparison
    inside our own framework.  Also tested: the union of the two spans, and
    the maximum of the two standardised statistics.

Everything is a function of the two rank vectors, so every statistic here is
exactly distribution-free under independence and is calibrated by simulation
rather than by asymptotics.
"""

import numpy as np
from numpy.random import default_rng

N = 100
PMAX = 6
R0 = 20000          # null replications, for exact critical values
R1 = 2000           # replications per (scheme, noise) cell
ALPHA = 0.05
NOISE = np.round(np.arange(0.2, 2.01, 0.2), 2)


# ----------------------------------------------------------------- bases ---
def gram_schmidt(V, n):
    """Orthonormalise columns against the uniform measure on the grid, twice."""
    Q = V.astype(float).copy()
    for _ in range(2):
        for j in range(Q.shape[1]):
            for i in range(j):
                Q[:, j] -= (Q[:, i] @ Q[:, j]) / n * Q[:, i]
            nrm = np.sqrt(Q[:, j] @ Q[:, j] / n)
            if nrm < 1e-9:
                Q[:, j] = 0.0
            else:
                Q[:, j] /= nrm
    return Q


def poly_basis(n, p):
    u = (np.arange(1, n + 1)) / (n + 1.0)
    z = 2 * u - 1
    V = np.empty((n, p + 1))
    V[:, 0] = 1.0
    V[:, 1] = z
    for k in range(2, p + 1):
        V[:, k] = ((2 * k - 1) * z * V[:, k - 1] - (k - 1) * V[:, k - 2]) / k
    return gram_schmidt(V, n)[:, 1:]


def haar_basis(n, levels=2):
    """Constant dropped; 2^levels - 1 functions, i.e. 3 at two levels."""
    u = (np.arange(1, n + 1)) / (n + 1.0)
    cols = [np.ones(n)]
    for L in range(levels):
        m = 2 ** L
        for b in range(m):
            lo, hi = b / m, (b + 1) / m
            mid = (lo + hi) / 2
            h = np.zeros(n)
            h[(u > lo) & (u <= mid)] = -1.0
            h[(u > mid) & (u <= hi)] = 1.0
            cols.append(h)
    return gram_schmidt(np.column_stack(cols), n)[:, 1:]


def union_basis(n, p_poly=3, levels=2):
    u = np.column_stack([np.ones(n), poly_basis(n, p_poly), haar_basis(n, levels)])
    Q = gram_schmidt(u, n)[:, 1:]
    keep = np.abs(Q).sum(0) > 1e-8
    return Q[:, keep]


# ------------------------------------------------------------- statistic ---
def cross(G, rx, ry, n):
    """Whitened cross-moment array, shape (R, p, q).  M_ij ~ N(0,1) under H0."""
    U = G[rx]                      # (R, n, p)
    V = G[ry]
    return np.einsum('rni,rnj->rij', U, V) / np.sqrt(n)


def ranks(X):
    return np.argsort(np.argsort(X, axis=1, kind="stable"), axis=1)


def stats_from_M(M, n, pmax):
    """All the statistics, from one cross-moment array."""
    S = np.cumsum(np.cumsum(M ** 2, axis=1), axis=2)      # S[:,p-1,q-1] = P(p,q)
    pq = np.outer(np.arange(1, pmax + 1), np.arange(1, pmax + 1))
    bic = S - pq[None] * np.log(n)
    flat = bic.reshape(len(M), -1).argmax(1)
    sel = S.reshape(len(M), -1)[np.arange(len(M)), flat]
    sel_p = flat // pmax + 1
    sel_q = flat % pmax + 1
    return S, sel, sel_p, sel_q


# ------------------------------------------------------------- schemes -----
def make(scheme, n, nu, rng, R):
    U = rng.random
    if scheme == "linear":
        x = U((R, n)); y = x + nu * rng.standard_normal((R, n))
    elif scheme == "quadratic":
        x = U((R, n)); y = 4 * (x - .5) ** 2 + nu * rng.standard_normal((R, n))
    elif scheme == "cubic":
        x = U((R, n)); z = x - 1 / 3
        y = 128 * z ** 3 - 48 * z ** 2 - 12 * z + nu * 3 * rng.standard_normal((R, n))
    elif scheme == "sine":
        x = U((R, n)); y = np.sin(4 * np.pi * x) + nu * 2 * rng.standard_normal((R, n))
    elif scheme == "sine 8pi":
        x = U((R, n)); y = np.sin(8 * np.pi * x) + nu * 2 * rng.standard_normal((R, n))
    elif scheme == "fourth root":
        x = U((R, n)); y = x ** .25 + nu * .5 * rng.standard_normal((R, n))
    elif scheme == "circle":
        x = U((R, n)); r = rng.integers(0, 2, (R, n))
        y = (2 * r - 1) * np.sqrt(np.maximum(0, 1 - (2 * x - 1) ** 2)) \
            + nu * .5 * rng.standard_normal((R, n))
    elif scheme == "two curves":
        x = U((R, n)); r = rng.integers(0, 2, (R, n))
        y = 2 * r * x + .5 * (1 - r) * np.sqrt(x) + nu * rng.standard_normal((R, n))
    elif scheme == "X":
        x = U((R, n)); r = rng.integers(0, 2, (R, n))
        y = r * x + (1 - r) * (1 - x) + nu * .2 * rng.standard_normal((R, n))
    elif scheme == "diamond":
        x = U((R, n)); w = U((R, n))
        lo = np.where(x < .5, .5 - x, x - .5)
        hi = np.where(x < .5, .5 + x, 1.5 - x)
        y = lo + w * (hi - lo) + nu * .1 * rng.standard_normal((R, n))
    elif scheme == "heteroscedastic":
        x = np.abs(rng.standard_normal((R, n)))
        y = rng.standard_normal((R, n)) * x + nu * rng.standard_normal((R, n))
    elif scheme == "lattice":
        i = rng.integers(1, 5, (R, n))
        inner = (i == 1) | (i == 4)
        j = np.where(inner, rng.choice([2, 3], (R, n)), rng.choice([1, 4], (R, n)))
        x = (i - U((R, n))) / 4
        y = (j - U((R, n))) / 4 + nu * .15 * rng.standard_normal((R, n))
    elif scheme == "Clayton":
        th = 1.0
        u = U((R, n)); w = U((R, n))
        v = ((w ** (-th / (1 + th)) - 1) * u ** (-th) + 1) ** (-1 / th)
        x = u; y = v + nu * .3 * rng.standard_normal((R, n))
    else:
        raise ValueError(scheme)
    return x, y


SCHEMES = ["linear", "quadratic", "cubic", "sine", "fourth root", "circle",
           "two curves", "X", "diamond", "heteroscedastic",
           "lattice", "sine 8pi", "Clayton"]

G_POLY = poly_basis(N, PMAX)
G_HAAR = haar_basis(N, 2)
G_UNION = union_basis(N, 3, 2)
print(f"n = {N}; basis dimensions: poly {G_POLY.shape[1]}, "
      f"haar {G_HAAR.shape[1]}, union {G_UNION.shape[1]}")
RULE = 3        # rule (12) at n = 100


def all_stats(x, y, n):
    rx, ry = ranks(x), ranks(y)
    Mp = cross(G_POLY, rx, ry, n)
    S, sel, sp, sq = stats_from_M(Mp, n, PMAX)
    Mh = cross(G_HAAR, rx, ry, n)
    Mu = cross(G_UNION, rx, ry, n)
    poly_rule = S[:, RULE - 1, RULE - 1]
    haar = (Mh ** 2).sum((1, 2))
    union = (Mu ** 2).sum((1, 2))
    poly6 = S[:, PMAX - 1, PMAX - 1]
    return dict(poly_rule=poly_rule, poly_bic=sel, poly6=poly6,
                haar=haar, union=union), (sp, sq)


# ------------------------------------------------------------- the null ----
rng = default_rng(20261011)
xn = rng.random((R0, N)); yn = rng.random((R0, N))
null, (sp, sq) = all_stats(xn, yn, N)
crit = {k: np.quantile(v, 1 - ALPHA) for k, v in null.items()}

from scipy.stats import chi2
print("\n--- (2) the BIC-selected statistic under the null -------------------")
print(f"selected (p,q) = (1,1) in {np.mean((sp == 1) & (sq == 1))*100:.1f}% of "
      f"{R0} null samples;  mean pq = {np.mean(sp*sq):.2f}")
q = [0.5, 0.9, 0.95, 0.99]
print(f"{'quantile':>10}{'simulated':>12}{'chi2_1':>10}")
for a in q:
    print(f"{a:>10.2f}{np.quantile(null['poly_bic'], a):>12.3f}{chi2.ppf(a,1):>10.3f}")
sz = np.mean(null['poly_bic'] > chi2.ppf(1 - ALPHA, 1))
print(f"size of the chi^2_1 test at nominal 5%: {sz:.4f}   "
      f"-> {'usable' if abs(sz-.05)<.01 else 'NOT usable, needs simulated critical values'}")
print(f"exact 95% point, simulated: {crit['poly_bic']:.3f}  "
      f"(chi2_1 gives {chi2.ppf(.95,1):.3f})")

# ------------------------------------------------------------- power -------
print("\n--- power, size-adjusted with the exact simulated critical values ---")
keys = ["poly_rule", "poly_bic", "poly6", "haar", "union"]
res = {k: {} for k in keys}
for s in SCHEMES:
    acc = {k: [] for k in keys}
    for nu in NOISE:
        x, y = make(s, N, nu, rng, R1)
        st, _ = all_stats(x, y, N)
        for k in keys:
            acc[k].append(np.mean(st[k] > crit[k]))
    for k in keys:
        res[k][s] = float(np.mean(acc[k]))

hdr = f"{'scheme':<17}" + "".join(f"{k:>11}" for k in keys)
print(hdr); print("-" * len(hdr))
for s in SCHEMES:
    best = max(res[k][s] for k in keys)
    print(f"{s:<17}" + "".join(
        (f"{res[k][s]:>10.3f}" + ("*" if abs(res[k][s] - best) < 1e-9 else " "))
        for k in keys))
print("-" * len(hdr))
print(f"{'MEAN':<17}" + "".join(f"{np.mean(list(res[k].values())):>10.3f} " for k in keys))
PY = None
