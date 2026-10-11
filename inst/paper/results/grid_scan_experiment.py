"""
Scanning the (p,q) grid instead of fixing or selecting one order.

Everything here is a function of the two rank vectors, so every rule is
exactly distribution-free and every critical value below is simulated, not
asymptotic.  The grid costs nothing to compute: P_n(p,q) for every (p,q) is a
nested partial sum of the squared whitened cross-moments, so one array at
(pmax, qmax) yields the whole surface.

Rules compared:
  rule        p = q = 3, rule (12) at n = 100
  bic         argmax over the grid of P_n(p,q) - pq log n
  maxz        max over the grid of (P_n(p,q) - mean) / sd, nulls simulated
  minp        min over the grid of the EXACT null p-value of P_n(p,q)
  minp_diag   the same, restricted to the diagonal p = q
  seq_geom    fixed sequence down the diagonal, alpha_k proportional to 2^-k,
              i.e. a weighted Bonferroni that spends most of the level at the
              smallest order -- the "sequential design" with its correction
  oracle      best single (p,q) for that scheme, an upper bound no test attains

minp is the principled scan: converting to exact p-values first puts the
cells on a common footing, which max-of-standardised does not, since the
cells have different degrees of freedom and different skewness.
"""

import numpy as np
from numpy.random import default_rng
import exp as E

N, PMAX, R0, R1, AL = E.N, 6, 20000, 2000, 0.05
GRID = [(p, q) for p in range(1, PMAX + 1) for q in range(1, PMAX + 1)]
DIAG = [(k, k) for k in range(1, PMAX + 1)]


def surface(x, y):
    """P_n(p,q) for every (p,q); shape (R, PMAX, PMAX)."""
    M = E.cross(E.G_POLY, E.ranks(x), E.ranks(y), N)
    return np.cumsum(np.cumsum(M ** 2, 1), 2)


rng = default_rng(20261011)
# two independent null batches: one to build the per-cell null laws, one to
# calibrate the scan statistics on top of them.  Using one batch for both
# would make the scan look better than it is.
Sa = surface(rng.random((R0, N)), rng.random((R0, N)))
Sb = surface(rng.random((R0, N)), rng.random((R0, N)))

srt = {pq: np.sort(Sa[:, pq[0] - 1, pq[1] - 1]) for pq in GRID}
mu = {pq: Sa[:, pq[0] - 1, pq[1] - 1].mean() for pq in GRID}
sd = {pq: Sa[:, pq[0] - 1, pq[1] - 1].std() for pq in GRID}


def pval(S, pq):
    v = S[:, pq[0] - 1, pq[1] - 1]
    return 1.0 - np.searchsorted(srt[pq], v, side="left") / len(srt[pq])


def rules(S):
    out = {}
    out["rule"] = S[:, 2, 2]
    pq = np.outer(np.arange(1, PMAX + 1), np.arange(1, PMAX + 1))
    bic = S - pq[None] * np.log(N)
    f = bic.reshape(len(S), -1).argmax(1)
    out["bic"] = S.reshape(len(S), -1)[np.arange(len(S)), f]
    out["maxz"] = np.max([(S[:, p - 1, q - 1] - mu[(p, q)]) / sd[(p, q)]
                          for p, q in GRID], axis=0)
    out["minp"] = -np.min([pval(S, pq) for pq in GRID], axis=0)
    out["minp_diag"] = -np.min([pval(S, pq) for pq in DIAG], axis=0)
    # weighted Bonferroni down the diagonal: alpha_k proportional to 2^-k
    w = np.array([2.0 ** -(k + 1) for k in range(len(DIAG))]); w /= w.sum()
    out["seq_geom"] = np.max([-pval(S, pq) / wk for pq, wk in zip(DIAG, w)], axis=0)
    return out


null = rules(Sb)
crit = {k: np.quantile(v, 1 - AL) for k, v in null.items()}

print(f"n = {N}, grid {PMAX}x{PMAX} = {len(GRID)} cells, "
      f"{R0} null draws for the cell laws and {R0} more for calibration")
print("\nexact 5% critical values (simulated):")
for k in crit:
    print(f"  {k:<11}{crit[k]:>12.4f}")

keys = ["rule", "bic", "maxz", "minp", "minp_diag", "seq_geom"]
res = {k: {} for k in keys + ["oracle"]}
for s in E.SCHEMES:
    acc = {k: [] for k in keys}
    orac = []
    for nu in E.NOISE:
        x, y = E.make(s, N, nu, rng, R1)
        S = surface(x, y)
        r = rules(S)
        for k in keys:
            acc[k].append(np.mean(r[k] > crit[k]))
        # oracle: best single cell, calibrated at its own exact 5% point
        orac.append(max(np.mean(S[:, p - 1, q - 1] > np.quantile(srt[(p, q)], 1 - AL))
                        for p, q in GRID))
    for k in keys:
        res[k][s] = float(np.mean(acc[k]))
    res["oracle"][s] = float(np.mean(orac))

cols = keys + ["oracle"]
hdr = f"{'scheme':<17}" + "".join(f"{k:>11}" for k in cols)
print("\n" + hdr); print("-" * len(hdr))
for s in E.SCHEMES:
    best = max(res[k][s] for k in keys)          # oracle excluded from the star
    print(f"{s:<17}" + "".join(
        (f"{res[k][s]:>10.3f}" + ("*" if k in keys and abs(res[k][s] - best) < 1e-9 else " "))
        for k in cols))
print("-" * len(hdr))
print(f"{'MEAN':<17}" + "".join(f"{np.mean(list(res[k].values())):>10.3f} " for k in cols))
print("\nloss against the oracle (mean over schemes):")
o = np.mean(list(res['oracle'].values()))
for k in keys:
    print(f"  {k:<11}{o - np.mean(list(res[k].values())):>8.3f}")
