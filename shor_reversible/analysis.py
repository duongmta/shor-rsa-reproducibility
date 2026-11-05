"""Exact finite-QPE distribution, factor-recovery rule and sampling statistics."""
from fractions import Fraction
from math import gcd, sqrt
import numpy as np


def order(a, N):
    k, x = 1, a % N
    while x != 1:
        x = x * a % N
        k += 1
    return k


def identity_blocks(N, a, t):
    return sum(1 for j in range(t) if pow(a, 2 ** j, N) == 1)


def qpe_distribution(r, Q):
    """P(y) = (1/r) sum_s |(1/Q) sum_k exp(2 pi i k (s/r - y/Q))|^2 ."""
    k = np.arange(Q)
    P = np.zeros(Q)
    for s in range(r):
        amp = np.fft.fft(np.exp(2j * np.pi * k * s / r)) / Q
        P += np.abs(amp) ** 2
    return P / r


def convergents(fr):
    p, q = fr.numerator, fr.denominator
    a = []
    while q:
        a.append(p // q)
        p, q = q, p - (p // q) * q
    h0, h1, k0, k1 = 1, a[0], 0, 1
    out = [Fraction(h1, k1)]
    for ai in a[1:]:
        h0, h1 = h1, ai * h1 + h0
        k0, k1 = k1, ai * k1 + k0
        out.append(Fraction(h1, k1))
    return out


def recovers_factor(y, Q, a, N):
    """Paper's rule: some convergent denominator d<N is even, a^d = 1 mod N and
    both gcd(a^(d/2) -+ 1, N) are non-trivial."""
    for cv in convergents(Fraction(y, Q)):
        d = cv.denominator
        if d < N and d % 2 == 0 and pow(a, d, N) == 1:
            x = pow(a, d // 2, N)
            g1, g2 = gcd(x - 1, N), gcd(x + 1, N)
            if 1 < g1 < N and 1 < g2 < N:
                return True
    return False


def success_probability(N, a, t):
    r = order(a, N)
    Q = 2 ** t
    P = qpe_distribution(r, Q)
    acc = np.array([P[y] > 1e-15 and recovers_factor(y, Q, a, N) for y in range(Q)])
    return float(P[acc].sum()), P, acc


def wilson(k, n, z=1.96):
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (c - h, c + h)


def compare_counts(counts_y, P, shots, min_expected=5.0):
    """counts_y: dict y->count.  Returns TVD, off-support counts, pooled chi-square."""
    Q = len(P)
    obs = np.zeros(Q)
    for y, c in counts_y.items():
        obs[y] += c
    tvd = 0.5 * np.abs(obs / shots - P).sum()
    off_support = int(obs[P < 1e-12].sum())
    exp = P * shots
    keep = exp >= min_expected
    o = list(obs[keep]); e = list(exp[keep])
    if (~keep).any():
        o.append(obs[~keep].sum()); e.append(exp[~keep].sum())
    o = np.array(o); e = np.array(e)
    good = e > 1e-6
    chi2 = float((((o - e) ** 2)[good] / e[good]).sum())
    dof = int(good.sum() - 1)
    try:
        from scipy.stats import chi2 as _chi2
        pval = float(_chi2.sf(chi2, dof)) if dof > 0 else float("nan")
    except Exception:
        pval = float("nan")
    return {"tvd": float(tvd), "off_support": off_support,
            "chi2": chi2, "dof": dof, "chi2_p": pval}
