"""Exhaustive functional verification of controlled-U_k for all families."""
import sys, time, numpy as np
from math import gcd
from revcirc import Circuit, run_bits, run_statevector
from families import controlled_u, cuccaro_add

def enc(L, ctrl, x):
    s = ctrl << L["ctrl"]
    for i, q in enumerate(L["x"]):
        if x >> i & 1: s |= 1 << q
    return s

def verify(family, N, k):
    c, L = controlled_u(family, N, k)
    n = len(L["x"])
    inputs = [(cb, x) for cb in (0, 1) for x in range(N)]
    expect = [enc(L, cb, (k * x) % N if cb else x) for cb, x in inputs]
    if c.is_permutation():
        out = [run_bits(c, enc(L, cb, x)) for cb, x in inputs]
        ok = out == expect
        err = 0.0 if ok else 1.0
    else:
        D = 2 ** c.nq
        psi = np.zeros((D, len(inputs)), dtype=complex)
        for j, (cb, x) in enumerate(inputs):
            psi[enc(L, cb, x), j] = 1
        psi = run_statevector(c, psi)
        err = max(1 - abs(psi[e, j]) ** 2 for j, e in enumerate(expect))
        ok = err < 1e-9
    return ok, err, c

if __name__ == "__main__":
    # Cuccaro adder sanity (all 4-bit inputs)
    for n in (3, 4):
        c = Circuit(2 * n + 2)
        a = list(range(n)); b = list(range(n, 2 * n)); z = 2 * n; c0 = 2 * n + 1
        cuccaro_add(c, a, b, z, c0)
        for A in range(2 ** n):
            for B in range(2 ** n):
                s = run_bits(c, A | B << n)
                assert s & (2 ** n - 1) == A and (s >> n) & (2 ** n - 1) == (A + B) % 2 ** n
                assert (s >> z & 1) == ((A + B) >> n) and not (s >> c0 & 1)
    print("cuccaro adder OK")
    fams = sys.argv[1:] or ["cuccaro", "beauregard", "montgomery"]
    for fam in fams:
        t0 = time.time(); bad = 0; tot = 0; maxerr = 0
        for N in (15, 21, 35):
            for k in range(1, N):
                if gcd(k, N) != 1: continue
                ok, err, c = verify(fam, N, k)
                tot += 1; bad += (not ok); maxerr = max(maxerr, err)
                if not ok: print("FAIL", fam, N, k, err)
        print(f"{fam}: {tot-bad}/{tot} blocks correct, max error {maxerr:.2e}, {time.time()-t0:.1f}s")
