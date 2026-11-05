"""Independent gate-by-gate NumPy statevector runs of full circuits, every base.
Records max |P_sim - P_exact| and a seeded 1,024-shot sample from P_sim."""
import sys, json, time, numpy as np
from math import gcd
from qpe import shor_circuit
from revcirc import run_statevector
from analysis import success_probability, recovers_factor, compare_counts, wilson, order, identity_blocks
fam, N = sys.argv[1], int(sys.argv[2])
for a in range(2, N):
    if gcd(a, N) != 1: continue
    c, meas, info = shor_circuit(fam, N, a, remove_identity=True)
    t0 = time.time()
    psi = np.zeros((2 ** c.nq, 1), complex); psi[0, 0] = 1
    psi = run_statevector(c, psi)[:, 0]
    prob = np.abs(psi) ** 2; del psi
    idx = np.arange(2 ** c.nq); y = np.zeros(2 ** c.nq, np.int64)
    for q, cb in meas: y |= ((idx >> q) & 1) << cb
    t = info["t"]; Q = 2 ** t
    Py = np.zeros(Q); np.add.at(Py, y, prob); del prob, idx, y
    p_exact, P, acc = success_probability(N, a, t)
    rng = np.random.default_rng(20260929 + 1000 * N + a + 100000 * {"beauregard": 0, "cuccaro": 1, "montgomery": 2}[fam])
    samp = rng.choice(Q, size=1024, p=Py / Py.sum())
    counts = {}
    for v in samp: counts[int(v)] = counts.get(int(v), 0) + 1
    k = sum(cn for yy, cn in counts.items() if recovers_factor(yy, Q, a, N))
    st = compare_counts(counts, P, 1024)
    rec = {"family": fam, "N": N, "a": a, "qubits": c.nq, "gates": len(c.ops), "order": order(a, N),
           "identity_blocks": identity_blocks(N, a, t), "max_dev": float(np.abs(Py - P).max()),
           "p_exact": p_exact, "p_sim_exact": float(Py[acc].sum()), "success_count": k,
           "p_hat": k / 1024, "wilson95": wilson(k, 1024), **st, "wall_s": round(time.time() - t0, 1)}
    print(json.dumps(rec), flush=True)
