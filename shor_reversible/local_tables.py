"""Logical (pre-transpilation) gate counts from the IR for every (family, N, a),
before/after identity-block removal, plus exact order / identity blocks / p(N,a).
These equal Qiskit's pre-transpile count_ops (1:1 adapter, plus t measurements)."""
import csv
from math import gcd
from qpe import shor_circuit
from analysis import order, identity_blocks, success_probability
rows = []
for N in (15, 21, 35):
    t = 2 * N.bit_length()
    for a in range(2, N):
        if gcd(a, N) != 1: continue
        p, _, _ = success_probability(N, a, t)
        for fam in ("beauregard", "cuccaro", "montgomery"):
            rec = {"family": fam, "N": N, "a": a, "order": order(a, N),
                   "identity_blocks": identity_blocks(N, a, t), "p_exact": round(p, 4)}
            for rem in (False, True):
                c, m, info = shor_circuit(fam, N, a, remove_identity=rem)
                co = c.count_ops(); st = "after" if rem else "before"
                rec["qubits"] = c.nq
                rec[f"{st}_gates"] = len(c.ops)
                rec[f"{st}_2q+"] = sum(v for k, v in co.items() if k not in ("h", "x", "p"))
                rec[f"{st}_depth"] = c.depth()
            rec["gate_reduction_%"] = round(100 * (1 - rec["after_gates"] / rec["before_gates"]), 1)
            rows.append(rec)
with open("logical_counts.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
# compact summary
from collections import defaultdict
g = defaultdict(list)
for r in rows: g[(r["family"], r["N"], r["order"], r["identity_blocks"])].append(r)
print(f"{'family':11s} {'N':>3s} {'r':>3s} {'id':>3s} {'#a':>3s} {'q':>3s} {'before(min-max)':>17s} {'after(min-max)':>17s} {'red%':>6s}")
for k in sorted(g):
    L = g[k]; b = [r["before_gates"] for r in L]; a_ = [r["after_gates"] for r in L]
    print(f"{k[0]:11s} {k[1]:3d} {k[2]:3d} {k[3]:3d} {len(L):3d} {L[0]['qubits']:3d} {min(b):8d}-{max(b):<8d} {min(a_):8d}-{max(a_):<8d} {L[0]['gate_reduction_%']:6.1f}")
