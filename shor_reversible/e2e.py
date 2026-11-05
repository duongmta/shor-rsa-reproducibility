import sys, time, numpy as np
from qpe import shor_circuit
from revcirc import run_statevector
from analysis import qpe_distribution, order
fam, N, a, rem = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4] == "after"
c, meas, info = shor_circuit(fam, N, a, remove_identity=rem)
t0 = time.time()
psi = np.zeros((2 ** c.nq, 1), complex); psi[0, 0] = 1
psi = run_statevector(c, psi)[:, 0]
prob = np.abs(psi) ** 2; idx = np.arange(2 ** c.nq); y = np.zeros(2 ** c.nq, np.int64)
for q, cb in meas: y |= ((idx >> q) & 1) << cb
Py = np.zeros(2 ** info["t"]); np.add.at(Py, y, prob)
# work/ancilla registers must end in |a^k mod N> only; check ancillas clean
P = qpe_distribution(order(a, N), 2 ** info["t"])
print(fam, N, a, "after" if rem else "before", "qubits", c.nq, "gates", len(c.ops),
      "max|P_sim-P_exact| =", float(np.abs(Py - P).max()), f"{time.time()-t0:.0f}s", flush=True)
