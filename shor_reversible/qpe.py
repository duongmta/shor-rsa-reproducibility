"""Full-register Shor order-finding circuit built from one arithmetic family.

Layout: phase register P[0..t-1] = qubits 0..t-1 (P[j] controls U^(2^j)),
work registers of the family after it; work register x starts in |1>.
Readout: inverse QFT (no swaps) on reg[q] = P[t-1-q]; clbit q <- P[t-1-q],
so the classical register directly holds the integer y.
"""
from revcirc import Circuit
from families import FAMILIES, iqft


def shor_circuit(family, N, a, t=None, remove_identity=False):
    n = N.bit_length()
    t = 2 * n if t is None else t
    F = FAMILIES[family]
    L = F.layout(n, base=t)
    nq = t + L["nq"]
    c = Circuit(nq)
    P = list(range(t))
    for q in P:
        c.h(q)
    c.x(L["x"][0])                          # |x> = |1>
    blocks = []
    for j in range(t):
        k = pow(a, 2 ** j, N)
        if remove_identity and k == 1:
            blocks.append((j, k, "removed"))
            continue
        start = len(c.ops)
        F.cu(c, L, k, N, P[j])
        blocks.append((j, k, len(c.ops) - start))
    reg = [P[t - 1 - q] for q in range(t)]
    iqft(c, reg)
    measure = [(reg[q], q) for q in range(t)]   # (qubit, clbit)
    info = {"family": family, "N": N, "a": a, "n": n, "t": t, "nq": nq,
            "remove_identity": remove_identity, "blocks": blocks, "layout": L}
    return c, measure, info
