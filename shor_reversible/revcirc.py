"""Minimal, dependency-free gate IR for reversible Shor circuits.

Gates (name, qubits, param):
  h(q) x(q) cx(c,t) ccx(c1,c2,t) p(q;th) cp(c,t;th) ccp(c1,c2,t;th)
  swap(a,b) cswap(c,a,b)
All builders emit these gates only, so the same list can be
  * simulated here (bit-level for permutation gates, numpy statevector otherwise),
  * converted 1:1 to a Qiskit QuantumCircuit (see qiskit_adapter.py).
"""
from math import pi

SELF_INVERSE = {"h", "x", "cx", "ccx", "swap", "cswap"}
PHASE = {"p", "cp", "ccp"}
ARITY = {"h": 1, "x": 1, "p": 1, "cx": 2, "cp": 2, "swap": 2,
         "ccx": 3, "ccp": 3, "cswap": 3}


def _norm(th):
    th = th % (2 * pi)
    if th > pi:
        th -= 2 * pi
    return th


class Circuit:
    def __init__(self, nq=0):
        self.nq = nq
        self.ops = []          # list of (name, tuple(qubits), param or None)

    # --- gate emitters -------------------------------------------------
    def _add(self, name, qs, th=None):
        assert len(qs) == ARITY[name], (name, qs)
        assert len(set(qs)) == len(qs), (name, qs)
        if name in PHASE:
            th = _norm(th)
            if abs(th) < 1e-12:      # exact identity: never emitted
                return
        self.ops.append((name, tuple(qs), th))

    def h(self, q): self._add("h", (q,))
    def x(self, q): self._add("x", (q,))
    def cx(self, c, t): self._add("cx", (c, t))
    def ccx(self, a, b, t): self._add("ccx", (a, b, t))
    def p(self, q, th): self._add("p", (q,), th)
    def cp(self, c, t, th): self._add("cp", (c, t), th)
    def ccp(self, a, b, t, th): self._add("ccp", (a, b, t), th)
    def swap(self, a, b): self._add("swap", (a, b))
    def cswap(self, c, a, b): self._add("cswap", (c, a, b))

    def phase(self, controls, t, th):
        """Phase on target conditioned on 0,1 or 2 extra controls."""
        controls = [c for c in controls if c is not None]
        if len(controls) == 0:
            self.p(t, th)
        elif len(controls) == 1:
            self.cp(controls[0], t, th)
        elif len(controls) == 2:
            self.ccp(controls[0], controls[1], t, th)
        else:
            raise ValueError("at most 2 controls")

    def extend(self, other):
        self.ops.extend(other.ops)
        return self

    def inverse(self):
        inv = Circuit(self.nq)
        for name, qs, th in reversed(self.ops):
            inv.ops.append((name, qs, -th if name in PHASE else None))
        return inv

    def remap(self, mapping):
        c = Circuit(self.nq)
        for name, qs, th in self.ops:
            c.ops.append((name, tuple(mapping[q] for q in qs), th))
        return c

    # --- statistics ----------------------------------------------------
    def count_ops(self):
        d = {}
        for name, _, _ in self.ops:
            d[name] = d.get(name, 0) + 1
        return d

    def depth(self):
        level = [0] * self.nq
        for _, qs, _ in self.ops:
            L = max(level[q] for q in qs) + 1
            for q in qs:
                level[q] = L
        return max(level) if level else 0

    def is_permutation(self):
        return all(name in SELF_INVERSE and name != "h" for name, _, _ in self.ops)


# ---------------------------------------------------------------------------
# Bit-level simulator (permutation gates only). State = python int bit mask.
# ---------------------------------------------------------------------------
def run_bits(circ, state):
    for name, qs, _ in circ.ops:
        if name == "x":
            state ^= 1 << qs[0]
        elif name == "cx":
            if state >> qs[0] & 1:
                state ^= 1 << qs[1]
        elif name == "ccx":
            if (state >> qs[0] & 1) and (state >> qs[1] & 1):
                state ^= 1 << qs[2]
        elif name == "swap":
            a, b = qs
            if (state >> a & 1) != (state >> b & 1):
                state ^= (1 << a) | (1 << b)
        elif name == "cswap":
            c, a, b = qs
            if (state >> c & 1) and (state >> a & 1) != (state >> b & 1):
                state ^= (1 << a) | (1 << b)
        else:
            raise ValueError("non-permutation gate " + name)
    return state


# ---------------------------------------------------------------------------
# Batched numpy statevector simulator: psi has shape (2**nq, B), qubit q is
# bit q of the row index (little-endian, same convention as Qiskit).
# ---------------------------------------------------------------------------
def run_statevector(circ, psi):
    import numpy as np
    nq = circ.nq
    B = psi.shape[1]
    psi = psi.reshape((2,) * nq + (B,))

    def ax(q):  # numpy axis of qubit q
        return nq - 1 - q

    s2 = 1 / np.sqrt(2)
    for name, qs, th in circ.ops:
        if name == "h":
            a = ax(qs[0])
            s0 = [slice(None)] * (nq + 1); s0[a] = 0
            s1 = [slice(None)] * (nq + 1); s1[a] = 1
            v0 = psi[tuple(s0)].copy(); v1 = psi[tuple(s1)].copy()
            psi[tuple(s0)] = (v0 + v1) * s2
            psi[tuple(s1)] = (v0 - v1) * s2
            continue
        if name in PHASE:
            idx = [slice(None)] * (nq + 1)
            for q in qs:
                idx[ax(q)] = 1
            psi[tuple(idx)] *= np.exp(1j * th)
            continue
        # permutation gates: controls + action on targets
        if name == "x":
            ctrl, tg = (), qs
        elif name == "cx":
            ctrl, tg = qs[:1], qs[1:]
        elif name == "ccx":
            ctrl, tg = qs[:2], qs[2:]
        elif name == "swap":
            ctrl, tg = (), qs
        elif name == "cswap":
            ctrl, tg = qs[:1], qs[1:]
        idx = [slice(None)] * (nq + 1)
        for c in ctrl:
            idx[ax(c)] = 1
        if name in ("x", "cx", "ccx"):
            i0 = list(idx); i0[ax(tg[0])] = 0
            i1 = list(idx); i1[ax(tg[0])] = 1
            v0 = psi[tuple(i0)].copy()
            psi[tuple(i0)] = psi[tuple(i1)]
            psi[tuple(i1)] = v0
        else:  # swap / cswap on (a,b): exchange |01> and |10>
            a, b = tg
            i01 = list(idx); i01[ax(a)] = 1; i01[ax(b)] = 0
            i10 = list(idx); i10[ax(a)] = 0; i10[ax(b)] = 1
            v = psi[tuple(i01)].copy()
            psi[tuple(i01)] = psi[tuple(i10)]
            psi[tuple(i10)] = v
    return psi.reshape(2 ** nq, B)
