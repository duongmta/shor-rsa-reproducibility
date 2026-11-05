"""Exact hybrid simulator for basis-state inputs: qubits in `qset` are held in a
statevector, all other qubits are classical bits.  A classical qubit may only be
written from a quantum qubit when that qubit is in a definite basis state (checked).
Used to verify blocks whose Fourier-basis arithmetic lives on a small register."""
import numpy as np
from revcirc import PHASE

class NotDefinite(Exception):
    pass

def run_hybrid(circ, qset, bits, tol=1e-9):
    qset = list(qset); loc = {q: i for i, q in enumerate(qset)}
    m = len(qset)
    psi = np.zeros(2 ** m, dtype=complex)
    idx0 = sum(1 << loc[q] for q in qset if bits >> q & 1)
    psi[idx0] = 1
    bits &= ~sum(1 << q for q in qset)
    ar = np.arange(2 ** m)
    def mask(q): return (ar >> loc[q]) & 1
    def definite(q):
        p1 = float(np.sum(np.abs(psi[mask(q) == 1]) ** 2))
        if p1 < tol: return 0
        if p1 > 1 - tol: return 1
        raise NotDefinite(q)
    gphase = 1.0 + 0j
    for name, qs, th in circ.ops:
        Q = [q for q in qs if q in loc]
        if name in PHASE:
            if any(not (bits >> q & 1) for q in qs if q not in loc):
                continue
            if not Q:
                gphase *= np.exp(1j * th); continue
            sel = np.ones(2 ** m, bool)
            for q in Q: sel &= mask(q) == 1
            psi[sel] *= np.exp(1j * th); continue
        if name == "h":
            q = qs[0]; assert q in loc, "h on classical qubit"
            b = loc[q]; i0 = ar[mask(q) == 0]; i1 = i0 | (1 << b)
            v0 = psi[i0].copy(); v1 = psi[i1].copy()
            psi[i0] = (v0 + v1) / np.sqrt(2); psi[i1] = (v0 - v1) / np.sqrt(2); continue
        if name in ("x", "cx", "ccx"):
            ctrl, tg = qs[:-1], qs[-1]
        elif name == "swap":
            ctrl, tg = (), qs
        else:
            ctrl, tg = qs[:1], qs[1:]
        # resolve controls
        cl_ok = all(bits >> c & 1 for c in ctrl if c not in loc)
        if not cl_ok: continue
        qctrl = [c for c in ctrl if c in loc]
        tgt_q = [t for t in (tg if isinstance(tg, tuple) else (tg,)) if t in loc]
        tgt_all = tg if isinstance(tg, tuple) else (tg,)
        if len(tgt_q) == len(tgt_all):          # all targets quantum -> unitary on psi
            sel = np.ones(2 ** m, bool)
            for c in qctrl: sel &= mask(c) == 1
            if name in ("x", "cx", "ccx"):
                b = loc[tgt_all[0]]
                i0 = ar[sel & (mask(tgt_all[0]) == 0)]; i1 = i0 | (1 << b)
                v = psi[i0].copy(); psi[i0] = psi[i1]; psi[i1] = v
            else:
                a, bq = tgt_all; ba, bb = loc[a], loc[bq]
                i01 = ar[sel & (mask(a) == 1) & (mask(bq) == 0)]; i10 = i01 ^ (1 << ba) ^ (1 << bb)
                v = psi[i01].copy(); psi[i01] = psi[i10]; psi[i10] = v
            continue
        # some target classical: all quantum controls must be definite
        if any(definite(c) == 0 for c in qctrl): continue
        if name in ("x", "cx", "ccx"):
            bits ^= 1 << tgt_all[0]; continue
        a, bq = tgt_all
        va = definite(a) if a in loc else bits >> a & 1
        vb = definite(bq) if bq in loc else bits >> bq & 1
        if va == vb: continue
        for q in (a, bq):
            if q in loc:
                b = loc[q]; i0 = ar[mask(q) == 0]; i1 = i0 | (1 << b)
                v = psi[i0].copy(); psi[i0] = psi[i1]; psi[i1] = v
            else:
                bits ^= 1 << q
    # final: quantum part must be a basis state
    k = int(np.argmax(np.abs(psi)))
    prob = abs(psi[k]) ** 2
    out = bits | sum(1 << q for q in qset if k >> loc[q] & 1)
    return out, prob
