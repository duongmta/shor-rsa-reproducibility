"""Three generic reversible constructions of the controlled modular multiplier
    CU_k : |c>|x>|0..0>  ->  |c>|k^c * x mod N>|0..0>      (0 <= x < N, gcd(k,N)=1)
used inside Shor order finding.

  A  "beauregard" : Draper QFT-basis constant adder, Beauregard modular adder,
                    CMULT(k) ; CSWAP ; CMULT(k^-1)^dagger          (2n+3 qubits incl. control)
  B  "cuccaro"    : Toffoli ripple-carry adder (Cuccaro-Draper-Kutin-Moulton),
                    VBE-style modular adder with a shared constant register
                                                                      (3n+4 qubits incl. control)
  C  "montgomery" : out-of-place Montgomery multiplier (parity-based REDC with
                    Draper additions), y ^= k*x mod N, CSWAP, y ^= k^-1*x' mod N
                                                                      (4n+4 qubits incl. control)
All three are generated for arbitrary odd N and k from the same gate set, with no
knowledge of the order of k.
"""
from math import pi
from revcirc import Circuit

# --------------------------------------------------------------------------
# QFT without final swaps.  F|b> puts phase exp(2*pi*i*b/2^(q+1)) on reg[q].
# --------------------------------------------------------------------------
def qft(c, reg):
    m = len(reg)
    for q in range(m - 1, -1, -1):
        c.h(reg[q])
        for k in range(q - 1, -1, -1):
            c.cp(reg[k], reg[q], pi / 2 ** (q - k))


def iqft(c, reg):
    tmp = Circuit(c.nq)
    qft(tmp, reg)
    c.extend(tmp.inverse())


def phi_add(c, reg, k, controls=()):
    """Draper constant adder in the Fourier frame of qft(): reg <- reg + k mod 2^m."""
    for q in range(len(reg)):
        c.phase(list(controls), reg[q], 2 * pi * k / 2 ** (q + 1))


# --------------------------------------------------------------------------
# Family A: Beauregard (quant-ph/0205095)
# --------------------------------------------------------------------------
class Beauregard:
    name = "beauregard"

    @staticmethod
    def layout(n, base=0):
        L = {"x": list(range(base, base + n))}
        L["b"] = list(range(base + n, base + 2 * n + 1))    # n+1 qubits
        L["anc"] = base + 2 * n + 1
        L["nq"] = 2 * n + 2                                  # + external control
        return L

    @staticmethod
    def _cc_mod_add(c, L, k, N, c1, c2):
        b, anc = L["b"], L["anc"]
        msb = b[-1]
        phi_add(c, b, k, (c1, c2))
        phi_add(c, b, -N)
        iqft(c, b); c.cx(msb, anc); qft(c, b)
        phi_add(c, b, N, (anc,))
        phi_add(c, b, -k, (c1, c2))
        iqft(c, b); c.x(msb); c.cx(msb, anc); c.x(msb); qft(c, b)
        phi_add(c, b, k, (c1, c2))

    @classmethod
    def _cmult(cls, c, L, k, N, ctrl):
        n = len(L["x"])
        qft(c, L["b"])
        for i in range(n):
            cls._cc_mod_add(c, L, (k << i) % N, N, ctrl, L["x"][i])
        iqft(c, L["b"])

    @classmethod
    def cu(cls, c, L, k, N, ctrl):
        kinv = pow(k, -1, N)
        cls._cmult(c, L, k, N, ctrl)
        for i in range(len(L["x"])):
            c.cswap(ctrl, L["x"][i], L["b"][i])
        tmp = Circuit(c.nq)
        cls._cmult(tmp, L, kinv, N, ctrl)
        c.extend(tmp.inverse())


# --------------------------------------------------------------------------
# Family B: Toffoli ripple-carry (Cuccaro et al., quant-ph/0410184) + VBE-style
# modular addition of a classical constant loaded into a shared register R.
# --------------------------------------------------------------------------
def _maj(c, x, y, z):
    c.cx(z, y); c.cx(z, x); c.ccx(x, y, z)


def _uma(c, x, y, z):
    c.ccx(x, y, z); c.cx(z, x); c.cx(x, y)


def cuccaro_add(c, a, b, z, c0):
    """b <- a + b mod 2^n, z ^= carry-out.  a,b: n-qubit lists (LSB first)."""
    n = len(a)
    _maj(c, c0, b[0], a[0])
    for i in range(1, n):
        _maj(c, a[i - 1], b[i], a[i])
    c.cx(a[n - 1], z)
    for i in range(n - 1, 0, -1):
        _uma(c, a[i - 1], b[i], a[i])
    _uma(c, c0, b[0], a[0])


class Cuccaro:
    name = "cuccaro"

    @staticmethod
    def layout(n, base=0):
        L = {"x": list(range(base, base + n))}
        L["b"] = list(range(base + n, base + 2 * n + 1))        # n+1 (b[n] = sign/carry)
        L["r"] = list(range(base + 2 * n + 1, base + 3 * n + 1))  # constant register
        L["c0"] = base + 3 * n + 1
        L["t"] = base + 3 * n + 2
        L["nq"] = 3 * n + 3
        return L

    @staticmethod
    def _add(c, L):
        cuccaro_add(c, L["r"], L["b"][:-1], L["b"][-1], L["c0"])

    @staticmethod
    def _sub(c, L):
        tmp = Circuit(c.nq)
        Cuccaro._add(tmp, L)
        c.extend(tmp.inverse())

    @staticmethod
    def _load(c, L, k, c1, c2):
        for j, q in enumerate(L["r"]):
            if k >> j & 1:
                c.ccx(c1, c2, q)

    @classmethod
    def _cc_mod_add(cls, c, L, k, N, c1, c2):
        r, t, msb = L["r"], L["t"], L["b"][-1]
        nbits = [q for j, q in enumerate(r) if N >> j & 1]
        cls._load(c, L, k, c1, c2); cls._add(c, L); cls._load(c, L, k, c1, c2)
        for q in nbits: c.x(q)                 # R = N
        cls._sub(c, L)
        c.cx(msb, t)                           # t = 1  <=>  b + k - N < 0
        c.x(t)
        for q in nbits: c.cx(t, q)             # R = t*N
        c.x(t)
        cls._add(c, L)
        for q in nbits: c.cx(t, q)             # R = 0
        cls._load(c, L, k, c1, c2); cls._sub(c, L)
        c.x(msb); c.cx(msb, t); c.x(msb)       # uncompute t
        cls._add(c, L); cls._load(c, L, k, c1, c2)

    @classmethod
    def _cmult(cls, c, L, k, N, ctrl):
        for i in range(len(L["x"])):
            cls._cc_mod_add(c, L, (k << i) % N, N, ctrl, L["x"][i])

    @classmethod
    def cu(cls, c, L, k, N, ctrl):
        kinv = pow(k, -1, N)
        cls._cmult(c, L, k, N, ctrl)
        for i in range(len(L["x"])):
            c.cswap(ctrl, L["x"][i], L["b"][i])
        tmp = Circuit(c.nq)
        cls._cmult(tmp, L, kinv, N, ctrl)
        c.extend(tmp.inverse())


# --------------------------------------------------------------------------
# Family C: Montgomery multiplication (REDC with R = 2^n), Draper additions.
#   MM_k: |c>|x>|y>|0> -> |c>|x>|y XOR (c*k*x mod N)>|0>
# --------------------------------------------------------------------------
class Montgomery:
    name = "montgomery"

    @staticmethod
    def layout(n, base=0):
        L = {"x": list(range(base, base + n))}
        L["y"] = list(range(base + n, base + 2 * n))
        L["T"] = list(range(base + 2 * n, base + 3 * n + 2))    # n+2 accumulator
        L["q"] = list(range(base + 3 * n + 2, base + 4 * n + 2))  # REDC quotient bits
        L["f"] = base + 4 * n + 2
        L["nq"] = 4 * n + 3
        return L

    @staticmethod
    def _redc_mult(c, L, k, N, ctrl):
        """Compute T = c*k*x mod N (garbage in q, f). Returns final logical T order."""
        n = len(L["x"])
        kp = (k << n) % N              # Montgomery form k*R mod N
        T = list(L["T"])
        for i in range(n):
            qft(c, T); phi_add(c, T, kp, (ctrl, L["x"][i])); iqft(c, T)
            c.cx(T[0], L["q"][i])      # quotient bit = parity
            qft(c, T); phi_add(c, T, N, (L["q"][i],)); iqft(c, T)
            T = T[1:] + T[:1]          # exact division by 2 = relabelling
        qft(c, T); phi_add(c, T, -N); iqft(c, T)
        c.cx(T[-1], L["f"])            # f = 1 <=> T < N
        qft(c, T); phi_add(c, T, N, (L["f"],)); iqft(c, T)
        return T

    @classmethod
    def _mm(cls, c, L, k, N, ctrl):
        tmp = Circuit(c.nq)
        T = cls._redc_mult(tmp, L, k, N, ctrl)
        c.extend(tmp)
        for j in range(len(L["y"])):
            c.cx(T[j], L["y"][j])
        c.extend(tmp.inverse())

    @classmethod
    def cu(cls, c, L, k, N, ctrl):
        kinv = pow(k, -1, N)
        cls._mm(c, L, k, N, ctrl)                  # y = k x
        for i in range(len(L["x"])):
            c.cswap(ctrl, L["x"][i], L["y"][i])    # x <-> y
        cls._mm(c, L, kinv, N, ctrl)               # y ^= k^-1 (k x) = x  -> 0


FAMILIES = {F.name: F for F in (Beauregard, Cuccaro, Montgomery)}


def work_layout(family, n):
    """Layout with the external control at qubit 0 and work registers after it."""
    F = FAMILIES[family]
    L = F.layout(n, base=1)
    L["ctrl"] = 0
    L["nq"] += 1
    return L


def controlled_u(family, N, k):
    """Stand-alone controlled-U_k block (control = qubit 0) for verification."""
    n = N.bit_length()
    L = work_layout(family, n)
    c = Circuit(L["nq"])
    FAMILIES[family].cu(c, L, k, N, 0)
    return c, L
