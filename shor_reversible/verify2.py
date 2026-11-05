import sys, time
from math import gcd
from families import controlled_u
from hybrid import run_hybrid
from verify import enc
fams = sys.argv[1:] or ["beauregard", "montgomery"]
for fam in fams:
    t0=time.time(); bad=tot=0; worst=1.0
    for N in (15, 21, 35):
        for k in range(1, N):
            if gcd(k, N) != 1: continue
            c, L = controlled_u(fam, N, k)
            qset = L["b"] if fam == "beauregard" else L["T"]
            ok = True
            for cb in (0, 1):
                for x in range(N):
                    out, p = run_hybrid(c, qset, enc(L, cb, x))
                    worst = min(worst, p)
                    if out != enc(L, cb, (k*x) % N if cb else x) or p < 1-1e-9: ok = False
            tot += 1; bad += not ok
            if not ok: print("FAIL", fam, N, k)
    print(f"{fam}: {tot-bad}/{tot} blocks correct (k=1 included), min output prob {worst:.12f}, {time.time()-t0:.0f}s")
