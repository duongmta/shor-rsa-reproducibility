import numpy as np, random
from revcirc import Circuit, run_statevector, run_bits
from math import pi
def dense(c):
    n=c.nq; D=2**n; U=np.eye(D,dtype=complex)
    for name,qs,th in c.ops:
        G=np.zeros((D,D),dtype=complex)
        for i in range(D):
            b=lambda q:(i>>q)&1
            if name=="h":
                q=qs[0]; j=i^(1<<q); s=1/np.sqrt(2)
                G[i,i]+= s*(-1 if b(q) else 1); G[j,i]+=s
            elif name in("p","cp","ccp"):
                G[i,i]=np.exp(1j*th) if all(b(q) for q in qs) else 1
            else:
                j=i
                if name=="x": j=i^(1<<qs[0])
                if name=="cx" and b(qs[0]): j=i^(1<<qs[1])
                if name=="ccx" and b(qs[0]) and b(qs[1]): j=i^(1<<qs[2])
                if name in("swap","cswap"):
                    ok = name=="swap" or b(qs[0]); a,bb=qs[-2:]
                    if ok and b(a)!=b(bb): j=i^(1<<a)^(1<<bb)
                G[j,i]=1
        U=G@U
    return U
random.seed(1)
for trial in range(30):
    c=Circuit(4)
    for _ in range(25):
        g=random.choice(["h","x","cx","ccx","p","cp","ccp","swap","cswap"])
        qs=random.sample(range(4),{"h":1,"x":1,"p":1,"cx":2,"cp":2,"swap":2,"ccx":3,"ccp":3,"cswap":3}[g])
        th=random.uniform(-3,3)
        getattr(c,g)(*qs,th) if g in("p","cp","ccp") else getattr(c,g)(*qs)
    U=dense(c); S=run_statevector(c,np.eye(16,dtype=complex))
    assert np.allclose(U,S), trial
    assert np.allclose(dense(c.inverse())@U,np.eye(16))
print("statevector sim OK")
