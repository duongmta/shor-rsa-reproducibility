import sys,json,math
from pathlib import Path
import numpy as np
from qiskit import transpile
from qiskit_aer import AerSimulator
ROOT=Path(__file__).resolve().parent;sys.path.insert(0,str(ROOT.parent/'shor_reversible'))
from qpe import shor_circuit
from run_qiskit import to_qiskit
from analysis import success_probability
out=[]
for cutoff in [None,4,2]:
 c,meas,info=shor_circuit('beauregard',21,2,remove_identity=True)
 total=len(c.ops);removed=0
 if cutoff is not None:
  new=[]
  for name,qs,th in c.ops:
   if name=='cp' and max(qs)<info['t'] and abs(th)<math.pi/2**cutoff-1e-14:removed+=1
   else:new.append((name,qs,th))
  c.ops=new
 qc=to_qiskit(c,[])
 tqc=transpile(qc,basis_gates=['rz','sx','x','cx'],optimization_level=1,seed_transpiler=20260929)
 tqc.save_probabilities([q for q,cb in sorted(meas,key=lambda z:z[1])])
 res=AerSimulator(method='matrix_product_state',matrix_product_state_truncation_threshold=1e-16,max_parallel_threads=1).run(tqc).result()
 P=np.array(res.data(0)['probabilities']);pex,E,acc=success_probability(21,2,10)
 out.append(dict(cutoff=cutoff,removed_cp=removed,logical_gates=len(c.ops),cx=tqc.count_ops().get('cx',0),p=float(P[acc].sum()),tvd=float(abs(P-E).sum()/2),max_error=float(max(abs(P-E)))))
 print(out[-1],flush=True)
(ROOT/'aqft_results.json').write_text(json.dumps(out,indent=2))
