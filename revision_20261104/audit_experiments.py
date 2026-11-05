import sys,json,csv,math,collections,platform,time
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parent
SRC=ROOT.parent/'shor_reversible'
sys.path.insert(0,str(SRC))
from qpe import shor_circuit
from analysis import order,convergents,wilson,recovers_factor,qpe_distribution
from fractions import Fraction
ROOT.mkdir(exist_ok=True)
def dump(name,obj):
 (ROOT/name).write_text(json.dumps(obj,indent=2,ensure_ascii=False,default=str),encoding='utf-8')
def accept(y,Q,a,N,L=1):
 for cv in convergents(Fraction(int(y),Q)):
  d=cv.denominator
  # Include d=1 explicitly: the classical budget can itself find short orders.
  for j in range(1,L+1):
   D=j*d
   if D>=N:break
   if D%2==0 and pow(a,D,N)==1:
    x=pow(a,D//2,N)
    if 1<math.gcd(x-1,N)<N and 1<math.gcd(x+1,N)<N:return True
 return False
def prob_chunk(r,Q,start,stop):
 y=np.arange(start,stop,dtype=np.int64)
 d=((r*y)%Q).astype(float)/Q
 d[d>.5]-=1
 m,l=divmod(Q,r)
 return (l*(m+1)**2*np.sinc((m+1)*d)**2+(r-l)*m*m*np.sinc(m*d)**2)/(Q*Q*np.sinc(d)**2)
def analytic():
 for r in [2,3,4,6,12]:
  assert np.max(abs(prob_chunk(r,1024,0,1024)-qpe_distribution(r,1024)))<1e-12
 small=[]
 for N in [15,21,35]:
  for a in range(2,N):
   if math.gcd(a,N)!=1:continue
   for t in [2*N.bit_length()-2,2*N.bit_length(),2*N.bit_length()+2]:
    Q=2**t;P=prob_chunk(order(a,N),Q,0,Q)
    for L in [1,8]:
     ac=np.array([accept(y,Q,a,N,L) for y in range(Q)])
     small.append(dict(N=N,a=a,t=t,L=L,p=float(P[ac].sum()),p_rand=float(ac.mean())))
 dump('postprocessing_small.json',small)
 large=[]
 for N in [899,1763,3233]:
  for a in [2,3,5]:
   if math.gcd(a,N)!=1:continue
   t=2*N.bit_length();Q=2**t;r=order(a,N);shots=20000
   rng=np.random.default_rng(20261104+N*10+a)
   u=np.sort(rng.random(shots));ys=np.zeros(shots,dtype=np.int64);cum=0.;pos=0
   for start in range(0,Q,262144):
    p=prob_chunk(r,Q,start,min(Q,start+262144));cs=np.cumsum(p)+cum
    end=int(np.searchsorted(u,cs[-1],side='right'))
    if end>pos:ys[pos:end]=start+np.searchsorted(cs,u[pos:end]);pos=end
    cum=float(cs[-1])
   assert abs(cum-1)<1e-7 and pos==shots,(N,a,cum,pos)
   yr=rng.integers(0,Q,size=shots)
   for L in [1,8]:
    k=sum(accept(y,Q,a,N,L) for y in ys);kr=sum(accept(y,Q,a,N,L) for y in yr)
    large.append(dict(N=N,a=a,n=N.bit_length(),t=t,r=r,L=L,shots=shots,k=k,p_hat=k/shots,ci=wilson(k,shots),random_k=kr,p_rand_hat=kr/shots,random_ci=wilson(kr,shots),normalization=cum))
   print('analytic large',N,a,r,flush=True)
 dump('postprocessing_large.json',large)
 rows=list(csv.DictReader(open(SRC/'logical_counts.csv',encoding='utf-8')))
 groups=collections.defaultdict(list)
 for row in rows:groups[(row['family'],row['N'],row['order'],row['identity_blocks'])].append(row)
 spreads={}
 for stage in ['before','after']:
  values=[]
  for key,g in groups.items():
   vals=[int(r[stage+'_gates']) for r in g]
   values.append(((max(vals)-min(vals))/np.mean(vals)*100,key,vals))
  spreads[stage]=max(values)
 dump('spreads.json',spreads)
def quantum():
 import qiskit,qiskit_aer
 from qiskit import transpile
 from qiskit_aer import AerSimulator
 from run_qiskit import to_qiskit,metrics
 dump('new_environment.json',dict(python=sys.version,qiskit=qiskit.__version__,aer=qiskit_aer.__version__,numpy=np.__version__,platform=platform.platform()))
 diag=[]
 for trunc in [True,False]:
  for rem in [False,True]:
   c,meas,inf=shor_circuit('beauregard',15,7,remove_identity=rem)
   qc=to_qiskit(c,meas);tqc=transpile(qc,basis_gates=['rz','sx','x','cx'],optimization_level=1,seed_transpiler=20260929)
   sim=AerSimulator(method='statevector',seed_simulator=20275936,enable_truncation=trunc,max_parallel_threads=1)
   res=sim.run(tqc,shots=1024).result()
   diag.append(dict(stage='after' if rem else 'before',truncation=trunc,counts=res.get_counts(),metadata=res.results[0].metadata))
 dump('sampling_diagnostic.json',diag)
 print('sampling diagnostic done',flush=True)
 out=[]
 for fam in ['beauregard','cuccaro','montgomery']:
  for rem in [False,True]:
   c,meas,inf=shor_circuit(fam,15,7,remove_identity=rem);qc=to_qiskit(c,meas)
   for level in [1,3]:
    tqc=transpile(qc,basis_gates=['rz','sx','x','cx'],optimization_level=level,seed_transpiler=20260929)
    out.append(dict(family=fam,stage='after' if rem else 'before',level=level,logical=metrics(qc),compiled=metrics(tqc)))
 dump('transpiler_sensitivity.json',out)
 # Explicit Clifford+T arithmetic-only accounting with QFT removed.
 from families import Cuccaro
 from revcirc import Circuit
 tstats=[]
 for rem in [False,True]:
  L=Cuccaro.layout(4,base=8);c=Circuit(8+L['nq'])
  for j in range(8):
   k=pow(7,2**j,15)
   if rem and k==1:continue
   Cuccaro.cu(c,L,k,15,j)
  qc=to_qiskit(c,[])
  tc=transpile(qc,basis_gates=['h','s','sdg','t','tdg','x','cx'],optimization_level=0)
  layers=[0]*tc.num_qubits
  for inst in tc.data:
   qs=[tc.find_bit(q).index for q in inst.qubits]
   lev=max([layers[q] for q in qs],default=0)+int(inst.operation.name in ['t','tdg'])
   for q in qs:layers[q]=lev
  tstats.append(dict(stage='after' if rem else 'before',ccx=c.count_ops().get('ccx',0),cswap=c.count_ops().get('cswap',0),T_count=tc.count_ops().get('t',0)+tc.count_ops().get('tdg',0),T_depth=max(layers)))
 dump('arithmetic_T_metrics.json',tstats)
 print('quantum metrics done',flush=True)
if __name__=='__main__':
 globals()[sys.argv[1]]()
