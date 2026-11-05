"""Validate release files, archived evidence and freshly generated circuits.

python check_release.py                 # no Aer simulations, no file changes
python check_release.py --compiled      # also recompile all 246 circuits
python check_release.py --report PATH   # optionally save the audit result
"""
from pathlib import Path
from math import gcd, log, ceil, isclose
from collections import Counter
import argparse, ast, csv, hashlib, json, re, sys

ROOT = Path(__file__).resolve().parent
SRC = ROOT/'shor_reversible'
REV = ROOT/'revision_20261104'
sys.path.insert(0,str(SRC))

def read_json(p): return json.loads(p.read_text(encoding='utf-8'))
def records(p): return [json.loads(s) for s in p.read_text(encoding='utf-8').splitlines() if s.strip()]
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def require(value, message):
    if not value: raise RuntimeError(message)
def number(s): return float(s.replace(',','').strip())

def audit(compiled=False, development=False):
    import numpy as np
    from analysis import order, identity_blocks, success_probability, recovers_factor, wilson, compare_counts
    from qpe import shor_circuit
    from families import FAMILIES
    report={'status':'RUNNING','fresh_compiled_check':compiled,'checks':{}}
    checks=report['checks']
    if not development:
        manifest=read_json(ROOT/'MANIFEST_SHA256.json')
        for name,digest in manifest.items():
            require((ROOT/name).is_file() and sha(ROOT/name)==digest,'File integrity: '+name)
        checks['manifest_files']=len(manifest)
    for row in read_json(ROOT/'SOURCE_PROVENANCE.json'):
        require(sha(ROOT/row['file'])==row['sha256'],'Original-source identity: '+row['file'])
    checks['original_files_unchanged']=len(read_json(ROOT/'SOURCE_PROVENANCE.json'))
    for p in list(SRC.glob('*.py'))+list(REV.glob('*.py'))+list(ROOT.glob('*.py')):
        ast.parse(p.read_text(encoding='utf-8'),filename=p.name)
    bases={(N,a) for N in [15,21,35] for a in range(2,N) if gcd(a,N)==1}
    counts=records(SRC/'results/counts.jsonl')
    samples=records(SRC/'results/simulate.jsonl')
    families=['beauregard','cuccaro','montgomery']
    expected={f'{f}|{N}|{a}|{s}' for N,a in bases for f in families for s in ['before','after']}
    require(len(counts)==246 and {r['key'] for r in counts}==expected,'246 compiled configurations')
    sample_keys={f'{f}|{N}|{a}|{s}' for N,a in bases for f in families for s in (['before','after'] if N==15 else ['after'])}
    require(len(samples)==144 and {r['key'] for r in samples}==sample_keys,'144 sampled configurations')
    checks['bases_by_modulus']={N:sorted(a for NN,a in bases if NN==N) for N in [15,21,35]}
    checks['compiled_record_count']=len(counts);checks['aer_record_count']=len(samples)
    logical=list(csv.DictReader((SRC/'logical_counts.csv').open(encoding='utf-8')))
    logical_map={(r['family'],int(r['N']),int(r['a'])):r for r in logical}
    require(len(logical_map)==len(logical)==123,'123 paired logical records')
    circuit_digests={}
    if compiled:
        from run_qiskit import to_qiskit, metrics, transpile_fixed
        import qiskit,qiskit_aer
        report['environment']={'python':sys.version.split()[0],'numpy':np.__version__,
                               'qiskit':qiskit.__version__,'qiskit_aer':qiskit_aer.__version__}
    for i,r in enumerate(counts):
        c,m,info=shor_circuit(r['family'],r['N'],r['a'],remove_identity=r['stage']=='after')
        original=logical_map[(r['family'],r['N'],r['a'])]
        require(c.nq==r['qubits']==int(original['qubits']),r['key']+' qubits')
        require(len(c.ops)==int(original[r['stage']+'_gates']),r['key']+' logical gates')
        require(c.depth()==int(original[r['stage']+'_depth']),r['key']+' logical depth')
        require(len(c.ops)+len(m)==r['logical']['size'],r['key']+' measurements')
        require(dict(Counter('mcphase' if x[0]=='ccp' else x[0] for x in c.ops)+Counter({'measure':len(m)}))==r['logical']['ops'],r['key']+' operations')
        require(order(r['a'],r['N'])==r['order'],r['key']+' order')
        require(identity_blocks(r['N'],r['a'],info['t'])==r['identity_blocks'],r['key']+' identities')
        circuit_digests[r['key']]=hashlib.sha256(json.dumps([c.nq,c.ops,m],separators=(',',':')).encode()).hexdigest()
        if compiled:
            qc=to_qiskit(c,m)
            require(metrics(qc)==r['logical'],r['key']+' Qiskit logical metrics')
            actual=metrics(transpile_fixed(qc))
            require(actual==r['transpiled'],r['key']+' compiled metrics')
        if (i+1)%30==0:print('Checked',i+1,'/ 246 circuits',flush=True)
    checks['fresh_logical_circuits']=246
    if compiled:checks['fresh_compiled_circuits']=246
    report['circuit_ir_sha256']=circuit_digests
    analytic={}
    for N,a in sorted(bases):
        t=2*N.bit_length();Q=2**t
        p,P,acc=success_probability(N,a,t)
        pr=sum(recovers_factor(y,Q,a,N) for y in range(Q))/Q
        analytic[(N,a)]=(p,pr,P)
    for N,a,num in [(15,4,123),(21,8,523),(35,6,2337)]:
        require(analytic[(N,a)][1]*2**(2*N.bit_length())==num,'Order-two acceptance numerator')
    for r in samples:
        p,pr,P=analytic[(r['N'],r['a'])]
        shots=r['shots'];Q=len(P);hist={int(k):v for k,v in r['counts'].items()}
        require(sum(hist.values())==shots==1024,r['key']+' shots')
        k=sum(v for y,v in hist.items() if recovers_factor(y,Q,r['a'],r['N']))
        require(k==r['success_count'] and k/shots==r['p_hat'],r['key']+' recovery')
        require(np.allclose(wilson(k,shots),r['wilson95'],atol=1e-14,rtol=0),r['key']+' Wilson interval')
        for name,value in compare_counts(hist,P,shots).items():
            require(isclose(value,r[name],abs_tol=1e-12,rel_tol=1e-10) or (np.isnan(value) and np.isnan(r[name])),r['key']+' '+name)
    checks['fresh_analytic_bases']=41;checks['aer_histograms_recalculated']=144
    context=read_json(ROOT/'MANUSCRIPT_CONTEXT.json');tables=context['tables']
    for row,f in zip(tables[0][1:],families):
        require([int(x) for x in row[3:]]==[2*N.bit_length()+FAMILIES[f].layout(N.bit_length())['nq'] for N in [15,21,35]],'Table 1')
    for row in tables[1][1:]:
        N=int(row[0])
        for a in map(int,row[1].split(',')):
            p,pr,_=analytic[(N,a)]
            require(order(a,N)==int(row[2]) and identity_blocks(N,a,2*N.bit_length())==int(row[3]),'Table 2 r/i')
            require(abs(p-number(row[4]))<=0.00005 and abs(pr-number(row[5]))<=0.00005,'Table 2 probabilities')
            if pr:
                require(abs(p/pr-number(row[6]))<=0.05,'Table 2 G')
                require(ceil(log(.05)/log(1-p))==int(row[7]),'Table 2 m95')
    for row in tables[2][1:]:
        N,r,ids,bcount=map(int,row[:4])
        for f,cell,red in zip(families,row[4:7],row[7].split('/')):
            group=[x for x in logical if x['family']==f and int(x['N'])==N and int(x['order'])==r and int(x['identity_blocks'])==ids]
            require(len(group)==bcount,'Table 3 bases')
            before=np.mean([int(x['before_gates']) for x in group]);after=np.mean([int(x['after_gates']) for x in group])
            parts=cell.split('→')
            require(number(parts[0])==round(before) and number(parts[-1])==round(after),'Table 3 counts')
            require(abs(100*(1-after/before)-number(red))<=0.05,'Table 3 reduction')
    cmap={r['key']:r for r in counts};smap={r['key']:r for r in samples};fammap=dict(zip('ABC',families))
    for row in tables[3][1:]:
        r=cmap[f'{fammap[row[0]]}|15|7|{row[1]}'];L=r['logical'];T=r['transpiled']
        require(list(map(number,row[2:]))==[r['qubits'],L['size'],T['size'],T['cx'],T['depth'],T['rz_t_like'],T['rz_arbitrary']],'Table 4')
    for row in tables[4][1:]:
        r=smap[f'{fammap[row[0]]}|15|7|{row[1]}']
        require([int(x) for x in row[3].split('/')]==[r['counts'].get(str(y),0) for y in [0,64,128,192]],'Table 5 histogram')
        require(row[2]=={'statevector':'SV','matrix_product_state':'MPS'}[r['method']],'Table 5 method')
        displayed=list(map(float,re.findall(r'\d+\.\d+',row[4])))
        require(len(displayed)==3 and abs(displayed[0]-r['p_hat'])<=.00005,'Table 5 acceptance')
        require(all(abs(a-b)<=.0005 for a,b in zip(displayed[1:],r['wilson95'])),'Table 5 interval')
        require(abs(number(row[5])-r['tvd'])<=.00005+1e-12,'Table 5 TVD')
        chi=list(map(float,re.findall(r'\d+\.\d+',row[6])))
        require(len(chi)==2 and abs(chi[0]-r['chi2'])<=.005+1e-12 and abs(chi[1]-r['chi2_p'])<=.0005+1e-12,'Table 5 chi-square')
        require(int(row[7])==r['off_support'],'Table 5 off-support')
        require(abs(number(row[-1])-r['wall_s'])<=0.050001,'Table 5 historical time')
    large=read_json(REV/'postprocessing_large.json');large_map={(r['N'],r['a'],r['L']):r for r in large}
    require(len(large)==len(large_map)==18,'18 large-modulus rule records')
    for row in tables[5][1:]:
        N,a,r=map(int,row[:3]);l1=large_map[(N,a,1)];l8=large_map[(N,a,8)]
        require(r==order(a,N)==l1['r']==l8['r'],'Table 6 order')
        require(list(map(number,row[3:]))==[l1['p_hat'],l8['p_hat'],l1['p_rand_hat'],l8['p_rand_hat']],'Table 6 values')
    for r in large:
        require(r['shots']==20000,'Large N sample size')
        for k,p,ci in [('k','p_hat','ci'),('random_k','p_rand_hat','random_ci')]:
            require(r[p]==r[k]/r['shots'],'Large N sampled rate')
            require(np.allclose(wilson(r[k],r['shots']),r[ci],atol=1e-14,rtol=0),'Large N Wilson interval')
            require((r[ci][1]-r[ci][0])/2<=.007,'Wilson half-width bound')
    checks['tables_matched_current_manuscript']=[1,2,3,4,5,6]
    small=read_json(REV/'postprocessing_small.json')
    require(len(small)==246,'Small sensitivity records')
    for r in small:
        if r['L']==1 and r['t']==2*r['N'].bit_length():
            p,pr,_=analytic[(r['N'],r['a'])]
            require(abs(p-r['p'])<1e-12 and abs(pr-r['p_rand'])<1e-12,'Small sensitivity original-rule match')
    independent=[]
    for name in ['local_A15.jsonl','local_A21.jsonl','local_B15.jsonl']:independent+=records(SRC/'results_independent'/name)
    require(len(independent)==len({(r['family'],r['N'],r['a']) for r in independent})==25,'25 independent after configurations')
    require(max(r['max_dev'] for r in independent)<1.5e-13,'Independent marginal errors')
    checks['independent_distinct_after_configurations']=25
    checks['independent_before_check']='Preserved in logs/05_e2e_A15_a7_before.txt; not rerun by this script'
    checks['full_aer_or_independent_statevector_sweep_rerun']=False
    checks['large_analytic_sampling_rerun']=False
    for name in ['transpiler_sensitivity.json','sampling_diagnostic.json','arithmetic_T_metrics.json','aqft_results.json']:
        require(bool(read_json(REV/name)),name)
    report['status']='PASS'
    return report

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--compiled',action='store_true')
    ap.add_argument('--report',type=Path)
    ap.add_argument('--development',action='store_true',help=argparse.SUPPRESS)
    args=ap.parse_args()
    report=audit(args.compiled,args.development)
    if args.report:
        args.report.parent.mkdir(parents=True,exist_ok=True)
        args.report.write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k!='circuit_ir_sha256'},indent=2))
