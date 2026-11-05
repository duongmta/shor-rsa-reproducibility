"""Run original scientific programs in a NEW working directory.

Examples: python reproduce.py tables
          python reproduce.py analytic
          python reproduce.py counts
Use --out runs/another-name for a later rerun. Existing directories are refused.
Archived evidence is never overwritten. See README.md for coverage and runtime.
"""
from pathlib import Path
import argparse, json, os, shutil, subprocess, sys

ROOT=Path(__file__).resolve().parent
MODES=['tables','analytic','quantum','aqft','counts','aer','functional','independent']
def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('mode',choices=MODES)
    ap.add_argument('--out',type=Path)
    args=ap.parse_args()
    out=(args.out or ROOT/'runs'/args.mode).resolve()
    if out.exists():raise SystemExit('Refusing to reuse existing output directory: '+str(out))
    # Prevent accidental output inside immutable evidence folders.
    for name in ['shor_reversible','revision_20261104','figures_manuscript']:
        if out.is_relative_to(ROOT/name):raise SystemExit('Choose a fresh directory under runs/ or outside the repository')
    out.mkdir(parents=True)
    for name in ['shor_reversible','revision_20261104']:
        (out/name).mkdir()
        for p in (ROOT/name).glob('*.py'):shutil.copy2(p,out/name/p.name)
    core=out/'shor_reversible';rev=out/'revision_20261104';logs=out/'logs';logs.mkdir()
    shutil.copy2(ROOT/'shor_reversible/logical_counts.csv',core/'logical_counts.csv')
    env=dict(os.environ,PYTHONIOENCODING='utf-8',PYTHONDONTWRITEBYTECODE='1')
    executions=[]
    def run(folder,arguments,logname):
        print('Running',arguments[0],*arguments[1:],flush=True)
        path=logs/logname
        with path.open('w',encoding='utf-8') as f:
            process=subprocess.run([sys.executable,'-B',*arguments],cwd=folder,env=env,stdout=f,stderr=subprocess.STDOUT)
        executions.append({'arguments':arguments,'log':str(path.relative_to(out)),'exit_code':process.returncode})
        if process.returncode:
            raise RuntimeError('Execution failed. See '+str(path))
        return path
    if args.mode=='tables':
        for folder in ['results','results_independent']:
            shutil.copytree(ROOT/'shor_reversible'/folder,core/folder)
        run(core,['paper_tables.py'],'tables.txt')
        for p in (core/'tables').glob('*.csv'):
            import csv
            old=list(csv.DictReader((ROOT/'shor_reversible/tables'/p.name).open(encoding='utf-8')))
            new=list(csv.DictReader(p.open(encoding='utf-8')))
            if old!=new:raise RuntimeError('Regenerated table differs: '+p.name)
    elif args.mode in ['analytic','quantum']:
        run(rev,['audit_experiments.py',args.mode],args.mode+'.txt')
    elif args.mode=='aqft':
        run(rev,['aqft_check.py'],'aqft.txt')
    elif args.mode=='counts':
        run(core,['run_qiskit.py','env'],'environment.txt')
        run(core,['run_qiskit.py','counts'],'counts.txt')
        old={r['key']:r for r in map(json.loads,(ROOT/'shor_reversible/results/counts.jsonl').read_text().splitlines())}
        new={r['key']:r for r in map(json.loads,(core/'results/counts.jsonl').read_text().splitlines())}
        if old.keys()!=new.keys():raise RuntimeError('Incomplete compiled configuration set')
        for k in old:
            for field in ['logical','transpiled','qubits','order','identity_blocks']:
                if old[k][field]!=new[k][field]:raise RuntimeError(k+' differs in '+field)
    elif args.mode=='aer':
        run(core,['run_qiskit.py','env'],'environment.txt')
        # Preserve the original six Table-5 method choices before resuming with MPS.
        run(core,['run_qiskit.py','simulate','--N','15','--bases','7','--which','both','--method','auto'],'table5.txt')
        run(core,['run_qiskit.py','simulate','--N','15','--which','both','--method','matrix_product_state'],'aer15.txt')
        run(core,['run_qiskit.py','simulate','--N','21','35','--method','matrix_product_state'],'aer21_35.txt')
    elif args.mode=='functional':
        run(core,['test_sim.py'],'simulator.txt')
        p=run(core,['verify.py','cuccaro'],'cuccaro.txt')
        if '44/44 blocks correct' not in p.read_text() or 'FAIL' in p.read_text():raise RuntimeError('Cuccaro verifier failed')
        p=run(core,['verify2.py','beauregard','montgomery'],'fourier_montgomery.txt')
        s=p.read_text()
        if s.count('44/44 blocks correct')!=2 or 'FAIL' in s:raise RuntimeError('Fourier/Montgomery verifier failed')
    elif args.mode=='independent':
        (core/'results_independent').mkdir()
        run(core,['e2e.py','beauregard','15','7','before'],'before_A15_a7.txt')
        for fam,N,label in [('beauregard',15,'A15'),('cuccaro',15,'B15'),('beauregard',21,'A21')]:
            p=run(core,['e2e_sweep.py',fam,str(N)],'local_'+label+'.jsonl')
            shutil.copy2(p,core/'results_independent'/p.name)
    report={'mode':args.mode,'executions':executions,'python':sys.version,
            'archived_results_overwritten':False,
            'note':'Successful execution alone does not imply identical sampled outcomes. Compare new outputs with archived evidence.'}
    (out/'execution_report.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print('Completed. New results:',out)

if __name__=='__main__':main()
