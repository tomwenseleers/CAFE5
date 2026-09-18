#!/usr/bin/env python3
"""Resumable full-family N11 model screening; predictive checks are not branch tests."""
import argparse,concurrent.futures,csv,hashlib,json,subprocess,sys,time
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__);p.add_argument('binary',type=Path);p.add_argument('--output',type=Path,required=True);p.add_argument('--counts',type=Path,required=True);p.add_argument('--tree',type=Path,required=True);p.add_argument('--workers',type=int,default=4);p.add_argument('--threads',type=int,default=2);p.add_argument('--models',nargs='*');p.add_argument('--iterations',type=int,default=1600);a=p.parse_args();a.binary=a.binary.resolve();a.output.mkdir(parents=True,exist_ok=True)
models={
 'baseline':(6,[]),
 'gamma12':(12,[]), 'gamma24':(24,[]),
 'asymmetric':(6,['--estimate-mu','--initial-mu','.02']),
 'zero_error0':(6,['--epsilon-zero','0']),
 'zero_error001':(6,['--epsilon-zero','.001']),
 'zero_error_free':(6,['--estimate-epsilon-zero','--epsilon-zero','.001']),
 'hurdle_poisson':(6,['--root-family','hurdle-poisson','--root-zero','.4','--estimate-root-zero']),
 'hurdle_nb':(6,['--root-family','hurdle-nb','--root-zero','.4','--estimate-root-zero','--root-shape','1','--estimate-root-shape']),
 'combined_poisson_zero0':(6,['--estimate-mu','--initial-mu','.02','--epsilon-zero','0','--root-family','hurdle-poisson','--root-zero','.4','--estimate-root-zero']),
 'combined_nb_zero_free':(6,['--estimate-mu','--initial-mu','.02','--estimate-epsilon-zero','--epsilon-zero','.001','--root-family','hurdle-nb','--root-zero','.4','--estimate-root-zero','--root-shape','1','--estimate-root-shape'])}
def read(prefix):
 with Path(str(prefix)+'_results.tsv').open() as h:return {r['parameter']:r['value'] for r in csv.DictReader(h,delimiter='\t')}
def one(item):
 name,(cats,extra)=item;folder=a.output/name;folder.mkdir(exist_ok=True);prefix=folder/'fit';done=folder/'audit.json'
 cmd=[str(a.binary),'--innovation','-t',str(a.tree),'-i',str(a.counts),'--root-mean','.1' if '--root-family' in extra else '.443237','--estimate-root-mean','--gamma-cats',str(cats),'--estimate-epsilon','--max-count','180','--threads',str(a.threads),'--starts','1','--iterations',str(a.iterations),'--skip-boundary-fits','--likelihood-only','--initial-lambda','.016782','--initial-nu','.0004609','--initial-alpha','.068416','--initial-epsilon','.015966']+extra+['-o',str(prefix)]
 manifest={'command':cmd,'binary_sha256':hashlib.sha256(a.binary.read_bytes()).hexdigest(),'counts_sha256':hashlib.sha256(a.counts.read_bytes()).hexdigest(),'tree_sha256':hashlib.sha256(a.tree.read_bytes()).hexdigest()}
 if done.exists():
  old=json.loads(done.read_text())
  if old['manifest']!=manifest:raise ValueError('Resume manifest differs for '+name)
  return name,old
 begun=time.time();print('START '+name,flush=True)
 with (folder/'fit.log').open('w') as h:
  h.write(' '.join(cmd)+'\n');h.flush();r=subprocess.run(cmd,stdout=h,stderr=h)
 result=read(prefix) if Path(str(prefix)+'_results.tsv').exists() else {}
 audit={'manifest':manifest,'returncode':r.returncode,'seconds':time.time()-begun,'fit':result,'scope':'Single-start screening with no boundary audit; not a final model comparison.'}
 if r.returncode==0:
  cmd=[sys.executable,str(Path(__file__).with_name('predictive_counts.py')),str(a.binary),'--tree',str(a.tree),'--fit',str(prefix),'--counts',str(a.counts),'--output',str(folder/'predictive'),'--replicates','4','--seed','930100']
  with (folder/'predictive.log').open('w') as h:rr=subprocess.run(cmd,stdout=h,stderr=h)
  audit['predictive_returncode']=rr.returncode
  if rr.returncode==0:audit['predictive']=json.loads((folder/'predictive/predictive_checks.json').read_text())
 done.write_text(json.dumps(audit,indent=2)+'\n');print('DONE '+name+' '+str(result.get('negative_log_likelihood'))+' status='+str(r.returncode),flush=True);return name,audit
selected={k:v for k,v in models.items() if not a.models or k in a.models}
if a.models and set(a.models)-set(models):raise ValueError('Unknown model')
summary={}
with concurrent.futures.ThreadPoolExecutor(max_workers=a.workers) as pool:
 for f in concurrent.futures.as_completed([pool.submit(one,x) for x in selected.items()]):
  name,result=f.result();summary[name]=result;(a.output/'screening_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
