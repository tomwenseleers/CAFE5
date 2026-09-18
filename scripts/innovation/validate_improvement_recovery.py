#!/usr/bin/env python3
"""Conditional parameter recovery for new fit dimensions, with independent likelihood replay."""
import argparse,csv,json,subprocess
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('binary',type=Path);p.add_argument('output',type=Path);a=p.parse_args();a.binary=a.binary.resolve();a.output.mkdir(parents=True,exist_ok=True)
def run(args):
 r=subprocess.run([str(a.binary),'--innovation']+list(map(str,args)),capture_output=True,text=True)
 if r.returncode:raise RuntimeError(r.stderr)
def read(prefix):
 with open(str(prefix)+'_results.tsv') as h:return {r['parameter']:r['value'] for r in csv.DictReader(h,delimiter='\t')}
tree=a.output/'tree.nwk';tree.write_text('((A:0.3,B:0.3):0.4,C:0.7);\n')
base={'--lambda':.2,'--mu':.4,'--nu':.3,'--epsilon':.07,'--epsilon-zero':.02,'--root-family':'hurdle-nb','--root-mean':.8,'--root-zero':.3,'--root-shape':1.7}
common=['-t',tree,'--max-count',40,'--threads',2]
def args(d):return [v for k,x in d.items() for v in [k,x]]
sample=a.output/'sample';run(common+args(base)+['--simulate',16000,'--seed',20260918,'-o',sample]);data=str(sample)+'_simulated.tsv'
checks={}
for name,flag,estimate,start,tolerance in [('mu','--mu','--estimate-mu',.2,.035),('epsilon_zero','--epsilon-zero','--estimate-epsilon-zero',.008,.01),('root_zero','--root-zero','--estimate-root-zero',.5,.04),('root_shape','--root-shape','--estimate-root-shape',.8,.35)]:
 d=base.copy()
 if name=='mu':del d[flag];d['--initial-mu']=start
 else:d[flag]=start
 prefix=a.output/name;run(common+args(d)+[estimate,'-i',data,'--starts',1,'--iterations',1200,'--skip-boundary-fits','--likelihood-only','-o',prefix]);r=read(prefix);value=float(r[name]);truth=base[flag]
 assert abs(value-truth)<tolerance,(name,value,truth)
 d=base.copy();d[flag]=value;replay=a.output/(name+'_replay');run(common+args(d)+['-i',data,'--starts',1,'--likelihood-only','-o',replay]);rr=read(replay)
 delta=float(rr['negative_log_likelihood'])-float(r['negative_log_likelihood']);assert abs(delta)<1e-8
 checks[name]={'truth':truth,'estimate':value,'replay_nll_difference':delta,'truncation_pass':r['truncation_pass']}
checks['scope']='One parameter estimated at a time, others known; validates parameter plumbing, not joint identifiability.'
(a.output/'validation.json').write_text(json.dumps(checks,indent=2)+'\n');print(json.dumps(checks,indent=2))
