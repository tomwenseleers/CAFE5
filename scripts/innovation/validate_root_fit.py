#!/usr/bin/env python3
"""Check estimated root law against simulation and independent fixed-law evaluation."""
import argparse,csv,json,subprocess
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('binary',type=Path);p.add_argument('output',type=Path);a=p.parse_args();a.output.mkdir(exist_ok=True,parents=True);a.binary=a.binary.resolve()
def run(args):
 result=subprocess.run(list(map(str,[a.binary,'--innovation']+args)),capture_output=True,text=True)
 if result.returncode:raise RuntimeError(result.stdout+result.stderr)
def results(prefix):
 with Path(str(prefix)+'_results.tsv').open() as f:return {r['parameter']:r['value'] for r in csv.DictReader(f,delimiter='\t')}
tree=a.output/'tree.nwk';tree.write_text('((A:0.3,B:0.3):0.4,C:0.7);\n')
common=['-t',tree,'--lambda',.2,'--nu',.3,'--gamma-cats',3,'--alpha',1.3,'--epsilon',.08,'--max-count',35,'--threads',2]
prefix=a.output/'sample';run(common+['--root-mean',1.4,'--simulate',4000,'--seed',61841,'-o',prefix]);data=str(prefix)+'_simulated.tsv'
fit=a.output/'fitted';run(common+['-i',data,'--root-mean',.5,'--estimate-root-mean','--starts',2,'--iterations',500,'--likelihood-only','-o',fit]);r=results(fit);mean=float(r['root_mean'])
assert abs(mean-1.4)<.12,mean
fixed=a.output/'fixed';run(common+['-i',data,'--root-mean',mean,'--likelihood-only','-o',fixed]);rr=results(fixed)
assert abs(float(r['negative_log_likelihood'])-float(rr['negative_log_likelihood']))<1e-9
assert r['truncation_pass']=='1' and r['root_mean_estimated']=='1'
summary={'generating_root_mean':1.4,'estimated_root_mean':mean,'families':4000,'fixed_law_likelihood_difference':float(rr['negative_log_likelihood'])-float(r['negative_log_likelihood']),'limitations':'Root-only recovery with other parameters fixed; does not establish five-parameter identifiability.'}
(a.output/'validation.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
