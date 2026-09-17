#!/usr/bin/env python3
"""Empirical operating characteristics of branch-tail tests under gamma + error."""
import argparse,csv,json,subprocess,sys
from pathlib import Path
import numpy as np
p=argparse.ArgumentParser();p.add_argument('binary',type=Path);p.add_argument('output',type=Path);p.add_argument('--outer-replicates',type=int,default=4);a=p.parse_args();a.binary=a.binary.resolve();a.output.mkdir(parents=True,exist_ok=True)
driver=Path(__file__).with_name('branch_bootstrap.py')
tree=a.output/'tree.nwk';tree.write_text('((A:.3,B:.3):.5,(C:.5,D:.5):.3);\n')
fixed=['--lambda',.2,'--nu',.3,'--gamma-cats',2,'--alpha',1.3,'--epsilon',.08]
common=['-t',tree,'--root-mean',1,'--max-count',45]
def run(cmd,name):
 r=subprocess.run(list(map(str,cmd)),capture_output=True,text=True)
 (a.output/f'{name}.log').write_text(r.stdout+r.stderr)
 if r.returncode:raise RuntimeError(f'{name} failed; inspect its log')
def test_summary(directory):
 with (directory/'branch_tests.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
 branches=list(dict.fromkeys(r['Node'] for r in rows))
 rates={b:float(np.mean([float(r['p_two_sided'])<=.05 for r in rows if r['Node']==b])) for b in branches}
 return {'fraction_p_le_005_by_branch':rates,'mean_fraction':float(np.mean(list(rates.values()))),'family_branch_tests':len(rows)}
# Large independent sample under known generating parameters.
prefix=a.output/'known_sample';run([a.binary,'--innovation']+common+fixed+['--simulate',2000,'--seed',8871,'-o',prefix],'known_simulate')
known=a.output/'known_tests';run([sys.executable,driver,a.binary]+common+fixed+['-i',str(prefix)+'_simulated.tsv','-o',known,'--fixed-parameters','--replicates',3999,'--seed',5517],'known_test')
known_summary=test_summary(known);assert max(known_summary['fraction_p_le_005_by_branch'].values())<.08,known_summary
print('Known-parameter calibration:',known_summary,flush=True)
# Full refitting pipeline, estimating lambda, nu, alpha and epsilon in every
# bootstrap data set. Small outer count is explicitly reported, never disguised
# as thousands of independent calibration data sets.
outer=[]
for j in range(a.outer_replicates):
 prefix=a.output/f'outer_{j}';run([a.binary,'--innovation']+common+fixed+['--simulate',300,'--seed',4400+j,'-o',prefix],f'outer_simulate_{j}')
 out=a.output/f'outer_tests_{j}'
 run([sys.executable,driver,a.binary]+common+['--gamma-cats',2,'--estimate-epsilon','--starts',1,'--iterations',1500,'-i',str(prefix)+'_simulated.tsv','-o',out,'--replicates',39,'--seed',19000+100*j],f'outer_test_{j}')
 summary=test_summary(out);summary['outer_replicate']=j;outer.append(summary);print(summary,flush=True)
# Conservative smoke criterion only: no claim of precise nominal 5% control.
mean=float(np.mean([r['mean_fraction'] for r in outer]));assert mean<.09,mean
summary={'known_parameters':known_summary,'refitted_outer_replicates':outer,'refitted_mean_fraction_p_le_005':mean,'bootstrap_datasets_per_outer_replicate':39,'families_per_outer_dataset':300,'estimated_parameters':['lambda','nu','alpha','epsilon'],'limitation':'Limited calibration at one generating regime. Only four independent outer datasets by default; cannot certify type-I error across all parameters or the empirical wasp analysis. Coarse Monte Carlo resolution 2/(39+1)=0.05.'}
(a.output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
