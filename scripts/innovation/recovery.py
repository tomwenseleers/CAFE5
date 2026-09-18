#!/usr/bin/env python3
"""Seeded parameter recovery and fixed-parameter tail calibration."""
import argparse,csv,json,subprocess
from pathlib import Path
import numpy as np
p=argparse.ArgumentParser();p.add_argument('binary',type=Path);p.add_argument('output',type=Path);p.add_argument('--replicates',type=int,default=12);a=p.parse_args();a.binary=a.binary.resolve();a.output.mkdir(parents=True,exist_ok=True)
def run(args):
 r=subprocess.run([str(a.binary),'--innovation']+list(map(str,args)),capture_output=True,text=True)
 if r.returncode:raise RuntimeError((args,r.returncode,r.stdout,r.stderr))
def results(prefix):
 with Path(str(prefix)+'_results.tsv').open() as f:return {x['parameter']:x['value'] for x in csv.DictReader(f,delimiter='\t')}
tree=a.output/'recovery.nwk';tree.write_text('(((A:0.4,B:0.4):0.5,(C:0.6,D:0.6):0.3):0.6,(E:0.7,F:0.7):0.8);\n')
common=['-t',tree,'--root-mean',1,'--max-count',45]
records=[]
for i in range(a.replicates):
 prefix=a.output/f'rep{i:02d}'
 run(common+['--lambda',.2,'--nu',.3,'--simulate',2500,'--seed',7100+i,'-o',prefix])
 run(common+['-i',str(prefix)+'_simulated.tsv','--starts',2,'--threads',1,'-o',prefix])
 r=results(prefix);record={'replicate':i,'lambda':float(r['lambda']),'nu':float(r['nu']),'truncation_delta':float(r['truncation_nll_difference'])};records.append(record)
 print(json.dumps(record),flush=True)
with (a.output/'recovery.tsv').open('w') as f:w=csv.DictWriter(f,fieldnames=list(records[0]),delimiter='\t');w.writeheader();w.writerows(records)
means={x:float(np.mean([r[x] for r in records])) for x in ['lambda','nu']}
assert abs(means['lambda']-.2)<.025 and abs(means['nu']-.3)<.04,means
# Fixed-parameter calibration on an independent observed-family sample.
prefix=a.output/'calibration';run(common+['--lambda',.2,'--nu',.3,'--simulate',1500,'--seed',81001,'-o',prefix])
run(common+['-i',str(prefix)+'_simulated.tsv','--lambda',.2,'--nu',.3,'--bootstrap',5000,'--seed',92001,'-o',prefix])
with Path(str(prefix)+'_experimental_pvalues.tsv').open() as f:pvals=np.array([float(r['fixed_parameter_MC_p']) for r in csv.DictReader(f,delimiter='\t')])
rate=float(np.mean(pvals<=.05));assert .025<rate<.075,rate
summary={'replicates':a.replicates,'families_per_replicate':2500,'true_lambda':.2,'true_nu':.3,'mean_estimates':means,'fixed_parameter_null_fraction_p_le_005':rate,'calibration_families':len(pvals),'calibration_null_simulations':5000,'limitation':'One generating regime; fixed parameters only. Does not establish calibration after fitting parameters to the tested data or branch-specific significance.'}
(a.output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
