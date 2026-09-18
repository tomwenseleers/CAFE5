#!/usr/bin/env python3
"""Independent transition-score check and null operating characteristics at 1% and 5%."""
import argparse,csv,json,subprocess
from pathlib import Path
import numpy as np
from scipy.linalg import expm
p=argparse.ArgumentParser();p.add_argument('binary',type=Path);p.add_argument('output',type=Path);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True);a.binary=a.binary.resolve()
def read(path):
 with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def run(args):subprocess.run(list(map(str,args)),check=True,capture_output=True,text=True)
tree=a.output/'tree.nwk';tree.write_text('((A:0.3,B:0.3):0.4,C:0.7);\n');counts=a.output/'counts.tsv';counts.write_text('Desc\tFamily ID\tA\tB\tC\nx\tf1\t0\t1\t2\nx\tf2\t1\t0\t0\n')
common=[a.binary,'--innovation','-t',tree,'--root-mean',1,'--lambda',.2,'--nu',.3,'--gamma-cats',3,'--alpha',1.3,'--epsilon',.12,'--max-count',35,'--threads',2]
prefix=a.output/'check';run(common+['-i',counts,'-o',prefix]);rates=[float(r['lambda_multiplier']) for r in read(Path(str(prefix)+'_categories.tsv'))];families={r['Family ID']:r for r in read(Path(str(prefix)+'_families.tsv'))}
Ps=[]
for g in rates:
 Q=np.zeros((161,161))
 for n in range(161):
  if n<160:Q[n,n+1]=.2*g*n+.3
  if n:Q[n,n-1]=.2*g*n
  Q[n,n]=-Q[n].sum()
 Ps.append({b:expm(Q*t) for b,t in [('Node1',.4),('A',.3),('B',.3),('C',.7)]})
errors=[]
for r in read(Path(str(prefix)+'_branch_statistics.tsv')):
 i,j=int(r['parent_MAP_count']),int(r['child_MAP_count']);row=families[r['Family ID']]
 lo=sum(float(row[f'P_category_{k}'])*Ps[k][r['Node']][i,:j+1].sum() for k in range(3));hi=sum(float(row[f'P_category_{k}'])*Ps[k][r['Node']][i,j:].sum() for k in range(3))
 expected=-np.log(min(1,2*min(lo,hi)));errors.append(abs(expected-float(r['transition_tail_score'])))
assert max(errors)<2e-12,errors
obs=a.output/'observed';run(common+['--simulate',2000,'--seed',1701,'-o',obs]);run(common+['-i',str(obs)+'_simulated.tsv','-o',obs])
out=a.output/'reference';driver=Path(__file__).with_name('refitted_bootstrap.py')
run(['python3',driver,a.binary,'--tree',tree,'--observed',obs,'--output',out,'--replicates',8,'--workers',2,'--threads',1,'--seed',1741,'--nodes','Node1','C','--minimum-cap',35,'--fixed'])
r=read(out/'branch_tests.tsv');summary={'independent_transition_score_max_abs_error':max(errors),'heldout_null_families':2000,'reference_families':16000,'branch_threshold':.01,'branch_rejection_rates':{b:sum(float(x['transition_branch_p'])<.01 for x in r if x['Node']==b)/2000 for b in ['Node1','C']},'family_rejection_rate':sum(float(x['family_p'])<.05 for x in r if x['Node']=='Node1')/2000,'limitations':'Known-parameter check at one regime; not certification of estimated-parameter calibration or model adequacy.'}
summary['branch_rejection_rates_at_05']={b:sum(float(x['transition_branch_p'])<.05 for x in r if x['Node']==b)/2000 for b in ['Node1','C']}
assert all(v<.025 for v in summary['branch_rejection_rates'].values()),summary
assert all(v<.075 for v in summary['branch_rejection_rates_at_05'].values()),summary
assert summary['family_rejection_rate']<.075,summary
(a.output/'validation.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
