#!/usr/bin/env python3
"""Model-based branch tail tests. Default: whole-dataset parametric bootstrap with refitting.

Statistic: posterior E[N_child - N_parent | all observed counts], marginalized
across root counts, gamma categories and observation error. The null is the
specified global BDI population of observed families, NOT 'no change on this edge'.
"""
import argparse,csv,json,subprocess
from pathlib import Path
import numpy as np

def read_results(prefix):
 with Path(str(prefix)+'_results.tsv').open() as f:return {r['parameter']:r['value'] for r in csv.DictReader(f,delimiter='\t')}
def read_stats(prefix):
 with Path(str(prefix)+'_branch_statistics.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
 families=list(dict.fromkeys(r['Family ID'] for r in rows));branches=list(dict.fromkeys(r['Node'] for r in rows))
 values=np.array([float(r['posterior_mean_change']) for r in rows]).reshape(len(families),len(branches))
 assert all(r['Family ID']==families[i//len(branches)] and r['Node']==branches[i%len(branches)] for i,r in enumerate(rows))
 return rows,values

def adjust(p):
 shape=p.shape;p=p.ravel();m=len(p);order=np.argsort(p);ranked=p[order]
 harmonic=np.sum(1/np.arange(1,m+1))
 by=np.minimum(1,np.minimum.accumulate((ranked*m*harmonic/np.arange(1,m+1))[::-1])[::-1])
 holm=np.minimum(1,np.maximum.accumulate(ranked*(m-np.arange(m))))
 out_by=np.empty(m);out_holm=np.empty(m);out_by[order]=by;out_holm[order]=holm
 return out_by.reshape(shape),out_holm.reshape(shape)

def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('binary',type=Path);p.add_argument('-t','--tree',type=Path,required=True);p.add_argument('-i','--infile',type=Path,required=True);p.add_argument('-o','--output',type=Path,required=True)
 root=p.add_mutually_exclusive_group(required=True);root.add_argument('--root-mean',type=float);root.add_argument('--root-prior',type=Path)
 p.add_argument('--gamma-cats',type=int,default=1);p.add_argument('--alpha',type=float)
 error=p.add_mutually_exclusive_group();error.add_argument('--epsilon',type=float);error.add_argument('--estimate-epsilon',action='store_true');error.add_argument('--error-model',type=Path)
 p.add_argument('--lambda',dest='lam',type=float);p.add_argument('--nu',type=float)
 p.add_argument('--max-count',type=int,default=80);p.add_argument('--starts',type=int,default=2);p.add_argument('--iterations',type=int,default=800);p.add_argument('--threads',type=int,default=1)
 p.add_argument('--replicates',type=int,default=199);p.add_argument('--seed',type=int,default=7142)
 p.add_argument('--fixed-parameters',action='store_true',help='Known/fixed-parameter tests: requires every model parameter fixed explicitly; replicates counts simulated families, no refitting.')
 p.add_argument('--unconditioned',action='store_true')
 a=p.parse_args();a.binary=a.binary.resolve();a.output.mkdir(parents=True,exist_ok=True)
 if a.replicates<19:raise ValueError('At least 19 bootstrap replicates required; use >=999 for final testing.')
 if a.fixed_parameters and (a.lam is None or a.nu is None or (a.gamma_cats>1 and a.alpha is None) or a.estimate_epsilon):
  raise ValueError('Fixed-parameter tests require fixed lambda, nu, alpha (if gamma), and an unfitted error law.')
 common=['-t',a.tree,'--gamma-cats',a.gamma_cats,'--max-count',a.max_count,'--threads',a.threads,'--iterations',a.iterations,'--starts',a.starts]
 common+=['--root-mean',a.root_mean] if a.root_mean is not None else ['--root-prior',a.root_prior]
 if a.unconditioned:common+=['--unconditioned']
 original=[]
 if a.lam is not None:original+=['--lambda',a.lam]
 if a.nu is not None:original+=['--nu',a.nu]
 if a.alpha is not None:original+=['--alpha',a.alpha]
 if a.epsilon is not None:original+=['--epsilon',a.epsilon]
 if a.estimate_epsilon:original+=['--estimate-epsilon']
 if a.error_model is not None:common+=['--error-model',a.error_model]
 def run(args,name):
  command=[str(a.binary),'--innovation']+list(map(str,args))
  result=subprocess.run(command,capture_output=True,text=True)
  (a.output/f'{name}.log').write_text(result.stdout+result.stderr)
  if result.returncode:raise RuntimeError(f'{name} failed with status {result.returncode}; no failed replicate is discarded. See its log.')
 observed=a.output/'observed';run(common+original+['-i',a.infile,'-o',observed],'observed')
 r=read_results(observed);rows,values=read_stats(observed);n,branches=values.shape
 fixed=['--lambda',r['lambda'],'--nu',r['nu']]
 if a.gamma_cats>1:fixed+=['--alpha',r['alpha']]
 if a.error_model is None:fixed+=['--epsilon',r['epsilon']]
 upper=np.zeros_like(values,dtype=np.int64);lower=upper.copy();fit_records=[]
 if a.fixed_parameters:
  prefix=a.output/'null';run(common+fixed+['--simulate',a.replicates,'--seed',a.seed,'-o',prefix],'null_simulation')
  run(common+fixed+['-i',str(prefix)+'_simulated.tsv','-o',prefix],'null_reconstruction')
  _,null=read_stats(prefix)
  for b in range(branches):
   ordered=np.sort(null[:,b]);upper[:,b]=len(ordered)-np.searchsorted(ordered,values[:,b],side='left');lower[:,b]=np.searchsorted(ordered,values[:,b],side='right')
 else:
  for b in range(a.replicates):
   prefix=a.output/f'replicate_{b:05d}'
   run(common+fixed+['--simulate',n,'--seed',a.seed+b,'-o',prefix],f'simulate_{b:05d}')
   run(common+original+['-i',str(prefix)+'_simulated.tsv','-o',prefix],f'refit_{b:05d}')
   _,null=read_stats(prefix)
   if null.shape!=values.shape:raise ValueError('Bootstrap shape mismatch')
   # Corresponding exchangeable family indices give B independent replicate
   # statistics per hypothesis, not falsely B*n independent refitted datasets.
   upper+=null>=values;lower+=null<=values
   rr=read_results(prefix);fit_records.append({'replicate':b,**{k:rr[k] for k in ['lambda','nu','alpha','epsilon','optimizer_converged','truncation_pass']}})
   print(f'Completed refit {b+1}/{a.replicates}',flush=True)
 p_upper=(1+upper)/(a.replicates+1);p_lower=(1+lower)/(a.replicates+1);p_two=np.minimum(1,2*np.minimum(p_upper,p_lower))
 by,holm=adjust(p_two)
 fields=list(rows[0])+['p_upper','p_lower','p_two_sided','q_BY_all_family_branches','p_Holm_all_family_branches','bootstrap_replicates','method']
 method='known_parameter_MC' if a.fixed_parameters else 'refitted_parametric_bootstrap'
 with (a.output/'branch_tests.tsv').open('w') as f:
  writer=csv.DictWriter(f,fieldnames=fields,delimiter='\t');writer.writeheader()
  for i,row in enumerate(rows):
   j,k=divmod(i,branches);writer.writerow({**row,'p_upper':p_upper[j,k],'p_lower':p_lower[j,k],'p_two_sided':p_two[j,k],'q_BY_all_family_branches':by[j,k],'p_Holm_all_family_branches':holm[j,k],'bootstrap_replicates':a.replicates,'method':method})
 if fit_records:
  with (a.output/'bootstrap_fits.tsv').open('w') as f:w=csv.DictWriter(f,fieldnames=list(fit_records[0]),delimiter='\t');w.writeheader();w.writerows(fit_records)
 metadata={'method':method,'families':n,'branches':branches,'replicates':a.replicates,'seed':a.seed,'minimum_two_sided_p':2/(a.replicates+1),'null':'An exchangeable observed family from the fitted global BDI+gamma+error population on this branch. Not absence of any gene-count change.','statistic':'Posterior mean child count minus posterior mean parent count.','calibration':'Known-parameter rank tests are finite-sample conservative under the specified generative model. Refitted plug-in bootstrap tests are approximate and require empirical operating-characteristic checks.','multiple_testing':'BY FDR and Holm FWER across all supplied family-branch hypotheses; small B can make discovery impossible.','parameters':r}
 (a.output/'method.json').write_text(json.dumps(metadata,indent=2)+'\n')
if __name__=='__main__':main()
