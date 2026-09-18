#!/usr/bin/env python3
"""Summarize whole-dataset bootstrap fits and transparent predictive checks."""
import argparse,csv,json,collections
from pathlib import Path
import numpy as np
p=argparse.ArgumentParser();p.add_argument('bootstrap',type=Path);p.add_argument('observed_counts',type=Path);a=p.parse_args()
def features(path):
 with path.open() as f:
  r=csv.reader(f,delimiter='\t');header=next(r);x=np.array([[int(v) for v in row[2:]] for row in r])
 positive=(x>0).sum(axis=1)
 return {'mean_count':float(x.mean()),'mean_species_present':float(positive.mean()),'single_species_fraction':float((positive==1).mean()),'all_species_present_fraction':float((positive==x.shape[1]).mean()),'max_count':int(x.max()),'differential_gt20_fraction':float(((x.max(axis=1)-x.min(axis=1))>20).mean()),'mean_counts_by_species':dict(zip(header[2:],map(float,x.mean(axis=0))))}
observed=features(a.observed_counts);fits=[];predictive=[]
for directory in sorted(a.bootstrap.glob('replicate_*')):
 if not (directory/'audit.json').exists():continue
 audit=json.loads((directory/'audit.json').read_text());fit=audit['fit']
 fits.append({'replicate':audit['replicate'],'seconds':audit['seconds'],'max_observed':audit['observed_max'],'cap':audit['cap'],**{key:fit.get(key,'') for key in ['lambda','nu','alpha','epsilon','root_mean','alpha_at_bound','alpha_identifiable','boundary_fits','negative_log_likelihood','optimizer_converged','truncation_nll_difference','truncation_pass','root_mean_estimated']}})
 predictive.append(features(directory/'null_simulated.tsv'))
if not fits:raise RuntimeError('No completed fits')
with (a.bootstrap/'bootstrap_fit_summary.tsv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(fits[0]),delimiter='\t');w.writeheader();w.writerows(fits)
summary={'completed_datasets':len(fits),'observed':observed,'predictive':{},'parameter_bootstrap_quantiles':{},'interpretation':'Predictive intervals describe replicated datasets under the fitted model. Parameter quantiles are parametric-bootstrap percentiles, not a full model-uncertainty analysis.'}
for k,v in observed.items():
 if isinstance(v,dict):continue
 vals=np.array([r[k] for r in predictive]);summary['predictive'][k]={'observed':v,'replicate_mean':float(vals.mean()),'replicate_2.5_percentile':float(np.quantile(vals,.025)),'replicate_97.5_percentile':float(np.quantile(vals,.975)),'observed_outside_interval':bool(v<np.quantile(vals,.025) or v>np.quantile(vals,.975))}
for k in ['lambda','nu','alpha','epsilon','root_mean']:
 vals=np.array([float(r[k]) for r in fits]);summary['parameter_bootstrap_quantiles'][k]=dict(zip(['q025','median','q975'],map(float,np.quantile(vals,[.025,.5,.975]))))
(a.bootstrap/'diagnostics.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
