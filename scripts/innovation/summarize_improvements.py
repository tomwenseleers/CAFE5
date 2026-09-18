#!/usr/bin/env python3
"""Summarize completed screening fits without silently dropping failures."""
import argparse,csv,json
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('screen',type=Path);p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
keys=['mean_count','mean_species_present','single_species_fraction','one_gene_total_fraction','all_species_present_fraction','one_copy_every_species_fraction','differential_gt20_fraction']
rows=[];pred=[]
for folder in sorted(a.screen.iterdir()):
 if not folder.is_dir():continue
 path=folder/'audit.json'
 if not path.exists():rows.append({'model':folder.name,'status':'running_or_incomplete'});continue
 audit=json.loads(path.read_text());r=audit['fit'];row={'model':folder.name,'status':'passed_numerical_checks' if audit['returncode']==0 else 'failed_numerical_checks','seconds':audit['seconds']}
 for k in ['negative_log_likelihood','lambda','mu','nu','alpha','epsilon','epsilon_zero','root_family','root_mean','root_zero','root_shape','P_root_zero','optimizer_converged','truncation_pass','gamma_categories']:row[k]=r.get(k,'')
 predinfo=audit.get('predictive',{}).get('comparison',{})
 for k in keys:
  if k not in predinfo:continue
  d=predinfo[k];row[k+'_observed']=d['observed'];row[k+'_predicted']=d['simulated_mean'];pred.append({'model':folder.name,'metric':k,**d})
 rows.append(row)
for name,data in [('fit_comparison.tsv',rows),('predictive_comparison.tsv',pred)]:
 if not data:continue
 fields=list(dict.fromkeys(k for row in data for k in row))
 with (a.output/name).open('w') as h:w=csv.DictWriter(h,fields,delimiter='\t');w.writeheader();w.writerows(data)
print(json.dumps(rows,indent=2))
