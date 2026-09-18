#!/usr/bin/env python3
"""Assemble the N11 model comparison, retaining incomplete/failed fits explicitly."""
import argparse,csv,json,math,shutil
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('analysis',type=Path);p.add_argument('output',type=Path);a=p.parse_args();b=a.analysis;a.output.mkdir(parents=True,exist_ok=True)
def read(prefix):
 with Path(str(prefix)+'_results.tsv').open() as h:return {r['parameter']:r['value'] for r in csv.DictReader(h,delimiter='\t')}
def write(name,rows):
 if not rows:return
 fields=list(dict.fromkeys(k for r in rows for k in r))
 with (a.output/name).open('w') as h:w=csv.DictWriter(h,fields,delimiter='\t');w.writeheader();w.writerows(rows)
models={};root=b.parent/'focal_bootstrap/N11_adequacy';old=json.loads((root/'predictive_a/predictive_checks.json').read_text())
models['baseline']={'fit':old['fit'],'predictive':old,'status':'converged_archived_baseline'}
for group in ['screen','root01_ablations']:
 for folder in sorted((b/group).iterdir()):
  if not folder.is_dir():continue
  audit=folder/'audit.json'
  if audit.exists():
   d=json.loads(audit.read_text());models[folder.name]={'fit':d['fit'],'predictive':d.get('predictive',{}),'status':'converged' if d['returncode']==0 else 'failed_numerical_checks'}
  else:models[folder.name]={'fit':{},'predictive':{},'status':'incomplete'}
if (b/'root01/fit_results.tsv').exists():
 r=read(b/'root01/fit');pred=b/'root01/predictive/predictive_checks.json';models['root01']={'fit':r,'predictive':json.loads(pred.read_text()) if pred.exists() else {},'status':'converged' if r['optimizer_converged']=='1' and r['truncation_pass']=='1' else 'failed_numerical_checks'}
for name,folder in [('asymmetric_second_start','asymmetric'),('combined_second_start','combined')]:
 prefix=b/'refinement'/folder/'second_start'
 if Path(str(prefix)+'_results.tsv').exists():
  r=read(prefix);pred=prefix.parent/'predictive/predictive_checks.json'
  models[name]={'fit':r,'predictive':json.loads(pred.read_text()) if pred.exists() else {},'status':'converged_boundary_caution' if r['optimizer_converged']=='1' and r['truncation_pass']=='1' else 'failed_numerical_checks'}
metrics=['mean_count','mean_species_present','single_species_fraction','one_gene_total_fraction','all_species_present_fraction','one_copy_every_species_fraction','differential_gt20_fraction']
rows=[];predictions=[];species=[]
for name,d in models.items():
 r=d['fit'];row={'model':name,'status':d['status']}
 for k in ['negative_log_likelihood','families','lambda','mu','nu','alpha','epsilon','epsilon_zero','root_family','root_mean','root_zero','root_shape','P_root_zero','root_prior_omitted_mass','gamma_categories','optimizer_converged','truncation_pass']:row[k]=r.get(k,'')
 if r:
  row['mu']=r.get('mu',r['lambda']);row['loss_duplication_ratio']=float(row['mu'])/float(r['lambda']) if float(r['lambda']) else math.nan
  row['epsilon_zero']=r.get('epsilon_zero',r['epsilon']);row['root_family']=r.get('root_family','poisson');row['P_root_zero']=r.get('P_root_zero',math.exp(-float(r['root_mean'])))
 for k,v in d['predictive'].get('comparison',{}).items():
  target=predictions if k in metrics else species
  target.append({'model':name,'metric':k,**v})
  if k in metrics:row[k+'_observed']=v['observed'];row[k+'_predicted']=v['simulated_mean']
 rows.append(row)
write('fit_comparison.tsv',rows);write('predictive_comparison.tsv',predictions);write('per_species_predictive_checks.tsv',species)
heldout=[]
for group in ['heldout','root01_heldout']:
 path=b/group/'heldout_scores.tsv'
 if path.exists():
  with path.open() as h:heldout+=list(csv.DictReader(h,delimiter='\t'))
write('heldout_scores.tsv',heldout)
for name in ['ultrametric_audit.json','gene_membership_audit.json']:
 if (b/name).exists():shutil.copyfile(b/name,a.output/name)
lines=['# N11 model-improvement comparison','',
 'All 13,836 observed families are retained on the matching 15-tip Vespidae tree. '
 'These are model-development and predictive checks, not new calibrated branch significance results. '
 'Convergence alone is not evidence of model adequacy. Failed or unfinished candidates remain visible.','',
 '| Model | Status | NLL | Loss / duplication | Single-species % predicted | One copy in every species % predicted |',
 '|---|---|---:|---:|---:|---:|']
for r in rows:
 def val(k,scale=1):
  x=r.get(k,'');return f'{float(x)*scale:.3f}' if x!='' else '—'
 lines.append(f"| {r['model']} | {r['status']} | {val('negative_log_likelihood')} | {val('loss_duplication_ratio')} | {val('single_species_fraction_predicted',100)} | {val('one_copy_every_species_fraction_predicted',100)} |")
lines+=['','Observed single-species fraction: **12.764%**. Observed exactly-one-copy-in-every-species fraction: **29.604%**. '
 'For hurdle roots, root_mean is the mean number of copies above one conditional on a positive root, '
 'not the unconditional mean. Root01 fixes that excess to zero and estimates the root-zero probability.','',
 'The predictive tables include mean occupancy, mean count, count tails and per-species checks. '
 'Simulation ranges are descriptive replicated-dataset ranges, not confidence intervals for evolutionary changes. '
 'Extreme-tail disagreement alone need not invalidate a background model: genuine family-specific changes can be outliers. '
 'The widespread occupancy and conserved-copy discrepancies are the main adequacy concern.','',
 '## Held-out scores','',
 'Models are re-fitted to 11,166 training families and scored without parameter re-estimation on 2,670 test families. '
 'Original OG identities define the split, keeping related N11 subfamilies together. '
 'These validation scores guide model selection; they are not an unbiased final performance estimate after selecting a model. '
 'The root01 ablations were added after the flexible root reached that boundary and are exploratory.','',
 '| Model | Status | Test NLL per family |','|---|---|---:|']
for r in heldout:lines.append(f"| {r['model']} | {r['status']} | {r.get('test_NLL_per_family','—')} |")
lines+=['','## Structural and numerical limits','',
 'All root-to-tip lengths are 110, up to 8e-9 rounding. A shared homogeneous process and shared observation law '
 'therefore imply identical marginal tip distributions. Conditioning on at least one observed tip does not remove '
 'this restriction. Real per-species mean counts range from 0.655 to 0.979. No gene identifier is repeated across N11 HOGs. '
 'Persistent interspecies discrepancies could reflect lineage-specific evolution, systematic observation differences or ascertainment; '
 'these causes are not identified by the count table alone.','',
 'See [methods](../../model_improvement_methods.md) for transition equations, root laws, the root-zero identifiability limit, '
 'and the distinction between these predictive simulations and a calibrated branch bootstrap. '
 'All comparisons use the same family universe; no TE, singleton, large-family or inferred-root filter is applied.','']
if (b/'interpretation.md').exists():lines+=['', (b/'interpretation.md').read_text()]
(a.output/'comparison.md').write_text('\n'.join(lines))
print(a.output/'comparison.md')
