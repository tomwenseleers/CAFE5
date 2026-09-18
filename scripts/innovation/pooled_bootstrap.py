#!/usr/bin/env python3
"""Refit whole simulated datasets; pool exchangeable families, with cluster MC uncertainty.

Pooling is a variance-reduction approximation to the global-family bootstrap CDF,
not B*N independent refits. Student-t intervals use B independent dataset clusters.
These are Monte Carlo intervals, not confidence intervals for biological effects.
"""
import argparse,csv,json,subprocess,time,concurrent.futures,os,hashlib
from pathlib import Path
import numpy as np
from scipy.stats import t
from branch_bootstrap import read_results,adjust

def table(p):
 with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def stats(prefix,nodes):
 fam=table(Path(str(prefix)+'_families.tsv'));ids=[r['Family ID'] for r in fam]
 values=np.array([float(r['conditional_log_likelihood']) for r in fam]);index={h:i for i,h in enumerate(ids)}
 changes=np.full((len(ids),len(nodes)),np.nan);surprise=changes.copy();mapchange=changes.copy()
 for row in table(Path(str(prefix)+'_branch_statistics.tsv')):
  if row['Node'] in nodes:
   i,j=index[row['Family ID']],nodes.index(row['Node'])
   changes[i,j]=float(row['posterior_mean_change']);surprise[i,j]=float(row['transition_tail_score']);mapchange[i,j]=int(row['MAP_change'])
 if not all(np.isfinite(x).all() for x in [values,changes,surprise,mapchange]):raise ValueError('Missing/nonfinite focal statistics')
 return ids,values,changes,surprise,mapchange

def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('binary',type=Path);p.add_argument('--tree',type=Path,required=True);p.add_argument('--observed',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
 p.add_argument('--replicates',type=int,default=24);p.add_argument('--workers',type=int,default=3);p.add_argument('--threads',type=int,default=4);p.add_argument('--seed',type=int,default=918100);p.add_argument('--nodes',nargs='+',default=['Node2','Node12'])
 p.add_argument('--interior-refits',action='store_true');p.add_argument('--boundary-audits',type=int,default=2)
 p.add_argument('--minimum-cap',type=int,default=180)
 p.add_argument('--fixed',action='store_true',help='Diagnostic plug-in reference simulation without parameter refitting')
 a=p.parse_args();a.binary=a.binary.resolve();a.tree=a.tree.resolve();a.output.mkdir(parents=True,exist_ok=True)
 r=read_results(a.observed);ids,ll,change,surprise,mapchange=stats(a.observed,a.nodes);n=len(ids)
 if r['truncation_pass']!='1' or r['optimizer_converged']!='1':raise ValueError('Observed fit failed diagnostics')
 if r.get('root_prior_file') or r.get('error_model_file'):raise ValueError('This pooled driver currently supports a Poisson root law and scalar epsilon only')
 theta={k:float(r[k]) for k in ['lambda','nu','alpha','epsilon']}
 common=['--innovation','-t',str(a.tree),'--root-mean',r['root_mean'],'--gamma-cats',r['gamma_categories'],'--threads',str(a.threads),'--iterations','1200','--starts','1']
 fixed=sum(([f'--{k}',str(v)] for k,v in theta.items() if k!='alpha' or int(r['gamma_categories'])>1),[])
 start=sum(([f'--initial-{k}',str(v)] for k,v in theta.items()),[])
 if r.get('root_mean_estimated')=='1':start+=['--estimate-root-mean']
 manifest={'method':'pooled_fixed_plugin' if a.fixed else 'pooled_refitted_parametric_bootstrap','replicates':a.replicates,'interior_refits':a.interior_refits,'boundary_audits':a.boundary_audits,'families_per_dataset':n,'seed':a.seed,'nodes':a.nodes,'parameters':theta,'root_mean':r['root_mean'],'root_mean_estimated':r.get('root_mean_estimated','0'),'null':'Exchangeable observed family from fitted global BDI+gamma+error model. Not a branch-specific no-change null.','pooling':'CDF averages over all families per full simulated dataset. Dataset is the independent unit for MC standard errors.','family_statistic':'Conditional log probability of complete observed tip pattern; lower tail.','branch_statistic':'Primary: upper tail of posterior-category-weighted transition surprise at marginal-MAP parent/child counts. Sensitivity: doubled minimum tail of posterior mean change.','selection':'All nonzero observed families; no size or root-state filtering.','thresholds':{'family':.05,'branch_two_sided':.01},'MC_intervals':'Approximate pointwise 95% Student-t intervals across independent simulated datasets; not simultaneous and not biological confidence intervals.'}
 manifest['tree_sha256']=hashlib.sha256(a.tree.read_bytes()).hexdigest()
 manifest['observed_family_results_sha256']=hashlib.sha256(Path(str(a.observed)+'_families.tsv').read_bytes()).hexdigest()
 manifest['binary_sha256']=hashlib.sha256(a.binary.read_bytes()).hexdigest()
 previous=a.output/'method.json'
 if previous.exists():
  old=json.loads(previous.read_text())
  for key in ['method','families_per_dataset','seed','nodes','parameters','branch_statistic','root_mean','root_mean_estimated','tree_sha256','observed_family_results_sha256','binary_sha256','interior_refits','boundary_audits']:
   if key in old and old[key]!=manifest[key]:raise ValueError('Cannot resume with changed '+key)
 previous.write_text(json.dumps(manifest,indent=2)+'\n')
 def run(args,log):
  with log.open('w') as f:
   f.write(' '.join(map(str,[a.binary]+args))+'\n');f.flush()
   return subprocess.run([str(a.binary)]+args,stdout=f,stderr=f).returncode
 def replicate(b):
  folder=a.output/f'replicate_{b:04d}';folder.mkdir(exist_ok=True);prefix=folder/'fit';saved=folder/'tails.npz'
  if saved.exists():return b,str(saved)
  begun=time.time();sim=folder/'null'
  if not Path(str(sim)+'_simulated.tsv').exists():
   if run(common+fixed+['--simulate',str(n),'--seed',str(a.seed+b),'--max-count','180','-o',str(sim)],folder/'simulate.log'):raise RuntimeError(f'Simulation {b} failed')
  infile=Path(str(sim)+'_simulated.tsv')
  with infile.open() as f:
   cr=csv.reader(f,delimiter='\t');next(cr);maximum=max(int(x) for row in cr for x in row[2:])
  cap=max(a.minimum_cap,maximum+10)
  while True:
   interior=a.interior_refits and b>=a.boundary_audits and not a.fixed
   args=common+(fixed if a.fixed else start+['--estimate-epsilon'])+(['--skip-boundary-fits'] if interior else [])+['-i',str(infile),'--max-count',str(cap),'-o',str(prefix)]
   status=run(args,folder/f'fit_cap{cap}.log')
   result=read_results(prefix) if Path(str(prefix)+'_results.tsv').exists() else {}
   if status==0 and result.get('truncation_pass')=='1' and result.get('optimizer_converged')=='1':
    if interior and any(theta[k]>0 and float(result[k])<theta[k]/10 for k in ['lambda','nu','epsilon']):
     args.remove('--skip-boundary-fits');status=run(args,folder/f'boundary_fallback_cap{cap}.log');result=read_results(prefix)
     if status!=0:raise RuntimeError(f'Boundary fallback {b} failed; see logs')
    if result.get('truncation_pass')=='1' and result.get('optimizer_converged')=='1':break
   if result.get('optimizer_converged')=='1' and result.get('truncation_pass')=='0':
    cap*=2
    if cap>2400:raise RuntimeError(f'Replicate {b} requires cap above 2400; retained for explicit numerical investigation')
    continue
   raise RuntimeError(f'Replicate {b} failed; retained, never silently omitted. See {folder}')
  _,null_ll,null_change,null_surprise,_=stats(prefix,a.nodes)
  if len(null_ll)!=n:raise ValueError('A simulated family was lost; bootstrap cannot continue')
  family=(1+np.searchsorted(np.sort(null_ll),ll,side='right'))/(n+1)
  upper=np.empty_like(change);lower=upper.copy();transition=upper.copy()
  for j in range(len(a.nodes)):
   ordered=np.sort(null_change[:,j]);upper[:,j]=(1+n-np.searchsorted(ordered,change[:,j],side='left'))/(n+1);lower[:,j]=(1+np.searchsorted(ordered,change[:,j],side='right'))/(n+1)
   ordered=np.sort(null_surprise[:,j]);transition[:,j]=(1+n-np.searchsorted(ordered,surprise[:,j],side='left'))/(n+1)
  np.savez_compressed(saved,family=family,upper=upper,lower=lower,transition=transition)
  (folder/'audit.json').write_text(json.dumps({'replicate':b,'seed':a.seed+b,'observed_max':maximum,'cap':cap,'seconds':time.time()-begun,'fit':result},indent=2)+'\n')
  return b,str(saved)
 with concurrent.futures.ThreadPoolExecutor(max_workers=a.workers) as pool:
  futures=[pool.submit(replicate,b) for b in range(a.replicates)]
  for future in concurrent.futures.as_completed(futures):
   try:b,path=future.result()
   except Exception:
    for pending in futures:pending.cancel()
    raise
   print(f'Completed dataset {b+1}/{a.replicates}: {path}',flush=True)
 family=[];upper=[];lower=[];transition=[]
 for b in range(a.replicates):
  d=np.load(a.output/f'replicate_{b:04d}/tails.npz');family.append(d['family']);upper.append(d['upper']);lower.append(d['lower']);transition.append(d['transition'])
 def estimate(x):
  x=np.array(x);mean=x.mean(axis=0);se=x.std(axis=0,ddof=1)/np.sqrt(len(x));half=t.ppf(.975,len(x)-1)*se
  return mean,np.maximum(0,mean-half),np.minimum(1,mean+half),se
 fp,fl,fh,fs=estimate(family);up,ul,uh,us=estimate(upper);lp,llow,lh,ls=estimate(lower)
 bp=np.minimum(1,2*np.minimum(up,lp));bl=np.minimum(1,2*np.minimum(ul,llow));bh=np.minimum(1,2*np.minimum(uh,lh))
 tp,tl,th,ts=estimate(transition)
 by,holm=adjust(tp);mean_by,mean_holm=adjust(bp);rows=[]
 for i,h in enumerate(ids):
  for j,node in enumerate(a.nodes):
   selected=bool(fp[i]<.05 and tp[i,j]<.01 and mapchange[i,j]!=0)
   uncertain=(fl[i]<.05<=fh[i] and tl[i,j]<.01) or (tl[i,j]<.01<=th[i,j] and fl[i]<.05)
   rows.append({'Family ID':h,'Node':node,'posterior_mean_change':change[i,j],'MAP_change':int(mapchange[i,j]),'transition_tail_score':surprise[i,j],'direction':'expansion' if mapchange[i,j]>0 else 'contraction' if mapchange[i,j]<0 else 'unchanged','family_p':fp[i],'family_MC_low':fl[i],'family_MC_high':fh[i],'branch_p':tp[i,j],'branch_MC_low':tl[i,j],'branch_MC_high':th[i,j],'mean_change_branch_p':bp[i,j],'mean_change_MC_low':bl[i,j],'mean_change_MC_high':bh[i,j],'mean_change_q_BY_focal_branches':mean_by[i,j],'mean_change_p_Holm_focal_branches':mean_holm[i,j],'selected_mean_change_test':bool(fp[i]<.05 and bp[i,j]<.01 and change[i,j]!=0),'p_upper':up[i,j],'p_lower':lp[i,j],'q_BY_focal_branches':by[i,j],'p_Holm_focal_branches':holm[i,j],'selected_raw_thresholds':selected,'MC_threshold_uncertain':bool(uncertain),'independent_datasets':a.replicates,'simulated_families':a.replicates*n})
 with (a.output/'focal_tests.tsv').open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
 print(json.dumps({node:{d:sum(r['selected_raw_thresholds'] and r['Node']==node and r['direction']==d for r in rows) for d in ['expansion','contraction']} for node in a.nodes}),flush=True)
if __name__=='__main__':main()
