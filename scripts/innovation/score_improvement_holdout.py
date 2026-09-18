#!/usr/bin/env python3
"""Replay training-fitted parameters on held-out families, without refitting."""
import argparse,csv,hashlib,json,subprocess
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('binary',type=Path);p.add_argument('--training',type=Path,required=True);p.add_argument('--test',type=Path,required=True);p.add_argument('--tree',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.binary=a.binary.resolve();a.output.mkdir(parents=True,exist_ok=True)
rows=[]
for folder in sorted(a.training.iterdir()):
 if not folder.is_dir():continue
 auditfile=folder/'audit.json'
 if not auditfile.exists():rows.append({'model':folder.name,'status':'training_incomplete'});continue
 audit=json.loads(auditfile.read_text());r=audit['fit']
 if audit['returncode']!=0:rows.append({'model':folder.name,'status':'training_failed_numerical_checks'});continue
 out=a.output/folder.name;out.mkdir(exist_ok=True);prefix=out/'heldout'
 command=[str(a.binary),'--innovation','-t',str(a.tree),'-i',str(a.test),'--lambda',r['lambda'],'--nu',r['nu'],'--root-mean',r['root_mean'],'--gamma-cats',r['gamma_categories'],'--epsilon',r['epsilon'],'--max-count',r['max_count'],'--threads','2','--starts','1','--likelihood-only','-o',str(prefix)]
 if int(r['gamma_categories'])>1:command+=['--alpha',r['alpha']]
 if r.get('model')=='BDI_separate_birth_death':command+=['--mu',r['mu']]
 if r.get('epsilon_zero_separate')=='1':command+=['--epsilon-zero',r['epsilon_zero']]
 if r.get('root_family','poisson')!='poisson':command+=['--root-family',r['root_family'],'--root-zero',r['root_zero'],'--root-shape',r['root_shape']]
 manifest={'command':command,'training_fit_sha256':hashlib.sha256((folder/'fit_results.tsv').read_bytes()).hexdigest(),'test_sha256':hashlib.sha256(a.test.read_bytes()).hexdigest(),'binary_sha256':hashlib.sha256(a.binary.read_bytes()).hexdigest()}
 manifest_path=out/'manifest.json'
 if not manifest_path.exists() or json.loads(manifest_path.read_text())!=manifest or not Path(str(prefix)+'_results.tsv').exists():
  with (out/'score.log').open('w') as h:subprocess.run(command,stdout=h,stderr=h,check=True)
  manifest_path.write_text(json.dumps(manifest,indent=2)+'\n')
 with Path(str(prefix)+'_results.tsv').open() as h:rr={x['parameter']:x['value'] for x in csv.DictReader(h,delimiter='\t')}
 assert rr['truncation_pass']=='1' and rr['optimizer_converged']=='1'
 rows.append({'model':folder.name,'status':'scored','training_NLL':r['negative_log_likelihood'],'training_families':r['families'],'test_NLL':rr['negative_log_likelihood'],'test_families':rr['families'],'test_NLL_per_family':float(rr['negative_log_likelihood'])/int(rr['families'])})
fields=list(dict.fromkeys(k for r in rows for k in r))
with (a.output/'heldout_scores.tsv').open('w') as h:w=csv.DictWriter(h,fields,delimiter='\t');w.writeheader();w.writerows(rows)
print(json.dumps(rows,indent=2))
