#!/usr/bin/env python3
"""Independent gamma/error likelihood, posterior, simulation and fitting checks."""
import argparse,csv,json,subprocess
from pathlib import Path
import numpy as np
from scipy.linalg import expm
from scipy.stats import gamma,poisson
p=argparse.ArgumentParser();p.add_argument('binary',type=Path);p.add_argument('output',type=Path);a=p.parse_args();a.binary=a.binary.resolve();a.output.mkdir(parents=True,exist_ok=True)
def run(args):
 r=subprocess.run([str(a.binary),'--innovation']+list(map(str,args)),capture_output=True,text=True)
 if r.returncode:raise RuntimeError((args,r.returncode,r.stdout,r.stderr))
 return r
def result(prefix):
 with Path(str(prefix)+'_results.tsv').open() as f:return {r['parameter']:r['value'] for r in csv.DictReader(f,delimiter='\t')}
tree=a.output/'small.nwk';tree.write_text('((A:0.3,B:0.3):0.4,C:0.7);\n')
counts=a.output/'counts.tsv';counts.write_text('Desc\tFamily ID\tA\tB\tC\nx\tf1\t0\t1\t2\nx\tf2\t1\t0\t0\n')
K=35;alpha=1.3;eps=.12;lam=.2;nu=.3;cats=3
prefix=a.output/'mixture';common=['-t',tree,'--root-mean',1,'--gamma-cats',cats,'--alpha',alpha,'--epsilon',eps,'--max-count',K]
run(common+['-i',counts,'--lambda',lam,'--nu',nu,'-o',prefix])
with Path(str(prefix)+'_categories.tsv').open() as f:r=list(csv.DictReader(f,delimiter='\t'));rates=np.array([float(x['lambda_multiplier']) for x in r]);weights=np.array([float(x['weight']) for x in r])
ref_rates=gamma.ppf((np.arange(cats)+.5)/cats,alpha,scale=1/alpha);ref_rates/=ref_rates.mean()
# Legacy CAFE's chi-square quantile algorithm is an approximation.
assert np.max(abs(ref_rates-rates))<1e-5
pi=poisson.pmf(np.arange(K+1),1);pi/=pi.sum()
def emission(y):
 e=np.zeros(K+1)
 for n in range(K+1):
  if n==y:e[n]=1-eps if n==0 else 1-2*eps
  elif y==n+1 or (n>0 and y==n-1):e[n]=eps
 return e
P=[]
for rate in rates:
 Q=np.zeros((161,161))
 for n in range(161):
  if n<160:Q[n,n+1]=lam*rate*n+nu
  if n:Q[n,n-1]=lam*rate*n
  Q[n,n]=-Q[n].sum()
 P.append({t:expm(Q*t)[:K+1,:K+1] for t in [.3,.4,.7]})
def roots(y):return [pi*(pp[.4]@((pp[.3]@emission(y[0]))*(pp[.3]@emission(y[1]))))*(pp[.7]@emission(y[2])) for pp in P]
def likelihood(y):return sum(w*x.sum() for w,x in zip(weights,roots(y)))
inc=1-likelihood([0,0,0]);ys=[[0,1,2],[1,0,0]];expected=np.array([np.log(likelihood(y)/inc) for y in ys])
with Path(str(prefix)+'_families.tsv').open() as f:actual=np.array([float(r['conditional_log_likelihood']) for r in csv.DictReader(f,delimiter='\t')])
err=float(np.max(abs(actual-expected)));assert err<2e-12,err
assert abs(float(result(prefix)['negative_log_likelihood'])+sum(expected))<2e-12
root=sum(w*x for w,x in zip(weights,roots(ys[0])));root/=root.sum()
with Path(str(prefix)+'_ancestral.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
r=next(r for r in rows if r['Family ID']=='f1' and r['parent']=='NA');assert abs(float(r['P_zero'])-root[0])<2e-12
with Path(str(prefix)+'_families.tsv').open() as f:first=next(csv.DictReader(f,delimiter='\t'))
category=np.array([w*x.sum() for w,x in zip(weights,roots(ys[0]))]);category/=category.sum()
assert max(abs(float(first[f'P_category_{k}'])-category[k]) for k in range(cats))<2e-12
# Fixed-file and parametric epsilon error laws must agree.
error=a.output/'error.txt';error.write_text('maxcnt: 100\ncntdiff: -1 0 1\n0 0 .88 .12\n1 .12 .76 .12\n')
fileprefix=a.output/'file_error';run(['-t',tree,'-i',counts,'--root-mean',1,'--gamma-cats',cats,'--alpha',alpha,'--error-model',error,'--lambda',lam,'--nu',nu,'--max-count',K,'-o',fileprefix])
assert abs(float(result(fileprefix)['negative_log_likelihood'])+sum(expected))<2e-12
# Mixture simulation must condition observed counts AFTER mixing categories.
sim=a.output/'simulation';run(common+['--lambda',lam,'--nu',nu,'--simulate',40000,'--seed',521,'-o',sim])
with Path(str(sim)+'_simulated.tsv').open() as f:r=csv.reader(f,delimiter='\t');next(r);samples=np.array([list(map(int,x[2:])) for x in r])
assert np.all(samples.sum(axis=1)>0)
patterns=[[0,0,1],[1,1,1],[1,0,0],[0,1,2]];simulation_checks=[]
for y in patterns:
 probability=likelihood(y)/inc;observed=np.mean(np.all(samples==y,axis=1));se=np.sqrt(probability*(1-probability)/len(samples))
 assert abs(observed-probability)<5*se+1/len(samples),(y,observed,probability)
 simulation_checks.append({'pattern':y,'expected':probability,'observed':observed})
# Simultaneous estimation of all four parameters with a richer tree.
rich=a.output/'rich.nwk';rich.write_text('(((A:.1,B:.1):.3,(C:.2,D:.2):.2):.4,(E:.35,F:.35):.45);\n')
records=[]
for seed in [691,692,693]:
 prefix=a.output/f'recovery_{seed}';c=['-t',rich,'--root-mean',1,'--gamma-cats',3,'--max-count',55]
 run(c+['--lambda',.2,'--nu',.3,'--alpha',1.3,'--epsilon',.08,'--simulate',4000,'--seed',seed,'-o',prefix])
 r=run(c+['-i',str(prefix)+'_simulated.tsv','--estimate-epsilon','--iterations',1000,'--starts',2,'--threads',2,'--likelihood-only','-o',prefix])
 rr=result(prefix);record={k:float(rr[k]) for k in ['lambda','nu','alpha','epsilon']};record['seed']=seed;records.append(record);print(record,flush=True)
mean={k:float(np.mean([r[k] for r in records])) for k in ['lambda','nu','alpha','epsilon']}
assert abs(mean['lambda']-.2)<.04 and abs(mean['nu']-.3)<.06 and abs(mean['epsilon']-.08)<.025 and .6<mean['alpha']<2.6,mean
# Verify the gamma-degenerate lambda=0 boundary reduces to one category.
pure=a.output/'pure_immigration';run(common+['-i',counts,'--lambda',0,'--nu',nu,'-o',pure])
flat=a.output/'pure_immigration_flat';run(['-t',tree,'-i',counts,'--root-mean',1,'--epsilon',eps,'--max-count',K,'--lambda',0,'--nu',nu,'-o',flat])
assert abs(float(result(pure)['negative_log_likelihood'])-float(result(flat)['negative_log_likelihood']))<2e-12
summary={'independent_mixture_loglik_max_error':err,'gamma_multiplier_approximation_max_error':float(np.max(abs(ref_rates-rates))),'fixed_error_file_equivalence':'pass','mixture_observed_simulation':simulation_checks,'joint_recovery':records,'joint_recovery_means':mean,'limitations':'Three joint-recovery replicates at one generating regime; not general identifiability or interval coverage.'}
(a.output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
