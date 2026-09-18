#!/usr/bin/env python3
"""Independent asymmetric-kernel, root/emission and simulation checks."""
import argparse,csv,json,subprocess
from pathlib import Path
import numpy as np
from scipy.linalg import expm
from scipy.stats import nbinom,poisson
p=argparse.ArgumentParser();p.add_argument('binary',type=Path);p.add_argument('output',type=Path);a=p.parse_args();a.binary=a.binary.resolve();a.output.mkdir(parents=True,exist_ok=True)
checks={}
def run(args,ok=(0,)):
 r=subprocess.run([str(a.binary),'--innovation']+list(map(str,args)),capture_output=True,text=True)
 assert r.returncode in ok,(args,r.returncode,r.stderr)
 return r
def matrix(l,m,n,t,K=40):
 f=a.output/'matrix.tsv';run(['--lambda',l,'--mu',m,'--nu',n,'--matrix-time',t,'--max-count',K,'--matrix-output',f]);return np.loadtxt(f)
def generator(l,m,n,K):
 q=np.zeros((K+1,K+1))
 for i in range(K+1):
  if i<K:q[i,i+1]=l*i+n
  if i:q[i,i-1]=m*i
  q[i,i]=-q[i].sum()
 return q
cases=[(1e-20,5e-21,.003,10),(5e-21,1e-20,.003,10),(.2,.4,.3,.7),(.4,.2,.3,.7),(.2,.2,.3,.7),(.2,.2+1e-12,.3,.7),(.2,.2-1e-12,.3,.7),(0,.3,.4,2),(.3,0,.4,2),(0,0,.4,2),(.3,.5,0,1),(1e-14,.3,.2,1),(.3,1e-14,.2,1),(.2,.3,.4,0)]
for l,m,n,t in cases:
 b=matrix(l,m,n,t);ref=expm(generator(l,m,n,180)*t)[:41,:41];err=float(np.max(abs(b-ref)))
 assert err<2e-11,(l,m,n,t,err)
 assert b.min()>=0 and b.sum(1).max()<1+1e-12
 checks[f'kernel_{l}_{m}_{n}_{t}']=err
assert np.max(abs((matrix(.2,.4,.3,.3,100)@matrix(.2,.4,.3,.4,100)-matrix(.2,.4,.3,.7,100))[:10,:10]))<1e-12
checks['semigroup']='pass'
tree=a.output/'tree.nwk';tree.write_text('((A:0.3,B:0.3):0.4,C:0.7);\n')
counts=a.output/'counts.tsv';counts.write_text('Desc\tFamily ID\tA\tB\tC\nx\tf1\t0\t1\t2\nx\tf2\t1\t0\t0\n')
l,m,n=.2,.4,.3;K=40;eps=.07;e0=.002
P={t:expm(generator(l,m,n,180)*t)[:K+1,:K+1] for t in [.3,.4,.7]}
E=np.zeros((K+1,K+1))
E[0,0]=1-e0;E[0,1]=e0
for j in range(1,K+1):
 E[j,j]=1-2*eps;E[j,j-1]=eps
 if j<K:E[j,j+1]=eps
common=['-t',tree,'-i',counts,'--lambda',l,'--mu',m,'--nu',n,'--epsilon',eps,'--epsilon-zero',e0,'--max-count',K,'--starts',1,'--likelihood-only']
for root in ['poisson','hurdle-poisson','hurdle-nb']:
 mean=.8;zero=.3;shape=1.7
 prior=poisson.pmf(np.arange(K+1),mean) if root=='poisson' else np.r_[zero,(1-zero)*(poisson.pmf(np.arange(K),mean) if root=='hurdle-poisson' else nbinom.pmf(np.arange(K),shape,shape/(shape+mean)))]
 prior/=prior.sum()
 def likelihood(y):return prior@(P[.4]@((P[.3]@E[:,y[0]])*(P[.3]@E[:,y[1]]))*(P[.7]@E[:,y[2]]))
 inc=1-likelihood([0,0,0]);expected=-sum(np.log(likelihood(y)/inc) for y in [[0,1,2],[1,0,0]])
 prefix=a.output/root;run(common+['--root-family',root,'--root-mean',mean,'--root-zero',zero,'--root-shape',shape,'-o',prefix])
 r={x['parameter']:x['value'] for x in csv.DictReader(open(str(prefix)+'_results.tsv'),delimiter='\t')}
 err=abs(float(r['negative_log_likelihood'])-expected);assert err<1e-10,(root,err)
 checks['likelihood_'+root]=err
# Empirical full-tree frequencies compared with independent pruning probabilities.
root='hurdle-nb';zero=.3;mean=.8;shape=1.7
prior=np.r_[zero,(1-zero)*nbinom.pmf(np.arange(K),shape,shape/(shape+mean))];prior/=prior.sum()
prefix=a.output/'simulation';run(['-t',tree,'--lambda',l,'--mu',m,'--nu',n,'--epsilon',eps,'--epsilon-zero',e0,'--root-family',root,'--root-mean',mean,'--root-zero',zero,'--root-shape',shape,'--max-count',K,'--simulate',50000,'--seed',314159,'-o',prefix])
with open(str(prefix)+'_simulated.tsv') as h: rr=csv.reader(h,delimiter='\t');next(rr);x=np.array([[int(v) for v in row[2:]] for row in rr])
inc=1-likelihood([0,0,0])
for pattern in [[1,1,1],[1,0,0],[0,1,2],[2,2,2]]:
 expected=likelihood(pattern)/inc;observed=np.all(x==pattern,axis=1).mean();se=np.sqrt(expected*(1-expected)/len(x))
 assert abs(expected-observed)<6*se+1/len(x),(pattern,expected,observed)
 checks['simulation_'+str(pattern)]={'expected':float(expected),'observed':float(observed)}
# At zero innovation and zero false occurrence, root-zero mass cancels exactly.
zero_scores=[]
for p0 in [.1,.8]:
 prefix=a.output/('zero_mass_'+str(p0));run(['-t',tree,'-i',counts,'--lambda',.2,'--mu',.4,'--nu',0,'--epsilon',.07,'--epsilon-zero',0,'--root-family','hurdle-nb','--root-mean',.8,'--root-zero',p0,'--root-shape',1.7,'--max-count',40,'--starts',1,'--likelihood-only','-o',prefix])
 result={r['parameter']:r['value'] for r in csv.DictReader(open(str(prefix)+'_results.tsv'),delimiter='\t')};zero_scores.append(float(result['negative_log_likelihood']))
assert abs(zero_scores[0]-zero_scores[1])<1e-11
checks['root_zero_nonidentifiability_at_zero_gain']=abs(zero_scores[0]-zero_scores[1])
# Supercritical simulation moments test the rescaled kernel branch independently.
mt=a.output/'moment_tree.nwk';mt.write_text('(A:4,B:4);\n');priorfile=a.output/'root_three.txt';priorfile.write_text('3 1\n');prefix=a.output/'supercritical'
run(['-t',mt,'--root-prior',priorfile,'--lambda',.4,'--mu',.1,'--nu',.3,'--unconditioned','--simulate',30000,'--seed',71825,'-o',prefix])
with open(str(prefix)+'_simulated.tsv') as h:
 reader=csv.reader(h,delimiter='\t');next(reader);z=np.array([[int(v) for v in row[2:]] for row in reader])
r=.3;t=4;growth=np.exp(r*t);expected_mean=3*growth+.3*np.expm1(r*t)/r
expected_var=3*(.4+.1)/r*growth*np.expm1(r*t)+.3/r*np.expm1(r*t)+.3*.4/r**2*np.expm1(r*t)**2
assert np.max(abs(z.mean(0)-expected_mean))<6*np.sqrt(expected_var/len(z))
checks['supercritical_simulation_mean']={'expected':float(expected_mean),'observed':z.mean(0).tolist()}
# Root/zero-error boundaries and CLI conflicts must fail explicitly.
for args in [['--root-family','bad'],['--root-zero',1],['--epsilon-zero',1],['--estimate-mu','--mu',.2],['--estimate-root-zero','--root-family','poisson']]:run(args,ok=(1,))
checks['invalid_options']='pass'
(a.output/'validation.json').write_text(json.dumps(checks,indent=2)+'\n');print(json.dumps(checks,indent=2))
