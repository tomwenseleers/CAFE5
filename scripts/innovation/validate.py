#!/usr/bin/env python3
"""Independent numerical/integration validation; requires numpy and scipy."""
import argparse,csv,json,subprocess
from pathlib import Path
import numpy as np
from scipy.linalg import expm
from scipy.stats import poisson
p=argparse.ArgumentParser();p.add_argument('binary',type=Path);p.add_argument('output',type=Path);a=p.parse_args();a.binary=a.binary.resolve();a.output.mkdir(parents=True,exist_ok=True)
checks={}
def run(args,allowed=(0,)):
 r=subprocess.run([str(a.binary),'--innovation']+list(map(str,args)),capture_output=True,text=True)
 if r.returncode not in allowed:raise AssertionError((args,r.returncode,r.stdout,r.stderr))
 return r

def matrix(l,n,t,K=35):
 f=a.output/'matrix.tsv';run(['--lambda',l,'--nu',n,'--matrix-time',t,'--max-count',K,'--matrix-output',f]);return np.loadtxt(f)
def generator(l,n,K):
 q=np.zeros((K+1,K+1))
 for i in range(K+1):
  if i<K:q[i,i+1]=l*i+n
  if i:q[i,i-1]=l*i
  q[i,i]=-q[i].sum()
 return q
for l,n,t in [(.2,.3,.7),(.2,0,.7),(0,.3,.7),(0,0,1),(.2,.3,0),(1e-14,.3,.7),(.5,.2,5)]:
 m=matrix(l,n,t);ref=expm(generator(l,n,160)*t)[:36,:36]
 err=float(np.max(np.abs(m-ref)));assert err<2e-12,(l,n,t,err)
 assert np.min(m)>=0 and np.max(m.sum(axis=1))<=1+1e-12
 checks[f'kernel_{l}_{n}_{t}']=err
m=matrix(.2,.3,.7,150);j=np.arange(151)
assert np.max(np.abs(m[:10]@j-(np.arange(10)+.3*.7)))<1e-10
checks['moment_max_error']=float(np.max(np.abs(m[:10]@j-(np.arange(10)+.3*.7))))
m=matrix(.2,.3,.7);assert abs(m[0,0]-(1+.2*.7)**(-.3/.2))<1e-14
assert np.max(np.abs(matrix(.2,.3,.3,100)@matrix(.2,.3,.4,100)-matrix(.2,.3,.7,100))[:10,:10])<1e-12
checks['semigroup']='pass'
large=matrix(0,900,1,1000)
assert abs(large[0,900]-poisson.pmf(900,900))<1e-12
checks['large_immigration_no_underflow']='pass'
tree=a.output/'small.nwk';tree.write_text('((A:0.3,B:0.3):0.4,C:0.7);\n')
counts=a.output/'counts.tsv';counts.write_text('Desc\tFamily ID\tA\tB\tC\nx\tf1\t0\t1\t2\nx\tf2\t1\t0\t0\n')
K=35;l=.2;n=.3;prior=poisson.pmf(np.arange(K+1),1);prior/=prior.sum()
# Independent transition construction uses a large finite generator, never the BDI kernel.
P={t:expm(generator(l,n,160)*t)[:K+1,:K+1] for t in [.3,.4,.7]}
def likelihood(y):return prior@(P[.4]@(P[.3][:,y[0]]*P[.3][:,y[1]])*P[.7][:,y[2]])
inc=1-likelihood([0,0,0]);expected=[np.log(likelihood(y)/inc) for y in [[0,1,2],[1,0,0]]]
prefix=a.output/'small';run(['-t',tree,'-i',counts,'--root-mean',1,'--lambda',l,'--nu',n,'--max-count',K,'-o',prefix])
with Path(str(prefix)+'_families.tsv').open() as f:actual=[float(r['conditional_log_likelihood']) for r in csv.DictReader(f,delimiter='\t')]
assert np.max(np.abs(np.array(actual)-expected))<1e-12
with Path(str(prefix)+'_results.tsv').open() as f:result={r['parameter']:r['value'] for r in csv.DictReader(f,delimiter='\t')}
assert abs(float(result['negative_log_likelihood'])+sum(actual))<1e-12
checks['independent_likelihood_error']=float(np.max(np.abs(np.array(actual)-expected)))
with Path(str(prefix)+'_ancestral.tsv').open() as f:anc=list(csv.DictReader(f,delimiter='\t'))
y=[0,1,2];root=prior*(P[.4]@(P[.3][:,0]*P[.3][:,1]))*P[.7][:,2];root/=root.sum()
row=next(r for r in anc if r['Family ID']=='f1' and r['parent']=='NA')
assert abs(float(row['P_zero'])-root[0])<1e-12
# Independent internal-node marginal by explicit root sum.
internal=(prior*P[.7][:,2])@P[.4]*(P[.3][:,0]*P[.3][:,1]);internal/=internal.sum()
row=next(r for r in anc if r['Family ID']=='f1' and r['Node']=='Node1')
assert abs(float(row['P_zero'])-internal[0])<1e-12
checks['ancestral_root_and_internal']='pass'
# Pure innovation from a deterministic zero root must generate nonzero leaves.
zero=a.output/'zero';run(['-t',tree,'--root-mean',0,'--lambda',0,'--nu',.5,'--simulate',3000,'--seed',91,'-o',zero])
with Path(str(zero)+'_simulated.tsv').open() as f:r=csv.reader(f,delimiter='\t');next(r);ys=np.array([list(map(int,x[2:])) for x in r])
assert np.all(ys.sum(axis=1)>0)
checks['zero_root_simulation_nonzero_families']=len(ys)
# Seed reproducibility.
other=a.output/'zero_repeat';run(['-t',tree,'--root-mean',0,'--lambda',0,'--nu',.5,'--simulate',3000,'--seed',91,'-o',other])
assert Path(str(zero)+'_simulated.tsv').read_bytes()==Path(str(other)+'_simulated.tsv').read_bytes()
# Explicit rejection of unsupported mode/input rather than silent fallback.
for args in [['--gamma','2'],['--lambda','-1'],['--nu','nan'],['--max-count','1.5'],['--lambda','1e300','--nu','1','--matrix-time','1e300','--matrix-output',a.output/'invalid_matrix.tsv']]:run(args,allowed=(1,))
checks['invalid_arguments']='pass'
# Simulation moment checks for a star tree with deterministic root 2, unconditioned.
star=a.output/'star.nwk';star.write_text('(A:0.7,B:0.7);\n');rp=a.output/'root.txt';rp.write_text('2 1\n')
sp=a.output/'moments';run(['-t',star,'--root-prior',rp,'--lambda',l,'--nu',n,'--simulate',30000,'--unconditioned','--seed',191,'-o',sp])
with Path(str(sp)+'_simulated.tsv').open() as f:r=csv.reader(f,delimiter='\t');next(r);ys=np.array([list(map(int,x[2:])) for x in r])
mean=2+n*.7;var=2*2*l*.7+n*.7+n*l*.7**2
assert np.max(np.abs(ys.mean(axis=0)-mean))<.035
assert np.max(np.abs(ys.var(axis=0)-var))<.07
checks['simulation_mean']=ys.mean(axis=0).tolist();checks['simulation_variance']=ys.var(axis=0).tolist()
(a.output/'checks.json').write_text(json.dumps(checks,indent=2)+'\n');print(json.dumps(checks,indent=2))
