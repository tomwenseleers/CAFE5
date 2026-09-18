#!/usr/bin/env python3
"""Independent joint parent/child posterior diagnostics for nominated focal HOGs."""
import argparse,csv,json,re
from pathlib import Path
import numpy as np
from scipy.signal import lfilter
from scipy.special import gammaln
p=argparse.ArgumentParser();p.add_argument('--tree',type=Path,required=True);p.add_argument('--counts',type=Path,required=True);p.add_argument('--fit',type=Path,required=True);p.add_argument('--hog-list',type=Path,required=True);p.add_argument('--output',type=Path,required=True);p.add_argument('--nodes',nargs='+',default=['Node3','Node13']);a=p.parse_args()
def read(path):
 with path.open() as f:return list(csv.DictReader(f,delimiter='\t'))
r={x['parameter']:x['value'] for x in read(Path(str(a.fit)+'_results.tsv'))};K=int(r['max_count']);n=K+1;lam,nu,eps,mean=[float(r[k]) for k in ['lambda','nu','epsilon','root_mean']]
if r.get('model')!='BDI_equal_birth_death' or r.get('root_family','poisson')!='poisson' or r.get('epsilon_zero_separate','0')=='1':raise ValueError('This independent posterior check supports the original equal-rate Poisson-root model only')
categories=read(Path(str(a.fit)+'_categories.tsv'));rates=[float(x['lambda_multiplier']) for x in categories];wanted=set(a.hog_list.read_text().split());counts={x['Family ID']:x for x in read(a.counts) if x['Family ID'] in wanted};family={x['Family ID']:x for x in read(Path(str(a.fit)+'_families.tsv')) if x['Family ID'] in wanted}
tokens=re.findall(r'[(),:;]|[^(),:;\s]+',a.tree.read_text());pos=0;nodes=[]
def parse(parent=-1):
 global pos
 i=len(nodes);nodes.append({'name':'Node'+str(i),'parent':parent,'children':[],'time':0.})
 if tokens[pos]=='(':
  pos+=1;nodes[i]['children'].append(parse(i))
  while tokens[pos]==',':pos+=1;nodes[i]['children'].append(parse(i))
  assert tokens[pos]==')';pos+=1
 else:nodes[i]['name']=tokens[pos];pos+=1
 if tokens[pos] not in [':',',',')',';']:pos+=1
 if tokens[pos]==':':pos+=1;nodes[i]['time']=float(tokens[pos]);pos+=1
 return i
parse();index={x['name']:i for i,x in enumerate(nodes)}
pi=np.exp(-mean+np.arange(n)*np.log(mean)-gammaln(np.arange(n)+1));pi/=pi.sum()
def kernel(l,t):
 P=np.zeros((n,n));x=l*t;q=x/(1+x);survival=1/(1+x)
 if l==0:
  mass=np.exp(-nu*t+np.arange(n)*np.log(nu*t)-gammaln(np.arange(n)+1)) if nu*t>0 else np.r_[1.,np.zeros(K)]
  for i in range(n):P[i,i:]=mass[:n-i]
  return P
 if nu==0:P[0,0]=1
 else:
  logs=np.empty(n);logs[0]=-nu/l*np.log1p(x)
  for j in range(1,n):logs[j]=logs[j-1]+np.log((nu+l*(j-1))*t)-np.log(j)-np.log1p(x)
  P[0]=np.exp(logs)
 for i in range(1,n):P[i]=q*P[i-1]+survival**2*lfilter([0,1],[1,-q],P[i-1])
 return P
matrices=[{t:kernel(lam*g,t) for t in {x['time'] for x in nodes[1:]}} for g in rates]
rows=[];difference=np.arange(n)[None,:]-np.arange(n)[:,None]
for h,obs in counts.items():
 joints={name:np.zeros((n,n)) for name in a.nodes}
 for k,Ps in enumerate(matrices):
  weight=float(family[h][f'P_category_{k}'])
  if weight==0:continue
  inside=[np.ones(n) for _ in nodes]
  for i in reversed(range(len(nodes))):
   node=nodes[i]
   if not node['children']:
    y=int(obs[node['name']]);inside[i]=np.zeros(n)
    for truth in [y-1,y,y+1]:
     if 0<=truth<n:inside[i][truth]=(1-eps if truth==0 else 1-2*eps) if truth==y else eps
   else:
    for c in node['children']:
     inside[i]*=Ps[nodes[c]['time']]@inside[c];inside[i]/=inside[i].max()
  outside=[np.zeros(n) for _ in nodes];outside[0]=pi.copy()
  for i,node in enumerate(nodes):
   for c in node['children']:
    context=outside[i].copy()
    for sibling in node['children']:
     if sibling!=c:context*=Ps[nodes[sibling]['time']]@inside[sibling];context/=context.max()
    P=Ps[nodes[c]['time']];outside[c]=context@P;outside[c]/=outside[c].max()
    if nodes[c]['name'] in joints:
     joint=context[:,None]*P*inside[c][None,:];joint/=joint.sum();joints[nodes[c]['name']]+=weight*joint
 for node,joint in joints.items():
  hist=np.bincount((difference+K).ravel(),weights=joint.ravel(),minlength=2*K+1);cdf=hist.cumsum();i,j=np.unravel_index(joint.argmax(),joint.shape)
  rows.append({'HOG':h,'Node':node,'P_expansion':float(joint[difference>0].sum()),'P_contraction':float(joint[difference<0].sum()),'P_no_change':float(np.trace(joint)),'posterior_mean_change':float((difference*joint).sum()),'change_q025':int(np.searchsorted(cdf,.025)-K),'change_median':int(np.searchsorted(cdf,.5)-K),'change_q975':int(np.searchsorted(cdf,.975)-K),'joint_MAP_parent':int(i),'joint_MAP_child':int(j),'joint_MAP_change':int(j-i)})
reference={(x['Family ID'],x['Node']):float(x['posterior_mean_change']) for x in read(Path(str(a.fit)+'_branch_statistics.tsv')) if x['Family ID'] in wanted and x['Node'] in a.nodes}
error=max(abs(x['posterior_mean_change']-reference[(x['HOG'],x['Node'])]) for x in rows);assert error<1e-8,error
a.output.parent.mkdir(exist_ok=True,parents=True)
with a.output.open('w') as f:w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
print(json.dumps({'HOGs':len(counts),'max_mean_difference_from_Cpp':error,'note':'Posterior direction probabilities and intervals are conditional on fitted parameters; they are not branch p-values.'}))
