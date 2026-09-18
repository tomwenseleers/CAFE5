#!/usr/bin/env python3
"""Check observed-family proposals against analytic star-tree probabilities."""
import argparse, csv, json, math, subprocess
from pathlib import Path
import numpy as np

p=argparse.ArgumentParser();p.add_argument('binary',type=Path);p.add_argument('output',type=Path);a=p.parse_args()
a.output.mkdir(parents=True,exist_ok=True);tree=a.output/'star.nwk';tree.write_text('(A:1,B:1,C:1);\n')
checks={}
for case,root,nu,error in [('rare_immigration',0,1e-12,0),('rare_error',0,0,1e-12),('rare_root_and_gain',1e-10,1e-11,0),('multiple_arrivals',0,.7,0),('mixed_error',0,.2,.08)]:
 prefix=a.output/case
 subprocess.run([str(a.binary.resolve()),'--innovation','-t',str(tree),'--lambda','0','--nu',str(nu),'--root-mean',str(root),'--epsilon',str(error),'--simulate','60000','--seed','914117','-o',str(prefix)],check=True,capture_output=True,timeout=120)
 with Path(str(prefix)+'_simulated.tsv').open() as h:
  r=csv.reader(h,delimiter='\t');next(r);x=np.array([[int(v) for v in row[2:]] for row in r])
 assert (x.sum(1)>0).all()
 # At root zero, independent Poisson immigration followed by +/-1 error.
 def emit(k,y):
  if k==0:return (1-error if y==0 else error if y==1 else 0)
  return (1-2*error if y==k else error if abs(y-k)==1 else 0)
 def tip(y):return sum(math.exp(-nu)*nu**k/math.factorial(k)*emit(k,y) for k in range(25))
 if root:
  # Rare-event limit; neglected terms are O(root + nu), far below MC error.
  expected={(1,1,1):root/(root+3*nu),(1,0,0):nu/(root+3*nu)}
 else:
  # Stable probability of at least one nonzero observed tip.
  # P(tip0) = exp(-nu)*(1-error+nu*error).
  log_p0=-nu+math.log1p(-error+nu*error)
  inc=-math.expm1(3*log_p0)
  expected={y:math.prod(tip(z) for z in y)/inc for y in [(1,0,0),(0,1,0),(1,1,1),(2,0,0)]}
 result={}
 for y,prob in expected.items():
  observed=float(np.all(x==y,axis=1).mean());se=math.sqrt(prob*(1-prob)/len(x))
  assert abs(observed-prob)<6*se+1/len(x),(case,y,prob,observed)
  result[str(y)]={'expected':prob,'observed':observed}
 checks[case]=result
(a.output/'validation.json').write_text(json.dumps(checks,indent=2)+'\n')
print(json.dumps(checks,indent=2))
