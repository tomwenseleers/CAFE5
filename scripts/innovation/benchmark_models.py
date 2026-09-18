#!/usr/bin/env python3
"""Compare CAFE-like process restrictions under one normalized likelihood.

This does not compare raw legacy CAFE objectives with BDI AIC. All fits share
root count one, observed-family conditioning and the same complete input table.
"""
import argparse
import csv
import json
from pathlib import Path
import subprocess
from model_replay import read_results

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('binary',type=Path);p.add_argument('counts',type=Path);p.add_argument('tree',type=Path)
p.add_argument('-o','--output',type=Path,required=True);p.add_argument('--threads',type=int,default=2)
a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
with a.counts.open() as f:
    rows=csv.reader(f,delimiter='\t');next(rows);maximum=max(int(x) for row in rows for x in row[2:])
summary=[]
for name,k,free,extra in [('critical_base',1,2,['--nu','0']),('critical_gamma',13,3,['--nu','0']),('bdi_gamma_error',13,6,['--estimate-mu','--estimate-epsilon-zero'])]:
    cap=max(180,maximum+20);prefix=a.output/name
    initial=['--initial-lambda','.02','--initial-epsilon','.01']+(['--initial-alpha','.2'] if k>1 else [])
    if name=='bdi_gamma_error':initial=['--initial-lambda','.01893405','--initial-mu','.03325355','--initial-nu','.00026007','--initial-alpha','.17720017','--initial-epsilon','.01045103','--epsilon-zero','.05291919']
    while cap<=2400:
        cmd=[str(a.binary.resolve()),'--innovation','-i',str(a.counts.resolve()),'-t',str(a.tree.resolve()),'-o',str(prefix),
             '--root-family','hurdle-poisson','--root-mean','0','--root-zero','0','--gamma-cats',str(k),'--estimate-epsilon',
             '--max-count',str(cap),'--iterations','2400','--starts','2','--threads',str(a.threads),'--skip-boundary-fits','--likelihood-only']+extra+initial
        with (a.output/f'{name}_cap{cap}.log').open('w') as log:
            status=subprocess.run(cmd,stdout=log,stderr=log).returncode
        r=read_results(prefix)
        (a.output/f'{name}_cap{cap}.json').write_text(json.dumps({'command':cmd,'exit_status':status,'fit':r},indent=2)+'\n')
        if status==0 and r['optimizer_converged']=='1' and r['truncation_pass']=='1':break
        if r['truncation_pass']=='0':
            cap*=2
            for flag,key in [('initial-lambda','lambda'),('initial-alpha','alpha'),('initial-epsilon','epsilon')]:
                if '--'+flag in initial:initial[initial.index('--'+flag)+1]=r[key]
            continue
        raise RuntimeError(f'{name} did not converge; inspect retained logs')
    else:raise RuntimeError(f'{name} exceeded cap limit; no AIC accepted')
    summary.append({'model':name,'K':k,'parameters':free,'NLL':float(r['negative_log_likelihood']),'AIC':2*float(r['negative_log_likelihood'])+2*free,'cap':cap,'cap_difference':float(r['truncation_nll_difference'])})
    with (a.output/'comparison.tsv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=summary[0],delimiter='\t');w.writeheader();w.writerows(summary)
