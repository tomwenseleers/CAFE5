#!/usr/bin/env python3
"""End-to-end release smoke test, including resume and provenance rejection."""
import csv
import json
import os
from pathlib import Path
import subprocess
import sys

binary=Path(sys.argv[1]).resolve(); out=Path(sys.argv[2]).resolve();out.mkdir(parents=True,exist_ok=True)
tree=out/'tree.nwk';tree.write_text('((A:1,B:1):1,C:2);\n')
subprocess.run([str(binary),'--innovation','-t',str(tree),'--root-family','hurdle-poisson','--root-mean','0','--root-zero','0','--lambda','.12','--mu','.2','--nu','.1','--epsilon','.02','--epsilon-zero','.04','--simulate','200','--seed','908','--max-count','40','-o',str(out/'null')],check=True)
cmd=[sys.executable,str(Path(__file__).with_name('analyze.py')),str(out/'null_simulated.tsv'),str(tree),'--binary',str(binary),'--gamma-cats','1','--replicates','2','--threads','1','-o',str(out/'result')]
subprocess.run(cmd,check=True)
path=out/'result/bootstrap/branch_tests.tsv'
with path.open() as f: rows=list(csv.DictReader(f,delimiter='\t'))
with (out/'result/bootstrap/family_tests.tsv').open() as f: families=list(csv.DictReader(f,delimiter='\t'))
assert len(families)==200 and all(float(r['family_p'])==float(families[i//4]['family_p']) for i,r in enumerate(rows))
assert len(rows)==800 and len({r['Family ID'] for r in rows})==200
assert not any('Holm' in k or 'q_' in k for k in rows[0])
for r in rows:
    p=min(1.,2*min(float(r['p_upper']),float(r['p_lower'])))
    assert abs(float(r['branch_p'])-p)<1e-14
    assert r['selected_raw_thresholds']==str(float(r['family_p'])<.05 and p<.01 and float(r['posterior_mean_change'])!=0)
    assert r['direction']==('expansion' if float(r['posterior_mean_change'])>0 else 'contraction' if float(r['posterior_mean_change'])<0 else 'unchanged')
tails=out/'result/bootstrap/replicate_0000/tails.npz';stamp=tails.stat().st_mtime_ns
subprocess.run(cmd,check=True)
assert tails.stat().st_mtime_ns==stamp
original=tree.read_text();tree.write_text(original.replace('C:2','C:3'))
try:
    failed=subprocess.run(cmd,capture_output=True,text=True)
    assert failed.returncode!=0 and 'different inputs/model/software' in failed.stderr
finally:tree.write_text(original)
print(json.dumps({'families':200,'branches':4,'datasets':2,'nominal_selection':'pass','resume':'pass','changed_input_rejection':'pass'}))
