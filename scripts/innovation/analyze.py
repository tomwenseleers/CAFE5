#!/usr/bin/env python3
"""Fit the BDI+gamma+error model and compute nominal refitted bootstrap tests."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
from model_replay import read_results


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('counts', type=Path, help='CAFE tab-separated counts: description, family ID, species columns')
    p.add_argument('tree', type=Path, help='Rooted Newick tree with branch lengths')
    p.add_argument('-o', '--output', type=Path, default=Path('innovation_results'))
    p.add_argument('--binary', type=Path, default=Path(__file__).resolve().parents[2] / 'build/cafe5')
    p.add_argument('--gamma-cats', type=int, default=13, help='Discrete gamma approximation (default 13; not an automatically selected optimum)')
    p.add_argument('--root', choices=['one', 'poisson'], default='one', help='Root fixed at one (clade-origin HOGs), or estimated Poisson including zero')
    p.add_argument('--replicates', type=int, default=96, help='Whole simulated datasets, each refitted')
    p.add_argument('--workers', type=int, default=1)
    p.add_argument('--threads', type=int, default=2)
    p.add_argument('--seed', type=int, default=918100)
    p.add_argument('--nodes', nargs='+', help='Optional node names; default tests every branch')
    p.add_argument('--fit-only', action='store_true')
    a = p.parse_args()
    if not 1 <= a.gamma_cats <= 128 or min(a.workers, a.threads) < 1 or a.replicates < 2:
        p.error('Require 1..128 categories, positive resources and >=2 bootstrap datasets')
    a.binary = a.binary.resolve(); a.counts = a.counts.resolve(); a.tree = a.tree.resolve(); a.output = a.output.resolve()
    with a.counts.open() as f:
        rows = csv.reader(f, delimiter='\t'); header = next(rows); maximum = 0; n = 0
        for row in rows:
            if len(row) != len(header): raise ValueError('Inconsistent count-table row width')
            values = [int(x) for x in row[2:]]
            if not values or min(values) < 0 or max(values) == 0: raise ValueError('Require nonnegative counts and at least one observed copy per family')
            maximum = max(maximum, max(values)); n += 1
    if not n: raise ValueError('Empty count table')
    a.output.mkdir(parents=True, exist_ok=True)
    digest = lambda path: hashlib.sha256(path.read_bytes()).hexdigest()
    manifest = {'counts_sha256': digest(a.counts), 'tree_sha256': digest(a.tree), 'binary_sha256': digest(a.binary),
                'root': a.root, 'gamma_categories': a.gamma_cats, 'seed': a.seed,
                'workflow_sha256': digest(Path(__file__)), 'families': n}
    provenance = a.output / 'analysis.json'
    if provenance.exists() and json.loads(provenance.read_text()) != manifest:
        raise ValueError('Output belongs to different inputs/model/software; use a new output directory')
    provenance.write_text(json.dumps(manifest, indent=2) + '\n')
    prefix = a.output / 'fit'
    cap = max(40, maximum + 20)
    root = ['--root-family', 'hurdle-poisson', '--root-mean', '0', '--root-zero', '0'] if a.root == 'one' else ['--root-mean', '1', '--estimate-root-mean']
    if not (a.output / 'fit_complete.json').exists():
        while cap <= 2400:
            command = [str(a.binary), '--innovation', '-i', str(a.counts), '-t', str(a.tree), '-o', str(prefix),
                       '--gamma-cats', str(a.gamma_cats), '--estimate-mu', '--estimate-epsilon', '--estimate-epsilon-zero',
                       '--max-count', str(cap), '--iterations', '3600', '--starts', '3', '--threads', str(a.threads), '--seed', str(a.seed)] + root
            with (a.output / f'fit_cap{cap}.log').open('w') as log:
                log.write(' '.join(command) + '\n'); log.flush()
                status = subprocess.run(command, stdout=log, stderr=log).returncode
            result = read_results(prefix) if Path(str(prefix) + '_results.tsv').exists() else {}
            if status == 0 and result.get('optimizer_converged') == '1' and result.get('truncation_pass') == '1': break
            if result.get('truncation_pass') == '0': cap *= 2; continue
            raise RuntimeError('Fit failed or did not converge; inspect fit logs. No families have been discarded.')
        else: raise RuntimeError('Required numerical cap exceeds workflow limit 2400; use the low-level interface for investigation')
        free = 5 + int(a.gamma_cats > 1) + int(a.root == 'poisson')
        (a.output / 'fit_complete.json').write_text(json.dumps({'parameters': free, 'AIC': 2*float(result['negative_log_likelihood'])+2*free, 'cap': cap}, indent=2)+'\n')
    else:
        cap = json.loads((a.output / 'fit_complete.json').read_text())['cap']
        result = read_results(prefix)
        if result.get('optimizer_converged') != '1' or result.get('truncation_pass') != '1': raise ValueError('Saved fit no longer passes checks')
    if a.fit_only: return
    command = [sys.executable, str(Path(__file__).with_name('refitted_bootstrap.py')), str(a.binary), '--tree', str(a.tree),
               '--observed', str(prefix), '--output', str(a.output/'bootstrap'), '--replicates', str(a.replicates),
               '--workers', str(a.workers), '--threads', str(a.threads), '--seed', str(a.seed), '--minimum-cap', str(cap)]
    if a.nodes: command += ['--nodes'] + a.nodes
    subprocess.run(command, check=True)

if __name__ == '__main__':
    main()
