#!/usr/bin/env python3
"""Descriptive whole-dataset predictive checks, including likelihood-only fits."""
import argparse
import csv
import json
import subprocess
from pathlib import Path
import numpy as np


def features(path):
    with path.open() as handle:
        reader = csv.reader(handle, delimiter='\t')
        header = next(reader)
        x = np.array([[int(v) for v in row[2:]] for row in reader])
    occupied = (x > 0).sum(axis=1)
    return {'families': len(x), 'species': len(header)-2, 'mean_count': float(x.mean()),
            'mean_species_present': float(occupied.mean()),
            'single_species_fraction': float((occupied == 1).mean()),
            'one_gene_total_fraction': float((x.sum(axis=1) == 1).mean()),
            'all_species_present_fraction': float((occupied == x.shape[1]).mean()),
            'one_copy_every_species_fraction': float((x == 1).all(axis=1).mean()),
            'differential_gt20_fraction': float(((x.max(axis=1)-x.min(axis=1)) > 20).mean()),
            **{'mean_count_' + taxon: float(x[:,i].mean()) for i,taxon in enumerate(header[2:])},
            **{'zero_fraction_' + taxon: float((x[:,i]==0).mean()) for i,taxon in enumerate(header[2:])}}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('binary', type=Path)
    p.add_argument('--tree', type=Path, required=True)
    p.add_argument('--fit', type=Path, required=True)
    p.add_argument('--counts', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--replicates', type=int, default=4)
    p.add_argument('--seed', type=int, default=927100)
    a = p.parse_args()
    a.output.mkdir(parents=True, exist_ok=True)
    with Path(str(a.fit)+'_results.tsv').open() as handle:
        r = {x['parameter']: x['value'] for x in csv.DictReader(handle, delimiter='\t')}
    if r['optimizer_converged'] != '1' or r['truncation_pass'] != '1':
        raise ValueError('Input fit failed diagnostics')
    if r.get('root_prior_file') or r['conditioned_on_observed'] != '1':
        raise ValueError('This diagnostic expects an observed-family parametric-root fit')
    observed = features(a.counts)
    command = [str(a.binary.resolve()), '--innovation', '-t', str(a.tree),
               '--root-mean', r['root_mean'], '--gamma-cats', r['gamma_categories'],
               '--lambda', r['lambda'], '--nu', r['nu'], '--max-count', r['max_count']]
    if r.get('model') == 'BDI_separate_birth_death':
        command += ['--mu', r['mu']]
    if r.get('epsilon_zero_separate') == '1':
        command += ['--epsilon-zero', r['epsilon_zero']]
    if r.get('root_family', 'poisson') != 'poisson':
        command += ['--root-family', r['root_family'], '--root-zero', r['root_zero'], '--root-shape', r['root_shape']]
    if int(r['gamma_categories']) > 1:
        command += ['--alpha', r['alpha']]
    command += ['--error-model', r['error_model_file']] if r.get('error_model_file') else ['--epsilon', r['epsilon']]
    predicted = []
    for i in range(a.replicates):
        prefix = a.output / f'replicate_{i:04d}'
        args = command + ['--simulate', str(observed['families']), '--seed', str(a.seed+i), '-o', str(prefix)]
        with Path(str(prefix)+'.log').open('w') as handle:
            handle.write(' '.join(args)+'\n')
            handle.flush()
            subprocess.run(args, stdout=handle, stderr=handle, check=True)
        predicted.append(features(Path(str(prefix)+'_simulated.tsv')))
    comparison = {}
    for key, value in observed.items():
        if key in ['families', 'species']:
            continue
        values = [x[key] for x in predicted]
        comparison[key] = {'observed': value, 'simulated_mean': float(np.mean(values)),
                           'simulated_min': min(values), 'simulated_max': max(values)}
    result = {'observed': observed, 'fit': r, 'replicates': a.replicates, 'seed': a.seed,
              'comparison': comparison,
              'interpretation': 'Descriptive replicated-dataset checks at fitted parameters. Not calibrated significance tests, biological confidence intervals or proof of a global likelihood optimum.'}
    (a.output / 'predictive_checks.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
