#!/usr/bin/env python3
"""Paired diagnostic of reference-statistic refitting, holding simulations fixed."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np


def read(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('observed', type=Path)
    p.add_argument('refitted', type=Path)
    p.add_argument('fixed', type=Path)
    p.add_argument('manuscript_events', type=Path)
    p.add_argument('output', type=Path)
    p.add_argument('--replicates', type=int, default=4)
    a = p.parse_args()
    if a.replicates < 2:
        raise ValueError('At least two independent dataset pairs are required')
    a.output.mkdir(parents=True, exist_ok=True)
    first = json.loads((a.refitted / 'method.json').read_text())
    second = json.loads((a.fixed / 'method.json').read_text())
    for key in ['parameters', 'root_mean', 'nodes', 'observed_family_results_sha256', 'binary_sha256']:
        if first[key] != second[key]:
            raise ValueError('Compared procedures differ in '+key)
    if first['nodes'] != ['Node3', 'Node13']:
        raise ValueError('Unexpected focal branch order')
    arrays, hashes = {}, []
    for mode, directory in [('refitted', a.refitted), ('fixed', a.fixed)]:
        arrays[mode] = {key: [] for key in ['family', 'transition']}
        for i in range(a.replicates):
            folder = directory / f'replicate_{i:04d}'
            digest = hashlib.sha256((folder / 'null_simulated.tsv').read_bytes()).hexdigest()
            if mode == 'refitted':
                hashes.append(digest)
            elif digest != hashes[i]:
                raise ValueError('The compared simulations are not identical')
            data = np.load(folder / 'tails.npz')
            for key in arrays[mode]:
                arrays[mode][key].append(data[key])
        arrays[mode] = {k: np.array(v) for k, v in arrays[mode].items()}
    family_ids = [r['Family ID'] for r in read(Path(str(a.observed)+'_families.tsv'))]
    index = {h: i for i, h in enumerate(family_ids)}
    nodes = ['Node3', 'Node13']
    changes = {(r['Family ID'], r['Node']): int(r['MAP_change'])
               for r in read(Path(str(a.observed)+'_branch_statistics.tsv')) if r['Node'] in nodes}
    means = {mode: {k: v.mean(axis=0) for k, v in data.items()} for mode, data in arrays.items()}
    rows = []
    for old in read(a.manuscript_events):
        h = old['HOG']
        node = {'20': 'Node3', '25': 'Node13'}[old['node']]
        i, j = index[h], nodes.index(node)
        row = {'HOG': h, 'Node': node, 'manuscript_direction': old['direction'], 'MAP_change': changes[h, node]}
        for mode in ['fixed', 'refitted']:
            fp, bp = means[mode]['family'][i], means[mode]['transition'][i, j]
            row[mode+'_family_p'] = fp
            row[mode+'_transition_p'] = bp
            for cutoff in [.01, .05]:
                row[f'{mode}_selected_branch_{cutoff}'] = bool(fp < .05 and bp < cutoff and changes[h, node] != 0)
        differences = arrays['refitted']['transition'][:, i, j] - arrays['fixed']['transition'][:, i, j]
        row['refitted_minus_fixed_branch_p'] = float(differences.mean())
        row['paired_difference_MC_SE'] = float(differences.std(ddof=1)/np.sqrt(a.replicates))
        rows.append(row)
    with (a.output / 'manuscript_events_refitting_diagnostic.tsv').open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)
    delta = np.array([r['refitted_minus_fixed_branch_p'] for r in rows])
    summary = {'paired_datasets': a.replicates, 'families_per_dataset': len(family_ids),
               'identical_simulation_sha256': hashes,
               'manuscript_event_branch_p_change': {'minimum': float(delta.min()), 'median': float(np.median(delta)), 'maximum': float(delta.max())},
               'manuscript_events_with_changed_raw_selection': {str(cutoff): sum(r[f'fixed_selected_branch_{cutoff}'] != r[f'refitted_selected_branch_{cutoff}'] for r in rows) for cutoff in [.01, .05]},
               'interpretation': 'Paired diagnostic of refitting parameters when scoring the same simulated datasets. It does not compare models with different parameter counts, measure power under alternatives, or replace the full refitted bootstrap. Monte Carlo uncertainty uses independent dataset pairs.'}
    (a.output / 'refitting_diagnostic.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
