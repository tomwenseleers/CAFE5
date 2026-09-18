#!/usr/bin/env python3
"""Audit retention of every input family and stability at an enlarged latent cap."""
import argparse
import csv
import json
from pathlib import Path


def read(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('counts', type=Path)
    p.add_argument('fit', type=Path)
    p.add_argument('larger_fit', type=Path)
    p.add_argument('output', type=Path)
    p.add_argument('--nodes', nargs='+', default=['Node3', 'Node13'])
    a = p.parse_args()
    a.output.mkdir(parents=True, exist_ok=True)
    input_rows = read(a.counts)
    families = read(Path(str(a.fit) + '_families.tsv'))
    expected = {r['Family ID'] for r in input_rows}
    actual = {r['Family ID'] for r in families}
    if len(expected) != len(input_rows) or len(actual) != len(families) or expected != actual:
        raise ValueError('Family identities or multiplicities changed')
    root = {r['Family ID']: r for r in read(Path(str(a.fit) + '_ancestral.tsv')) if r['parent'] == 'NA'}
    if set(root) != expected:
        raise ValueError('Missing root reconstruction')
    sizes = []
    for row in input_rows:
        values = [int(v) for k, v in row.items() if k not in ['Desc', 'Family ID']]
        sizes.append({'HOG': row['Family ID'], 'max_count': max(values),
                      'differential': max(values) - min(values),
                      'root_MAP_zero': int(root[row['Family ID']]['MAP_count']) == 0})
    large = [r for r in sizes if r['differential'] > 20]
    retention = {'families': len(actual), 'root_MAP_zero': sum(int(r['MAP_count']) == 0 for r in root.values()),
                 'root_Pzero_gt_0.95': sum(float(r['P_zero']) > .95 for r in root.values()),
                 'differential_gt20_retained': len(large), 'largest_family': max(sizes, key=lambda r: r['max_count']),
                 'zero_root_families_excluded': 0}
    (a.output / 'family_retention.json').write_text(json.dumps(retention, indent=2) + '\n')
    if large:
        with (a.output / 'large_family_audit.tsv').open('w') as handle:
            writer = csv.DictWriter(handle, fieldnames=list(large[0]), delimiter='\t')
            writer.writeheader()
            writer.writerows(large)
    one = {(r['Family ID'], r['Node']): r for r in read(Path(str(a.fit) + '_branch_statistics.tsv'))}
    two = {(r['Family ID'], r['Node']): r for r in read(Path(str(a.larger_fit) + '_branch_statistics.tsv'))}
    if set(one) != set(two):
        raise ValueError('Branch/family identities changed with cap')
    keys = [key for key in one if key[1] in a.nodes]
    mismatch = [key for key in keys if one[key]['MAP_change'] != two[key]['MAP_change']]
    caps = []
    for prefix in [a.fit, a.larger_fit]:
        results = {r['parameter']: r['value'] for r in read(Path(str(prefix) + '_results.tsv'))}
        if results['truncation_pass'] != '1':
            raise ValueError('Fit failed its internal truncation test')
        caps.append(int(results['max_count']))
    summary = {'families': len(actual), 'branches_compared': len(one), 'caps': caps,
               'focal_MAP_change_disagreements': mismatch,
               'focal_max_abs_mean_difference': max(abs(float(one[k]['posterior_mean_change']) - float(two[k]['posterior_mean_change'])) for k in keys),
               'focal_max_abs_transition_score_difference': max(abs(float(one[k]['transition_tail_score']) - float(two[k]['transition_tail_score'])) for k in keys)}
    (a.output / 'state_space_check.json').write_text(json.dumps(summary, indent=2) + '\n')
    if mismatch or summary['focal_max_abs_mean_difference'] > 1e-7 or summary['focal_max_abs_transition_score_difference'] > 1e-7:
        raise ValueError('Focal reconstruction is sensitive to latent cap')
    print(json.dumps({'retention': retention, 'state_space': summary}, indent=2))


if __name__ == '__main__':
    main()
