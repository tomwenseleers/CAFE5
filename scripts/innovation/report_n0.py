#!/usr/bin/env python3
"""Assemble the completed N0 comparison; refuse partial bootstrap results."""
import argparse
import csv
import hashlib
import json
import shutil
from pathlib import Path


def read(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def write(path, rows):
    if not rows:
        return
    with path.open('w') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('analysis', type=Path)
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    base, out = args.analysis, args.output
    boot = base / 'refitted'
    method = json.loads((boot / 'method.json').read_text())
    diagnostics = json.loads((boot / 'diagnostics.json').read_text())
    retention = json.loads((base / 'family_retention.json').read_text())
    if diagnostics['completed_datasets'] != method['replicates']:
        raise ValueError('The requested bootstrap is incomplete')
    for i in range(method['replicates']):
        audit = json.loads((boot / f'replicate_{i:04d}' / 'audit.json').read_text())
        if any(audit['fit'][k] != '1' for k in ['optimizer_converged', 'truncation_pass']):
            raise ValueError(f'Replicate {i} failed diagnostics')
    out.mkdir(parents=True, exist_ok=True)
    posterior = {(r['HOG'], r['Node']): r for r in read(base / 'branch_posterior_details.tsv')}
    summaries, old, selected = {}, {}, {}
    for statistic, folder in [('transition', 'comparison'), ('mean', 'mean_comparison'),
                              ('transition_p05', 'comparison_p05'), ('mean_p05', 'mean_comparison_p05')]:
        source = boot / folder
        summaries[statistic] = json.loads((source / 'summary.json').read_text())
        old[statistic] = read(source / 'manuscript_22_events_comparison.tsv')
        selected[statistic] = read(source / 'significant_HOGs_all.tsv')
        for kind, rows in [('selected', selected[statistic]), ('manuscript_22', old[statistic])]:
            merged = []
            for row in rows:
                key = (row['Family ID'], row['Node'])
                if key not in posterior:
                    raise ValueError(f'Missing joint posterior for {key}')
                extras = {f'joint_{k}': v for k, v in posterior[key].items()
                          if k not in ['HOG', 'Node']}
                merged.append({**row, **extras})
            write(out / f'{statistic}_{kind}_with_posterior.tsv', merged)
        for filename in ['summary.json', 'overlap_R.tsv', 'node_mapping.json', 'significant_HOG_member_evidence.tsv']:
            shutil.copyfile(source / filename, out / f'{statistic}_{filename}')
    paired = []
    mean_index = {(r['Family ID'], r['Node']): r for r in old['mean']}
    for row in old['transition']:
        key = (row['Family ID'], row['Node'])
        other = mean_index[key]
        pair = {k: row[k] for k in ['Family ID', 'Node', 'manuscript_node', 'manuscript_direction',
                'manuscript_change', 'manuscript_functional_class', 'annotation', 'family_p',
                'family_MC_low', 'family_MC_high', 'MAP_change', 'posterior_mean_change']}
        for label, source in [('transition', row), ('mean', other)]:
            pair.update({f'{label}_{k}': source[k] for k in ['branch_p', 'branch_MC_low',
                        'branch_MC_high', 'selected_raw_thresholds', 'same_direction_replication',
                        'MC_threshold_uncertain']})
        pair.update({k: posterior[key][k] for k in ['P_expansion', 'P_contraction', 'P_no_change',
                    'change_q025', 'change_q975', 'joint_MAP_change']})
        for statistic in ['transition_p05', 'mean_p05']:
            source = next(r for r in old[statistic] if (r['Family ID'], r['Node']) == key)
            pair[f'{statistic}_selected'] = source['selected_raw_thresholds']
            pair[f'{statistic}_same_direction_replication'] = source['same_direction_replication']
        paired.append(pair)
    write(out / 'manuscript_22_events_both_tests.tsv', paired)
    for filename in ['family_retention.json', 'manuscript_count_audit.json',
                     'state_space_check.json', 'large_family_audit.tsv', 'branch_posterior_details.tsv',
                     'estimated_root_gamma6_results.tsv', 'estimated_root_gamma6_categories.tsv',
                     'estimated_root_gamma6_optimization.tsv', 'estimated_root_gamma3_results.tsv',
                     'independent_start_results.tsv', 'runtime.json',
                     'boundary_optimization.tsv', 'boundary_results.tsv']:
        shutil.copyfile(base / filename, out / filename)
    audit = json.loads((base / 'input_audit.json').read_text())
    audit['sha256'] = {Path(key).name: value for key, value in audit['sha256'].items()}
    (out / 'input_audit.json').write_text(json.dumps(audit, indent=2) + '\n')
    for filename in ['method.json', 'diagnostics.json', 'bootstrap_fit_summary.tsv']:
        shutil.copyfile(boot / filename, out / filename)
    features = []
    for name, row in diagnostics['predictive'].items():
        features.append({'feature': name, **row})
    write(out / 'predictive_checks.tsv', features)
    fit = {r['parameter']: r['value'] for r in read(base / 'estimated_root_gamma6_results.tsv')}
    independent = {r['parameter']: r['value'] for r in read(base / 'independent_start_results.tsv')}
    nll_difference = float(independent['negative_log_likelihood']) - float(fit['negative_log_likelihood'])
    if nll_difference < -.01:
        raise ValueError('Independent start found a better empirical fit; investigate before reporting')
    names = {'Node3': 'Social Vespidae (manuscript 20)', 'Node13': 'Vespinae (manuscript 25)'}
    lines = ['# N0 / 17-species manuscript comparison', '',
        'The innovation extension fits all observed N0 families together, but this empirical model '
        'fails substantial predictive checks. Its significant HOG sets should therefore **not replace '
        'the manuscript results as validated biological conclusions**. Both branch statistics are '
        'reported because they test different aspects of reconstructed change.', '',
        '## Completed analysis', '',
        f"The analysis includes {method['families_per_dataset']:,} families and {method['replicates']} "
        f"complete simulated-dataset refits ({method['replicates'] * method['families_per_dataset']:,} "
        'simulated families). All five parameters were re-estimated in each dataset. The first two '
        'refits audit exact-zero parameter faces; later refits use the documented interior shortcut '
        'with boundary fallback. Every replicate passed convergence and count-cap checks. '
        'Monte Carlo intervals use datasets as independent clusters.', '',
        'The tree is exactly the archived manuscript 17-tip tree. All 11,025 originally prepared '
        'family count rows agree exactly. The input now also includes 1,550 single-species families '
        'and two formerly omitted multispecies families. All 56 families with count differential '
        f">20 and all {retention['root_MAP_zero']:,} families with root MAP zero are retained. "
        'Maximum observed count: 141.', '',
        'Equal per-copy duplication and loss are retained: q(n,n+1)=lambda*g*n+nu and '
        'q(n,n-1)=lambda*g*n. Six equiprobable gamma categories affect lambda only; innovation nu '
        'is global. A zero-inclusive Poisson root mean and symmetric one-copy error epsilon are '
        'estimated. The likelihood conditions on the family being observed in at least one tip. '
        'Innovation permits re-entry into an existing HOG from zero; it does not estimate the '
        'number of entirely unobserved HOG identities or distinguish all biological mechanisms '
        'that can yield an apparent gain.', '',
        '| Parameter | Estimate |', '|---|---:|']
    for key in ['lambda', 'nu', 'alpha', 'epsilon', 'root_mean']:
        lines.append(f'| {key} | {float(fit[key]):.9g} |')
    lines += ['', 'Rates are per tree time unit; lambda is per existing copy, nu per family. '
        'At a positive true count, epsilon applies to each one-copy error direction. '
        f"Conditional negative log likelihood: {float(fit['negative_log_likelihood']):.6f}. "
        f'Independent-start NLL minus primary NLL: {nll_difference:.6g}. '
        'The latent cap is independently doubled from 180 to 360; the accompanying '
        '`state_space_check.json` records reconstruction and score differences.', '',
        '## Requested significance comparison', '',
        'Calls use raw family p<0.05 and branch p<0.01 plus nonzero inferred direction. '
        'Counts below exclude TE-associated HOGs for comparison with the manuscript; TE families '
        'remain in every model fit and in the all-family tables. The manuscript branch quantity '
        'is a Viterbi probability, whereas the new values are simulated population-model tail '
        'probabilities. The numerical cutoffs are the same, but the tests are not equivalent.', '',
        '| Branch | Manuscript exp / con | Transition test exp / con | Mean-change test exp / con | Same-direction manuscript events recovered: transition / mean |',
        '|---|---:|---:|---:|---:|']
    for node, label in names.items():
        tr, mn = summaries['transition'][node], summaries['mean'][node]
        def count_string(counts):
            return f"{counts.get('expansion', 0)} / {counts.get('contraction', 0)}"
        lines.append(f"| {label} | {count_string(tr['manuscript_counts'])} | "
                     f"{count_string(tr['new_nonTE_counts'])} | {count_string(mn['new_nonTE_counts'])} | "
                     f"{tr['same_direction_replicated']} / {mn['same_direction_replicated']} |")
    lines += ['', '**Additional nominal branch cutoff requested by the user:** family p<0.05 '
        'and branch p<0.05. These use exactly the same fitted model, simulations and p-values; '
        'only the decision cutoff changes. The manuscript column remains its original '
        'family<0.05 / branch<0.01 set, allowing us to count recovery of those 22 events.', '',
        '| Branch | Transition test exp / con at branch<0.05 | Mean-change test exp / con at branch<0.05 | Same-direction manuscript events recovered: transition / mean |',
        '|---|---:|---:|---:|']
    for node, label in names.items():
        tr, mn = summaries['transition_p05'][node], summaries['mean_p05'][node]
        lines.append(f"| {label} | {count_string(tr['new_nonTE_counts'])} | "
                     f"{count_string(mn['new_nonTE_counts'])} | "
                     f"{tr['same_direction_replicated']} / {mn['same_direction_replicated']} |")
    lines += ['', 'No FDR or other multiplicity correction determines either table. Estimating '
        'additional parameters does not automatically penalize these bootstrap p-values by '
        'subtracting degrees of freedom. It can make some observed variation less unusual '
        'under the fitted model, while parameter refitting also changes the bootstrap '
        'distribution. Effects need not all have the same sign. The different branch statistic, '
        'root likelihood, family inclusion and innovation model also change inference; this '
        'comparison cannot attribute the differences solely to parameter count.', '']
    lines += ['', 'The transition test calibrates the surprise of a transition between marginal '
        'ancestral modes, weighting gamma categories by their posterior probabilities. The '
        'sensitivity test calibrates signed posterior mean change against the global family '
        'population. Neither is a test of literal absence of all turnover on a branch. '
        'Marginal modes need not form a likely joint ancestral state pair; the posterior '
        'direction columns are essential for interpreting transition flags.', '',
        '### HOG identities and functions', '']
    for statistic in ['transition', 'mean', 'transition_p05', 'mean_p05']:
        matches = [r for r in old[statistic] if r['same_direction_replication'] == 'True']
        label = statistic.replace('_p05', ' (branch p<0.05)').capitalize()
        lines.append(f'**{label} test: manuscript events retained in the same direction**')
        lines.append('')
        if not matches:
            lines.append('None.')
        for row in matches:
            lines.append(f"- {names[row['Node']]}: `{row['Family ID']}`, {row['direction']}; "
                         f"{row['annotation']}; family p={float(row['family_p']):.5g}, "
                         f"branch p={float(row['branch_p']):.5g}.")
        lines.append('')
    lines += ['**All non-TE transition-test flags**', '',
              '| Branch / HOG | Annotation | Family p | Branch p [95% MC interval] | P(expansion) / P(contraction) | 95% posterior change interval |',
              '|---|---|---:|---:|---:|---:|']
    for row in selected['transition']:
        if row['TE_related'] == 'True':
            continue
        post = posterior[(row['Family ID'], row['Node'])]
        lines.append(f"| {names[row['Node']]} / `{row['Family ID']}` | {row['annotation'].replace('|', '/')} | "
                     f"{float(row['family_p']):.5g} | {float(row['branch_p']):.5g} "
                     f"[{float(row['branch_MC_low']):.5g}, {float(row['branch_MC_high']):.5g}] | "
                     f"{float(post['P_expansion']):.3f} / {float(post['P_contraction']):.3f} | "
                     f"[{post['change_q025']}, {post['change_q975']}] |")
    lines += ['', 'Posterior probabilities and intervals condition on the fitted parameters. '
        'They are not p-values and do not include parameter uncertainty. Annotation correspondence '
        'is descriptive; no new GO enrichment claim is made. All 22 manuscript events, including '
        'those failing either new threshold, are in [the paired comparison](manuscript_22_events_both_tests.tsv). '
        'Each statistic has a complete selected-HOG table with direct annotation evidence, TE flags, '
        'Monte Carlo intervals, multiplicity adjustments, and joint-posterior diagnostics.', '',
        '## Model adequacy limits the biological inference', '',
        '| Diagnostic | Observed | Simulated dataset mean | Simulated 2.5–97.5 percentiles |',
        '|---|---:|---:|---:|']
    for row in features:
        lines.append(f"| {row['feature']} | {row['observed']:.6g} | {row['replicate_mean']:.6g} | "
                     f"{row['replicate_2.5_percentile']:.6g}–{row['replicate_97.5_percentile']:.6g} |")
    lines += ['', 'The fitted population generates too many single-species families and too few '
        'widely distributed families. Refitting simulated datasets calibrates the procedure under '
        'that population model; it cannot repair its mismatch to the real data. A more flexible '
        'root distribution and observation/ascertainment model are candidates for subsequent '
        'model development, evaluated by predictive performance rather than recovery of a desired '
        'HOG list. Equal birth/death rates remain a simplifying restriction. The present analysis '
        'therefore establishes computational feasibility and an auditable comparison, not a '
        'validated replacement biological analysis.', '',
        '## Precision, validation and reproducibility', '',
        'The 95% Monte Carlo intervals are approximate pointwise Student-t intervals over '
        f"{method['replicates']} independent datasets. They are not simultaneous intervals or "
        'confidence intervals for biological change. Calls whose intervals cross a decision '
        'threshold are explicitly marked. Pooling uses every exchangeable family in a dataset '
        'but does not treat them as independent parameter refits. No Gaussian approximation '
        'to the statistic tails is used.', '',
        'The source passes 206 legacy tests (582 assertions), independent transition and mixture '
        'likelihood checks, root-mean recovery, and small-tree null calibration. Independent '
        'joint-posterior means agree with the C++ reconstruction within numerical precision. '
        'The small-tree calibration experiments cover limited parameter settings and do not '
        'prove universal error control. See [methods](../../focal_bootstrap_methods.md) for '
        'the exact nulls, pooling argument, optimization shortcut and numerical safeguards.', '',
        'Input/result hashes, bootstrap seeds, fitted parameters and all replicate diagnostics '
        'are supplied in the accompanying JSON and TSV files. Full simulated datasets, logs '
        'and all-family reconstructions remain in the local validation directory; they are '
        'not copied into the source repository. This report does not change the manuscript.', '']
    (out / 'README.md').write_text('\n'.join(lines))
    manifest = {path.name: hashlib.sha256(path.read_bytes()).hexdigest()
                for path in sorted(out.iterdir()) if path.is_file() and path.name != 'SHA256.json'}
    (out / 'SHA256.json').write_text(json.dumps(manifest, indent=2) + '\n')
    print(out / 'README.md')


if __name__ == '__main__':
    main()
