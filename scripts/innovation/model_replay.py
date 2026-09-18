"""Explicitly replay generating and refitted model options for bootstrap datasets."""
import csv
from pathlib import Path


def read_results(prefix):
    with Path(str(prefix) + '_results.tsv').open() as handle:
        return {row['parameter']: row['value'] for row in csv.DictReader(handle, delimiter='\t')}


def arguments(r):
    common=['--root-family',r.get('root_family','poisson'),'--root-mean',r['root_mean'],'--gamma-cats',r['gamma_categories']]
    fixed=[];start=[]
    for name in ['lambda','nu','epsilon']:
        fixed += ['--'+name,r[name]]
        start += (['--initial-'+name,r[name]] if r.get(name+'_estimated','1')=='1' else ['--'+name,r[name]])
    if r.get('epsilon_estimated','1')=='1':start+=['--estimate-epsilon']
    if int(r['gamma_categories'])>1:
        fixed+=['--alpha',r['alpha']];start+=(['--initial-alpha',r['alpha']] if r.get('alpha_estimated','1')=='1' else ['--alpha',r['alpha']])
    if r.get('model')=='BDI_separate_birth_death':
        fixed+=['--mu',r['mu']]
        start+=['--estimate-mu','--initial-mu',r['mu']] if r.get('mu_estimated')=='1' else ['--mu',r['mu']]
    if r.get('epsilon_zero_separate')=='1':
        fixed+=['--epsilon-zero',r['epsilon_zero']];start+=['--epsilon-zero',r['epsilon_zero']]
        if r.get('epsilon_zero_estimated')=='1':start+=['--estimate-epsilon-zero']
    if r.get('root_mean_estimated')=='1':start+=['--estimate-root-mean']
    if r.get('root_family','poisson')!='poisson':
        common+=['--root-zero',r['root_zero'],'--root-shape',r['root_shape']]
        if r.get('root_zero_estimated')=='1':start+=['--estimate-root-zero']
        if r.get('root_shape_estimated')=='1':start+=['--estimate-root-shape']
    return common,fixed,start
