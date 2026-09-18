#!/usr/bin/env python3
"""Direct manuscript comparison for N0 families on the 17-species dated tree."""
import argparse,csv,collections,json,re,hashlib
from pathlib import Path
from compare_manuscript import read,write

def main():
 p=argparse.ArgumentParser();p.add_argument('repository',type=Path);p.add_argument('observed',type=Path);p.add_argument('tests',type=Path);p.add_argument('output',type=Path);p.add_argument('--statistic',choices=['transition','mean'],default='transition');p.add_argument('--family-threshold',type=float,default=.05);p.add_argument('--branch-threshold',type=float,default=.01);a=p.parse_args();a.output.mkdir(exist_ok=True,parents=True)
 if not 0<a.family_threshold<1 or not 0<a.branch_threshold<1:raise ValueError('Thresholds must be between zero and one')
 source=a.repository/'nextflow_runs/2_EXCON/3_EXCON_CAFE_run/results_EXCON/cafe/base/N0.tsv'
 hogs={r['HOG']:r for r in read(source)}
 old=read(a.repository/'output/focal_node_significant_changes_annotated.tsv')
 old_index={(r['HOG'],r['node']):r for r in old}
 oldclasses={r['HOG']:r['Fig. 3 class'] for r in read(a.repository/'output/supplemental_tables/TableS12.tsv')}
 annotations={r['HOG']:r for r in read(a.repository/'input_CAFE/annotation/HOG_node_changes_and_annotations.tsv')}
 for filename in ['root_zero_filtered_HOGs_annotated.tsv','tested_TE_HOGs_excluded_from_GO_and_interpretation.tsv','focal_significant_HOG_annotations_unique.tsv']:
  annotations.update({r['HOG']:r for r in read(a.repository/'output'/filename)})
 archived_te=set((a.repository/'output/wasp/backgrounds/tested_TE_HOGs_excluded.txt').read_text().split())
 genes={}
 for sp,short in [('Vespula_vulgaris','vv'),('Polistes_dominula','pd')]:
  with (a.repository/f'input_CAFE/annotation/{short}_genetable.csv').open() as f:
   for r in csv.DictReader(f):genes[(sp,r[f'gene_{short}'])]=(r[f'symbol_{short}'],r[f'description_{short}'],'NCBI gene table')
 for r in read(a.repository/'output/focal_significant_HOG_member_annotation_evidence.tsv'):
  genes.setdefault((r['species'],r['gene_id']),(r['gene_symbol'],r['gene_name'],r['annotation_source']))
 # Define branches by their descendants, never by the old CAFE numbering.
 nodes={};first=None
 for r in read(Path(str(a.observed)+'_ancestral.tsv')):
  if first is None:first=r['Family ID']
  if r['Family ID']!=first:break
  nodes[r['Node']]=r['parent']
 def tips(n):
  cs=[k for k,v in nodes.items() if v==n]
  return set().union(*(tips(c) for c in cs)) if cs else {n}
 all_tips=tips('Node0');social={s for s in all_tips if s.startswith(('Polistes_','Mischocyttarus_','Vespula_','Vespa_','Dolichovespula_'))};vespine={s for s in social if s.startswith(('Vespula_','Vespa_','Dolichovespula_'))}
 mapping={}
 for oldnode,desc in [('20',social),('25',vespine)]:
  matches=[n for n in nodes if tips(n)==desc]
  if len(matches)!=1:raise ValueError('Focal clade is not unique')
  mapping[matches[0]]=oldnode
 (a.output/'node_mapping.json').write_text(json.dumps({k:{'manuscript_node':v,'descendants':sorted(tips(k)),'parent':nodes[k]} for k,v in mapping.items()},indent=2)+'\n')
 te=re.compile(r'transposon|retrotransposon|retroviral|reverse transcriptase|RNA.directed DNA polymerase|transposase|piggybac|tigger|mariner|helitron|gypsy|copia|gag.pol',re.I)
 evidence=[];annotated={}
 for h,r in hogs.items():
  names=set();symbols=set()
  for sp in ['Vespula_vulgaris','Polistes_dominula','Ancistrocerus_nigricornis']:
   for g in r[sp].split(','):
    g=g.strip()
    if (sp,g) in genes:
     sym,name,source_name=genes[(sp,g)]
     if name and name!='NA':
      names.add(name);symbols.add(sym);evidence.append(dict(HOG=h,species=sp,gene_id=g,gene_symbol=sym,gene_name=name,source=source_name))
  ann=annotations.get(h,{});direct='; '.join(sorted(names));consensus=ann.get('consensus_gene_name','')
  if consensus.lower() in ['unannotated','unannotated hog','na'] and direct:consensus=direct
  annotated[h]=dict(manuscript_functional_class=oldclasses.get(h,''),annotation=consensus or direct or 'unannotated',member_gene_names=direct,member_gene_symbols='; '.join(sorted(symbols)),TE_in_manuscript_exclusion_list=h in archived_te,TE_related=h in archived_te or ann.get('TE_related_HOG_flag','').upper()=='TRUE' or bool(te.search(direct+' '+consensus)),annotation_source='Existing same-N0 annotation plus direct member gene tables')
 tests=read(a.tests);out=[]
 if a.statistic=='mean':
  for r in tests:
   r['transition_branch_p']=r['branch_p'];r['branch_p']=r['mean_change_branch_p'];r['branch_MC_low']=r['mean_change_MC_low'];r['branch_MC_high']=r['mean_change_MC_high'];r['selected_raw_thresholds']=r['selected_mean_change_test']
   r['direction']='expansion' if float(r['posterior_mean_change'])>0 else 'contraction' if float(r['posterior_mean_change'])<0 else 'unchanged'
   r['q_BY_focal_branches']=r.get('mean_change_q_BY_focal_branches','NA');r['p_Holm_focal_branches']=r.get('mean_change_p_Holm_focal_branches','NA')
 for r in tests:
  r['branch_test_statistic']=a.statistic;r['family_threshold']=a.family_threshold;r['branch_threshold']=a.branch_threshold
  r['selected_raw_thresholds']=str(float(r['family_p'])<a.family_threshold and float(r['branch_p'])<a.branch_threshold and r['direction']!='unchanged')
  r['MC_threshold_uncertain']=str((float(r['family_MC_low'])<a.family_threshold<=float(r['family_MC_high']) and float(r['branch_MC_low'])<a.branch_threshold) or (float(r['branch_MC_low'])<a.branch_threshold<=float(r['branch_MC_high']) and float(r['family_MC_low'])<a.family_threshold))
 for r in tests:
  if r['Node'] not in mapping:continue
  h=r['Family ID'];oldnode=mapping[r['Node']];prior=old_index.get((h,oldnode),{})
  out.append({**r,'manuscript_node':oldnode,**annotated[h],'manuscript_significant':bool(prior),'manuscript_direction':prior.get('direction',''),'manuscript_change':prior.get('change',''),'manuscript_family_p':prior.get('family_p',''),'manuscript_branch_p':prior.get('branch_p',''),'manuscript_track':prior.get('result_source',''),'same_direction_replication':bool(prior) and r['selected_raw_thresholds']=='True' and r['direction']==prior['direction']})
 write(a.output/'all_focal_tests_annotated.tsv',out)
 selected=[r for r in out if r['selected_raw_thresholds']=='True'];write(a.output/'significant_HOGs_all.tsv',selected)
 nonte=[r for r in selected if not r['TE_related']];write(a.output/'significant_HOGs_nonTE.tsv',nonte)
 write(a.output/'manuscript_22_events_comparison.tsv',[r for r in out if r['manuscript_significant']])
 sel={r['Family ID'] for r in selected};write(a.output/'significant_HOG_member_evidence.tsv',[r for r in evidence if r['HOG'] in sel])
 summary={}
 for node in mapping:
  sub=[r for r in out if r['Node']==node];oldsub=[r for r in sub if r['manuscript_significant']]
  summary[node]={'manuscript_node':mapping[node],'manuscript_counts':dict(collections.Counter(r['manuscript_direction'] for r in oldsub)),'new_all_counts':dict(collections.Counter(r['direction'] for r in selected if r['Node']==node)),'new_nonTE_counts':dict(collections.Counter(r['direction'] for r in nonte if r['Node']==node)),'same_direction_replicated':sum(r['same_direction_replication'] for r in oldsub),'MC_uncertain_selected':sum(r['selected_raw_thresholds']=='True' and r['MC_threshold_uncertain']=='True' for r in sub)}
 (a.output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
