#!/usr/bin/env python3
"""Attach direct member annotations and explicit N0 overlap to N11 bootstrap calls."""
import argparse,csv,collections,json,re
from pathlib import Path
from compare_manuscript import read,write

def main():
 p=argparse.ArgumentParser();p.add_argument('repository',type=Path);p.add_argument('comparison',type=Path);p.add_argument('tests',type=Path);p.add_argument('output',type=Path);a=p.parse_args();a.output.mkdir(exist_ok=True,parents=True)
 base=a.repository/'nextflow_runs/2_EXCON/2_EXCON_orthofinder_eggnogmapper_run/results_EXCON/orthofinder/Phylogenetic_Hierarchical_Orthogroups'
 n11={r['HOG']:r for r in read(base/'N11.tsv')};lookup={}
 for sp,short in [('Vespula_vulgaris','vv'),('Polistes_dominula','pd')]:
  with (a.repository/f'input_CAFE/annotation/{short}_genetable.csv').open() as f:
   for r in csv.DictReader(f):lookup[(sp,r[f'gene_{short}'])]=(r[f'symbol_{short}'],r[f'description_{short}'],'direct NCBI gene table')
 for r in read(a.repository/'output/focal_significant_HOG_member_annotation_evidence.tsv'):
  lookup.setdefault((r['species'],r['gene_id']),(r['gene_symbol'],r['gene_name'],r['annotation_source']))
 oldann={r['HOG']:r for r in read(a.repository/'input_CAFE/annotation/HOG_node_changes_and_annotations.tsv')}
 oldann.update({r['HOG']:r for r in read(a.repository/'output/focal_significant_HOG_annotations_unique.tsv')})
 cross=collections.defaultdict(list)
 for r in read(a.comparison/'N0_N11_membership_crosswalk.tsv'):cross[r['N11_HOG']].append(r)
 te=re.compile(r'transposon|retrotransposon|retroviral|reverse transcriptase|RNA.directed DNA polymerase|transposase|piggybac|tigger|mariner|helitron|gypsy|copia|gag.pol',re.I)
 def category(text):
  groups=[('Chemoperception',r'odorant|olfactor|gustatory|chemosensor'),('Lipid metabolism',r'fatty.acid|elongation of very long|elongase|lipase|acyl|lipid'),('Digestion/proteolysis',r'aminopeptidase|trypsin|chymotrypsin|carboxypeptidase|protease'),('Endocrine-associated',r'ecdysone|juvenile.hormone|farnesol|vitellogenin'),('Chromatin',r'histone|chromatin|ESCO'),('Other metabolism',r'dehydrogenase|cytochrome P450|ornithine|antizyme|carbonic anhydrase'),('Neural signalling',r'glutamate receptor')]
  return '; '.join(g for g,pattern in groups if re.search(pattern,text,re.I)) or ('Unannotated' if not text else 'Other/uncategorised')
 annotations={};evidence=[]
 for h,row in n11.items():
  names=set();symbols=set();sources=set();total=0;named=0
  for sp in ['Vespula_vulgaris','Polistes_dominula','Ancistrocerus_nigricornis']:
   for gene in row[sp+'.clean'].split(','):
    gene=gene.strip()
    if not gene:continue
    total+=1
    if (sp,gene) in lookup:
     sym,name,source=lookup[(sp,gene)]
     if name and name!='NA':
      names.add(name);symbols.add(sym);sources.add(source);named+=1
      evidence.append(dict(N11_HOG=h,species=sp,gene_id=gene,gene_symbol=sym,gene_name=name,source=source))
  overlaps=sorted(cross[h],key=lambda x:-int(x['shared_genes']))
  inherited=sorted({oldann[x['N0_HOG']].get('consensus_gene_name','') for x in overlaps if x['N0_HOG'] in oldann})
  oldte=any(oldann.get(x['N0_HOG'],{}).get('TE_related_HOG_flag','').upper()=='TRUE' for x in overlaps)
  direct='; '.join(sorted(names));context='; '.join(x for x in inherited if x)
  annotations[h]=dict(member_gene_names=direct,member_gene_symbols='; '.join(sorted(symbols)),annotation_sources='; '.join(sorted(sources)),annotated_reference_species_members=named,reference_species_members=total,N0_overlap_HOGs='; '.join(x['N0_HOG'] for x in overlaps),N0_overlap_annotation_context=context,TE_member_evidence=bool(te.search(direct)),TE_any_N0_overlap=oldte,functional_categories=category(direct),annotation_note='N0 overlap labels are contextual; not guaranteed annotations of every N11 subfamily.')
 tests=read(a.tests);merged=[{**r,**annotations[r['Family ID']]} for r in tests]
 write(a.output/'all_focal_tests_annotated.tsv',merged)
 selected=[r for r in merged if r['selected_raw_thresholds']=='True'];write(a.output/'significant_focal_HOGs_all.tsv',selected)
 nonte=[r for r in selected if not r['TE_member_evidence'] and not r['TE_any_N0_overlap']];write(a.output/'significant_focal_HOGs_nonTE_conservative.tsv',nonte)
 selh={r['Family ID'] for r in selected};write(a.output/'significant_HOG_member_evidence.tsv',[r for r in evidence if r['N11_HOG'] in selh])
 idx={(r['Family ID'],r['Node']):r for r in merged};comparison=[]
 for old in read(a.comparison/'manuscript_focal_HOG_crosswalk.tsv'):
  new=idx.get((old['N11_HOG'],old['new_node']),{})
  comparison.append({**old,**{f'new_{k}':new.get(k,'') for k in ['posterior_mean_change','direction','family_p','branch_p','selected_raw_thresholds','MC_threshold_uncertain','member_gene_names','functional_categories']},'same_direction_significant':new.get('selected_raw_thresholds')=='True' and new.get('direction')==old['old_direction']})
 write(a.output/'manuscript_comparison_by_membership.tsv',comparison)
 event_groups=collections.defaultdict(list)
 for r in comparison:event_groups[(r['old_node'],r['old_HOG'])].append(r)
 events=[]
 for key,rs in event_groups.items():
  events.append(dict(old_node=key[0],old_HOG=key[1],old_direction=rs[0]['old_direction'],annotation=rs[0]['annotation'],N11_overlapping_families=len(rs),any_same_direction_significant=any(r['same_direction_significant'] for r in rs),N11_significant_matches='; '.join(r['N11_HOG'] for r in rs if r['same_direction_significant']),comparison_caution='Overlap-based correspondence, not a one-to-one replication test.'))
 write(a.output/'manuscript_event_summary.tsv',events)
 summary={node:{'all_selected':dict(collections.Counter(r['direction'] for r in selected if r['Node']==node)),'nonTE_conservative':dict(collections.Counter(r['direction'] for r in nonte if r['Node']==node)),'MC_uncertain_selected':sum(r['MC_threshold_uncertain']=='True' for r in selected if r['Node']==node)} for node in sorted({r['Node'] for r in tests})}
 summary['old_events_with_same_direction_significant_overlap']=sum(r['any_same_direction_significant'] for r in events)
 summary['old_events']=len(events)
 (a.output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
if __name__=='__main__':main()
