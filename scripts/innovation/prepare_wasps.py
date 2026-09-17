#!/usr/bin/env python3
"""Audit HOG scope and produce exact copy counts; never infer outgroup zeros from N13."""
import argparse,csv,hashlib,json,re
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('repository',type=Path);p.add_argument('output',type=Path);a=p.parse_args()
a.output.mkdir(parents=True,exist_ok=True)

def parse_tree(text):
 tokens=re.findall(r"[(),:;]|[^(),:;\s]+",text);pos=0
 def parse():
  nonlocal pos
  children=[]
  if tokens[pos]=='(':
   pos+=1;children.append(parse())
   while tokens[pos]==',':pos+=1;children.append(parse())
   assert tokens[pos]==')';pos+=1
  label=''
  if tokens[pos] not in [':',',',')',';']:label=tokens[pos];pos+=1
  if tokens[pos]==':':pos+=2
  return {'label':label.removesuffix('.clean'),'children':children}
 return parse()
def leaves(node):
 return set().union(*(leaves(c) for c in node['children'])) if node['children'] else {node['label']}
def find(node,label):
 if node['label']==label:return node
 for c in node['children']:
  result=find(c,label)
  if result:return result
def splits(node,selected):
 result=set();ingroup=leaves(node)&selected
 if 1<len(ingroup)<len(selected):result.add(tuple(sorted(ingroup)))
 for c in node['children']:result|=splits(c,selected)
 return result
base=a.repository/'nextflow_runs/2_EXCON/2_EXCON_orthofinder_eggnogmapper_run/results_EXCON/orthofinder'
tree=a.repository/'input_CAFE/tree_dating/Vespidae_with_outgroups_dated_primary.nwk'
source_tree=base/'Species_Tree/SpeciesTree_rooted_node_labels.txt'
source_nodes=parse_tree(source_tree.read_text())
nw=tree.read_text().strip();taxa=re.findall(r'(?<=[(,])([^():,;]+):',nw)
assert len(taxa)==17
(a.output/tree.name).write_text(nw+'\n')
report={'orthofinder_tree_source':str(source_tree.resolve()),'orthofinder_tree_sha256':hashlib.sha256(source_tree.read_bytes()).hexdigest(),'tree_sha256':hashlib.sha256(tree.read_bytes()).hexdigest(),'tree_source':str(tree.resolve()),'taxa':taxa,'tables':{}}
for node in ['N10','N11','N12','N13']:
 source=base/'Phylogenetic_Hierarchical_Orthogroups'/f'{node}.tsv'
 with source.open() as f: rows=list(csv.DictReader(f,delimiter='\t'))
 counts=[];nonempty={t:0 for t in taxa};maximum=0
 for row in rows:
  y=[]
  for t in taxa:
   genes=[x.strip() for x in row[t+'.clean'].split(',') if x.strip()]
   assert len(genes)==len(set(genes)),(row['HOG'],t)
   y.append(len(genes));nonempty[t]+=bool(genes)
  if any(y):counts.append([node,row['HOG']]+y);maximum=max(maximum,max(y))
 info={'source':str(source.resolve()),'sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'source_rows':len(rows),'observed_rows_in_selected_taxa':len(counts),'max_count':maximum,'nonempty_families_by_taxon':nonempty}
 report['tables'][node]=info
 if node=='N10':
  source_node=find(source_nodes,node)
  assert set(taxa)<=leaves(source_node)
  assert splits(source_node,set(taxa))==splits(parse_tree(nw),set(taxa))
  info['rooted_topology_matches']=True
  with (a.output/'N10_Vespidae_with_outgroups_counts.tsv').open('w') as f:
   w=csv.writer(f,delimiter='\t');w.writerow(['Desc','Family ID']+taxa);w.writerows(counts)
 if node=='N11':
  vesp_tree=a.repository/'input_CAFE/tree_dating/Vespidae_dated_primary.nwk'
  vnw=vesp_tree.read_text().strip();vesp_taxa=re.findall(r'(?<=[(,])([^():,;]+):',vnw)
  assert leaves(find(source_nodes,node))==set(vesp_taxa)
  assert splits(find(source_nodes,node),set(vesp_taxa))==splits(parse_tree(vnw),set(vesp_taxa))
  info['rooted_topology_matches']=True
  (a.output/vesp_tree.name).write_text(vnw+'\n')
  info['tree_source']=str(vesp_tree.resolve());info['tree_sha256']=hashlib.sha256(vesp_tree.read_bytes()).hexdigest()
  with (a.output/'N11_Vespidae_counts.tsv').open('w') as f:
   w=csv.writer(f,delimiter='\t');w.writerow(['Desc','Family ID']+vesp_taxa)
   for row in counts:w.writerow(row[:2]+[row[2+taxa.index(t)] for t in vesp_taxa])
 if node=='N13':
  social=[t for t in taxa if nonempty[t]]
  # Vespidae's second root child is precisely the social-wasp crown.
  # Preserve its dated branch lengths, remove only the stem to its old root.
  vnw=(a.repository/'input_CAFE/tree_dating/Vespidae_dated_primary.nwk').read_text().strip()
  start=vnw.index(',')+1;depth=0;end=None
  for j in range(start,len(vnw)):
   if vnw[j]=='(':depth+=1
   elif vnw[j]==')':
    depth-=1
    if depth==0:end=j+1;break
  social_nw=vnw[start:end]+';'
  assert set(re.findall(r'(?<=[(,])([^():,;]+):',social_nw))==set(social)
  assert leaves(find(source_nodes,node))==set(social)
  assert splits(find(source_nodes,node),set(social))==splits(parse_tree(social_nw),set(social))
  info['rooted_topology_matches']=True
  (a.output/'social_wasps_dated_from_Vespidae.nwk').write_text(social_nw+'\n')
  with (a.output/'N13_social_wasps_counts.tsv').open('w') as f:
   w=csv.writer(f,delimiter='\t');w.writerow(['Desc','Family ID']+social)
   for row in counts:w.writerow(row[:2]+[row[2+taxa.index(t)] for t in social])
report['interpretation']='N10 is the MRCA of Vespidae and the requested outgroups in the source OrthoFinder tree. N13 outgroup entries are outside HOG scope, not observed absences. No TE or large-family filtering applied; empirical fits are methodological tests.'
(a.output/'input_audit.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps({k:{x:v[x] for x in ['source_rows','observed_rows_in_selected_taxa','max_count']} for k,v in report['tables'].items()},indent=2))
