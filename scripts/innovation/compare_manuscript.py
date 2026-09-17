#!/usr/bin/env python3
"""Map manuscript N0 HOGs to N11 by species-qualified gene membership."""
import argparse,csv,collections,json,hashlib
from pathlib import Path

def read(p):
 with p.open() as f:return list(csv.DictReader(f,delimiter='\t'))
def write(p,rows):
 if not rows:return
 with p.open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def main():
 p=argparse.ArgumentParser();p.add_argument('repository',type=Path);p.add_argument('output',type=Path);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
 base=a.repository/'nextflow_runs/2_EXCON/2_EXCON_orthofinder_eggnogmapper_run/results_EXCON/orthofinder/Phylogenetic_Hierarchical_Orthogroups'
 n0path=a.repository/'nextflow_runs/2_EXCON/3_EXCON_CAFE_run/results_EXCON/cafe/base/N0.tsv'
 n11=read(base/'N11.tsv');n0=read(n0path)
 n0=[{(k+'.clean' if k not in ['HOG','OG','Gene Tree Parent Clade'] and not k.endswith('.clean') else k):v for k,v in r.items()} for r in n0]
 species=[s for s in n11[0] if s.endswith('.clean') and any(r[s].strip() for r in n11)]
 def genes(r):return {(s,g.strip()) for s in species for g in r[s].split(',') if g.strip()}
 memberships={r['HOG']:genes(r) for r in n11}
 lookup=collections.defaultdict(set)
 for h,gs in memberships.items():
  for g in gs:lookup[g].add(h)
 oldgenes={r['HOG']:genes(r) for r in n0}
 mapping=[]
 for r in n0:
  gs=oldgenes[r['HOG']];hits=collections.Counter(h for g in gs for h in lookup[g])
  for h,overlap in hits.items():
   mapping.append(dict(N0_HOG=r['HOG'],N11_HOG=h,shared_genes=overlap,N0_genes_in_Vespidae=len(gs),N11_genes=len(memberships[h]),fraction_N0=overlap/len(gs),fraction_N11=overlap/len(memberships[h]),exact_Vespidae_membership=gs==memberships[h]))
 write(a.output/'N0_N11_membership_crosswalk.tsv',mapping)
 byold=collections.defaultdict(list)
 for r in mapping:byold[r['N0_HOG']].append(r)
 original=read(a.repository/'output/focal_node_significant_changes_annotated.tsv')
 comparisons=[]
 for r in original:
  for m in byold[r['HOG']] or [dict(N11_HOG='UNMAPPED',shared_genes=0,fraction_N0=0,fraction_N11=0,exact_Vespidae_membership=False)]:
   comparisons.append(dict(old_node=r['node'],new_node={'20':'Node2','25':'Node12'}[r['node']],old_HOG=r['HOG'],N11_HOG=m['N11_HOG'],old_direction=r['direction'],old_change=r['change'],old_family_p=r['family_p'],old_branch_p=r['branch_p'],annotation=r['consensus_gene_name'],TE_related_HOG_flag=r['TE_related_HOG_flag'],old_result_source=r['result_source'],shared_genes=m['shared_genes'],fraction_old=m['fraction_N0'],fraction_new=m['fraction_N11'],exact_membership=m['exact_Vespidae_membership']))
 write(a.output/'manuscript_focal_HOG_crosswalk.tsv',comparisons)
 annotations=read(a.repository/'input_CAFE/annotation/HOG_node_changes_and_annotations.tsv')
 print('Annotation columns:',list(annotations[0]))
 manifest={'species':species,'N11_families':len(n11),'N0_families':len(n0),'old_focal_events':len(original),'mapped_focal_rows':len(comparisons),'source_sha256':{str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in [n0path,base/'N11.tsv',a.repository/'output/focal_node_significant_changes_annotated.tsv']}}
 (a.output/'comparison_input_audit.json').write_text(json.dumps(manifest,indent=2)+'\n')
 print(json.dumps({k:v for k,v in manifest.items() if k not in ['source_sha256','species']}))
if __name__=='__main__':main()
