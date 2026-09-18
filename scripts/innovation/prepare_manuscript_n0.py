#!/usr/bin/env python3
"""Prepare all observed N0 families on the exact 17-tip manuscript tree."""
import argparse,csv,hashlib,json,re
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('repository',type=Path);p.add_argument('output',type=Path);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
base=a.repository/'nextflow_runs/2_EXCON/3_EXCON_CAFE_run/results_EXCON/cafe/base';source=base/'N0.tsv';tree=a.repository/'input_CAFE/tree_dating/Vespidae_with_outgroups_dated_primary.nwk'
nw=tree.read_text().strip();taxa=re.findall(r'(?<=[(,])([^():,;]+):',nw);assert len(taxa)==len(set(taxa))==17
if nw!=(base/'cafe_input_tree.txt').read_text().strip():raise ValueError('Dated tree differs from archived manuscript input; resolve provenance before an exact comparison')
with source.open() as f:raw=list(csv.DictReader(f,delimiter='\t'))
rows=[];empty=[];index={}
for r in raw:
 counts=[]
 for sp in taxa:
  genes=[g.strip() for g in r[sp].split(',') if g.strip()]
  if len(set(genes))!=len(genes):raise ValueError(f'Duplicate gene identifier within {r["HOG"]}, {sp}')
  counts.append(len(genes))
 if any(counts):rows.append(['N0',r['HOG']]+counts);index[r['HOG']]=dict(zip(taxa,counts))
 else:empty.append(r['HOG'])
with (a.output/'N0_counts.tsv').open('w') as f:
 w=csv.writer(f,delimiter='\t');w.writerow(['Desc','Family ID']+taxa);w.writerows(rows)
(a.output/'tree.nwk').write_text(nw+'\n')
original=set();discrepancies=[]
for name in ['hog_gene_counts.tsv','hog_gene_counts_large.tsv']:
 with (base/name).open() as f:
  for r in csv.DictReader(f,delimiter='\t'):
   original.add(r['HOG'])
   for sp in taxa:
    if int(r[sp])!=index[r['HOG']][sp]:discrepancies.append([r['HOG'],sp,r[sp],index[r['HOG']][sp]])
if discrepancies:raise ValueError('Counts differ from original preparation: '+str(discrepancies[:5]))
audit={'source_rows':len(raw),'observed_families':len(rows),'omitted_all_zero_rows':len(empty),'single_species_families':sum(sum(x>0 for x in r[2:])==1 for r in rows),'differential_gt20_families':sum(max(r[2:])-min(r[2:])>20 for r in rows),'maximum_count':max(max(r[2:]) for r in rows),'original_prepared_families':len(original),'shared_count_discrepancies':discrepancies,'extra_multispecies_families':[h for h,x in index.items() if h not in original and sum(v>0 for v in x.values())>1],'policy':'All observed families enter one fit, including TE-associated, single-species, large, and potentially root-zero families. TE exclusion applies only to separately labelled biological summaries.','sha256':{str(x.resolve()):hashlib.sha256(x.read_bytes()).hexdigest() for x in [source,tree,a.output/'N0_counts.tsv']}}
(a.output/'input_audit.json').write_text(json.dumps(audit,indent=2)+'\n');print(json.dumps({k:v for k,v in audit.items() if k!='sha256'},indent=2))
# Audit the broader family-construction ascertainment, without pretending to know
# an exact inclusion rule for the complete orthology pipeline.
species=[k for k in raw[0] if k not in ['HOG','OG','Gene Tree Parent Clade']]
one_gene=one_species=tip_one=outside=0
for r in raw:
 c={sp:len([g for g in r[sp].split(',') if g.strip()]) for sp in species}
 one_gene+=sum(c.values())==1;one_species+=sum(v>0 for v in c.values())==1
 if sum(c[sp] for sp in taxa)==1:
  tip_one+=1;outside+=sum(c[sp] for sp in species if sp not in taxa)>0
ascertainment={'source_species':len(species),'source_HOGs':len(raw),'source_HOGs_with_exactly_one_gene':one_gene,'source_HOGs_occupied_in_one_species':one_species,'retained_HOGs_with_exactly_one_gene_in_17tips':tip_one,'of_these_also_present_outside_17tips':outside,'interpretation':'Conditioning on nonzero observations in 17 tips does not reproduce the larger HOG-construction ascertainment. This may contribute to predictive mismatch; it is not a demonstrated exclusive cause or an exact characterization of OrthoFinder inclusion rules.'}
(a.output/'HOG_ascertainment_audit.json').write_text(json.dumps(ascertainment,indent=2)+'\n')
