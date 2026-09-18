#!/usr/bin/env python3
"""Audit N11 membership and split by original OG, without filtering any family."""
import argparse,csv,hashlib,json
from collections import Counter
from pathlib import Path
p=argparse.ArgumentParser();p.add_argument('--membership',type=Path,required=True);p.add_argument('--counts',type=Path,required=True);p.add_argument('--output',type=Path,required=True);a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
with a.membership.open() as h: members=list(csv.DictReader(h,delimiter='\t'))
lookup={r['HOG']:r for r in members}
with a.counts.open() as h: r=csv.reader(h,delimiter='\t');header=next(r);rows=list(r)
source_species=[k for k in members[0] if k.endswith('.clean')];selected=[t+'.clean' for t in header[2:]]
train=[];test=[];assign=[];outside=0;single=0;ogs=Counter();parents=Counter()
for row in rows:
 m=lookup[row[1]];ogs[m['OG']]+=1;parents[m['Gene Tree Parent Clade']]+=1
 computed=[len([g for g in m[t].split(',') if g.strip()]) for t in selected]
 assert computed==list(map(int,row[2:])),row[1]
 single+=sum(computed)==1
 outside+=sum(len([g for g in m[t].split(',') if g.strip()]) for t in source_species if t not in selected)
 # Keep all descendants of the same original orthogroup in one partition.
 fold=int(hashlib.sha256(('N11-heldout-20260918:'+m['OG']).encode()).hexdigest()[:8],16)%5
 which='test' if fold==0 else 'train';(test if fold==0 else train).append(row);assign.append([row[1],m['OG'],which])
for name,data,columns in [('train.tsv',train,header),('test.tsv',test,header),('assignment.tsv',assign,['HOG','OG','partition'])]:
 with (a.output/name).open('w') as h:w=csv.writer(h,delimiter='\t');w.writerow(columns);w.writerows(data)
summary={'families':len(rows),'training':len(train),'test':len(test),'original_OGs':len(ogs),'OGs_split_into_multiple_N11_HOGs':sum(n>1 for n in ogs.values()),'largest_OG_HOG_count':max(ogs.values()),'single_gene_N11_HOGs':single,'genes_outside_focal_clade':outside,'gene_tree_parent_clade_counts':dict(parents),'membership_sha256':hashlib.sha256(a.membership.read_bytes()).hexdigest(),'counts_sha256':hashlib.sha256(a.counts.read_bytes()).hexdigest(),'selection':'All original N11 families retained. Holdout assignment uses OG identity only, not observed copy counts.','ascertainment_limit':'N11 membership is scoped to Vespidae. Empty outside columns are not biological absences. Presence of single-gene HOGs rules out an assumed universal two-gene minimum at N11. Exact gene-tree/HOG detection selection remains unmodelled.'}
(a.output/'audit.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps({k:v for k,v in summary.items() if k!='gene_tree_parent_clade_counts'},indent=2))
