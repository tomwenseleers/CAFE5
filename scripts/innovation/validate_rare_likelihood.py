"""Analytic observed-family likelihood checks far below machine epsilon inclusion."""
import csv,json,math,subprocess,sys
from pathlib import Path
binary=Path(sys.argv[1]).resolve();out=Path(sys.argv[2]);out.mkdir(parents=True,exist_ok=True)
tree=out/'star.nwk';tree.write_text('(A:1,B:1,C:1);\n');counts=out/'counts.tsv';counts.write_text('Desc\tFamily ID\tA\tB\tC\nx\ta\t1\t1\t1\nx\tb\t1\t0\t0\n')
results=[]
for r in [1e-2,.00526,.00527,1e-8,1e-16,1e-30,1e-120]:
 nu=.3*r;prefix=out/('mean_'+str(r));cmd=[str(binary),'--innovation','-t',str(tree),'-i',str(counts),'--root-mean',str(r),'--lambda','0','--mu','0','--nu',str(nu),'--epsilon','0','--max-count','35','--starts','1','--likelihood-only','-o',str(prefix)];subprocess.run(cmd,capture_output=True,check=True)
 fit={x['parameter']:x['value'] for x in csv.DictReader(Path(str(prefix)+'_results.tsv').open(),delimiter='\t')}
 inc=-math.expm1(-r-3*nu);expected=-math.log((r+nu**3)/inc)-math.log(nu/inc)+2*(r+3*nu)
 error=abs(float(fit['negative_log_likelihood'])-expected);assert error<1e-10,(r,error)
 assert abs(float(fit['log_inclusion_probability'])-math.log(inc))<1e-10
 results.append({'root_mean':r,'inclusion_probability':inc,'nll_expected':expected,'nll_actual':float(fit['negative_log_likelihood']),'absolute_error':error})
(out/'validation.json').write_text(json.dumps(results,indent=2)+'\n');print(json.dumps(results,indent=2))

# A broad critical birth/death tail must not require enumerating all possible
# observed counts merely to normalize inclusion. Root is exactly one.
tree=out/'broad_tree.nwk';tree.write_text('(A:110,B:110);\n')
counts=out/'broad_counts.tsv';counts.write_text('Desc\tFamily ID\tA\tB\nx\tlarge\t145\t1\n')
prefix=out/'broad';cmd=[str(binary),'--innovation','-t',str(tree),'-i',str(counts),'--root-family','hurdle-poisson','--root-mean','0','--root-zero','0','--lambda','.3','--nu','0','--epsilon','0','--max-count','160','--likelihood-only','-o',str(prefix)]
subprocess.run(cmd,capture_output=True,check=True)
fit={x['parameter']:x['value'] for x in csv.DictReader(Path(str(prefix)+'_results.tsv').open(),delimiter='\t')}
q=33/34;inc=1-q*q;expected=-math.log((1/34**2)**2*q**144/inc)
assert abs(float(fit['negative_log_likelihood'])-expected)<1e-10
assert fit['truncation_pass']=='1'
(out/'broad_tail_validation.json').write_text(json.dumps({'expected_nll':expected,'actual_nll':float(fit['negative_log_likelihood']),'inclusion_probability':inc,'cap':160,'truncation_difference':float(fit['truncation_nll_difference'])},indent=2)+'\n')
