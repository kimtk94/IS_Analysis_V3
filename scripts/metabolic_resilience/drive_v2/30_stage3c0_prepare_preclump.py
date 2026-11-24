#!/usr/bin/env python3
from pathlib import Path
import csv, gzip, json

ROOT=Path('/srv/is-analysis')
CIS=ROOT/'results/metabolic_resilience/stage3_full_pgwas/cis_marginal'
OUT=ROOT/'results/metabolic_resilience/stage3_full_pgwas/ld_clump'
AUDIT=ROOT/'results/metabolic_resilience/stage3_full_pgwas/audit'
OUT.mkdir(parents=True,exist_ok=True)
summary=[]
for path in sorted(CIS.glob('*_cis1Mb_marginal.tsv.gz')):
    primary=[]; stringent=[]
    with gzip.open(path,'rt',encoding='utf-8',newline='') as f:
        for r in csv.DictReader(f,delimiter='\t'):
            info=float(r['info']); F=float(r['f_stat']); logp=float(r['log10p'])
            if info>=0.8 and F>10 and logp>=7.301029995663981: primary.append(r)
            if info>=0.8 and F>10 and logp>=10.769551078621726: stringent.append(r)
    gene=(primary[0]['gene_symbol'] if primary else path.name.split('_')[0])
    fields=['protein_id','gene_symbol','chrom_hg19','pos_hg19','variant_id','allele0','allele1','effect_allele','other_allele','eaf','info','beta','se','f_stat','log10p','p_approx']
    for name,rows in [('primary',primary),('stringent',stringent)]:
        p=OUT/f'{gene}_{name}_preclump.tsv'
        with p.open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n',extrasaction='ignore'); w.writeheader(); w.writerows(rows)
    summary.append({'gene':gene,'primary_preclump':len(primary),'stringent_preclump':len(stringent)})

p=AUDIT/'STAGE3C0_PRECLUMP_COUNTS.tsv'
with p.open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=summary[0].keys(),delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(summary)
js={'targets':len(summary),'primary_total':sum(x['primary_preclump'] for x in summary),'stringent_total':sum(x['stringent_preclump'] for x in summary),
    'primary_definition':'P<5e-8, F>10, INFO>=0.8','stringent_definition':'P<1.7e-11, F>10, INFO>=0.8'}
(AUDIT/'STAGE3C0_PRECLUMP_COUNTS.json').write_text(json.dumps(js,indent=2)+'\n')
print(json.dumps(js,indent=2))
