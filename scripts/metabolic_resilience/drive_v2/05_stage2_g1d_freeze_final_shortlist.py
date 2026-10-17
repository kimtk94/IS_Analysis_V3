#!/usr/bin/env python3
from pathlib import Path
import csv, json

ROOT=Path('/srv/is-analysis')
SHORT=ROOT/'results/metabolic_resilience/stage2_gwas/glycemia_screen/STAGE2F1_METABOLIC_SHORTLIST_GE3.tsv'
T2D=ROOT/'results/metabolic_resilience/stage2_gwas/t2d_screen/STAGE2G1C_T2D_SHORTLIST_VALIDATION.tsv'
OUTDIR=ROOT/'results/metabolic_resilience/stage2_gwas/final_shortlist'
OUTDIR.mkdir(parents=True,exist_ok=True)
OUT=OUTDIR/'STAGE2_FINAL_SCREENING_SHORTLIST.tsv'
TARGETS=OUTDIR/'STAGE3_FULL_PGWAS_TARGETS.txt'
SUMMARY=OUTDIR/'STAGE2_FINAL_SUMMARY.json'

with SHORT.open() as f: srows=list(csv.DictReader(f,delimiter='\t'))
with T2D.open() as f: trows=list(csv.DictReader(f,delimiter='\t'))
by_pid={r['protein_id']:r for r in trows}

rows=[]
for s in srows:
    pid=s['protein_id']; t=by_pid.get(pid,{})
    r=dict(s)
    for k,v in t.items():
        if k not in r or r[k]=='': r[k]=v
    r['stage3_confirmatory']='1'
    rows.append(r)

# Operational ordering only; not a biological rank.
order={'T1_CORRECTED_VALIDATION':1,'T2_NOMINAL_SUPPORT':2,'T3_DIRECTION_ONLY':3,'T4_OPPOSITE_DIRECTION':4,'T0_UNAVAILABLE':5}
rows.sort(key=lambda r:(order.get(r.get('T2D_validation_class',''),99),-int(float(r.get('n_core_domains',0) or 0)),r.get('gene_symbol','')))

fields=[]
for r in rows:
    for k in r:
        if k not in fields: fields.append(k)
with OUT.open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
TARGETS.write_text('\n'.join(r['protein_id'] for r in rows)+'\n')
summary={
    'stage2_status':'LOCKED',
    'n_metabolic_shortlist':len(rows),
    'criterion':'support in >=3 of 4 metabolic domains',
    'n_ge2_core':sum(int(float(r.get('n_core_domains',0) or 0))>=2 for r in rows),
    'T2D_T1_corrected':sum(r.get('T2D_validation_class')=='T1_CORRECTED_VALIDATION' for r in rows),
    'T2D_T2_nominal':sum(r.get('T2D_validation_class')=='T2_NOMINAL_SUPPORT' for r in rows),
    'T2D_T3_direction':sum(r.get('T2D_validation_class')=='T3_DIRECTION_ONLY' for r in rows),
    'T2D_T4_opposite_nonsignificant':sum(r.get('T2D_validation_class')=='T4_OPPOSITE_DIRECTION' for r in rows),
    'next_stage':'candidate-only full UKB-PPP pGWAS; cis marginal variants; LD clumping; multi-SNP MR; colocalization'
}
SUMMARY.write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
