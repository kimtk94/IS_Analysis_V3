#!/usr/bin/env python3
from __future__ import annotations
import csv,json,os
from pathlib import Path
ROOT=Path(os.environ.get('IS_ANALYSIS_ROOT','/srv/is-analysis'))
SHORT=ROOT/'results/metabolic_resilience/stage2_gwas/final_shortlist/STAGE2_FINAL_SCREENING_SHORTLIST.tsv'
OUTDIR=ROOT/'results/metabolic_resilience/stage6_functional'; OUTDIR.mkdir(parents=True,exist_ok=True)
OUT=OUTDIR/'STAGE6A_FUNCTIONAL_ANNOTATION_REGISTRY.tsv'; SUMMARY=OUTDIR/'STAGE6A_FUNCTIONAL_ANNOTATION_REGISTRY.json'
if not SHORT.exists():
    p={'status':'HOLD_SHORTLIST_MISSING'}; SUMMARY.write_text(json.dumps(p,indent=2)+'\n'); print(json.dumps(p,indent=2)); raise SystemExit(2)
rows=[]
for r in csv.DictReader(SHORT.open('r',encoding='utf-8'),delimiter='\t'):
    gene=r['gene_symbol']; rows.append({
      'gene':gene,'protein_id':r['protein_id'],
      'HPA_kidney_cell_type':'','HPA_expression':'','kidney_scRNA_dataset':'','kidney_scRNA_cell_type':'','kidney_scRNA_specificity':'',
      'GTEx_kidney_expression':'','drug_target_evidence':'','known_pathway':'','annotation_source_versions':'','status':'REVIEW_ANNOTATION'
    })
fields=list(rows[0].keys()) if rows else ['gene']
with OUT.open('w',encoding='utf-8',newline='') as f:
    w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(rows)
p={'status':'HOLD_ANNOTATION_REQUIRED','genes':len(rows),'registry':str(OUT),'rule':'Functional annotation is supportive and must not override MR/coloc evidence.'}
SUMMARY.write_text(json.dumps(p,indent=2)+'\n'); print(json.dumps(p,indent=2))
