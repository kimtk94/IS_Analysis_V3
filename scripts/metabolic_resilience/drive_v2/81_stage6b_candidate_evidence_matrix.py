#!/usr/bin/env python3
"""Assemble a descriptive evidence matrix. Does not rank candidates."""
from __future__ import annotations
import csv,json,os
from pathlib import Path
ROOT=Path(os.environ.get('IS_ANALYSIS_ROOT','/srv/is-analysis'))
OUTDIR=ROOT/'results/metabolic_resilience/stage6_functional'; OUTDIR.mkdir(parents=True,exist_ok=True)
SHORT=ROOT/'results/metabolic_resilience/stage2_gwas/final_shortlist/STAGE2_FINAL_SCREENING_SHORTLIST.tsv'
MR=ROOT/'results/metabolic_resilience/stage3_confirmatory_mr/STAGE3D5_CANDIDATE_CONFIRMATORY_SUMMARY.tsv'
COLOC=ROOT/'results/metabolic_resilience/stage3_coloc/STAGE3E3_COLOC_SUMMARY.tsv'
EAS=ROOT/'results/metabolic_resilience/stage4_eas/STAGE4E_CROSS_ANCESTRY_SUMMARY.tsv'
FUNC=OUTDIR/'STAGE6A_FUNCTIONAL_ANNOTATION_REGISTRY.tsv'
OUT=OUTDIR/'STAGE6B_CANDIDATE_EVIDENCE_MATRIX.tsv'; SUMMARY=OUTDIR/'STAGE6B_CANDIDATE_EVIDENCE_MATRIX.json'
if not SHORT.exists():
    p={'status':'HOLD_SHORTLIST_MISSING'}; SUMMARY.write_text(json.dumps(p,indent=2)+'\n'); print(json.dumps(p,indent=2)); raise SystemExit(2)
short=list(csv.DictReader(SHORT.open('r',encoding='utf-8'),delimiter='\t'))
mr=list(csv.DictReader(MR.open('r',encoding='utf-8'),delimiter='\t')) if MR.exists() else []
co=list(csv.DictReader(COLOC.open('r',encoding='utf-8'),delimiter='\t')) if COLOC.exists() else []
eas=list(csv.DictReader(EAS.open('r',encoding='utf-8'),delimiter='\t')) if EAS.exists() else []
func=list(csv.DictReader(FUNC.open('r',encoding='utf-8'),delimiter='\t')) if FUNC.exists() else []
rows=[]
for s in short:
    g=s['gene_symbol']; m=[r for r in mr if r.get('gene')==g and r.get('mode')=='primary']; c=[r for r in co if r.get('gene')==g]; e=[r for r in eas if r.get('gene')==g]; f=[r for r in func if r.get('gene')==g]
    rows.append({
      'gene':g,'protein_id':s['protein_id'],'stage2_domains':s.get('n_domains_supported',''),'stage2_core_domains':s.get('n_core_domains',''),'T2D_validation_class':s.get('T2D_validation_class',''),
      'confirmatory_metabolic_fdr_lt_0p05':sum(int(float(r.get('fdr_lt_0p05',0))) for r in m if r.get('fdr_lt_0p05','')!=''),
      'coloc_strong_pairs':sum(str(r.get('evidence',''))=='STRONG_H4' for r in c),
      'eas_concordant_significant':sum(r.get('replication_class')=='CONCORDANT_SIGNIFICANT' for r in e),
      'eas_concordant_nonsignificant':sum(r.get('replication_class')=='CONCORDANT_NONSIGNIFICANT' for r in e),
      'functional_annotation_status':f[0].get('status','') if f else 'PENDING',
      'note':'Descriptive evidence matrix only; no automatic ranking or final causal verdict.'
    })
fields=list(rows[0]) if rows else ['gene']
with OUT.open('w',encoding='utf-8',newline='') as fh:
    w=csv.DictWriter(fh,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
p={'status':'PASS','candidates':len(rows),'rule':'No automatic ranking. Evidence dimensions remain separate.'}; SUMMARY.write_text(json.dumps(p,indent=2)+'\n'); print(json.dumps(p,indent=2))
