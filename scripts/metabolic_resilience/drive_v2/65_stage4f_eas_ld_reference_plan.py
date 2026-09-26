#!/usr/bin/env python3
from __future__ import annotations
import csv,json,os
from pathlib import Path
ROOT=Path(os.environ.get('IS_ANALYSIS_ROOT','/srv/is-analysis'))
REG=ROOT/'results/metabolic_resilience/stage4_eas/STAGE4D_EAS_RESOURCE_REGISTRY.tsv'
OUTDIR=ROOT/'results/metabolic_resilience/stage4_eas'; OUTDIR.mkdir(parents=True,exist_ok=True)
OUT=OUTDIR/'STAGE4F_EAS_LD_REFERENCE_PLAN.tsv'; SUMMARY=OUTDIR/'STAGE4F_EAS_LD_REFERENCE_PLAN.json'
if not REG.exists():
    p={'status':'HOLD_EAS_REGISTRY_MISSING'}; SUMMARY.write_text(json.dumps(p,indent=2)+'\n'); print(json.dumps(p,indent=2)); raise SystemExit(2)
rows=[]
for r in csv.DictReader(REG.open('r',encoding='utf-8'),delimiter='\t'):
    if r.get('role')!='outcome_GWAS': continue
    ancestry=r.get('ancestry','').upper(); resource=r.get('resource','')
    if 'KOREAN' in ancestry: ref='Korean-specific LD preferred; otherwise ancestry-matched EAS reference with sensitivity'
    elif 'JAPAN' in ancestry: ref='Japanese-specific LD preferred; otherwise ancestry-matched EAS reference with sensitivity'
    elif 'TAIWAN' in ancestry: ref='Taiwanese/Chinese-specific LD preferred; otherwise ancestry-matched EAS reference with sensitivity'
    else: ref='1000G EAS Phase 3 as fallback, with resource-specific LD preferred'
    rows.append({'resource':resource,'ancestry':r.get('ancestry'),'trait':r.get('trait'),'recommended_ld_reference':ref,'eur_ld_allowed_primary':0,'status':'REVIEW_REFERENCE_AVAILABILITY'})
fields=list(rows[0].keys()) if rows else ['resource']
with OUT.open('w',encoding='utf-8',newline='') as f:
    w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(rows)
p={'status':'PASS_PLAN','rows':len(rows),'locked_rule':'EUR LD must not be used as the primary clumping/coloc LD reference for EAS replication.'}
SUMMARY.write_text(json.dumps(p,indent=2)+'\n'); print(json.dumps(p,indent=2))
