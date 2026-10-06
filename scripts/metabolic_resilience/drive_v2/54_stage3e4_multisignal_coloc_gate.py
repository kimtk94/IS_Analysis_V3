#!/usr/bin/env python3
from __future__ import annotations
import csv,json,os
from pathlib import Path
ROOT=Path(os.environ.get('IS_ANALYSIS_ROOT','/srv/is-analysis'))
CLUMP=ROOT/'results/metabolic_resilience/stage3_full_pgwas/audit/STAGE3C1_LD_CLUMP_SUMMARY.tsv'
COLOC=ROOT/'results/metabolic_resilience/stage3_coloc/STAGE3E3_COLOC_SUMMARY.tsv'
OUTDIR=ROOT/'results/metabolic_resilience/stage3_coloc'; OUTDIR.mkdir(parents=True,exist_ok=True)
OUT=OUTDIR/'STAGE3E4_MULTISIGNAL_DECISION.tsv'; SUMMARY=OUTDIR/'STAGE3E4_MULTISIGNAL_DECISION.json'
if not CLUMP.exists() or not COLOC.exists():
    p={'status':'HOLD_MISSING_INPUT','clump':CLUMP.exists(),'coloc':COLOC.exists()}; SUMMARY.write_text(json.dumps(p,indent=2)+'\n'); print(json.dumps(p,indent=2)); raise SystemExit(2)
cl=list(csv.DictReader(CLUMP.open('r',encoding='utf-8'),delimiter='\t'))
co=list(csv.DictReader(COLOC.open('r',encoding='utf-8'),delimiter='\t'))
k_by_gene={r['gene']:int(float(r['independent_n'])) for r in cl if r.get('mode')=='primary' and r.get('independent_n')}
rows=[]
for r in co:
    gene=r.get('gene'); h4=float(r.get('PP.H4','nan')) if r.get('PP.H4') not in ('','NA',None) else None
    h3=float(r.get('PP.H3','nan')) if r.get('PP.H3') not in ('','NA',None) else None
    k=k_by_gene.get(gene,0)
    if k>=2 and (h4 is None or h4<0.8 or (h3 is not None and h3>=h4)):
        decision='RUN_MULTISIGNAL_COLOC'
    elif k>=2:
        decision='OPTIONAL_MULTISIGNAL_SENSITIVITY'
    else:
        decision='ABF_SINGLE_SIGNAL_ACCEPTABLE'
    rows.append({'gene':gene,'trait':r.get('trait'),'primary_independent_iv_n':k,'PP.H3':h3,'PP.H4':h4,'decision':decision})
fields=list(rows[0].keys()) if rows else ['gene']
with OUT.open('w',encoding='utf-8',newline='') as f:
    w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(rows)
p={'status':'PASS','pairs':len(rows),'run_multisignal':sum(r['decision']=='RUN_MULTISIGNAL_COLOC' for r in rows),'optional':sum(r['decision']=='OPTIONAL_MULTISIGNAL_SENSITIVITY' for r in rows)}
SUMMARY.write_text(json.dumps(p,indent=2)+'\n'); print(json.dumps(p,indent=2))
