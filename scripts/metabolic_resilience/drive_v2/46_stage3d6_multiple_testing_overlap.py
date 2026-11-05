#!/usr/bin/env python3
from __future__ import annotations
import csv,json,math,os
from pathlib import Path
ROOT=Path(os.environ.get('IS_ANALYSIS_ROOT','/srv/is-analysis'))
MR=ROOT/'results/metabolic_resilience/stage3_confirmatory_mr/STAGE3D5_PRIMARY_MR_SUMMARY.tsv'
REG=ROOT/'results/metabolic_resilience/stage3_confirmatory_mr/STAGE3D0_OUTCOME_REGISTRY.tsv'
OUTDIR=ROOT/'results/metabolic_resilience/stage3_confirmatory_mr'
OUTDIR.mkdir(parents=True,exist_ok=True)
OUT=OUTDIR/'STAGE3D6_MULTIPLE_TESTING_AND_OVERLAP.tsv'
SUMMARY=OUTDIR/'STAGE3D6_MULTIPLE_TESTING_AND_OVERLAP.json'

def fnum(x):
    try:
        v=float(x); return v if math.isfinite(v) else None
    except Exception:return None

def bh(vals):
    vals=[(i,p) for i,p in vals if p is not None and math.isfinite(p)]
    vals.sort(key=lambda x:x[1]); n=len(vals); out={}; prev=1.0
    for rank in range(n,0,-1):
        i,p=vals[rank-1]; q=min(prev,p*n/rank); out[i]=q; prev=q
    return out

if not MR.exists() or not REG.exists():
    payload={'status':'HOLD_MISSING_INPUT','mr_exists':MR.exists(),'registry_exists':REG.exists()}
    SUMMARY.write_text(json.dumps(payload,indent=2)+'\n'); print(json.dumps(payload,indent=2)); raise SystemExit(2)

mr=list(csv.DictReader(MR.open('r',encoding='utf-8'),delimiter='\t'))
reg=list(csv.DictReader(REG.open('r',encoding='utf-8'),delimiter='\t'))
reg_by_trait={r['trait']:r for r in reg if r.get('trait')}
rows=[]
for r in mr:
    if r.get('mode') not in ('','primary'): continue
    p=fnum(r.get('p')); trait=r.get('trait',''); domain=r.get('domain','')
    rr=dict(r); meta=reg_by_trait.get(trait,{})
    rr['sample_overlap_class']=meta.get('sample_overlap_class',meta.get('sample_overlap_with_UKBPPP','REVIEW'))
    rr['sample_overlap_note']=meta.get('sample_overlap_note','')
    rr['testing_family']='T2D_VALIDATION' if domain=='disease_validation' else domain.upper()
    rr['p_numeric']=p
    rows.append(rr)

for fam in sorted({r['testing_family'] for r in rows}):
    idx=[i for i,r in enumerate(rows) if r['testing_family']==fam]
    q=bh([(i,rows[i]['p_numeric']) for i in idx])
    for i in idx: rows[i]['FDR_family']=q.get(i)
met_idx=[i for i,r in enumerate(rows) if r['testing_family']!='T2D_VALIDATION']
qg=bh([(i,rows[i]['p_numeric']) for i in met_idx])
for i,r in enumerate(rows): r['FDR_metabolic_global']=qg.get(i) if i in qg else None

fields=['gene','trait','domain','testing_family','method','k','beta','se','p','FDR_family','FDR_metabolic_global','sample_overlap_class','sample_overlap_note']
with OUT.open('w',encoding='utf-8',newline='') as f:
    w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n',extrasaction='ignore'); w.writeheader(); w.writerows(rows)

payload={
 'status':'PASS','rows':len(rows),
 'families':sorted({r['testing_family'] for r in rows}),
 'rule':'Primary multiplicity: BH-FDR within prespecified metabolic domain; global metabolic FDR reported as sensitivity. T2D remains a separate disease-validation family.',
 'overlap_rule':'Outcome sample overlap is metadata and must not be silently inferred; REVIEW blocks high-confidence interpretation.'
}
SUMMARY.write_text(json.dumps(payload,indent=2)+'\n'); print(json.dumps(payload,indent=2))
