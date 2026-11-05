#!/usr/bin/env python3
from __future__ import annotations
import csv,json,math,os,gzip
from pathlib import Path
ROOT=Path(os.environ.get('IS_ANALYSIS_ROOT','/srv/is-analysis'))
H=ROOT/'results/metabolic_resilience/stage3_confirmatory_mr/harmonized/STAGE3D2_HARMONIZED.tsv.gz'
OUTDIR=ROOT/'results/metabolic_resilience/stage3_confirmatory_mr'; OUTDIR.mkdir(parents=True,exist_ok=True)
OUT=OUTDIR/'STAGE3D7_STEIGER_SENSITIVITY.tsv'; SUMMARY=OUTDIR/'STAGE3D7_STEIGER_SENSITIVITY.json'
if not H.exists():
    payload={'status':'HOLD_HARMONIZED_INPUT_MISSING'}; SUMMARY.write_text(json.dumps(payload,indent=2)+'\n'); print(json.dumps(payload,indent=2)); raise SystemExit(2)

def f(x):
    try:return float(x)
    except:return None

def r2_quant(beta,eaf,n=None):
    if beta is None or eaf is None:return None
    return 2*eaf*(1-eaf)*(beta**2)

rows=[]
with gzip.open(H,'rt',encoding='utf-8',newline='') as fh:
    for r in csv.DictReader(fh,delimiter='\t'):
        bx=f(r.get('beta_exposure')); by=f(r.get('beta_outcome')); eaf=f(r.get('exposure_eaf'))
        if None in (bx,by,eaf): continue
        if r.get('outcome_type')=='cc':
            status='HOLD_BINARY_LIABILITY_ASSUMPTIONS'; rx=ry=None
        else:
            rx=r2_quant(bx,eaf); ry=r2_quant(by,eaf); status='PASS' if rx is not None and ry is not None else 'REVIEW'
        rows.append({
            'gene':r.get('gene'),'mode':r.get('mode'),'trait':r.get('trait'),'variant_id':r.get('variant_id'),
            'r2_exposure_approx':rx,'r2_outcome_approx':ry,
            'steiger_support_exposure_to_outcome': int(rx>ry) if status=='PASS' else '',
            'status':status,
            'note':'Approximate continuous-trait Steiger sensitivity; not a substitute for exact liability-scale calculation for binary outcomes.'
        })
fields=['gene','mode','trait','variant_id','r2_exposure_approx','r2_outcome_approx','steiger_support_exposure_to_outcome','status','note']
with OUT.open('w',encoding='utf-8',newline='') as f2:
    w=csv.DictWriter(f2,fieldnames=fields,delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(rows)
payload={'status':'PASS','rows':len(rows),'continuous_pass':sum(r['status']=='PASS' for r in rows),'binary_hold':sum(r['status'].startswith('HOLD_BINARY') for r in rows)}
SUMMARY.write_text(json.dumps(payload,indent=2)+'\n'); print(json.dumps(payload,indent=2))
