#!/usr/bin/env python3
from __future__ import annotations
import csv,json,os
from pathlib import Path
ROOT=Path(os.environ.get('IS_ANALYSIS_ROOT','/srv/is-analysis'))
SRC=ROOT/'results/metabolic_resilience/stage5_koges/KOGES_ANALYSIS_TABLE.tsv'
OUTDIR=ROOT/'results/metabolic_resilience/stage5_koges'; OUTDIR.mkdir(parents=True,exist_ok=True)
OUT=OUTDIR/'STAGE5G_MISSINGNESS_QC.tsv'; SUMMARY=OUTDIR/'STAGE5G_MISSINGNESS_QC.json'
CORE=['participant_id','time_years','MBI','incident_mets','event_time','age','sex']
if not SRC.exists():
    p={'status':'HOLD_ANALYSIS_TABLE_MISSING','expected':str(SRC)}; SUMMARY.write_text(json.dumps(p,indent=2)+'\n'); print(json.dumps(p,indent=2)); raise SystemExit(2)
rows=list(csv.DictReader(SRC.open('r',encoding='utf-8'),delimiter='\t'))
if not rows: raise SystemExit('empty analysis table')
cols=list(rows[0])
summary=[]
for c in cols:
    miss=sum(str(r.get(c,'')).strip() in {'','NA','NaN','nan','None'} for r in rows)
    summary.append({'variable':c,'n':len(rows),'missing_n':miss,'missing_pct':round(100*miss/len(rows),3),'core_variable':int(c in CORE)})
with OUT.open('w',encoding='utf-8',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(summary[0]),delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(summary)
core_bad=[r for r in summary if r['core_variable']==1 and r['missing_pct']>20]
p={'status':'REVIEW' if core_bad else 'PASS','rows':len(rows),'variables':len(cols),'core_gt20pct_missing':[r['variable'] for r in core_bad],'rule':'Do not silently impute. Multiple imputation, complete-case, or variable-specific handling must be declared before final modeling.'}
SUMMARY.write_text(json.dumps(p,indent=2)+'\n'); print(json.dumps(p,indent=2))
