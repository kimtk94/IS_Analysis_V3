#!/usr/bin/env python3
"""Prepare variant-level artifact audit for cis-pQTL instruments."""
from __future__ import annotations
import csv,json,os
from pathlib import Path
ROOT=Path(os.environ.get('IS_ANALYSIS_ROOT','/srv/is-analysis'))
INST=ROOT/'results/metabolic_resilience/stage3_full_pgwas/instruments/STAGE3C2_FROZEN_INSTRUMENTS_ALL.tsv'
OUTDIR=ROOT/'results/metabolic_resilience/stage3_variant_artifact'; OUTDIR.mkdir(parents=True,exist_ok=True)
MAN=OUTDIR/'STAGE3F0_VARIANT_ANNOTATION_MANIFEST.tsv'; REVIEW=OUTDIR/'STAGE3F0_VARIANT_ARTIFACT_REVIEW.tsv'; SUMMARY=OUTDIR/'STAGE3F0_VARIANT_ARTIFACT_SUMMARY.json'
if not INST.exists():
    p={'status':'HOLD_INSTRUMENTS_MISSING'}; SUMMARY.write_text(json.dumps(p,indent=2)+'\n'); print(json.dumps(p,indent=2)); raise SystemExit(2)
inst=list(csv.DictReader(INST.open('r',encoding='utf-8'),delimiter='\t'))
uniq={}
for r in inst:
    key=r['variant_id']; uniq[key]={
      'variant_id':key,'rsid':r.get('rsid',''),'chrom_hg19':r.get('chrom_hg19',''),'pos_hg19':r.get('pos_hg19',''),
      'gene':r.get('gene_symbol',''),'ref':'','alt':'','most_severe_consequence':'','consequence_terms':'','impact':'',
      'protein_altering':'','missense':'','splice':'','coding':'','annotation_source':'','annotation_version':'','status':'REVIEW_ANNOTATION'
    }
if not MAN.exists():
    fields=list(next(iter(uniq.values())).keys()) if uniq else ['variant_id']
    with MAN.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(uniq.values())
    p={'status':'HOLD_ANNOTATION_REQUIRED','variants':len(uniq),'manifest':str(MAN),'instruction':'Annotate with reviewed VEP/Ensembl output; do not infer protein-altering status from position alone.'}
    SUMMARY.write_text(json.dumps(p,indent=2)+'\n'); print(json.dumps(p,indent=2)); raise SystemExit(2)
ann=list(csv.DictReader(MAN.open('r',encoding='utf-8'),delimiter='\t'))
rows=[]
for r in ann:
    def truth(x): return str(x).strip().lower() in {'1','true','yes','y'}
    flagged=truth(r.get('protein_altering')) or truth(r.get('missense')) or truth(r.get('splice')) or truth(r.get('coding'))
    rows.append({**r,'artifact_flag':int(flagged),'recommended_sensitivity':'EXCLUDE_AND_RERUN_MR_COLOC' if flagged else 'KEEP_PRIMARY'})
fields=list(rows[0].keys()) if rows else ['variant_id']
with REVIEW.open('w',encoding='utf-8',newline='') as f:
    w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(rows)
p={'status':'PASS','variants':len(rows),'artifact_flagged':sum(int(r['artifact_flag']) for r in rows),'rule':'Protein-altering/coding/splice cis-pQTLs are flagged for exclusion sensitivity because assay-binding or abundance interpretation may differ.'}
SUMMARY.write_text(json.dumps(p,indent=2)+'\n'); print(json.dumps(p,indent=2))
