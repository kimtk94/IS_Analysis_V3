#!/usr/bin/env python3
"""Eight independent bounded scientific input gates; no inference promotion."""
import argparse, csv, hashlib, io, json, os, pathlib, re, urllib.request, urllib.error
from datetime import datetime, timezone
P=pathlib.Path("results/master_20261110")
P.mkdir(parents=True,exist_ok=True)
TASKS={
"is":("GIGASTROKE accession identity","https://www.ebi.ac.uk/gwas/rest/api/studies/GCST90104540","GCST90104540","Full ancestry-matched GWAS checksum and allele harmonization"),
"glaucoma":("GSE313736 official sample metadata","https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE313736","GSE313736","Raw integer-count leave-pool-out QL rerun"),
"hcc":("PDC000198 source identity","https://pdc.cancer.gov/pdc/study/PDC000198","PDC000198","330 quantitative channels to official aliquot and RNA join"),
"radiomics":("IBSI benchmark provenance","https://ibsi.readthedocs.io/en/latest/","IBSI","NIfTI voxel spacing and clinical eligibility"),
"ocular":("ENA 46-run sample/experiment identity","https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA1243143&result=read_run&fields=run_accession,experiment_accession,sample_accession,fastq_md5,fastq_bytes&format=tsv","run_accession","FASTQ bytewise checksum and timepoint/strain map"),
"muscle":("GTEx v10 tissue reference","https://gtexportal.org/api/v2/reference/tissueSiteDetail?datasetId=gtex_v10","Muscle_Skeletal","20 loci build-aware harmonization and ancestry-matched LD"),
"oral":("GSE103322 cohort provenance","https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103322","GSE103322","Patient-resolved expression and clinical endpoint QC"),
"skin":("Bristol perceived-age allele schema","https://data.bris.ac.uk/datasets/21crwsnj4xwjm2g4qi8chathha/readme.txt","ALLELE1","Complete GWAS build, allele and instrument QC")
}
a=argparse.ArgumentParser();a.add_argument("--track",choices=TASKS,required=True);x=a.parse_args();name,url,token,next_gate=TASKS[x.track]
out={"date_kst":"2026-11-10","track":x.track,"task":name,"source":url,"next_gate":next_gate,"scientific_gate":"HOLD","source_gate":"HOLD","timestamp_utc":datetime.now(timezone.utc).isoformat(),"bytes_checked":0,"reason":None}
try:
 req=urllib.request.Request(url,headers={"User-Agent":"MasterDegree-CI-Research-Audit/1.0","Accept":"text/plain,text/html,application/json"})
 with urllib.request.urlopen(req,timeout=25) as f:
  data=f.read(262145)
  out["http_status"]=f.status
 if len(data)>262144: raise ValueError("Source exceeds bounded 256KiB inspection; full integrity not tested")
 out["bytes_checked"]=len(data);out["sha256_bounded"]=hashlib.sha256(data).hexdigest()
 s=data.decode("utf-8","replace")
 if token.casefold() not in s.casefold(): raise ValueError("Expected identifier absent from bounded response")
 if x.track=="ocular":
  rows=list(csv.DictReader(io.StringIO(s),delimiter="\t"))
  runs=[r.get("run_accession","") for r in rows]
  if len(rows)!=46 or len(set(runs))!=46: raise ValueError("ENA 46 unique run constraint failed")
  out["run_count"]=len(rows)
  out["unique_experiments"]=len(set(r.get("experiment_accession") for r in rows))
  out["unique_samples"]=len(set(r.get("sample_accession") for r in rows))
  if any(not r.get("fastq_md5") or not r.get("fastq_bytes") for r in rows): raise ValueError("Missing declared FASTQ checksum or byte count")
 out["source_gate"]="PASS";out["reason"]="Bounded official source identity verified; not raw-data integrity or biological inference"
except Exception as e:
 out["reason"]=type(e).__name__+": "+str(e)[:300]
(P/(x.track+".json")).write_text(json.dumps(out,indent=2,ensure_ascii=False)+"\n")
print(json.dumps(out,ensure_ascii=False))
# Expected source blocks are scientific HOLD, not CI failure. Unexpected code exceptions still fail.
