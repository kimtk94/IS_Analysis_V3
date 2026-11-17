#!/usr/bin/env python3
"""Fail-closed independent scientific gate audit, executed only on GitHub Actions."""
import json,os,sys,urllib.request,hashlib,datetime
from pathlib import Path
TASKS={
"is":("https://www.ebi.ac.uk/gwas/rest/api/studies/GCST90104540","GCST90104540","EUR/EAS allele-build harmonization with verified full-file MD5 and matched LD"),
"glaucoma":("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE313736","GSE313736","paired integer count donor design and leave-pool-out"),
"hcc":("https://pdc.cancer.gov/pdc/","PDC","official TMT channel-to-aliquot-to-case one-to-one mapping"),
"radiomics":("https://ibsi.readthedocs.io/en/latest/","IBSI","patient-level NIfTI/RNA/clinical mapping and spacing"),
"ocular":("https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA1243143&result=read_run&fields=run_accession,sample_accession,experiment_accession,fastq_md5","run_accession","actual FASTQ MD5 and timepoint/strain identity"),
"muscle":("https://gtexportal.org/api/v2/reference/tissueSiteDetail?datasetId=gtex_v10","Muscle","20 sentinel loci with ancestry-matched LD and cis-eQTL"),
"oral":("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103322","GSE103322","unique donor-level expression and clinical endpoint mapping"),
"skin":("https://www.ebi.ac.uk/gwas/rest/api/studies","_links","phenotype-specific summary-statistics build/allele/IV QC")}
def main():
 key=sys.argv[1];url,token,gate=TASKS[key]
 result={"research":key,"date":"2026-11-18","source_url":url,"scientific_status":"HOLD","blocking_gate":gate,"source_status":"SOURCE_UNAVAILABLE"}
 try:
  req=urllib.request.Request(url,headers={"User-Agent":"master-research-ci-gate/1.0","Accept":"*/*"})
  with urllib.request.urlopen(req,timeout=35) as response:
   body=response.read(2_000_001);result["http_status"]=response.status
  result["bytes_read"]=len(body);result["sha256_bounded_response"]=hashlib.sha256(body).hexdigest()
  result["source_status"]="SOURCE_REACHABLE" if token.lower().encode() in body.lower() else "IDENTITY_UNVERIFIED"
  if len(body)>2_000_000: result["source_status"]="RESPONSE_TRUNCATED"
  if key=="ocular" and result["source_status"]=="SOURCE_REACHABLE":
   lines=body.decode("utf-8",errors="replace").splitlines()
   header=lines[0].split("\t");rows=[dict(zip(header,line.split("\t"))) for line in lines[1:] if line.strip()]
   runs=[r.get("run_accession") for r in rows]
   result["run_count"]=len(runs);result["unique_run_count"]=len(set(runs))
   result["declared_md5_count"]=sum(bool(r.get("fastq_md5")) for r in rows)
   if len(runs)!=46 or len(set(runs))!=46:result["source_status"]="RUN_COUNT_MISMATCH"
 except Exception as exc:
  result["block_reason"]=type(exc).__name__+": "+str(exc)[:240]
 Path("results").mkdir(exist_ok=True)
 Path("results",key+".json").write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n")
 print(json.dumps(result,ensure_ascii=False))
if __name__=="__main__":main()
