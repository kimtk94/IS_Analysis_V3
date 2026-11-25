#!/usr/bin/env python3
"""CI-only, fail-closed 2026-11-11 scientific prerequisite gates.
PASS means the stated bounded metadata check only, NEVER scientific completion.
"""
import argparse,csv,hashlib,io,json,pathlib,re,urllib.request,datetime
TASKS={
"is":("GIGASTROKE EAS study accession and ancestry metadata","https://www.ebi.ac.uk/gwas/rest/api/studies/GCST90104540","GCST90104540","Full ancestry-specific GWAS and LD/harmonization"),
"glaucoma":("GSE313736 official supplementary/count file inventory","https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE313736","GSE313736","Raw integer counts and 9x leave-pool-out edgeR"),
"hcc":("PDC000198 official study identity prerequisite","https://pdc.cancer.gov/pdc/study/PDC000198","PDC000198","330 quantitative columns to official aliquots, RNA patient join"),
"radiomics":("IBSI benchmark provenance prerequisite","https://ibsi.readthedocs.io/en/latest/","IBSI","Source NIfTI spacing and clinical eligibility"),
"ocular":("ENA 46-run experiment/sample cardinality","https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA1243143&result=read_run&fields=run_accession,experiment_accession,sample_accession,fastq_md5,fastq_bytes&format=tsv","run_accession","Actual FASTQ checksum and timepoint/strain mapping"),
"muscle":("GTEx v10 skeletal-muscle reference identity","https://gtexportal.org/api/v2/reference/tissueSiteDetail?datasetId=gtex_v10","Muscle_Skeletal","20-locus harmonization, EAS LD and eQTL"),
"oral":("GSE103322 GEO patient/sample source prerequisite","https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103322","GSE103322","18-patient expression matrix and donor/endpoint QC"),
"skin":("Bristol perceived-age effect/allele schema","https://data.bris.ac.uk/datasets/21crwsnj4xwjm2g4qi8chathha/readme.txt","ALLELE1","Full sumstats build, instruments and pleiotropy")
}
p=argparse.ArgumentParser();p.add_argument("--track",choices=TASKS,required=True);a=p.parse_args()
task,url,token,nxt=TASKS[a.track]
out={"date":"2026-11-11","track":a.track,"task":task,"source":url,"status":"HOLD","scientific_status":"HOLD","next_gate":nxt,"utc":datetime.datetime.now(datetime.timezone.utc).isoformat()}
try:
 req=urllib.request.Request(url,headers={"User-Agent":"MasterDegree-ScientificGate-CI/1.0"})
 with urllib.request.urlopen(req,timeout=35) as resp:
  body=resp.read(524289);out["http_status"]=resp.status
 if len(body)>524288: raise ValueError("BOUNDED_RESPONSE_EXCEEDED_512K")
 out["bounded_bytes"]=len(body);out["bounded_sha256"]=hashlib.sha256(body).hexdigest()
 s=body.decode("utf-8","replace")
 if token.casefold() not in s.casefold(): raise ValueError("EXPECTED_IDENTITY_MISSING")
 if a.track=="ocular":
  rows=list(csv.DictReader(io.StringIO(s),delimiter="\t"))
  runs=[r.get("run_accession","") for r in rows]
  if len(rows)!=46 or len(set(runs))!=46: raise ValueError("RUN_CARDINALITY_NOT_46_UNIQUE")
  if any(not r.get("experiment_accession") or not r.get("sample_accession") for r in rows): raise ValueError("MISSING_SAMPLE_EXPERIMENT")
  out["run_count"]=len(rows);out["experiment_count"]=len(set(r["experiment_accession"] for r in rows));out["sample_count"]=len(set(r["sample_accession"] for r in rows))
  out["declared_fastq_md5_rows"]=sum(bool(r.get("fastq_md5")) for r in rows)
  if out["declared_fastq_md5_rows"]!=46: raise ValueError("MISSING_FASTQ_MD5_METADATA")
 if a.track=="skin":
  for field in ("ALLELE0","A1FREQ","LOG_OR","LOG_SE"):
   if field not in s: raise ValueError("MISSING_SCHEMA_"+field)
 out["status"]="PASS";out["scope"]="Bounded source/metadata only; no raw integrity or statistical inference"
except Exception as e:
 out["blocker"]=type(e).__name__+": "+str(e)[:250]
 out["ci_alternative"]="Retry bounded official metadata on Actions; do not bypass missing raw or mapping evidence"
path=pathlib.Path("results/master_20261111");path.mkdir(parents=True,exist_ok=True)
(path/(a.track+".json")).write_text(json.dumps(out,ensure_ascii=False,indent=2)+"\n")
print(json.dumps(out,ensure_ascii=False))
