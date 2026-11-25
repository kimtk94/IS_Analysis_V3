#!/usr/bin/env python3
"""Independent Nov 4 follow-up source checks; run in GitHub Actions only."""
import argparse,csv,datetime,hashlib,io,json,os,pathlib,re,urllib.request,urllib.parse,urllib.error
OUT=pathlib.Path("results/research_operations/20261104")
URLS={
"is":"https://www.ebi.ac.uk/gwas/rest/api/v2/studies/GCST90104545",
"glaucoma":"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE313736&targ=self&form=text&view=quick",
"hcc":"https://pdc.cancer.gov/graphql",
"radiomics":"https://zenodo.org/api/records/8179129",
"ocular":"https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA1243143&result=read_run&fields=run_accession,sample_accession,experiment_accession,library_layout,fastq_ftp,fastq_md5&format=tsv",
"muscle":"https://gtexportal.org/api/v2/dataset/variant?snpId=rs12916&datasetId=gtex_v10",
"oral":"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103322&targ=self&form=text&view=quick",
"skin":"https://data.bris.ac.uk/datasets/21crwsnj4xwjm2g4qi8chathha/readme.txt"}
NEXT={
"is":"Full summary-statistics checksum and EUR/EAS allele harmonization",
"glaucoma":"Donor-level count matrix and independent 9x leave-pool-out QL",
"hcc":"Official 330 quantitative column to aliquot mapping and RNA join",
"radiomics":"14 original NIfTI headers and patient clinical eligibility",
"ocular":"Run to exposure time and strain pairing, full FASTQ MD5",
"muscle":"All 20 source loci EAS LD and muscle eQTL colocalization",
"oral":"Patient-level expression/clinical metadata and donor QC",
"skin":"Full GWAS checksum, genome build and instrument validity"}
def fetch(url,data=None):
 req=urllib.request.Request(url,data=data,headers={"User-Agent":"MASTER-DEGREE-CI-20261104","Accept":"*/*",**({"Content-Type":"application/json"} if data else {})})
 with urllib.request.urlopen(req,timeout=28) as f:
  b=f.read(2000001)
  if len(b)>2000000: raise RuntimeError("SOURCE_SIZE_LIMIT")
  return b
def task(k):
 url=URLS[k]; evidence={}
 if k=="hcc":
  for study in ("PDC000198","PDC000199"):
   q='{ biospecimenPerStudy(pdc_study_id:"'+study+'",acceptDUA:true){aliquot_id sample_id case_id sample_type} }'
   obj=json.loads(fetch(url,json.dumps({"query":q}).encode()))
   rows=obj.get("data",{}).get("biospecimenPerStudy") or []
   if not rows: raise RuntimeError("NO_PDC_ROWS")
   parent={}
   for r in rows:
    a=r.get("aliquot_id");p=(r.get("sample_id"),r.get("case_id"))
    if not a or not all(p) or (a in parent and parent[a]!=p):raise RuntimeError("ALIQUOT_PARENT_CONFLICT")
    parent[a]=p
   evidence[study]={"rows":len(rows),"unique_aliquots":len(parent),"unique_cases":len(set(x[1] for x in parent.values())),"mapping_sha256":hashlib.sha256(json.dumps(sorted((a,*p) for a,p in parent.items())).encode()).hexdigest()}
  if evidence["PDC000198"]["mapping_sha256"]!=evidence["PDC000199"]["mapping_sha256"]:raise RuntimeError("CROSS_STUDY_PARENT_MISMATCH")
 elif k=="ocular":
  b=fetch(url);rows=list(csv.DictReader(io.StringIO(b.decode("utf-8-sig")),delimiter="\t"))
  if len(rows)!=46:raise RuntimeError("RUN_COUNT_DRIFT")
  sample_to_experiments={}
  for r in rows:sample_to_experiments.setdefault(r["sample_accession"],set()).add(r["experiment_accession"])
  evidence={"runs":len(rows),"unique_samples":len(sample_to_experiments),"multi_experiment_samples":sum(len(x)>1 for x in sample_to_experiments.values()),"manifest_sha256":hashlib.sha256(b).hexdigest()}
 elif k=="muscle":
  b=fetch(url);obj=json.loads(b);evidence={"api_data_count":len(obj.get("data",[])),"api_sha256":hashlib.sha256(b).hexdigest(),"scope":"single API contract only"}
 elif k in ("glaucoma","oral"):
  b=fetch(url);s=b.decode("utf-8","replace");accession="GSE313736" if k=="glaucoma" else "GSE103322"
  if accession not in s:raise RuntimeError("ACCESSION_MISMATCH")
  evidence={"accession":accession,"series_samples":len(re.findall(r"^!Series_sample_id",s,re.M)),"supplementary_links":len(re.findall(r"^!Series_supplementary_file",s,re.M)),"soft_sha256":hashlib.sha256(b).hexdigest()}
 else:
  b=fetch(url);evidence={"source_sha256":hashlib.sha256(b).hexdigest(),"bytes":len(b)}
  if k=="is":
   o=json.loads(b)
   if "GCST90104545" not in json.dumps(o):raise RuntimeError("ACCESSION_MISMATCH")
   evidence["study"]="GCST90104545"
  if k=="radiomics":
   o=json.loads(b)
   if str(o.get("id"))!="8179129":raise RuntimeError("ZENODO_ID_MISMATCH")
   evidence["listed_files"]=len(o.get("files",[]))
  if k=="skin":
   s=b.decode("utf-8","replace").upper()
   evidence["required_columns_present"]=all(x in s for x in ("ALLELE1","ALLELE0","LOG_OR","LOG_SE"))
   if not evidence["required_columns_present"]:raise RuntimeError("ALLELE_SCHEMA_MISSING")
 return evidence
def main():
 p=argparse.ArgumentParser();p.add_argument("--track",required=True,choices=sorted(URLS));a=p.parse_args()
 r={"date_kst":"2026-11-04","track":a.track,"selected_task":"Independent fresh official-source consistency check","scientific_gate":"HOLD","next_gate":NEXT[a.track],"execution":"GITHUB_ACTIONS_ONLY"}
 try:r["evidence"]=task(a.track);r["source_check"]="PASS"
 except Exception as e:r["source_check"]="HOLD";r["blocker"]=type(e).__name__+":"+str(e)[:180]
 OUT.mkdir(parents=True,exist_ok=True);(OUT/(a.track+".json")).write_text(json.dumps(r,ensure_ascii=False,indent=2)+"\n")
 print(json.dumps(r,ensure_ascii=False))
if __name__=="__main__":main()
