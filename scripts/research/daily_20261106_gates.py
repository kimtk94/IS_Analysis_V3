#!/usr/bin/env python3
"""Bounded official-source gates. Execute only in GitHub Actions; fail closed."""
import argparse,csv,datetime,hashlib,io,json,os,pathlib,re,urllib.request,urllib.parse
TRACKS=("is","glaucoma","hcc","radiomics","ocular","muscle","oral","skin")
TASKS={
"is":"Confirm GIGASTROKE v2 official accession and ancestry endpoint metadata",
"glaucoma":"Inspect GSE313736 official GEO sample and supplementary-file inventory",
"hcc":"Reconcile PDC000198 and PDC000199 official aliquot-parent mappings",
"radiomics":"Audit official Zenodo imaging file checksums and sizes",
"ocular":"Reconcile ENA run-level FASTQ/MD5/size triplets and timepoint coverage",
"muscle":"Check GTEx v10 official variant API availability and documented unresolved input",
"oral":"Audit GSE103322 official GEO sample links and supplementary inventory",
"skin":"Check Bristol perceived-age readme effect/allele schema and genome-build declaration"}
def get(url,limit=3000000):
 req=urllib.request.Request(url,headers={"User-Agent":"Research-CI-evidence-audit/2026-11-06"})
 with urllib.request.urlopen(req,timeout=25) as r:
  b=r.read(limit+1)
  if len(b)>limit: raise RuntimeError("SOURCE_EXCEEDS_BOUNDED_LIMIT")
  return b
def geo(acc):
 url="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?"+urllib.parse.urlencode({"acc":acc,"targ":"self","form":"text","view":"quick"})
 b=get(url);s=b.decode("utf8","replace")
 if acc not in s:raise ValueError("GEO_ACCESSION_MISSING")
 return {"accession":acc,"source_sha256":hashlib.sha256(b).hexdigest(),"sample_links":len(re.findall(r"^!Series_sample_id",s,re.M)),"supplementary_links":len(re.findall(r"^!Series_supplementary_file",s,re.M))}
def probe(t):
 if t=="is":
  b=get("https://www.ebi.ac.uk/gwas/rest/api/v2/studies/GCST90104540")
  if b"GCST90104540" not in b:raise ValueError("ACCESSION_MISMATCH")
  return {"accession":"GCST90104540","metadata_sha256":hashlib.sha256(b).hexdigest(),"full_sumstats_verified":False}
 if t=="glaucoma":return geo("GSE313736")
 if t=="oral":return geo("GSE103322")
 if t=="hcc":
  results=[]
  for study in ("PDC000198","PDC000199"):
   q='{ biospecimenPerStudy(pdc_study_id:"'+study+'",acceptDUA:true){aliquot_id sample_id case_id} }'
   req=urllib.request.Request("https://pdc.cancer.gov/graphql",data=json.dumps({"query":q}).encode(),headers={"Content-Type":"application/json","User-Agent":"Research-CI/1.0"})
   with urllib.request.urlopen(req,timeout=35) as resp:b=resp.read(5000001)
   if len(b)>5000000:raise RuntimeError("PDC_SIZE_LIMIT")
   obj=json.loads(b)
   if obj.get("errors"):raise RuntimeError("PDC_GRAPHQL_ERROR")
   rows=obj["data"]["biospecimenPerStudy"]
   triples=[(r["aliquot_id"],r["sample_id"],r["case_id"]) for r in rows]
   if len(triples)!=len(set(x[0] for x in triples)):raise ValueError("DUPLICATE_ALIQUOT")
   results.append((study,set(triples)))
  if results[0][1]!=results[1][1]:raise ValueError("CROSS_STUDY_PARENT_MISMATCH")
  return {"studies":[x[0] for x in results],"unique_aliquot_parents":len(results[0][1]),"quant_column_mapping_verified":False}
 if t=="radiomics":
  o=json.loads(get("https://zenodo.org/api/records/8179129"))
  if str(o.get("id"))!="8179129":raise ValueError("ZENODO_ID_MISMATCH")
  files=o.get("files",[])
  return {"record_id":8179129,"files":len(files),"bytes_declared":sum(int(f.get("size") or 0) for f in files),"files_with_checksum":sum(bool(f.get("checksum")) for f in files),"nifti_headers_verified":False}
 if t=="ocular":
  url="https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA1243143&result=read_run&format=tsv&fields=run_accession,sample_accession,fastq_ftp,fastq_md5,fastq_bytes,sample_title,experiment_title"
  rows=list(csv.DictReader(io.StringIO(get(url).decode("utf8-sig")),delimiter="\t"))
  if len(rows)!=46 or len(set(r["run_accession"] for r in rows))!=46:raise ValueError("ENA_RUN_COUNT_OR_DUPLICATE")
  triplets=0;malformed=0
  for r in rows:
   parts=[[x for x in (r.get(k) or "").split(";") if x] for k in ("fastq_ftp","fastq_md5","fastq_bytes")]
   if not all(len(x)==2 for x in parts):malformed+=1
   else:triplets+=1
  return {"runs":len(rows),"paired_fastq_triplets":triplets,"malformed_triplets":malformed,"raw_fastq_md5_verified":False,"timepoint_verified":False}
 if t=="muscle":
  # No substitute loci: API health check is not a 20-locus association audit.
  o=json.loads(get("https://gtexportal.org/api/v2/dataset/variant?snpId=rs429358&datasetId=gtex_v10&itemsPerPage=10"))
  return {"gtex_v10_api_json":isinstance(o,dict),"probe_rsid":"rs429358","target_20_loci_rerun":False,"unresolved_two_loci":True}
 if t=="skin":
  b=get("https://data.bris.ac.uk/datasets/21crwsnj4xwjm2g4qi8chathha/readme.txt",500000)
  s=b.decode("utf8","replace").upper()
  fields=("ALLELE1","ALLELE0","A1FREQ","LOG_OR","LOG_SE")
  if not all(x in s for x in fields):raise ValueError("ALLELE_EFFECT_SCHEMA_MISSING")
  return {"readme_sha256":hashlib.sha256(b).hexdigest(),"required_fields":list(fields),"genome_build_explicit":any(x in s for x in ("GRCH37","GRCH38","HG19","HG38")),"full_gwas_qc":False}
def main():
 p=argparse.ArgumentParser();p.add_argument("--track",required=True,choices=TRACKS);a=p.parse_args()
 if os.getenv("GITHUB_ACTIONS")!="true":raise SystemExit("GITHUB_ACTIONS_ONLY")
 d={"date_kst":"2026-11-06","track":a.track,"selected_task":TASKS[a.track],"scientific_gate":"HOLD","ckd":"HOLD_NO_EXECUTION","status":"HOLD","evidence":{},"blocker":None}
 try:
  d["evidence"]=probe(a.track);d["status"]="PASS_BOUNDED_SOURCE_GATE"
 except Exception as e:
  d["status"]="HOLD_SOURCE_OR_INTEGRITY_REVIEW";d["blocker"]=type(e).__name__+":"+str(e)[:200]
 out=pathlib.Path("results/research_operations/20261106");out.mkdir(parents=True,exist_ok=True)
 (out/(a.track+".json")).write_text(json.dumps(d,ensure_ascii=False,indent=2)+"\n")
 print(json.dumps({"track":a.track,"status":d["status"],"scientific_gate":d["scientific_gate"],"blocker":d["blocker"]}))
if __name__=="__main__":main()
