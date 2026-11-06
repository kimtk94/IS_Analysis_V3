#!/usr/bin/env python3
"""One bounded, fail-closed scientific prerequisite per track. GitHub Actions only."""
import argparse,csv,datetime,hashlib,io,json,os,pathlib,re,urllib.parse,urllib.request
TRACKS=("is","glaucoma","hcc","radiomics","ocular","muscle","oral","skin")
TASKS={
"is":"Check GWAS ancestry metadata for explicit EUR/EAS population provenance",
"glaucoma":"Check GSE313736 GEO supplementary count-file declarations",
"hcc":"Check official PDC case-level tumor-normal paired sample structure",
"radiomics":"Check Zenodo source file checksum algorithm and size metadata",
"ocular":"Check ENA 46-run unique sample-experiment-run parent consistency",
"muscle":"Check GTEx v10 rsID lookup for both sentinel variants and distinguish absent versus API failure",
"oral":"Check GSE103322 sample IDs uniqueness and declared sample count",
"skin":"Check Bristol GWAS readme effect scale and explicit genome build"}
def get(url,limit=3000000):
 req=urllib.request.Request(url,headers={"User-Agent":"MasterDegree-CI-gates/2026-11-07","Accept":"*/*"})
 with urllib.request.urlopen(req,timeout=30) as r:
  b=r.read(limit+1)
  if len(b)>limit:raise RuntimeError("BOUNDED_DOWNLOAD_LIMIT")
  return b
def geo(acc):
 u="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?"+urllib.parse.urlencode({"acc":acc,"targ":"self","form":"text","view":"quick"})
 b=get(u);s=b.decode("utf8","replace")
 if not re.search(r"\^SERIES\s*=\s*"+acc+r"\b",s,re.I):raise ValueError("GEO_ACCESSION_NOT_VERIFIED")
 return s,hashlib.sha256(b).hexdigest()
def probe(t):
 if t=="is":
  u="https://www.ebi.ac.uk/gwas/rest/api/v2/studies/GCST90104540/ancestries"
  b=get(u);o=json.loads(b)
  if not isinstance(o,dict):raise ValueError("INVALID_ANCESTRY_JSON")
  terms=sorted(set(re.findall(r"European|East Asian|EUR|EAS",json.dumps(o),re.I)))
  return {"ancestry_endpoint_sha256":hashlib.sha256(b).hexdigest(),"explicit_ancestry_terms":terms,"snp_level_ancestry_ld_verified":False}
 if t=="glaucoma":
  s,h=geo("GSE313736");links=re.findall(r"^!Series_supplementary_file\s*=\s*(.+)$",s,re.M)
  return {"geo_sha256":h,"supplementary_files":len(links),"count_like_file_names":sum(bool(re.search("count|matrix|mtx",x,re.I)) for x in links),"integer_count_matrix_verified":False}
 if t=="hcc":
  q='{ biospecimenPerStudy(pdc_study_id:"PDC000198",acceptDUA:true){aliquot_id sample_id case_id sample_type} }'
  req=urllib.request.Request("https://pdc.cancer.gov/graphql",data=json.dumps({"query":q}).encode(),headers={"Content-Type":"application/json","User-Agent":"MasterDegree-CI/1.0"})
  with urllib.request.urlopen(req,timeout=40) as r:b=r.read(5000001)
  if len(b)>5000000:raise RuntimeError("PDC_RESPONSE_LIMIT")
  o=json.loads(b)
  if o.get("errors"):raise RuntimeError("PDC_GRAPHQL_ERRORS")
  rows=o["data"]["biospecimenPerStudy"];cases={}
  for r in rows:
   c=r.get("case_id");a=r.get("aliquot_id");s=r.get("sample_id")
   if not all((c,a,s)):raise ValueError("MISSING_PARENT_ID")
   cases.setdefault(c,[]).append(str(r.get("sample_type") or ""))
  return {"case_count":len(cases),"aliquot_rows":len(rows),"sample_type_categories":{x:sum(x in v for v in cases.values()) for x in sorted({t for v in cases.values() for t in v})},"quant_column_mapping_verified":False,"rna_patient_join_verified":False}
 if t=="radiomics":
  o=json.loads(get("https://zenodo.org/api/records/8179129"))
  if str(o.get("id"))!="8179129":raise ValueError("RECORD_MISMATCH")
  f=o.get("files",[])
  if not f:raise ValueError("NO_FILES")
  return {"files":len(f),"size_positive":sum(int(x.get("size") or 0)>0 for x in f),"checksum_declared":sum(bool(x.get("checksum")) for x in f),"checksum_algorithms":sorted(set(str(x.get("checksum") or "").split(":")[0] for x in f if x.get("checksum"))),"actual_nifti_header_verified":False}
 if t=="ocular":
  u="https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA1243143&result=read_run&format=tsv&fields=run_accession,sample_accession,experiment_accession"
  rows=list(csv.DictReader(io.StringIO(get(u).decode("utf8-sig")),delimiter="\t"))
  if len(rows)!=46:raise ValueError("RUN_DENOMINATOR_DRIFT")
  runs={r["run_accession"] for r in rows};experiments={r["experiment_accession"] for r in rows}
  if len(runs)!=46:raise ValueError("DUPLICATE_RUN")
  if any(not r["sample_accession"] or not r["experiment_accession"] for r in rows):raise ValueError("MISSING_PARENT")
  return {"runs":len(runs),"experiments":len(experiments),"biosamples":len({r["sample_accession"] for r in rows}),"timepoint_donor_mapping_verified":False}
 if t=="muscle":
  probes={}
  for rsid in ("rs429358","rs7412"):
   u="https://gtexportal.org/api/v2/dataset/variant?"+urllib.parse.urlencode({"snpId":rsid,"datasetId":"gtex_v10","itemsPerPage":10})
   o=json.loads(get(u,1000000));v=o.get("data")
   if not isinstance(v,list):raise ValueError("INVALID_GTEX_RESPONSE")
   probes[rsid]=len(v)
  return {"api_sentinels":probes,"target_20_loci_association_verified":False,"note":"Sentinels are API checks, not the target 20 loci"}
 if t=="oral":
  s,h=geo("GSE103322");ids=re.findall(r"^!Series_sample_id\s*=\s*(GSM\d+)",s,re.M)
  if not ids:raise ValueError("NO_GEO_SAMPLE_IDS")
  if len(ids)!=len(set(ids)):raise ValueError("DUPLICATE_GEO_SAMPLE_IDS")
  return {"geo_sha256":h,"series_sample_links":len(ids),"unique_sample_links":len(set(ids)),"patient_to_cell_join_verified":False}
 if t=="skin":
  b=get("https://data.bris.ac.uk/datasets/21crwsnj4xwjm2g4qi8chathha/readme.txt",500000)
  s=b.decode("utf8","replace").upper()
  cols=("ALLELE1","ALLELE0","A1FREQ","LOG_OR","LOG_SE")
  missing=[c for c in cols if c not in s]
  if missing:raise ValueError("MISSING_SCHEMA_"+",".join(missing))
  return {"readme_sha256":hashlib.sha256(b).hexdigest(),"effect_scale":"LOG_OR","explicit_genome_build":sorted(set(re.findall("GRCH37|GRCH38|HG19|HG38",s))),"full_sumstats_qc_verified":False}
def main():
 p=argparse.ArgumentParser();p.add_argument("--track",choices=TRACKS,required=True);a=p.parse_args()
 if os.getenv("GITHUB_ACTIONS")!="true":raise SystemExit("CI_ONLY")
 out=pathlib.Path("results/research_operations/20261107");out.mkdir(parents=True,exist_ok=True)
 result={"track":a.track,"task":TASKS[a.track],"date_kst":"2026-11-07","scientific_gate":"HOLD","checked_utc":datetime.datetime.now(datetime.timezone.utc).isoformat()}
 try:result.update(probe_status="PASS_BOUNDED_PREREQUISITE",evidence=probe(a.track))
 except Exception as e:result.update(probe_status="HOLD_OR_FAIL_SOURCE",blocker=type(e).__name__+": "+str(e)[:180],evidence={})
 (out/(a.track+".json")).write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n")
 print(json.dumps(result,ensure_ascii=False))
if __name__=="__main__":main()
