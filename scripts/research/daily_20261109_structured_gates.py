#!/usr/bin/env python3
"""Bounded, fail-closed official metadata checks; never biological validation."""
import argparse,datetime,hashlib,json,pathlib,re,urllib.request,urllib.error,xml.etree.ElementTree as ET
P={
"is":("GIGASTROKE accession structured GWAS metadata","https://www.ebi.ac.uk/gwas/rest/api/studies/GCST90104540","GCST90104540"),
"glaucoma":("GSE313736 accession and series record","https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE313736&targ=self&form=text&view=quick","GSE313736"),
"hcc":("PDC000198 study identity","https://pdc.cancer.gov/pdc/study/PDC000198","PDC000198"),
"radiomics":("IBSI reference document availability","https://ibsi.readthedocs.io/en/latest/","IBSI"),
"ocular":("ENA PRJNA1243143 distinct run and sample counts","https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA1243143&result=read_run&fields=run_accession,experiment_accession,sample_accession&format=tsv","run_accession"),
"muscle":("GTEx v10 Muscle_Skeletal tissue identity","https://gtexportal.org/api/v2/dataset/tissueSiteDetail?itemsPerPage=100","Muscle_Skeletal"),
"oral":("GSE103322 accession and series record","https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103322&targ=self&form=text&view=quick","GSE103322"),
"skin":("Bristol perceived-age allele/effect field declarations","https://data.bris.ac.uk/datasets/21crwsnj4xwjm2g4qi8chathha/readme.txt","ALLELE1")
}
a=argparse.ArgumentParser();a.add_argument("--track",choices=P,required=True);a.add_argument("--output",required=True);z=a.parse_args()
task,url,token=P[z.track]
o={"date_kst":"2026-11-09","track":z.track,"task":task,"url":url,"scientific_gate":"HOLD","status":"HOLD","checked_utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),"scope":"official bounded metadata; not full raw-data integrity"}
try:
 req=urllib.request.Request(url,headers={"User-Agent":"ResearchMetadataAudit/1.0","Accept-Encoding":"identity"})
 with urllib.request.urlopen(req,timeout=45) as r:
  raw=r.read(1024*1024+1);o["http_status"]=r.status;o["bytes_read"]=len(raw);o["sha256"]=hashlib.sha256(raw).hexdigest()
  if len(raw)>1024*1024:raise ValueError("RESPONSE_EXCEEDS_1MIB_BOUND")
  s=raw.decode("utf-8","replace")
  if z.track=="ocular":
   lines=[x.split("\t") for x in s.strip().splitlines()]
   header=lines[0];rows=[dict(zip(header,x)) for x in lines[1:] if len(x)==len(header)]
   o["run_count"]=len(set(x["run_accession"] for x in rows));o["experiment_count"]=len(set(x["experiment_accession"] for x in rows));o["sample_count"]=len(set(x["sample_accession"] for x in rows))
   ok=o["run_count"]==46 and o["experiment_count"]>0 and o["sample_count"]>0
  elif z.track=="skin":
   o["required_fields"]={k:k in s for k in ("ALLELE1","ALLELE0","A1FREQ","LOG_OR","LOG_SE")}
   ok=all(o["required_fields"].values())
  else:
   ok=token.lower() in s.lower()
  o["identity_assertion"]=bool(ok)
  o["status"]="PASS_METADATA_ONLY" if ok else "HOLD"
  if not ok:o["blocker"]="expected official identifier or metadata missing; do not infer biological PASS"
except Exception as e:o["blocker"]=type(e).__name__+": "+str(e)[:300]
p=pathlib.Path(z.output);p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(o,ensure_ascii=False,indent=2)+"\n")
print(json.dumps(o,ensure_ascii=False))
