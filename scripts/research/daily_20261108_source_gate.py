#!/usr/bin/env python3
"""Eight bounded source gates; no biological PASS is inferred."""
import argparse, datetime, hashlib, json, pathlib, urllib.request, urllib.error, re
SOURCES={
"is":("GIGASTROKE GWAS accession response","https://www.ebi.ac.uk/gwas/rest/api/studies/GCST90104540"),
"glaucoma":("GSE313736 GEO record","https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE313736"),
"hcc":("PDC000198 official study record","https://pdc.cancer.gov/pdc/study/PDC000198"),
"radiomics":("IBSI reference provenance","https://ibsi.readthedocs.io/en/latest/"),
"ocular":("PRJNA1243143 ENA run metadata","https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA1243143&result=read_run&format=tsv"),
"muscle":("GTEx v10 official API availability","https://gtexportal.org/api/v2/dataset/tissueSiteDetail?itemsPerPage=1"),
"oral":("GSE103322 GEO record","https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103322"),
"skin":("Bristol perceived-age GWAS readme","https://data.bris.ac.uk/datasets/21crwsnj4xwjm2g4qi8chathha/readme.txt")
}
p=argparse.ArgumentParser();p.add_argument("--track",choices=SOURCES,required=True);p.add_argument("--output",required=True);a=p.parse_args()
task,url=SOURCES[a.track];out={"track":a.track,"task":task,"url":url,"retrieved_utc":datetime.datetime.now(datetime.timezone.utc).isoformat(),"scientific_gate":"HOLD","scope":"bounded source response only","status":"HOLD"}
try:
 req=urllib.request.Request(url,headers={"User-Agent":"MasterDegreeResearchAudit/1.0","Range":"bytes=0-65535","Accept-Encoding":"identity"})
 with urllib.request.urlopen(req,timeout=35) as r:
  b=r.read(65536);out.update(http_status=r.status,bytes_read=len(b),sha256_prefix=hashlib.sha256(b).hexdigest(),content_type=r.headers.get("Content-Type"))
  if r.status==200 or r.status==206:
   out["status"]="PASS_SOURCE_REACHABILITY_ONLY"
  else:out["blocker"]="unexpected HTTP response"
except Exception as e:out["blocker"]=type(e).__name__+": "+str(e)[:200]
pathlib.Path(a.output).parent.mkdir(parents=True,exist_ok=True)
pathlib.Path(a.output).write_text(json.dumps(out,ensure_ascii=False,indent=2)+"\n")
print(json.dumps(out,ensure_ascii=False))
