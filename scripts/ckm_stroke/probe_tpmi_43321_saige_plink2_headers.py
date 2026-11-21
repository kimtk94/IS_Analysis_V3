#!/usr/bin/env python3
"""Inspect authentic TPMI PheWeb PLINK2 vs SAIGE first-line schemas, not whole data.

The phenotype page includes summary statistics for both methods; reading the
first compressed MB does NOT substantiate per-variant case/control counts.
Require exactly each official PheWeb download path before probing.
No summary-stat rows exported or score files touched.
"""
import argparse,json,urllib.parse,urllib.request,zlib
from pathlib import Path

BASE="https://pheweb.ibms.sinica.edu.tw"
ENDPOINTS={"PLINK2":"/download/gwas_plink2/433.21",
           "SAIGE":"/download/saige_gwas/433.21"}
MAX_BYTES=4_194_304
def inspect(name,path):
 url=BASE+path
 req=urllib.request.Request(url,headers={"User-Agent":"IS-Analysis-V3-TPMI-43321-method-header/1.0",
                                         "Range":"bytes=0-1048575"})
 out={"provider":name,"official_link":url}
 try:
  with urllib.request.urlopen(req,timeout=75) as r:
   resolved=urllib.parse.urlsplit(r.geturl())
   if resolved.scheme!="https" or not resolved.hostname:
    raise ValueError("unapproved redirect")
   out.update({"http":r.status,"resolved_host":resolved.hostname,
               "content_type":r.headers.get("Content-Type"),
               "content_length":r.headers.get("Content-Length"),
               "content_range":r.headers.get("Content-Range"),
               "content_disposition":r.headers.get("Content-Disposition")})
   if r.status==206:
    data=r.read(1_048_578)
    if len(data)>1_048_576:raise ValueError("partial response over range size")
   else:
    data=r.read(MAX_BYTES+1)
    if len(data)>MAX_BYTES:raise ValueError("server did not honor range within bound")
  if data[:2]!=b"\x1f\x8b":
   out.update({"status":"HOLD_NOT_GZIP","first_32_hex":data[:32].hex()})
   return out
  unpack=zlib.decompressobj(zlib.MAX_WBITS|16)
  head=unpack.decompress(data,1_048_576)
  lines=head.decode("utf-8-sig",errors="replace").splitlines()[:3]
  columns=lines[0].split("\t") if lines else []
  out.update({"status":"PASS_ORIGINAL_HEADER_PROBED",
              "header_columns":columns[:60],"header_width":len(columns),
              "line_length_bounds":[len(s) for s in lines],
              "SAIGE_original_case_control_per_variant_in_schema":bool(
                  any("case" in c.lower() or "control" in c.lower() for c in columns)),
              "gzip_full_content_verified":False,
              "only_header_examined":True})
 except Exception as e:out.update({"status":"HOLD_DOWNLOAD_OR_SCHEMA",
                                   "error_type":type(e).__name__,"error":str(e)[:270]})
 return out
def main():
 p=argparse.ArgumentParser();p.add_argument("--out",required=True,type=Path);a=p.parse_args()
 result={"schema":"TPMI43321-SAIGE-vs-PLINK2-original-header-probe-v1",
         "methods":{k:inspect(k,v) for k,v in ENDPOINTS.items()},
         "phenotype_reported_cases":9249,"phenotype_reported_controls":304660,
         "pinned_PLINK2_OBS_CT":204099,
         "regression_specific_case_control_breakdown_confirmed":False,
         "TPMI_binary_Neff_approved":False,
         "scientific_status":"HOLD"}
 a.out.parent.mkdir(parents=True,exist_ok=True)
 a.out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n")
 print(json.dumps(result,ensure_ascii=False))
if __name__=="__main__":main()
