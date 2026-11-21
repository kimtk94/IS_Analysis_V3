#!/usr/bin/env python3
"""Conservative source-access audit for independent ischemic stroke outcomes.
No inferred coverage: report actual queried source and explicit access issues.
"""
import argparse, json, urllib.request, urllib.error, re
from pathlib import Path

TARGET="rs187759247"
CANDIDATES=[
 {"id":"bbj-a-129","url":"https://gwas.mrcieu.ac.uk/datasets/bbj-a-129/",
  "ancestry":"EAS","build":"GRCh37","cases":17671,"controls":192383,
  "source":"IEU OpenGWAS BBJ ischemic stroke"},
 {"id":"FinnGen_R12_I9_STR_EXH","url":"https://storage.googleapis.com/finngen-public-data-r12/summary_stats/finngen_R12_manifest.tsv",
  "ancestry":"FIN","build":"GRCh38","source":"FinnGen R12 public summary manifest"}
]
def fetch(c):
    req=urllib.request.Request(c["url"],headers={"User-Agent":"Stroke-GWAS-access-audit/1.0"})
    try:
        with urllib.request.urlopen(req,timeout=35) as r:
            raw=r.read(500001)
            body=raw.decode("utf-8",errors="replace")
            if c["id"]=="bbj-a-129":
                # Search public metadata only. Not the GWAS full summary file.
                positive=all(s in body.lower() for s in ("ischemic stroke","17,671","192,383"))
                return {"access":"PASS_HTTP","http":r.status,"metadata_text_consistent":positive,
                        "snp_coverage":"NOT_QUERIED","observed_bytes":len(raw)}
            marker="I9_STR_EXH"
            rows=[ln[:500] for ln in body.splitlines() if marker in ln]
            return {"access":"PASS_HTTP","http":r.status,"endpoint_manifest_rows_in_first_500kb":rows[:5],
                    "manifest_truncated":len(raw)>=500001,"snp_coverage":"NOT_QUERIED",
                    "observed_bytes":len(raw)}
    except Exception as e:
        return {"access":"HOLD_HTTP","error_type":type(e).__name__,"error":str(e)[:220],
                "snp_coverage":"NOT_QUERIED"}
def run():
    ap=argparse.ArgumentParser();ap.add_argument("--out",required=True,type=Path);args=ap.parse_args()
    rows=[{**c,"verification":fetch(c)} for c in CANDIDATES]
    data={"schema":"tpmi_43321_alternative_stroke_gwas_access_v1","target":TARGET,
          "tested_GRCh37":"19:15573379:T:C","candidate_datasets":rows,
          "scientific_status":"HOLD","independent_replication":False,
          "next_gate":"Verify rsID or genomic coordinate in authentic per-variant outcome association and dataset overlap before effect alignment."}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(data,ensure_ascii=False,indent=2)+"\n")
    print(json.dumps({"candidates":[{"id":x["id"],**x["verification"]} for x in rows],
                      "scientific_status":"HOLD"}))
if __name__=="__main__":run()
