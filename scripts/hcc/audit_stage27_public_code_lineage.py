#!/usr/bin/env python3
"""Stage27 source-code lineage audit pinned to public Suter HCC Git commit.

Only aggregate script filename/hash/query counts are exported. A lack of
generation code in this repository is not proof it never existed elsewhere.
"""
import argparse
import datetime as dt
import hashlib
import json
import pathlib
import re
import urllib.request

REPO="cbg-ethz/HCC"
COMMIT="438f52271d0a86a6701baa7a96e80c516ba575c2"
BASE=f"https://raw.githubusercontent.com/{REPO}/{COMMIT}/"
TREE=f"https://api.github.com/repos/{REPO}/git/trees/{COMMIT}?recursive=1"
SCRIPTS=("HCC_prior.R","HCCclutsering.R","HCCfullclustering.R",
         "helpresfns.R","Tutorial_bnClustOmics.R","README.txt")
MAX_BYTES=500000
def get(url):
    req=urllib.request.Request(url,headers={
        "User-Agent":"HCC-Stage27-public-method-lineage/1.0",
        "Accept":"application/json,text/plain"})
    with urllib.request.urlopen(req,timeout=25) as res:
        content=res.read(MAX_BYTES+1)
    if not 0<len(content)<=MAX_BYTES:
        raise ValueError("source unexpectedly large")
    return content
def summarize(code):
    # Exact file-string searches are proof of only what's present in this code.
    return {
        "HCCfullz_filename_mentions":len(re.findall(r"HCCfullz(?:\.rds)?",code,re.I)),
        "HCCfull_filename_mentions":len(re.findall(r"HCCfull\.rds",code,re.I)),
        "readRDS_call_mentions":len(re.findall(r"\breadRDS\s*\(",code)),
        "saveRDS_call_mentions":len(re.findall(r"\bsaveRDS\s*\(",code)),
        "HCCfullz_generation_assignment_mentions":len(re.findall(
            r"HCCfullz\s*<-\s*(?:saveRDS|scale|sweep|apply|lapply|normalize)",
            code,re.I)),
        "HCCfullz_file_writer_mentions":len(re.findall(
            r"(?:saveRDS|write\.table|write\.csv)\s*\([^;\n]*HCCfullz",code,re.I)),
        "feature_variance_screen_mentions":len(re.findall(r"\bFSbyVar\s*\(",code)),
        "HCCfullz_consumed_in_clustering_script":bool(re.search(
            r'HCCfullz\s*<-\s*readRDS\("HCCinputs/HCCfullz\.rds"\)',code)),
    }
def audit():
    tree=json.loads(get(TREE))
    if tree.get("truncated") or not isinstance(tree.get("tree"),list):
        raise ValueError("incomplete source repository tree")
    paths={e.get("path") for e in tree["tree"] if e.get("type")=="blob"}
    if not set(SCRIPTS).issubset(paths):
        raise ValueError("source code manifest missing")
    scripts={}
    for name in SCRIPTS:
        raw=get(BASE+name)
        src=raw.decode("utf-8")
        scripts[name]={
            "sha256":hashlib.sha256(raw).hexdigest(),
            "size_bytes":len(raw),
            **summarize(src)
        }
    n_src_fullz=sum(x["HCCfullz_filename_mentions"] for x in scripts.values())
    n_writes=sum(x["HCCfullz_file_writer_mentions"] for x in scripts.values())
    return {
      "schema":"hcc-stage27-source-code-provenance-v1",
      "source_repository":REPO,"source_commit":COMMIT,
      "checked_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
      "repository_truncated":False,
      "repository_blob_count":len(paths),
      "tracked_R_script_count":sum(p.lower().endswith(".r") for p in paths),
      "audited_named_scripts":scripts,
      "HCCfullz_source_filename_mentions_total":n_src_fullz,
      "HCCfullz_source_file_writer_mentions_total":n_writes,
      "HCCfullz_clustering_script_load_verified":
          scripts["HCCfullclustering.R"]["HCCfullz_consumed_in_clustering_script"],
      "HCCfullz_generation_code_verified_in_source_repository":False,
      "HCCfullz_generation_code_absent_elsewhere_proven":False,
      "original_Ng_2022_Fig2b_gene_universe_recreated":False,
      "source_patient_gene_or_expression_exported":False,
      "interpretation":"This pinned repository includes prebuilt RDS artifacts and downstream analysis code. No source generator established from these six scripts alone. Do not claim nonexistence outside repo."
    }
def main():
    p=argparse.ArgumentParser();p.add_argument("--out",required=True)
    a=p.parse_args()
    o=audit()
    f=pathlib.Path(a.out);f.parent.mkdir(parents=True,exist_ok=True)
    f.write_text(json.dumps(o,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(o,indent=2))
if __name__=="__main__":main()
