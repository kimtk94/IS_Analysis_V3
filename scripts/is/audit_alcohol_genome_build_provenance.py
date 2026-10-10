#!/usr/bin/env python3
"""Read-only provenance inventory. Never infers genome build from coordinates alone."""
import argparse, hashlib, json, re
from pathlib import Path

BUILD = re.compile(r"(?i)(?<![a-z0-9])(grch\s*3[78]|hg\s*1[89]|b3[78])(?![a-z0-9])")
TEXT_EXT={".json",".yaml",".yml",".md",".txt",".tsv",".csv",".log",".sh",".py",".r",".rds"}
MAX_BYTES=262144

def inventory(p):
    d={"path":str(p),"exists":p.is_file()}
    if not p.is_file(): return d
    d["size_bytes"]=p.stat().st_size
    if p.suffix.lower()==".rds":
        d["note"]="binary RDS not scanned; inspect source scripts and fit metadata"
        return d
    with p.open("rb") as f: raw=f.read(MAX_BYTES)
    text=raw.decode("utf-8",errors="replace")
    d["sample_sha256"]=hashlib.sha256(raw).hexdigest()
    d["sample_bytes"]=len(raw)
    hits=[]
    for i,line in enumerate(text.splitlines(),1):
        if BUILD.search(line):
            hits.append({"line":i,"text":line[:300]})
            if len(hits)>=20: break
    d["build_mentions"]=hits
    d["preview"]=text.splitlines()[:3] if p.suffix.lower() in TEXT_EXT else []
    return d

def main():
    a=argparse.ArgumentParser()
    a.add_argument("base",type=Path)
    a.add_argument("--repo",type=Path,default=Path.cwd())
    a.add_argument("--out",type=Path)
    args=a.parse_args()
    base=args.base.resolve()
    candidates=[
        base/"alcohol_conditional_susie_v1"/gene/"variants.tsv"
        for gene in ("ADH1B","ALDH2")
    ]
    candidates += [
        base/"alcohol_conditional_susie_v1"/gene/"ld_qc.json"
        for gene in ("ADH1B","ALDH2")
    ]
    repo=args.repo.resolve()
    for rel in (
        "scripts/is/run_alcohol_finite_ld_susie_sandbox.R",
        "scripts/is/export_alcohol_cs_variants.R",
        "scripts/is/prepare_alcohol_cs_annotation_manifest.py",
    ): candidates.append(repo/rel)
    report={"status":"PROVENANCE_INVENTORY_ONLY","build_resolution":"UNRESOLVED",
      "warnings":["Build mentions are clues, not verified provenance",
                  "No reference FASTA allele comparison performed",
                  "No rsID or VEP annotation performed"],
      "files":[inventory(p) for p in candidates]}
    out=args.out or base/"susie_rss_finite_ref_sandbox_v1"/"ALCOHOL_GENOME_BUILD_PROVENANCE_READONLY.json"
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,ensure_ascii=False,indent=2)+"\n")
    print("BUILD_PROVENANCE_INVENTORY_OK",out)
    for f in report["files"]:
        print("FILE",f["path"],"EXISTS",f["exists"],"MENTIONS",len(f.get("build_mentions",[])))
if __name__=="__main__":main()
