#!/usr/bin/env python3
"""Build a reproducible dataset registry for MASTER stages.

Manifest TSV (one row per logical dataset):
dataset_id,role,path,ancestry,genome_build,phenotype,platform,source_name
Optional:
sample_n,cases,controls,doi,url,notes

The registry validates local file existence, captures file size/SHA256, samples
the header without loading the full dataset, and records schema columns.
"""
from __future__ import annotations
import argparse,csv,gzip,bz2,lzma,hashlib,json
from pathlib import Path

def open_text(path:Path):
    n=path.name.lower()
    if n.endswith(".gz"): return gzip.open(path,"rt",encoding="utf-8",errors="replace",newline="")
    if n.endswith(".bz2"): return bz2.open(path,"rt",encoding="utf-8",errors="replace",newline="")
    if n.endswith(".xz"): return lzma.open(path,"rt",encoding="utf-8",errors="replace",newline="")
    return path.open("r",encoding="utf-8",errors="replace",newline="")

def sha256(path:Path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""):h.update(b)
    return h.hexdigest()

def read_manifest(path):
    with Path(path).open("r",encoding="utf-8",newline="") as f:
        yield from csv.DictReader(f,delimiter="	")

def header(path:Path):
    with open_text(path) as f:
        line=f.readline().rstrip("
")
    if not line:return []
    delim="	" if "	" in line else "," if "," in line else None
    return line.split(delim) if delim else line.split()

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--manifest",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    ap.add_argument("--summary",type=Path,required=True)
    ap.add_argument("--skip-hash",action="store_true")
    a=ap.parse_args()

    required={"dataset_id","role","path","ancestry","genome_build","phenotype","platform","source_name"}
    rows=[]
    seen=set()
    for m in read_manifest(a.manifest):
        miss=required-set(m)
        if miss:raise SystemExit("manifest missing columns: "+",".join(sorted(miss)))
        did=str(m["dataset_id"]).strip()
        if not did:raise SystemExit("empty dataset_id")
        if did in seen:raise SystemExit(f"duplicate dataset_id: {did}")
        seen.add(did)
        p=Path(m["path"]).expanduser()
        exists=p.is_file()
        cols=header(p) if exists else []
        row=dict(m)
        row.update({
          "exists":int(exists),
          "size_bytes":p.stat().st_size if exists else "",
          "sha256":"" if (a.skip_hash or not exists) else sha256(p),
          "column_count":len(cols),
          "columns":";".join(cols),
        })
        rows.append(row)

    fields=list(rows[0]) if rows else list(required)+["exists","size_bytes","sha256","column_count","columns"]
    a.output.parent.mkdir(parents=True,exist_ok=True)
    with a.output.open("w",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter="	",lineterminator="
",extrasaction="ignore")
        w.writeheader();w.writerows(rows)

    summary={
      "datasets":len(rows),
      "missing_files":[r["dataset_id"] for r in rows if r["exists"]==0],
      "roles":sorted({str(r.get("role","")) for r in rows}),
      "ancestries":sorted({str(r.get("ancestry","")) for r in rows}),
      "genome_builds":sorted({str(r.get("genome_build","")) for r in rows}),
      "platforms":sorted({str(r.get("platform","")) for r in rows}),
    }
    a.summary.parent.mkdir(parents=True,exist_ok=True)
    a.summary.write_text(json.dumps(summary,indent=2)+"
",encoding="utf-8")
    if summary["missing_files"]:
        raise SystemExit("dataset registry contains missing files: "+",".join(summary["missing_files"]))
    print(json.dumps(summary,indent=2))
    print("MASTER_DATASET_REGISTRY_PASS")

if __name__=="__main__":main()
