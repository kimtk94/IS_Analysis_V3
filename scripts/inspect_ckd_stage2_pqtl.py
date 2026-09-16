#!/usr/bin/env python3
"""Inspect Stage 2 UKB-PPP archive member names and headers without full extraction."""
from __future__ import annotations
import argparse, csv, gzip, io, tarfile
from pathlib import Path

def read_head_from_member(archive, member):
    raw=archive.extractfile(member)
    if raw is None:
        return ""
    stream=raw
    if member.name.lower().endswith(".gz"):
        stream=gzip.GzipFile(fileobj=raw,mode="rb")
    line=stream.readline()
    return line.decode("utf-8",errors="replace").rstrip("\r\n")

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args=ap.parse_args()
    rows=[]
    for path in sorted(args.root.rglob("*.tar")):
        ancestry=path.relative_to(args.root).parts[0] if len(path.relative_to(args.root).parts)>1 else ""
        gene=path.parent.name
        with tarfile.open(path,"r:*") as tar:
            members=[m for m in tar if m.isfile()]
            if not members:
                raise SystemExit(f"{path}: no regular archive member")
            for m in members[:5]:
                rows.append({
                    "ancestry":ancestry, "gene_symbol":gene, "archive":str(path),
                    "member":m.name, "member_size":m.size, "header":read_head_from_member(tar,m),
                })
    if not rows:
        raise SystemExit(f"no .tar archives under {args.root}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fields=list(rows[0])
    with args.output.open("w",encoding="utf-8",newline="") as fh:
        w=csv.DictWriter(fh,fieldnames=fields,delimiter="\t",lineterminator="\n")
        w.writeheader(); w.writerows(rows)
    print(f"CKD_STAGE2_SCHEMA_PASS archives={len(set(r['archive'] for r in rows))} output={args.output}")

if __name__=="__main__":
    main()
