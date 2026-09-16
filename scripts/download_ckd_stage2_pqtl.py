#!/usr/bin/env python3
"""Download only Stage 2 candidate UKB-PPP archives through the Synapse CLI."""
from __future__ import annotations
import argparse, csv, hashlib, json, shutil, subprocess
from pathlib import Path

def read_tsv(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))

def digest(path, kind):
    h = hashlib.md5(usedforsecurity=False) if kind == "md5" else hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def validate(row, path):
    size = row.get("expected_size_bytes", "").strip()
    if size and path.stat().st_size != int(size):
        return f"size mismatch: {path.stat().st_size} != {size}"
    for kind in ("md5", "sha256"):
        expected = row.get(kind, "").strip().lower()
        if expected and digest(path, kind) != expected:
            return f"{kind} mismatch"
    return None

def write_progress(path, rows):
    fields = ["gene_symbol","protein_id","ancestry","synapse_id","source_file","status","bytes","path"]
    with path.open("w", encoding="utf-8", newline="") as fh:
        w=csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(rows)

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--manifest", type=Path, required=True)
    ap.add_argument("--dest", type=Path, required=True)
    ap.add_argument("--ancestry", choices=["EUR","EAS","BOTH"], default="EUR")
    ap.add_argument("--progress", type=Path)
    args=ap.parse_args()

    if not shutil.which("synapse"):
        raise SystemExit("Synapse CLI not found. Install synapseclient in the CKD venv and authenticate with: synapse login")
    selected=[]
    for r in read_tsv(args.manifest):
        anc=r.get("ancestry","").upper()
        if args.ancestry != "BOTH" and anc != args.ancestry:
            continue
        if args.ancestry == "EUR" and str(r.get("download_now","")) != "1":
            continue
        selected.append(r)
    if not selected:
        raise SystemExit("no manifest rows selected")

    progress=[]
    for r in selected:
        gene, anc = r["gene_symbol"], r["ancestry"].upper()
        target_dir=args.dest/anc/gene
        target_dir.mkdir(parents=True, exist_ok=True)
        target=target_dir/r["source_file"]
        status="reused"
        if not target.is_file() or validate(r,target):
            if target.exists():
                target.unlink()
            subprocess.run(["synapse","get",r["synapse_id"],"--downloadLocation",str(target_dir)], check=True)
            downloaded=target_dir/r["source_file"]
            if not downloaded.is_file():
                raise SystemExit(f"Synapse download missing expected file: {downloaded}")
            target=downloaded
            status="downloaded"
        problem=validate(r,target)
        if problem:
            raise SystemExit(f"{gene} {anc}: {problem}")
        progress.append({
            "gene_symbol":gene, "protein_id":r["protein_id"], "ancestry":anc,
            "synapse_id":r["synapse_id"], "source_file":r["source_file"],
            "status":status, "bytes":target.stat().st_size, "path":str(target),
        })
        print(f"[{status}] {gene} {anc} {target.stat().st_size} bytes")

    p=args.progress or (args.dest/"stage2_pqtl_download_progress.tsv")
    p.parent.mkdir(parents=True, exist_ok=True)
    write_progress(p,progress)
    print(json.dumps({"downloads":len(progress),"ancestry":args.ancestry,"progress":str(p)},indent=2))
    print("CKD_STAGE2_PQTL_DOWNLOAD_PASS")

if __name__=="__main__":
    main()
