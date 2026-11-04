#!/usr/bin/env python3
"""Validate the presence and dimensions of an externally staged UMOD CI bundle.

Never claim scientific verification on synthetic, missing, or incomplete inputs.
Bundle layout under CKD_UMOD_INPUT_DIR, default: input/ckd/umod:
  UMOD.tsv.gz, UMOD.ld.meta.tsv, UMOD.ld.vars, UMOD.ld.bin
"""
import csv,gzip,json,os,sys
from pathlib import Path
p=Path(os.getenv("CKD_UMOD_INPUT_DIR","input/ckd/umod"))
req=("UMOD.tsv.gz","UMOD.ld.meta.tsv","UMOD.ld.vars","UMOD.ld.bin")
missing=[x for x in req if not (p/x).is_file()]
if missing:
 print(json.dumps({"status":"BLOCKED_MISSING_PRIMARY_DATA","input_dir":str(p),"missing":missing,
 "analysis_executed":False,"server_ssh_required":False},indent=2))
 sys.exit(0)  # gate reports blocked rather than failing a repository CI
with (p/"UMOD.ld.meta.tsv").open() as f:
 meta=list(csv.DictReader(f,delimiter="\t"))
with (p/"UMOD.ld.vars").open() as f:
 vars=[s.strip() for s in f if s.strip()]
n=len(meta);size=(p/"UMOD.ld.bin").stat().st_size
with gzip.open(p/"UMOD.tsv.gz","rt") as f:
 data={r["snp"] for r in csv.DictReader(f,delimiter="\t")}
assert n>0 and n==len(vars) and size==4*n*n and all(r["ref_id"]==v and r["snp"] in data for r,v in zip(meta,vars))
assert any(r["snp"]=="rs11860566" for r in meta)
print(json.dumps({"status":"PRIMARY_INPUT_QC_PASS","variants":n,"binary_bytes":size,"analysis_executed":False}))
