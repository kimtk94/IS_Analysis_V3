#!/usr/bin/env python3
"""Validate transcript-to-gene identity for the exact pilot GENCODE FASTA."""
import argparse
import collections
import csv
import gzip
import hashlib
import json
import math
import re
from pathlib import Path
from ocular_fixed_panel_score_preflight import read_panel, PIN

REFERENCE_MD5 = "a30369535731a2fdaa6add9b3530756a"

def parse_headers(lines):
    transcripts = {}
    gene_versions = collections.defaultdict(set)
    for line in lines:
        if not line.startswith(">"):
            continue
        fields = line[1:].strip().split("|")
        if len(fields) < 2:
            raise ValueError("unexpected GENCODE header")
        t, g = fields[:2]
        # GENCODE may distinguish chrY pseudoautosomal copies using _PAR_Y.
        # Keep those identities distinct rather than silently collapsing copies.
        if not re.fullmatch(r"ENST\d+\.\d+(?:_PAR_Y)?", t) or not re.fullmatch(r"ENSG\d+\.\d+(?:_PAR_Y)?", g):
            raise ValueError("unexpected versioned transcript or gene ID")
        if t in transcripts:
            raise ValueError("duplicate transcript")
        stable = re.sub(r"\.\d+(?=_PAR_Y$|$)", "", g)
        transcripts[t] = (g, stable)
        gene_versions[stable].add(g)
    if not transcripts or any(len(v)!=1 for v in gene_versions.values()):
        raise ValueError("empty mapping or ambiguous stable gene version")
    return transcripts, gene_versions

def audit(fasta, panel_path, output):
    panel = read_panel(panel_path)
    md5 = hashlib.md5()
    sha = hashlib.sha256()
    with fasta.open("rb") as f:
        for block in iter(lambda: f.read(1024*1024), b""):
            md5.update(block)
            sha.update(block)
    if md5.hexdigest() != REFERENCE_MD5:
        raise ValueError("reference checksum mismatch")
    with gzip.open(fasta, "rt", encoding="utf8") as f:
        tx, genes = parse_headers(f)
    groups = {}
    for name, ids in (("all241", list(panel)), ("ECM16", [g for g,(_,e) in panel.items() if e])):
        mapped = set(ids) & genes.keys()
        total = len(mapped)
        pos = sum(panel[g][0]==1 for g in mapped)
        neg = total-pos
        min_total = math.ceil(.75*len(ids))
        min_pos = math.ceil(.7*sum(panel[g][0]==1 for g in ids))
        min_neg = math.ceil(.7*sum(panel[g][0]==-1 for g in ids))
        groups[name] = {"frozen":len(ids),"mapped":total,"positive":pos,"negative":neg,"unmapped":sorted(set(ids)-mapped),"reference_coverage_gate":total>=min_total and pos>=min_pos and neg>=min_neg}
    output.mkdir(parents=True,exist_ok=True)
    path = output/"tx2gene.tsv"
    with path.open("w",newline="",encoding="utf8") as f:
        w=csv.writer(f,delimiter="\t",lineterminator="\n")
        w.writerow(["transcript_id","gene_id_versioned","gene_id"])
        for t,(g,s) in sorted(tx.items()):
            w.writerow([t,g,s])
    result={"status":"REFERENCE_MAPPING_VERIFIED" if all(x["reference_coverage_gate"] for x in groups.values()) else "HOLD_REFERENCE_COVERAGE","scientific_validation":False,"reference_md5":md5.hexdigest(),"reference_sha256":sha.hexdigest(),"panel_sha256":PIN,"transcripts":len(tx),"genes":len(genes),"tx2gene_sha256":hashlib.sha256(path.read_bytes()).hexdigest(),"groups":groups,"note":"Reference ID coverage only; no expression, detectable-gene or sample QC. Mapping keys retain versioned transcript ID, and PAR_Y copies are separate."}
    (output/"summary.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))
    if result["status"]!="REFERENCE_MAPPING_VERIFIED":
        raise SystemExit(1)

def selftest():
    t,g=parse_headers([">ENST0001.1|ENSG0001.2|x\n","ACGT\n",">ENST0002.3|ENSG0001.2|x\n",">ENST0001.1_PAR_Y|ENSG0001.2_PAR_Y|x\n"])
    assert len(t)==3 and len(g)==2 and t["ENST0001.1_PAR_Y"][1]=="ENSG0001_PAR_Y"
    for bad in ([">ENST0001.1|ENSG0001.2|x",">ENST0001.1|ENSG0001.2|x"],[">ENST0001.1|ENSG0001.2|x",">ENST0002.1|ENSG0001.3|x"],[">not_a_transcript"],["ACGT"]):
        try: parse_headers(bad)
        except ValueError: pass
        else: raise AssertionError("expected invalid identity rejection")
    print("TX2GENE_SYNTHETIC_PASS: 5 cases; no real expression data.")

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--selftest",action="store_true")
    p.add_argument("--fasta",type=Path)
    p.add_argument("--panel",type=Path)
    p.add_argument("--output",type=Path)
    a=p.parse_args()
    if a.selftest:selftest()
    else:
        if not all((a.fasta,a.panel,a.output)):p.error("fasta, panel and output required")
        audit(a.fasta,a.panel,a.output)
