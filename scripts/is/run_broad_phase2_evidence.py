#!/usr/bin/env python3
"""Reproduce v2 evidence / all-candidate task inputs without modifying canonical GWAS."""
import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1")
BASE=Path("/srv/is-analysis/results/is")
REGISTRY=Path("/srv/is-analysis/data/is/reference/gigastroke/metadata/GIGASTROKE_STUDY_MAP_V2.tsv")
DATA=Path("/srv/is-analysis/data/is")
def sha(path):
    digest=hashlib.sha256()
    with path.open("rb") as h:
        for block in iter(lambda:h.read(1024*1024),b""):
            digest.update(block)
    return digest.hexdigest()
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--root",type=Path,default=ROOT)
    ap.add_argument("--base",type=Path,default=BASE)
    ap.add_argument("--registry",type=Path,default=REGISTRY)
    ap.add_argument("--data",type=Path,default=DATA)
    a=ap.parse_args()
    steps=[
        ["build_locus_gene_evidence_v2.py","--root",str(a.root),"--base",str(a.base)],
        ["build_molecular_work_queue.py","--root",str(a.root)],
        ["audit_gigastroke_full_registry.py","--root",str(a.root),"--registry",str(a.registry),"--data",str(a.data)]
    ]
    for args in steps:
        p=subprocess.run([sys.executable,str(HERE/args[0])]+args[1:],capture_output=True,text=True)
        if p.returncode!=0:
            print(p.stdout,file=sys.stderr)
            print(p.stderr,file=sys.stderr)
            raise SystemExit(f"STOP: {args[0]} failed without changing canonical inputs")
    gene=json.loads((a.root/"IS_ALL_GENE_UNIVERSE_SUMMARY.json").read_text())
    evidence=json.loads((a.root/"IS_LOCUS_GENE_EVIDENCE_V2_SUMMARY.json").read_text())
    queue=json.loads((a.root/"IS_MOLECULAR_WORK_QUEUE_SUMMARY.json").read_text())
    if queue["gene_region_tasks"]!=gene["region_gene_associations"]:
        raise SystemExit("FAILED invariant: all positional gene-region pairs must be queued")
    if evidence["mapped_gene_locus_rows"]!=gene["region_gene_associations"]:
        raise SystemExit("FAILED invariant: evidence must account for every mapped pair")
    files=["IS_LOCUS_GENE_EVIDENCE_V2.tsv","IS_ALL_CANDIDATE_MOLECULAR_WORK_QUEUE.tsv",
           "IS_MOLECULAR_TASKS_BY_REGION.tsv","GIGASTROKE_FULL_STUDY_INVENTORY.tsv"]
    audit={"status":"PASS","inputs":"pre-existing gene universe/30-region manifest/legacy molecular screens",
       "region_count":gene["regions"],"gene_region_pairs":gene["region_gene_associations"],
       "locus_matched_coloc":evidence["locus_gene_matched_coloc"],
       "new_eqtl_tasks":queue["new_eqtl_candidates"],
       "sha256":{file:sha(a.root/file) for file in files},
       "guardrail":"Cannot infer independent loci or causal genes from these outputs"}
    (a.root/"IS_BROAD_PHASE2_REPRODUCIBILITY.json").write_text(json.dumps(audit,indent=2))
    print(json.dumps(audit,indent=2))
if __name__=="__main__":
    main()
