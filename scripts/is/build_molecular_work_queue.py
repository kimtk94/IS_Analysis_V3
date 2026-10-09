#!/usr/bin/env python3
"""Build full-universe molecular-QTL and LD execution queue, never a gene ranking."""
import argparse
import csv
import json
from collections import Counter,defaultdict
from pathlib import Path

DEFAULT=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1")
def read(p):
    with p.open(newline="") as f: return list(csv.DictReader(f,delimiter="\t"))

def build(root):
    genes=read(root/"IS_LOCUS_GENE_EVIDENCE_V2.tsv")
    regions=read(root/"IS_EXPANDED_REGION_EXECUTION_MANIFEST.tsv")
    allele=read(root/"LEAD_ALLELE_REFERENCE_AUDIT.tsv")
    if not genes or not regions or not allele:raise ValueError("Missing required inputs")
    regions_by={r["group_id"]:r for r in regions}
    allele_by={r["group_id"]:r for r in allele}
    queue=[]
    for row in genes:
        group=row["group_id"]
        if not group:
            continue
        region=regions_by[group]
        a=allele_by[group]
        measured=row["coloc_status"]=="LOCUS_AND_GENE_ID_MATCHED"
        # Existing coloc is legacy, not necessarily measured in every subtype.
        queue.append({
            "group_id":group,
            "gene_id":row["gene_id"],
            "gene_id_stable":row["gene_id_stable"],
            "gene_symbol":row["gene_symbol"],
            "biotype":row["biotype"],
            "gwas_phenotypes":row["group_phenotypes"],
            "gws_supported_region":row["group_has_gws"],
            "lead_p":row["group_best_p"],
            "window_chr":region["chr"],
            "window_start":region["window_start"],
            "window_end":region["window_end"],
            "lead_variant":region["lead_variant"],
            "existing_coloc":"LEGACY_MATCHED" if measured else "NOT_TESTED",
            "eqtl_task":"REASSESS_EXPANDED_GWAS" if measured else "NEW_QTL_SCREEN_REQUIRED",
            "sqtl_task":"NEW_QTL_SCREEN_REQUIRED",
            "pqtl_task":"NEW_QTL_SCREEN_REQUIRED",
            "celltype_qtl_task":"NEW_QTL_SCREEN_REQUIRED",
            "ld_reference_status":a["status"],
            "ld_finemap_gate":"BLOCKED_PENDING_ANCESTRY_AND_ALLELE_QC",
            "outcome":"NO_CAUSAL_GENE_VERDICT"
        })
    queue.sort(key=lambda r:(r["group_id"],r["gene_id"]))
    with (root/"IS_ALL_CANDIDATE_MOLECULAR_WORK_QUEUE.tsv").open("w",newline="") as f:
        w=csv.DictWriter(f,delimiter="\t",fieldnames=list(queue[0]))
        w.writeheader();w.writerows(queue)
    counts=defaultdict(lambda:Counter())
    for r in queue:
        c=counts[r["group_id"]]
        c["n_genes"]+=1
        c["new_qtl_screen"]+=int(r["eqtl_task"]=="NEW_QTL_SCREEN_REQUIRED")
        c["legacy_matched"]+=int(r["existing_coloc"]=="LEGACY_MATCHED")
    per=[]
    for region in regions:
        key=region["group_id"];c=counts[key]
        per.append({"group_id":key,"has_gws":region["has_gws"],
            "phenotypes":region["phenotypes"],"gene_pairs":c["n_genes"],
            "new_eqtl_needed":c["new_qtl_screen"],
            "legacy_coloc_to_reassess":c["legacy_matched"],
            "lead_allele_status":allele_by[key]["status"],
            "qtl_screen_status":"NEEDS_NEW_DATA_OR_RECOMPUTE",
            "ld_status":"BLOCKED_PENDING_LD"})
    with (root/"IS_MOLECULAR_TASKS_BY_REGION.tsv").open("w",newline="") as f:
        w=csv.DictWriter(f,delimiter="\t",fieldnames=list(per[0]))
        w.writeheader();w.writerows(per)
    output={"region_count":len(per),"gene_region_tasks":len(queue),
        "new_eqtl_candidates":sum(r["eqtl_task"]=="NEW_QTL_SCREEN_REQUIRED" for r in queue),
        "legacy_coloc_reassess":sum(r["existing_coloc"]=="LEGACY_MATCHED" for r in queue),
        "sqtl_screen_tasks":len(queue),"pqtl_screen_tasks":len(queue),
        "genotype_ready_regions":0,
        "policy":"No ancestry-matched LD plus allele harmonization => block fine-mapping; missing molecular data are not negative",
        "scope":"Gene-region work queue; no downloaded molecular QTL has been added"}
    (root/"IS_MOLECULAR_WORK_QUEUE_SUMMARY.json").write_text(json.dumps(output,indent=2))
    return output

if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--root",type=Path,default=DEFAULT)
    a=ap.parse_args()
    print(json.dumps(build(a.root),indent=2))
