#!/usr/bin/env python3
"""GWS per-study lead-allele and harmonized-input manifest.

Requires completed, already-audited reference-overlap rows. Produces a queue,
not SuSiE fine-mapping or independent locus verification.
"""
import argparse
import csv
import gzip
import json
from collections import Counter,defaultdict
from pathlib import Path

ROOT=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1")
FIELDS=["dataset","phenotype","group_id","chr","pos","variant_id","ref","alt",
        "effect_allele","other_allele","beta","se","p","eaf",
        "alt_effect_beta","alt_effect_eaf","qc_status","source_variant_id",
        "source_ref","source_alt","ancestry","build"]

def process(root):
    with (root/"IS_EXPANDED_REGION_EXECUTION_MANIFEST.tsv").open() as f:
        groups=[r for r in csv.DictReader(f,delimiter="\t") if r["has_gws"]=="1"]
    with (root/"IS_GWS_GWAS_REFERENCE_OVERLAP.tsv").open() as f:
        cover=list(csv.DictReader(f,delimiter="\t"))
    if len(groups)!=7 or not cover:
        raise ValueError("Expected completed seven-GWS reference overlap audit")
    expected={(r["dataset"],r["group_id"]):r for r in cover}
    observed=defaultdict(Counter)
    leads=defaultdict(list)
    inp=root/"IS_GWS_REFERENCE_MATCHED_VARIANTS.tsv.gz"
    with gzip.open(inp,"rt") as f:
        rows=csv.DictReader(f,delimiter="\t")
        if list(rows.fieldnames or []) != FIELDS:
            raise ValueError("Harmonization artifact schema mismatch")
        pos_by={g["group_id"]:int(g["lead_variant"].split(":")[1]) for g in groups}
        for r in rows:
            key=(r["dataset"],r["group_id"])
            observed[key][r["qc_status"]]+=1
            if int(r["pos"])==pos_by[r["group_id"]]:
                leads[key].append(r)
    output=[]
    for key,cov in sorted(expected.items()):
        dataset,gid=key
        c=observed[key]
        n=sum(c.values())
        if n!=int(cov["exported_matching_pairs"]):
            raise ValueError("Matched TSV vs overlap summary mismatch for "+str(key))
        lead=leads[key]
        group=next(g for g in groups if g["group_id"]==gid)
        wanted=group["lead_variant"].split(":")
        allelematches=[x for x in lead if set([x["ref"],x["alt"]])==set(wanted[2:4])]
        if len(allelematches)>1:raise ValueError("Ambiguous multiple GWAS lead alleles "+str(key))
        best=allelematches[0] if allelematches else None
        leadstatus=("LEAD_PRESENT_ORIENTED" if best and best["qc_status"].startswith("MATCH_")
            else "LEAD_PALINDROMIC_UNRESOLVED" if best else "LEAD_ABSENT_OR_ALLELE_UNMATCHED")
        output.append({
            "dataset":dataset,"group_id":gid,"chr":group["chr"],"phenotypes":group["phenotypes"],
            "lead_variant":group["lead_variant"],"source_window_n":int(cov["window_variants"]),
            "matched_nonpal_n":c["MATCH_REF_EFFECT"]+c["MATCH_ALT_EFFECT"],
            "palindromic_review_n":c["PALINDROMIC_REVIEW"],
            "matched_total_n":n,
            "n_lead_records":len(lead),"lead_allele_status":leadstatus,
            "lead_effect_alt_beta":best["alt_effect_beta"] if best else "",
            "lead_p":best["p"] if best else "",
            "study_input_ready":"AUDITED_VARIANTS_AVAILABLE_ONLY" if n else "NO_MATCHED_VARIANTS",
            "multi_signal_finemap_ready":"NO_LD_MATRIX_OR_SIGNAL_QC",
            "independent_replication":"UNKNOWN_COHORT_OVERLAP"
        })
    if sum(x["matched_total_n"] for x in output)!=sum(observed[k].total() for k in expected):
        raise ValueError("Study count consistency failure")
    dest=root/"IS_GWS_STUDY_INPUT_READINESS.tsv"
    with dest.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(output[0]),delimiter="\t")
        w.writeheader();w.writerows(output)
    summary={"group_count":len(groups),"study_region_pairs":len(output),
        "with_nonpal_variants":sum(x["matched_nonpal_n"]>0 for x in output),
        "with_oriented_lead":sum(x["lead_allele_status"]=="LEAD_PRESENT_ORIENTED" for x in output),
        "with_palindromic_lead":sum(x["lead_allele_status"]=="LEAD_PALINDROMIC_UNRESOLVED" for x in output),
        "no_lead_or_pair_mismatch":sum(x["lead_allele_status"]=="LEAD_ABSENT_OR_ALLELE_UNMATCHED" for x in output),
        "ready_for_multi_signal_finemap":0,
        "note":"Alleles are reference-matched in a limited regional GWAS subset. No ancestry-matched LD covariance matrix or locus independence verified."}
    (root/"IS_GWS_STUDY_INPUT_READINESS_SUMMARY.json").write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2))
    return summary
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--root",type=Path,default=ROOT)
    a=p.parse_args()
    process(a.root)
