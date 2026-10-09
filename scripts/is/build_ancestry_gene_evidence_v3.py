#!/usr/bin/env python3
"""Ancestry-safe gene/region evidence queue for 80 candidate discovery intervals.

Existing BBJ/EAS QTL/coloc evidence is attached *only* to its original EAS
group + stable Ensembl gene ID. EUR-specific evidence remains UNTESTED.
"""
import argparse,csv,json
from collections import defaultdict
from pathlib import Path

BASE=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1")
OUT=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v2_ancestry")
def rows(p):
    with p.open() as f:return list(csv.DictReader(f,delimiter="\t"))
def stable(x):return str(x or "").split(".")[0].upper()
def build(out,base):
    genes=rows(out/"IS_ALL_GENE_WINDOW_UNIVERSE.tsv")
    manifest=rows(out/"IS_EXPANDED_REGION_EXECUTION_MANIFEST.tsv")
    original=rows(base/"IS_LOCUS_GENE_EVIDENCE_V2.tsv")
    anchors=rows(base/"LEGACY_ANCHOR_CANDIDATES.tsv")
    if len(manifest)<30 or not genes:raise ValueError("Missing expanded ancestry manifest/genes")
    mg={r["group_id"]:r for r in manifest}
    ae={(r["group_id"],stable(r["gene_id"])):r for r in original if r["group_id"]}
    original_symbols={r["gene"].upper() for r in anchors}
    result=[]
    for r in genes:
        group=r["group_id"];locus=mg[group]
        eur=locus["ancestry"]=="EUR"
        baseline=ae.get((group,stable(r["gene_id"]))) if not eur else None
        coloc=baseline.get("coloc_status","NOT_TESTED_IN_THIS_GROUP") if baseline else "NOT_TESTED_IN_THIS_GROUP"
        if not baseline:coloc="NOT_TESTED_IN_THIS_GROUP"
        status="EUR_ANCESTRY_SPECIFIC_TEST_REQUIRED" if eur else "EAS_STUDY_SPECIFIC_TEST_REQUIRED"
        if baseline and coloc=="LOCUS_AND_GENE_ID_MATCHED":
            status="EAS_LEGACY_COLOC_REASSESS_FULL_GWAS"
        result.append({
            "group_id":group,"gene_id":r["gene_id"],"gene_id_stable":stable(r["gene_id"]),
            "gene_symbol":r["gene_symbol"],"biotype":r["biotype"],
            "ancestry":locus["ancestry"],"source":locus["source"],
            "region_gws":locus["has_gws"],"region_lead_variant":locus["lead_variant"],
            "region_best_p":locus["lead_p"],
            "gene_lead_distance_bp":r["lead_distance_bp"],
            "legacy_anchor":int(r["gene_symbol"].upper() in original_symbols),
            "legacy_coloc_status":coloc,
            "legacy_coloc_pp_h4":baseline.get("coloc_max_pp_h4","") if baseline else "",
            "legacy_original_bbj_locus":baseline.get("coloc_original_locus","") if baseline else "",
            "mol_qtl_task":status,
            "eqtl_status":"REASSESS_LEGACY" if baseline and coloc=="LOCUS_AND_GENE_ID_MATCHED" else "NOT_TESTED",
            "sqtl_status":"NOT_TESTED","pqtl_status":"NOT_TESTED",
            "ancestry_ld_state":"EUR_REFERENCE_NOT_AVAILABLE" if eur else "EAS_REFERENCE_AVAILABLE_ONLY_FOR_7_GWS_REGIONS",
            "fine_mapping_status":"BLOCKED_NO_ANCESTRY_MATCHED_VALIDATED_SIGNAL",
            "causal_gene_status":"UNRESOLVED"})
    present={r["gene_symbol"].upper() for r in result if r["legacy_anchor"]}
    for anchor in sorted(original_symbols-present):
        row={c:"" for c in result[0]}
        row.update({"gene_symbol":anchor,"legacy_anchor":1,"source":"LEGACY_ANCHOR_ONLY",
            "mol_qtl_task":"KEEP_AS_EXTERNAL_ANCHOR",
            "legacy_coloc_status":"NOT_EVALUATED_IN_DISCOVERY_GROUP",
            "eqtl_status":"NOT_TESTED","sqtl_status":"NOT_TESTED","pqtl_status":"NOT_TESTED",
            "ancestry_ld_state":"NOT_APPLICABLE","fine_mapping_status":"BLOCKED_NO_DISCOVERY_GROUP",
            "causal_gene_status":"UNRESOLVED"})
        result.append(row)
    with (out/"IS_ANCESTRY_GENE_EVIDENCE_V3.tsv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(result[0]),delimiter="\t")
        w.writeheader();w.writerows(result)
    counts={
        "source_group_rows":len(manifest),
        "gene_region_pairs":len(genes),
        "eas_gene_region_pairs":sum(r["ancestry"]!="EUR" for r in result if r["group_id"]),
        "eur_gene_region_pairs":sum(r["ancestry"]=="EUR" for r in result),
        "legacy_anchor_only_rows":len(result)-len(genes),
        "candidate_evidence_rows":len(result),
        "unique_gene_ids":len({stable(r["gene_id"]) for r in result if r["gene_id"]}),
        "legacy_coloc_region_pairs":sum(r["legacy_coloc_status"]=="LOCUS_AND_GENE_ID_MATCHED" for r in result),
        "eur_region_pairs_inappropriately_reusing_coloc":sum(r["ancestry"]=="EUR" and
            r["legacy_coloc_status"]=="LOCUS_AND_GENE_ID_MATCHED" for r in result),
        "unassessed_eqtl_region_pairs":sum(r["eqtl_status"]=="NOT_TESTED" for r in result if r["group_id"]),
        "unassessed_sqtl_region_pairs":len(genes),
        "unassessed_pqtl_region_pairs":len(genes),
        "warnings":"Legacy QTL cannot be transferred to EUR locus. Positional mapping does not imply causality."}
    if counts["eur_region_pairs_inappropriately_reusing_coloc"]:
        raise RuntimeError("Cross-ancestry evidence propagation violation")
    if counts["gene_region_pairs"]!=counts["eas_gene_region_pairs"]+counts["eur_gene_region_pairs"]:
        raise RuntimeError("Lost gene-region pair")
    (out/"IS_ANCESTRY_GENE_EVIDENCE_V3_SUMMARY.json").write_text(json.dumps(counts,indent=2))
    print(json.dumps(counts,indent=2))
    return counts

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,default=OUT)
    p.add_argument("--base",type=Path,default=BASE)
    a=p.parse_args()
    build(a.out,a.base)
