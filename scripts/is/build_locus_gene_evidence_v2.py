#!/usr/bin/env python3
"""IS locus-aware QTL evidence ledger. No positional causality assertion.

Requires explicit overlap between legacy BBJ GWAS interval and expanded region,
AND exact stable Ensembl gene ID and original legacy locus for QTL evidence.
Always use gene_id/group_id pairs; never propagate a coloc by symbol alone.
"""
import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path

DEFAULT_ROOT = Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1")
DEFAULT_BASE = Path("/srv/is-analysis/results/is")

def read(path):
    if not path.is_file():
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))

def stable_id(value):
    return str(value or "").split(".")[0].upper()

def number(value):
    try:
        n = float(value)
        return n if n == n and 0 <= n <= 1 else None
    except (ValueError, TypeError):
        return None

def crosswalk(groups, legacy):
    linked = defaultdict(list)
    for group in groups:
        chrom = group["chr"].removeprefix("chr")
        lo, hi = int(group["start"]), int(group["end"])
        for old in legacy:
            if old["chr"].removeprefix("chr") != chrom:
                continue
            if int(old["end"]) < lo or int(old["start"]) > hi:
                continue
            linked[group["group_id"]].append(old["locus_id"])
    return linked

def build(root, base):
    genes = read(root/"IS_ALL_GENE_WINDOW_UNIVERSE.tsv")
    groups = read(root/"CROSS_DATASET_INTERVAL_GROUPS.tsv")
    legacy = read(base/"stage1_gwas_qc/japan/bbj/BBJ_IS_PROVISIONAL_LOCI.tsv")
    coloc = read(base/"stage5_functional/phase9c_convergence/COLOC_ABF_MASTER_ANNOTATED_V2.tsv")
    human = read(base/"stage5_functional/phase9f_e_human_r3_4_import/HUMAN_TARGET_GENE_SUMMARY.tsv")
    mouse = read(base/"stage5_functional/phase9f_c_target_expression/CELLTYPE_TARGET_EVIDENCE_MASTER.tsv")
    anchors = read(root/"LEGACY_ANCHOR_CANDIDATES.tsv")
    if not genes or not groups or not legacy:
        raise RuntimeError("Required positional genes, groups, or BBJ loci are unavailable")
    by_group = {r["group_id"]: r for r in groups}
    old_group = crosswalk(groups, legacy)
    qtl = defaultdict(list)
    qc = {"raw_coloc_rows":len(coloc),"pass_rows_with_probability":0}
    for row in coloc:
        p = number(row.get("PP.H4"))
        if row.get("status","").upper() != "PASS" or p is None:
            continue
        gid = stable_id(row.get("gene_base"))
        old = row.get("locus","")
        if not gid or not old:
            continue
        qtl[(old, gid)].append((p, number(row.get("PP.H3")), row))
        qc["pass_rows_with_probability"] += 1
    hidx = {r.get("gene","").upper():r for r in human if r.get("gene")}
    midx = defaultdict(list)
    for row in mouse:
        if row.get("human_gene"):
            midx[row["human_gene"].upper()].append(row)
    anchor_names = {r["gene"].upper() for r in anchors}
    output=[]
    for r in genes:
        group = r["group_id"]
        stable = stable_id(r["gene_id"])
        symbol = r["gene_symbol"].upper()
        linked = [(old, entry) for old in old_group.get(group,[]) for entry in qtl.get((old,stable),[])]
        best = max(linked, key=lambda z:z[1][0]) if linked else None
        h = hidx.get(symbol)
        m = midx.get(symbol,[])
        expanded = by_group[group]
        output.append({
            "group_id":group,"gene_id":r["gene_id"],"gene_id_stable":stable,
            "gene_symbol":r["gene_symbol"],"biotype":r["biotype"],
            "gene_start":r["gene_start"],"gene_end":r["gene_end"],
            "lead_distance_bp":r["lead_distance_bp"],
            "group_has_gws":expanded["has_gws"],"group_best_p":expanded["lead_p"],
            "group_phenotypes":expanded["phenotypes"],
            "matched_legacy_bbj_loci":";".join(sorted(old_group.get(group,[]))),
            "legacy_anchor":int(symbol in anchor_names),
            "coloc_status":"LOCUS_AND_GENE_ID_MATCHED" if best else "NOT_TESTED_IN_THIS_GROUP",
            "coloc_max_pp_h4":best[1][0] if best else "",
            "coloc_pp_h3_at_best":best[1][1] if best and best[1][1] is not None else "",
            "coloc_original_locus":best[0] if best else "",
            "coloc_tissue":best[1][2].get("tissue_label","") if best else "",
            "coloc_dataset":best[1][2].get("dataset_key","") if best else "",
            "human_atlas_status":("PRESENT" if h.get("present","").upper()=="TRUE" else "NOT_DETECTED_IN_TARGET_ATLAS") if h else "NOT_TESTED",
            "human_top_celltype":h.get("top_celltype","") if h else "",
            "mouse_status":"TARGET_DESCRIPTIVE_ONLY" if m else "NOT_TESTED",
            "molecular_scope":"TARGETED_LEGACY" if (best or h or m) else "UNSCREENED",
            "causal_gene_status":"UNRESOLVED"})
    seen_anchors={r["gene_symbol"].upper() for r in output if r["legacy_anchor"]}
    # Keep anchor-only symbols without falsely declaring them part of a discovery region.
    for name in sorted(anchor_names-seen_anchors):
        empty={k:"" for k in output[0].keys()}
        h = hidx.get(name)
        m = midx.get(name, [])
        empty.update({"gene_symbol":name,"legacy_anchor":1,"coloc_status":"NOT_TESTED_IN_THIS_GROUP",
            "human_atlas_status":("PRESENT" if h.get("present","").upper()=="TRUE" else "NOT_DETECTED_IN_TARGET_ATLAS") if h else "NOT_TESTED",
            "human_top_celltype":h.get("top_celltype","") if h else "",
            "mouse_status":"TARGET_DESCRIPTIVE_ONLY" if m else "NOT_TESTED",
            "molecular_scope":"ANCHOR_OUTSIDE_CURRENT_WINDOWS","causal_gene_status":"UNRESOLVED"})
        output.append(empty)
    fields=list(output[0])
    with (root/"IS_LOCUS_GENE_EVIDENCE_V2.tsv").open("w",newline="") as handle:
        writer=csv.DictWriter(handle,fieldnames=fields,delimiter="\t")
        writer.writeheader();writer.writerows(output)
    qc.update({"total_gene_locus_rows":len(output),
       "mapped_gene_locus_rows":len(genes),
       "anchor_only_rows":len(output)-len(genes),
       "locus_gene_matched_coloc":sum(x["coloc_status"]=="LOCUS_AND_GENE_ID_MATCHED" for x in output),
       "human_target_atlas_assessed":sum(x["human_atlas_status"]!="NOT_TESTED" for x in output),
       "mouse_target_assessed":sum(x["mouse_status"]!="NOT_TESTED" for x in output),
       "legacy_loci_linked":{g:";".join(sorted(v)) for g,v in sorted(old_group.items())},
       "interpretation":"No cross-locus or symbol-only coloc propagation. Missing is not negative. Not a causal ranking."})
    (root/"IS_LOCUS_GENE_EVIDENCE_V2_SUMMARY.json").write_text(json.dumps(qc,indent=2,ensure_ascii=False))
    return qc

if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--root",type=Path,default=DEFAULT_ROOT)
    parser.add_argument("--base",type=Path,default=DEFAULT_BASE)
    args=parser.parse_args()
    print(json.dumps(build(args.root,args.base),indent=2,ensure_ascii=False))
