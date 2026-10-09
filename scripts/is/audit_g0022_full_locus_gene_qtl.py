#!/usr/bin/env python3
"""G0022 full-locus AIS credible-set positional annotation + *prior* BBJ QTL audit.

Never convert positional overlap or old BBJ GTEx coloc to AIS causal evidence.
"""
import argparse,csv,json,math
from collections import defaultdict
from pathlib import Path

ROOT=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1")
GENES=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/IS_ALL_GENE_WINDOW_UNIVERSE.tsv")
COLOC=Path("/srv/is-analysis/results/is/stage5_functional/phase9c_convergence/COLOC_ABF_MASTER_ANNOTATED_V2.tsv")
KEY=("ALDH2","PTPN11","IFT81","ATP2A2")

def rows(path):
    with path.open() as f:return list(csv.DictReader(f,delimiter="\t"))
def annotate(root,genes,prior):
    summary=json.loads((root/"G0022_FULL_AIS_SUSIE_SUMMARY.json").read_text())
    if summary.get("status")!="FULL_LOCUS_EXPLORATORY_ONLY" or not summary.get("converged"):
        raise ValueError("No converged full-locus exploration available")
    if summary.get("n_credible_sets",0)<1:raise ValueError("No local CS supported")
    pips={r["variant_id"]:r for r in rows(root/"G0022_FULL_AIS_PIP.tsv")}
    region=[x for x in rows(genes) if x["group_id"]=="IS_XDATA_G0022"]
    if not region:raise ValueError("Missing GENCODE G0022 gene annotation")
    ids=summary["credible_set_variant_ids"]
    if not isinstance(ids,dict):raise ValueError("Unexpected credible set JSON shape")
    records=[]
    for signal,variants in ids.items():
        for variant in variants:
            if variant not in pips:raise ValueError("CS variant missing from PIP table")
            pos=int(variant.split(":")[1])
            neighbours=[]
            for g in region:
                lo=int(g["gene_start"]);hi=int(g["gene_end"])
                dist=max(lo-pos,pos-hi,0)
                if dist<=500000:neighbours.append((dist,g))
            neighbours.sort(key=lambda v:(v[0],v[1]["gene_id"]))
            for rank,(dist,g) in enumerate(neighbours[:6],1):
                records.append({"signal":signal,"cs_snp":variant,
                    "cs_variant_pip":round(float(pips[variant]["pip"]),10),
                    "position_grch37":pos,"nearest_gene_rank":rank,
                    "gene_id":g["gene_id"],"gene_symbol":g["gene_symbol"],
                    "biotype":g["biotype"],"gene_distance_bp":dist,
                    "overlaps_gene_body":int(dist==0),
                    "status":"POSITIONAL_ONLY_NOT_ALLELE_SPECIFIC_QTL_EVIDENCE"})
    if not records:raise ValueError("No candidate annotations")
    with (root/"G0022_FULL_CS_GENCODE_NEAREST_GENES.tsv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(records[0]),delimiter="\t")
        w.writeheader();w.writerows(records)
    prior_rows=rows(prior)
    gene_by={x["gene_symbol"]:x for x in region}
    aud=[]
    for sym in KEY:
        if sym not in gene_by:raise ValueError("Missing key candidate "+sym)
        gene=gene_by[sym]
        baseid=gene["gene_id"].split(".")[0]
        matching=[x for x in prior_rows if
            x["gene_base"].split(".")[0]==baseid and x["locus"]=="BBJ_IS_L003"]
        valid=[]
        for x in matching:
            try:
                value=float(x["PP.H4"])
                if math.isfinite(value):valid.append((value,x))
            except (TypeError,ValueError):continue
        best=max(valid,key=lambda x:x[0]) if valid else None
        aud.append({"gene_symbol":sym,"gene_id":gene["gene_id"],
            "eas_ais_full_locus_tested_coloc":0,
            "prior_coloc_source":"BBJ_IS_L003",
            "prior_gtex_tissues_count":len(matching),
            "prior_gtex_best_abf_pp_h4":round(best[0],8) if best else "",
            "prior_gtex_best_tissue":best[1]["tissue_label"] if best else "",
            "prior_max_h4_supports_shared_signal":int(bool(best and best[0]>=0.8)),
            "new_qtl_status":"NOT_TESTED_IN_AIS_FULL_LOCUS",
            "association":"POSITIONAL_LINK_ONLY"})
    with (root/"G0022_CANDIDATE_BBJ_QTL_PROVENANCE.tsv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(aud[0]),delimiter="\t")
        w.writeheader();w.writerows(aud)
    result={"g0022_full_ais_credible_set_count":summary["n_credible_sets"],
        "credible_set_snp_total":sum(len(v) for v in ids.values()),
        "nearest_gene_rows":len(records),
        "key_candidates":len(aud),
        "prior_bbj_gtEx_coloc_rows":sum(x["prior_gtex_tissues_count"] for x in aud),
        "new_AIS_QTL_coloc_completed":0,
        "legacy_BBJ_coloc_h4_above_0_8":sum(x["prior_max_h4_supports_shared_signal"] for x in aud),
        "scientific_gate":"INTEGRATED_POSITION_QTL_PROVENANCE_ONLY"}
    (root/"G0022_FULL_LOCUS_POSITIONAL_QTL_AUDIT_SUMMARY.json").write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))
    return result

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--root",type=Path,default=ROOT)
    p.add_argument("--genes",type=Path,default=GENES)
    p.add_argument("--prior",type=Path,default=COLOC)
    a=p.parse_args()
    annotate(a.root,a.genes,a.prior)
