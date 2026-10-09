#!/usr/bin/env python3
"""Tissue-stratified GTEx v8 DAP-G SNP overlap; not a coloc test."""
import argparse,csv,json
from collections import defaultdict
from pathlib import Path
ROOT=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1/gtex_v8_dapg_liftover")
TISSUE_GROUP={
  "Artery_Aorta":"vascular_artery",
  "Artery_Coronary":"vascular_artery",
  "Artery_Tibial":"vascular_artery",
  "Brain_Cortex":"brain",
  "Brain_Caudate_basal_ganglia":"brain",
  "Brain_Nucleus_accumbens_basal_ganglia":"brain",
  "Brain_Cerebellum":"brain",
  "Whole_Blood":"blood",
}
def summarize(root):
    m=json.loads((root/"G0022_GTEX_V8_DAPG_OVERLAP_SUMMARY.json").read_text())
    inp=root/"G0022_GTEX_V8_DAPG_GWAS_VARIANT_OVERLAPS.tsv"
    acc=defaultdict(list)
    with inp.open() as f:
        for row in csv.DictReader(f,delimiter="\t"):
            if row["interpretation"]!="FINEMAP_VARIANT_OVERLAP_NOT_COLOCALIZATION":
                raise ValueError("Scientific status mismatch")
            if row["gwas_is_95pct_cs"]!="0":
                raise ValueError("Direct GWAS credible-set overlap found; requires separate scrutiny")
            acc[(row["gtex_gene_symbol"],row["tissue"])].append(row)
    total=sum(map(len,acc.values()))
    if total!=m["gwas_allele_exact_finemap_overlap_records"]:
        raise ValueError("Source rows changed")
    out=[]
    for (gene,tissue),records in sorted(acc.items()):
        qtlsets={(r["gtex_method"],r["gtex_finemap_set_id"]) for r in records}
        qtlpip=[float(r["gtex_finemap_pip"]) for r in records]
        gwaspip=[float(r["gwas_ais_full_locus_pip"]) for r in records]
        out.append({"gene_symbol":gene,"gtex_v8_tissue":tissue,
            "tissue_class":TISSUE_GROUP.get(tissue,"other"),
            "n_exact_allele_matching_records":len(records),
            "n_distinct_matched_variants":len({r["variant_grch37_ref_alt"] for r in records}),
            "n_distinct_qtl_finemap_sets":len(qtlsets),
            "max_gtex_qtl_pip":max(qtlpip),
            "max_gwas_full_pip_among_qtl_variants":max(gwaspip),
            "overlap_gwas_95pct_credible_set":0,
            "scientific_gate":"POSITIONAL_QTL_FINEMAP_OVERLAP_NO_COLOCALIZATION"})
    if not out:raise ValueError("No assay variant records")
    path=root/"G0022_GTEX_V8_DAPG_TISSUE_PRIORITIZATION.tsv"
    with path.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(out[0]),delimiter="\t")
        w.writeheader();w.writerows(out)
    summary={"gene_tissue_pairs_with_exact_allele_matched_qtl_variants":len(out),
        "vascular_gene_tissue_pairs":sum(x["tissue_class"]=="vascular_artery" for x in out),
        "brain_gene_tissue_pairs":sum(x["tissue_class"]=="brain" for x in out),
        "gtex_dapg_overlapping_gwas_95pct_cs":0,
        "new_AIS_eQTL_colocalizations":0,
        "status":"QTL_SIGNAL_PROXIMITY_SCREEN_ONLY"}
    (root/"G0022_GTEX_TISSUE_PRIORITIZATION_SUMMARY.json").write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2))
    for item in sorted((r for r in out if r["tissue_class"] in ("brain","vascular_artery")),
                       key=lambda x:-x["n_exact_allele_matching_records"])[:12]:
        print(item)
    return summary
if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--root",type=Path,default=ROOT)
    summarize(p.parse_args().root)
