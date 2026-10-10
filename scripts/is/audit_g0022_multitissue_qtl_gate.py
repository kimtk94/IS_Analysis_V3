#!/usr/bin/env python3
"""Fail-closed G0022 4-tissue molecular-QTL evidence ledger.

Retain original AIS CS, no post-hoc gene selection by diagnostic coloc PP.H4.
"""
import argparse,csv,json
from collections import Counter
from pathlib import Path
from prepare_g0022_eqtl_catalogue_multitissue import OUT,SOURCES,source_file
def csvrows(path):
    with path.open() as h:return list(csv.DictReader(h,delimiter="\t"))
def audit(root):
    abf=json.loads((root/"G0022_MULTITISSUE_ABF_SUMMARY.json").read_text())
    results=csvrows(root/"G0022_GTEX_MULTITISSUE_ABF_DIAGNOSTIC.tsv")
    source=json.loads((root/"G0022_FULL_AIS_FOUR_CS_QTL_SOURCE_COVERAGE_SUMMARY.json").read_text())
    if abf.get("analyzed_studies")!=4 or abf.get("gene_tissue_tests")!=32 or len(results)!=32:
        raise RuntimeError("Unexpected study/target-gene denominator")
    if abf.get("computed_posteriors")!=24 or abf.get("posterior_results_eligible_for_publication")!=0:
        raise RuntimeError("Unexpected diagnostic/validated coloc boundary")
    if source.get("dataset_variant_tests")!=16 or source.get("by_coverage_status")!={"NOT_TESTED_VARIANT_POSITION_ABSENT":16}:
        raise RuntimeError("Source 4 SNP x 4 tissue position coverage differs from recorded interpretation")
    seen={(x["dataset_id"],x["gene"]) for x in results}
    if len(seen)!=32:raise RuntimeError("Repeated gene-tissue result")
    if any(x["gate"]=="VALIDATED_COLOCALIZATION" or x["conclusion"]!="NOT_VALIDATED_COLOCALIZATION" for x in results):
        raise RuntimeError("Unexpected causal colocalization promotion")
    dataset_manifest=[];detail=[]
    for study,(tissue,n) in SOURCES.items():
        summary=json.loads((root/study/"G0022_MULTITISSUE_QTL_INPUT_SUMMARY.json").read_text())
        if summary["dataset_id"]!=study or summary["catalogue_sample_size"]!=n:
            raise RuntimeError("Dataset source sample provenance mismatch")
        if summary["GWAS_CS_SNPs_in_QTL_extract_any_gene"] or summary["genes_ready_for_diagnostic"]>0:
            raise RuntimeError("Four-CS input coverage gate unexpectedly passed")
        dataset_manifest.append({"dataset_id":study,"tissue":tissue,"catalogue_sample_size":n,
            "GRCh38_region":summary["region_grch38"],
            "QTL_source_url":f"https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/QTS000015/{study}/{study}.all.tsv.gz",
            "regional_rows":summary["qc"]["source_records"],"source_sha256":summary["raw_sha256"],
            "gene_snp_matches":summary["target_gene_snp_matches"],"readiness":"NO_SHARED_GWAS_CS_SNP"})
        for r in results:
            if r["dataset_id"]!=study:continue
            detail.append({"dataset_id":study,"tissue":tissue,"gene":r["gene"],
                "matched_snps":r["matched_snps"],
                "coloc_abf_diagnostic_PP_H4":r["PP.H4"],
                "co_loc_snp_coverage":"0_OF_4_GWAS_CS_VARIANTS_ASSAYED",
                "scientific_gate":"BLOCKED_INCOMPLETE_QTL_SNP_UNIVERSE",
                "causal_gene_status":"NOT_VALIDATED"})
    with (root/"G0022_MULTITISSUE_QTL_SOURCE_MANIFEST.tsv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(dataset_manifest[0]),delimiter="\t")
        w.writeheader();w.writerows(dataset_manifest)
    with (root/"G0022_MULTITISSUE_QTL_PUBLICATION_GATE.tsv").open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(detail[0]),delimiter="\t")
        w.writeheader();w.writerows(detail)
    verdict={"status":"PUBLICATION_GATE_FAIL_COVERAGE_NOT_BIOLOGICAL_NEGATIVE",
      "gtEx_gene_tissue_comparisons":len(detail),
      "distinct_QTL_datasets":len(dataset_manifest),
      "total_QTL_source_associations":sum(x["regional_rows"] for x in dataset_manifest),
      "total_gene_snp_matches":sum(x["gene_snp_matches"] for x in dataset_manifest),
      "four_credible_set_snps_source_positions_found_any_tissue":0,
      "diagnostic_ABF_posteriors_computed":abf["computed_posteriors"],
      "valid_colocalizations":0,"validated_causal_genes":0,
      "next_action":"Find ancestry-matched molecular QTL full SNP coverage or independent functional/clinical evidence; do not treat PP.H4 diagnostic as causal."}
    (root/"G0022_MULTITISSUE_QTL_PUBLICATION_GATE_SUMMARY.json").write_text(json.dumps(verdict,indent=2))
    print(json.dumps(verdict,indent=2))
    return verdict
if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--root",type=Path,default=OUT)
    audit(p.parse_args().root)
