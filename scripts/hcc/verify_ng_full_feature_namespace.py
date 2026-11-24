#!/usr/bin/env python3
"""Fail-closed source identity namespace checks; NO full-mapping overclaims."""
import argparse
import json
import pathlib
import re

COMMIT="438f52271d0a86a6701baa7a96e80c516ba575c2"
def check(r):
    if r.get("schema")!="hcc-ng-full-namespace-v1" or r.get("source_commit")!=COMMIT:
        raise ValueError("wrong full-feature source")
    for k in ("full_sha256","info_sha256"):
        if not re.fullmatch("[a-f0-9]{64}",str(r.get(k,""))):
            raise ValueError("original source SHA256 missing")
    if r.get("n_full_RNA")!=20356 or r.get("n_full_protein")!=5585:
        raise ValueError("source geometry changed")
    caps={"RNA_ENSG_like_n":20356,"protein_UniProt_like_n":5585,
          "protein_piped_UniProt_like_n":5585,
          "RNA_unique_feature_labels_n":20356,
          "protein_unique_feature_labels_n":5585,
          "raw_feature_names_direct_overlap_n":5585,
          "RNA_ENSG_version_suffix_n":20356,
          "selected_RNA_lookup_matches_full_n":188,
          "selected_protein_lookup_matches_full_n":116}
    for k,cap in caps.items():
        n=r.get(k)
        if type(n)!=int or not 0<=n<=cap:
            raise ValueError("invalid source ID metric: "+k)
    if r.get("RNA_ENSG_like_n")!=20356:
        raise ValueError("RNA ENSG namespace must match source")
    if r.get("Ng2022_paper_comparable_genes_benchmark")!=5490 or r.get(
        "Ng2022_paper_mRNA_protein_correlation_benchmark")!=5481:
        raise ValueError("incorrect source paper validation benchmarks")
    for k in ("full_ensembl_uniprot_mapping_verified",
              "same_biopsy_rna_protein_abundance_effect_tested",
              "feature_or_patient_ids_exported"):
        if r.get(k) is not False:
            raise ValueError("unreviewed mapping/biological claim")
    id_accessions=r["protein_UniProt_like_n"]+r["protein_piped_UniProt_like_n"]
    return {"source_status":"SOURCE_FULL_FEATURE_NAMESPACE_AUDITED",
            "RNA_full_Ensembl_labels":r["RNA_ENSG_like_n"],
            "protein_full_accession_like_labels":id_accessions,
            "source_P_selected_mapping_full_label_overlap":
                r["selected_protein_lookup_matches_full_n"],
            "original_paper_benchmark_comparable_n":5490,
            "original_paper_benchmark_correlation_n":5481,
            "full_gene_protein_join":"HOLD",
            "biological_inference":"NOT_PERMITTED",
            "no_individual_identifiers_or_values_in_output":True}
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--source",required=True)
    args=ap.parse_args()
    p=pathlib.Path(args.source)
    if not p.is_file() or p.stat().st_size>30000:
        raise SystemExit("source report unavailable")
    print(json.dumps(check(json.loads(p.read_text(encoding="utf-8"))),indent=2))
if __name__=="__main__":main()
