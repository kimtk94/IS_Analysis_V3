#!/usr/bin/env python3
"""Stage 23: fail-closed aggregate QC for Ng/Suter UniProt-to-Ensembl mapping.

Validate the *mapping evidence*, never declare original gene-set identity or
biological significance from a similar count. Input/output contain no IDs.
"""
import argparse
import datetime as dt
import json
import pathlib

INT_FIELDS = (
    "source_RNA_ENSG_n", "source_proteome_features_n",
    "source_protein_features_with_recognized_UniProt_accession_n",
    "source_protein_accessions_found_in_uniprot_reviewed_n",
    "protein_features_with_Ensembl_gene_in_source_RNA_n",
    "protein_features_with_unique_Ensembl_gene_n",
    "distinct_Ensembl_genes_corresponding_to_proteins_n",
    "one_gene_one_protein_unambiguous_genes_n",
    "multi_gene_protein_feature_count", "unmapped_or_noncanonical_source_protein_n",
    "source_paper_comparable_gene_benchmark_n",
    "source_paper_correlated_gene_benchmark_n",
    "absolute_deviation_from_paper_comparable_n",
    "uniprot_reviewed_entries_with_gene_ENSG_in_TSV_field_n",
    "uniprot_reviewed_entries_with_gene_ENSG_after_json_n",
    "canonical_accession_isoform_collapse_excess_features_n",
)

def evaluate(data):
    if data.get("schema") != "hcc-ng-full-uniprot-ensembl-live-audit-v1":
        raise ValueError("unrecognized source schema")
    for key in INT_FIELDS:
        val = data.get(key)
        if type(val) is not int or val < 0:
            raise ValueError("missing/invalid nonnegative integer: "+key)
    if data.get("reviewed_only_annotation_policy") is not True:
        raise ValueError("reviewed-only JSON policy not evidenced")
    if data.get("exact_original_paper_2022_mapping_reproduced") is not False:
        raise ValueError("cannot assert exact paper identity without paper gene list")
    if data.get("clinical_effect_or_pathway_significance_tested") is not False:
        raise ValueError("mapping is not biological inference")
    if data.get("patient_sample_gene_accession_or_expression_values_exported") is not False:
        raise ValueError("unexpected feature or patient-level export flag")
    a=data
    rna=a["source_RNA_ENSG_n"]
    total=a["source_proteome_features_n"]
    recognized=a["source_protein_features_with_recognized_UniProt_accession_n"]
    reviewed=a["source_protein_accessions_found_in_uniprot_reviewed_n"]
    mapped=a["protein_features_with_Ensembl_gene_in_source_RNA_n"]
    unique=a["protein_features_with_unique_Ensembl_gene_n"]
    multi=a["multi_gene_protein_feature_count"]
    distinct=a["distinct_Ensembl_genes_corresponding_to_proteins_n"]
    unambig=a["one_gene_one_protein_unambiguous_genes_n"]
    paper=a["source_paper_comparable_gene_benchmark_n"]
    if (rna,total,paper,a["source_paper_correlated_gene_benchmark_n"]) != (20356,5585,5490,5481):
        raise ValueError("source/paper cardinality benchmark drift")
    if not (0 <= reviewed <= recognized <= total and
            0 <= unique <= mapped <= recognized and
            unique+multi == mapped and
            0 <= unambig <= unique and unambig <= distinct <= rna and
            total-mapped == a["unmapped_or_noncanonical_source_protein_n"]):
        raise ValueError("crosswalk cardinality contradiction")
    if a["absolute_deviation_from_paper_comparable_n"] != abs(distinct-paper):
        raise ValueError("paper comparator deviation inconsistent")
    if a["uniprot_reviewed_entries_with_gene_ENSG_in_TSV_field_n"] > a["uniprot_reviewed_entries_with_gene_ENSG_after_json_n"]:
        raise ValueError("JSON addition cannot erase valid TSV references")
    if a["canonical_accession_isoform_collapse_excess_features_n"] > recognized:
        raise ValueError("invalid canonical accession collapse")
    fb=a.get("accession_specific_json_enrichment")
    if fb is not None:
        if type(fb.get("batches")) is not int or type(fb.get("batches_completed")) is not int or fb["batches"] != fb["batches_completed"]:
            raise ValueError("JSON batches incomplete")
    try:
        dt.datetime.fromisoformat(a["annotation_snapshot_utc"])
    except (ValueError,KeyError,TypeError) as exc:
        raise ValueError("annotation timestamp missing") from exc
    # Count matching is not gene identity matching, even if the count is exact.
    return {
        "schema": "hcc-ng-stage23-crosswalk-review-v1",
        "audit_status": "PASS_AGGREGATE_COUNTS_ONLY",
        "scientific_gate": "V03_PARTIAL_V07_NOT_STARTED",
        "quantitative_whole_proteome_mapping_engineering": "GO_EXPLORATORY",
        "original_paper_5490_gene_set_replicated": False,
        "original_paper_5481_correlation_test_replicated": False,
        "new_biological_significance": False,
        "feature_to_gene_mapping_n": mapped,
        "unique_feature_to_gene_mapping_n": unique,
        "multi_gene_feature_n": multi,
        "distinct_ensembl_gene_union_n": distinct,
        "strict_one_gene_one_protein_genes_n": unambig,
        "paper_comparable_count_n": paper,
        "paper_count_delta": distinct-paper,
        "paper_count_ratio_not_identity_recall": round(distinct/paper, 6),
        "unmapped_proteins_n": total-mapped,
        "isoform_collapse_excess_features_n": a["canonical_accession_isoform_collapse_excess_features_n"],
        "provenance": {
            "annotation_source": a["annotation_source"],
            "annotation_snapshot_utc": a["annotation_snapshot_utc"],
            "annotation_tsv_sha256": a["annotation_tsv_sha256"],
            "reviewed_only_json": True,
            "tsv_ensembl_gene_rows_n": a["uniprot_reviewed_entries_with_gene_ENSG_in_TSV_field_n"],
            "reviewed_ensembl_gene_rows_after_json_n": a["uniprot_reviewed_entries_with_gene_ENSG_after_json_n"]
        },
        "limitations": [
            "Count similarity cannot establish overlap of the original 5490 identities",
            "Multi-gene and isoform aliases need source-consistent resolution",
            "Current UniProt annotations need comparison with historical Ng GRCh37",
            "Preprocessed Suter RDS cannot be treated as raw RSEM expected counts",
            "50 biopsies represent 48 patients and biological analyses require patient-level dependence handling",
            "No paired RNA-protein correlation, clinical effect, or recurrence test has been performed"
        ]
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--input",required=True)
    ap.add_argument("--output",required=True)
    args=ap.parse_args()
    data=json.loads(pathlib.Path(args.input).read_text(encoding="utf-8"))
    result=evaluate(data)
    dest=pathlib.Path(args.output)
    dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(result,ensure_ascii=False,indent=2))

if __name__=="__main__":
    main()
