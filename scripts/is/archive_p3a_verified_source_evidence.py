#!/usr/bin/env python3
"""Pin all prior 36×2 P3A original-source aggregate evidence permanently in git.

Input files are exact, already source-verified GitHub Actions JSON artifacts,
NOT source SNP rows, donor metadata or genotype data. Original SHA-256 digests,
schemas, 36×2 cross-audit invariants and privacy are checked before copying.

This utility *never* pushes git commits itself. GitHub Actions may commit only
to the explicitly whitelisted research branch after all validations pass.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

COVERAGE = __import__(
    "scripts.is.audit_p3a_36x2_original_fullcis_gwas_coverage",
    fromlist=["INPUTS", "compute_36x2", "load_verified_artifact"],
)
PINNED_COVERAGE_SHA256 = (
    "6fc64e2068b44e45a1bf16e1162e20062e4913cca9fc1d3b2e9b702cd118e344"
)
ARCHIVE_SCHEMA = "is-p3a-36x2-persistent-github-pinned-source-evidence-v1"
ARCHIVE_FILES = {
    "reciprocal": COVERAGE.INPUTS["reciprocal"]["filename"],
    "cis": COVERAGE.INPUTS["cis"]["filename"],
    "qtl_n": COVERAGE.INPUTS["qtl_n"]["filename"],
    "coverage": "IS_P3A_36x2_original_cis_allele_coverage_HOLD.json",
}
SOURCE_RUNS = {
    **{k:v["run_id"] for k,v in COVERAGE.INPUTS.items()},
    "coverage":45189787019,
}
PRIVACY_DISALLOWED_EXACT_KEYS = {
    "beta","se","varbeta","pvalue","effect_allele","other_allele",
    "rsid","sample_id","donor_id","subject_id","individual_id",
    "genotype","genotypes","dosage","expression_values",
    "variant_id","snp_id","participant_id","pseudo_id",
}
MAX_ORIGINAL_AGGREGATE_BYTES=220_000


def require(p, m):
    if not p:raise ValueError(m)


def sha256(raw:bytes):
    return hashlib.sha256(raw).hexdigest()


def validate_privacy_and_shape(obj, kind):
    require(isinstance(obj,dict), "Only aggregate JSON documents are permitted")
    rows=obj.get("rows") or obj.get("groups") or []
    require(isinstance(rows,list) and len(rows)==(72 if kind=="coverage" else 36),
            "Per-gene grouped source count required; raw source rows prohibited")
    require(kind!="coverage" or all(
        "QTL_dataset" in item and "gene_symbol" in item and "ancestry" in item
        for item in rows), "Coverage row identity missing")
    forbidden=[]

    def walk(value, depth=0):
        require(depth<=10,"Unexpected deep nested original source schema")
        if isinstance(value,dict):
            for key,subvalue in value.items():
                require(isinstance(key,str), "Unexpected nonstring key")
                if key.lower() in PRIVACY_DISALLOWED_EXACT_KEYS:
                    forbidden.append(key)
                walk(subvalue,depth+1)
        elif isinstance(value,list):
            require(len(value)<=72,"Possible per-SNP vector not source aggregate")
            for subvalue in value:walk(subvalue,depth+1)
        elif isinstance(value,str):
            require(len(value)<=600,"Suspicious long genotype or sample data value")
        else:
            require(value is None or isinstance(value,(int,float,bool)),
                    "Unexpected nonserializable original source metadata")
    walk(obj)
    require(not forbidden,"Disallowed individual SNP effect or donor keys")
    require(obj["scientific_status"].startswith("HOLD"),
            "All original evidence must preserve scientific HOLD")
    for row in rows:
        require(row.get("QTL_dataset") in COVERAGE.DATASETS,
                "Unknown genuine original QTL study identity")
        require(isinstance(row.get("gene_id"),str) and
                row["gene_id"].startswith("ENSG"),
                "Unexpected real original gene identifier")
    return len(rows)


def check_and_archive(inputs, manifest, dest):
    require(set(inputs)==set(ARCHIVE_FILES),
            "Require exactly four pinned original evidence files")
    sources={}
    provenance={}
    raw_bytes={}
    for kind,name in ARCHIVE_FILES.items():
        path=Path(inputs[kind])
        raw=path.read_bytes()
        require(len(raw)<=MAX_ORIGINAL_AGGREGATE_BYTES,
                "Original evidence too large, may contain per-variant data")
        expected_sha=(PINNED_COVERAGE_SHA256 if kind=="coverage"
                      else COVERAGE.INPUTS[kind]["sha256"])
        require(sha256(raw)==expected_sha,
                "Previously audited original source aggregate SHA256 mismatch: "+kind)
        obj=json.loads(raw.decode("utf8"))
        require(validate_privacy_and_shape(obj,kind)==
                (72 if kind=="coverage" else 36),
                "Source row count drift")
        if kind=="coverage":
            require(obj["schema"]=="is-p3a-original-36x2-gwas-qtl-fullcis-nonpal-coverage-gate-v1"
                    and obj["technical_status"]==
                        "PASS_PINNED_3_ORIGINAL_SOURCE_AUDITS_AND_72_RECONCILED_GROUPS",
                    "Latest source-verified coverage schema/status changed")
        else:
            require(obj["schema"]==COVERAGE.INPUTS[kind]["schema"]
                    and obj["technical_status"]==COVERAGE.INPUTS[kind]["status"],
                    "Publisher aggregate source schema/status changed")
        sources[kind]=obj
        raw_bytes[kind]=raw
        provenance[kind]={
            "filename":name,"sha256":sha256(raw),"original_source_CI_run_id":SOURCE_RUNS[kind],
            "verified_source_GitHub_Actions_URL":
                "https://github.com/kimtk94/IS_Analysis_V3/actions/runs/"+
                str(SOURCE_RUNS[kind]),
            "original_artifact_group_level_only":True,
            "original_archive_file_bytes":len(raw),
        }
    rebuilt=COVERAGE.compute_36x2(manifest,sources["reciprocal"],
                                   sources["cis"],sources["qtl_n"])
    actual=sources["coverage"]
    for key in ("rows","aggregate_by_ancestry","EUR_only_eight_discordant_gene_assay_rows",
                "n_assay_ancestry_rows","n_permutation_vs_nominal_variant_count_exact_groups"):
        require(rebuilt[key]==actual[key],
                "Original 72-row independent recomputation differs: "+key)
    require(len({(r["QTL_dataset"],r["gene_id"],r["ancestry"])
                for r in rebuilt["rows"]})==72,
            "72 unique ancestry-gene-tissue identifiers not recovered")
    require(rebuilt["aggregate_by_ancestry"]["EAS"]["matched"]==89597 and
            rebuilt["aggregate_by_ancestry"]["EUR"]["matched"]==128282 and
            rebuilt["aggregate_by_ancestry"]["EUR"]["discordant"]==8 and
            rebuilt["n_permutation_vs_nominal_variant_count_exact_groups"]==36,
            "Original 36×2 scientific source count baseline changed")
    require(rebuilt["formal_coloc_abf_eligible"]==0
            and rebuilt["formal_coloc_abf_executed"]==0
            and rebuilt["formal_coloc_susie_executed"]==0,
            "Must never auto-promote Coloc on archived aggregation alone")
    dst=Path(dest)
    require(dst.as_posix().startswith("docs/is/evidence/verified_p3a_20261124"),
            "Evidence archive destination must remain controlled repository subtree")
    dst.mkdir(parents=True,exist_ok=True)
    for kind,name in ARCHIVE_FILES.items():
        (dst/name).write_bytes(raw_bytes[kind])
    index={
        "schema":ARCHIVE_SCHEMA,
        "source_audit_asof":"2026-11-23",
        "original_genuine_gene_tissue_assays":36,
        "original_assay_x_ancestry_aggregate_rows":72,
        "original_distinct_genes":20,
        "original_independent_GWAS_loci":1,
        "original_all_permutation_tested_cis_gene_assay_variants":202400,
        "original_strict_nonpal_QTL_GWAS_candidate_gene_assay_SNPs":148344,
        "original_signed_match_gene_assay_records":{
            "EAS":89597,"EUR":128282},
        "original_unresolved_allele_discordant_assay_records":{
            "EAS":0,"EUR":8},
        "public_aggregate_evidence_permanently_in_git":True,
        "source_evidence_file_sha256_and_GitHub_run":provenance,
        "source_inputs_are_not_original_per_variant_effect_data":True,
        "QTL_AN_over_2_is_not_independently_confirmed_exact_gene_regression_N":True,
        "ROSMAP_2021_QC_pass_samples":569,
        "ROSMAP_R7_published_nominal_N":560,
        "ROSMAP_2021_QC_vs_R7_9_sample_gap_explained":False,
        "official_source_QTL_normalized_phenotype_sdY_directly_measured":False,
        "GWAS_cc_beta_SE_path_requires_case_fraction_s":False,
        "scientific_status":"HOLD_QTL_MODEL_N_ORIGINAL_SD_Y_AND_FULLCIS_PRIOR_SENSITIVITY",
        "formal_coloc_abf_eligible":0,
        "formal_coloc_abf_executed":0,
        "formal_coloc_susie_executed":0,
        "original_aggregate_evidence_donor_or_variant_effect_rows_archived":False,
    }
    idx=dst/"INDEX.json"
    idx.write_text(json.dumps(index,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    return index


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--manifest",type=Path,required=True)
    p.add_argument("--reciprocal",type=Path,required=True)
    p.add_argument("--cis",type=Path,required=True)
    p.add_argument("--qtl-n",type=Path,required=True)
    p.add_argument("--coverage",type=Path,required=True)
    p.add_argument("--output-dir",type=Path,required=True)
    a=p.parse_args()
    manifest=json.loads(a.manifest.read_text(encoding="utf8"))
    idx=check_and_archive({"reciprocal":a.reciprocal,"cis":a.cis,
                           "qtl_n":a.qtl_n,"coverage":a.coverage},
                          manifest,a.output_dir)
    print(json.dumps({
        "technical_status":"PASS_FOUR_SOURCE_SHA256_AND_DURABLE_AGGREGATE_ONLY_ARCHIVE",
        "scientific_status":idx["scientific_status"],
        "n_original_files":len(ARCHIVE_FILES),
        "n_ancestry_group_rows":idx["original_assay_x_ancestry_aggregate_rows"],
        "n_formal_coloc_executed":0,
    }),flush=True)


if __name__=="__main__":
    main()
