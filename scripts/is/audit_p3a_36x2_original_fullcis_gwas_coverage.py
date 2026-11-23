#!/usr/bin/env python3
"""Fail-closed original IS P3A 36 QTL assays x EAS/EUR GWAS coverage accounting.

Reconcile three immutable source-verified GitHub Actions artifacts:
- 36-assay reciprocal original GRCh38↔GRCh37 FASTA/allele mapping
- 36-assay original QTL nominal vs permutation cis tested variant census
- 36-assay original QTL genotype AN/2 provenance
Do NOT conflate all-cis variants with tested canonical nonpal SNPs; do NOT
interpret gene×tissue SNP records as genomic independent loci or coloc results.
No raw SNP beta/se/MAF or participant records are ingested in this job.
"""
from __future__ import annotations
import argparse
import collections
import hashlib
import json
import math
from pathlib import Path

INPUTS = {
    "reciprocal": {
        "run_id": 45031834139,
        "sha256": "53b6d71229ed98ab79f404194fa80e94e95bd838fa119ba2961ef234f1308ba8",
        "filename": "P3A_36_real_source_2way_ref_SNP_EAS_EUR_GWAS_HOLD.json",
        "schema": "is-p3a-36-real-source-reciprocal-SNP-GWAS-harmonization-v1",
        "status": "PASS_QTL_GWAS_SOURCE_SNP_HARMONIZATION_AGGREGATES_ONLY",
    },
    "cis": {
        "run_id": 45036525615,
        "sha256": "4a65a2e4ae9f260c9c154cb8f7a47af96ca95c1492185cd0a9739dbd178f7481",
        "filename": "IS_P3A_36_actual_source_permuted_cis_tested_nvariants_HOLD.json",
        "schema": "is-p3a-36-original-eQTL-permuted-n_variants-nominal-completeness-v1",
        "status": "PASS_3_ORIGINAL_PERMUTATION_AND_NOMINAL_REGION_SOURCES_READ",
    },
    "qtl_n": {
        "run_id": 45038801844,
        "sha256": "023750fde4a2a57c1604793dfa2b78ba65a29dbe06ff849e17382a559906eb2a",
        "filename": "IS_P3A_36_original_QTL_SE_MAF_AN_residual_proxy_N_model_HOLD.json",
        "schema": "is-p3a-36-primary-qtlsdY-proxy-GWAS-case-control-effect-size-path-v1",
        "status": "PASS_3_REAL_ORIGINAL_QTL_SD_Y_PROXY_INPUT_DIAGNOSTICS",
    },
}
DATASETS = {
    "QTD000131": {"label":"GTEx Artery Aorta", "n":10,
                  "all_cis":55808,"strict":44693},
    "QTD000171": {"label":"GTEx Brain Cortex", "n":11,
                  "all_cis":56748,"strict":45570},
    "QTD000434": {"label":"ROSMAP DLPFC", "n":15,
                  "all_cis":89844,"strict":58081},
}
ANCESTRIES = {"EAS":"GCST90104545","EUR":"GCST90104540"}
EXCLUSIVE_GWAS_STATES = (
    "MATCH_ALT_EFFECT", "MATCH_REF_EFFECT_SIGN_FLIPPED",
    "NO_GWAS_POSITION", "HOLD_ALLELE_PAIR_MISMATCH",
    "HOLD_COMPLEMENT_ONLY_NO_ORIENTATION_APPROVAL",
    "HOLD_PALINDROMIC_STRAND",
    "HOLD_INVALID_GWAS_BETA_OR_EAF",
    "HOLD_NONCANONICAL_ALLELES",
)
EXPECTED_TOTAL = {
    "cis": 202400,
    "strict_nonpal": 148344,
    "EAS_signed":89597,
    "EUR_signed":128282,
    "EAS_absent":58747,
    "EUR_absent":20054,
    "EUR_conflict":8,
}
PREDECLARED_NONPAL_COVERAGE_REVIEW_THRESHOLD = 0.8


def require(condition, message):
    if not condition:
        raise ValueError(message)


def load_verified_artifact(path, kind):
    spec=INPUTS[kind]
    raw=Path(path).read_bytes()
    actual_sha=hashlib.sha256(raw).hexdigest()
    require(actual_sha==spec["sha256"], "Pinned "+kind+" source artifact SHA256 mismatch")
    obj=json.loads(raw.decode("utf8"))
    require(obj.get("schema")==spec["schema"]
            and obj.get("technical_status")==spec["status"]
            and obj.get("scientific_status")=="HOLD",
            "Original audited "+kind+" source evidence status/schema changed")
    return obj, {
        "source_run_id":spec["run_id"],
        "source_original_report_name":spec["filename"],
        "sha256":actual_sha,
        "verified_against_immutable_previously_audited_artifact":True,
        "original_SNP_level_data_redownloaded":False,
    }


def indexed_groups(rows):
    require(isinstance(rows,list) and len(rows)==36,
            "Wrong number of original group rows")
    result={}
    for row in rows:
        key=(row["QTL_dataset"],row["gene_id"])
        require(key not in result and row["QTL_dataset"] in DATASETS,
                "Duplicate or unknown original gene-tissue source group")
        result[key]=row
    return result


def rank_summary(groups):
    # Technical prioritization only: not causal gene prioritization or a
    # threshold for coloc eligibility. Prefer high common-SNP denominator.
    ordered=sorted(groups,key=lambda r:(-r["fraction_QTL_strict_nonpal_with_GWAS_match"],
                                        -r["matched_nonpal_SNP_assay_records"],
                                        r["QTL_dataset"],r["gene_symbol"]))
    return [{"gene_symbol":row["gene_symbol"],"QTL_dataset":row["QTL_dataset"],
             "matched_nonpal_SNP_assay_records":row["matched_nonpal_SNP_assay_records"],
             "fraction_QTL_strict_nonpal_with_GWAS_match":
                  row["fraction_QTL_strict_nonpal_with_GWAS_match"],
             "this_is_NOT_causal_priority":True}
            for row in ordered[:5]]


def compute_36x2(manifest, source_snp, source_cis, source_n):
    orig=manifest.get("groups",[])
    require(manifest.get("schema")=="is-p3a-36-real-group-keys-from-legacy-393-v1"
            and manifest.get("scientific_status")=="HOLD"
            and len(orig)==36
            and len({(x["QTL_dataset"],x["gene_id"]) for x in orig})==36
            and len({x["gene_id"] for x in orig})==20,
            "Original 36 source assay manifest identities drift")
    keys={(x["QTL_dataset"],x["gene_id"]) for x in orig}
    tables={
        "reciprocal":indexed_groups(source_snp["rows"]),
        "cis":indexed_groups(source_cis["groups"]),
        "qtl_n":indexed_groups(source_n["groups"])
    }
    require(all(set(t)==keys for t in tables.values()),
            "Original 36 assay source per-gene ID mismatch")
    require(source_snp.get("n_real_assay_groups")==36
            and source_cis.get("n_gene_tissue_groups_matching_full_cis_tested_variant_COUNT")==36
            and source_n.get("n_real_gene_tissue_assays")==36,
            "Original source 36-assay completeness missing")
    require(source_snp.get("n_formal_coloc_abf_executed")==0
            and source_cis.get("n_formal_coloc_abf_executed")==0
            and source_n.get("n_formal_coloc_abf_executed")==0,
            "Do not silently switch into a formal Coloc workflow")
    for ancestry, accession in ANCESTRIES.items():
        x=source_snp["GWAS_ancestries"][ancestry]
        require(x["original_publisher_accession"]==accession
                and x["original_publisher_full_file_MD5_SHA256_bytes_pass"] is True
                and x["GWAS_target_invalid_numeric_or_allele_rows"]==0
                and x["GWAS_target_duplicate_locus_rows"]==0,
                "Original GIGASTROKE uploaded GWAS source failure")
    full_per_tissue=collections.Counter()
    strict_per_tissue=collections.Counter()
    by_ancestry=collections.defaultdict(lambda:collections.Counter())
    results=[]
    mismatches=[]
    for record in orig:
        key=(record["QTL_dataset"],record["gene_id"])
        s=tables["reciprocal"][key]
        c=tables["cis"][key]
        n=tables["qtl_n"][key]
        require(record["gene_symbol"]==s["gene_symbol"]==c["gene_symbol"]==n["gene_symbol"]
                and s["QTL_source_indexed_nonpal_valid"]==s["QTL_uniquely_reciprocally_lifted_REF_valid"]
                and not s["QTL_liftover_blockers"],
                "Source gene phenotype identity / reciprocal FASTA liftover changed")
        count=c["nominal_distinct_all_variant_records"]
        strict=s["QTL_source_indexed_nonpal_valid"]
        classes=c["original_variant_class_counts"]
        require(c["published_permutation_vs_recovered_nominal_unique_agree"] is True
                and c["source_original_cis_n_variants_from_permutation"]==count
                and c["nominal_strict_nonpal_SNP_source_records"]==strict
                and sum(classes.values())==count
                and classes.get("INDEL",0)>=0 and classes.get("SNP",0)>=strict
                and count>=strict>0,
                "Real QTL original complete cis and SNP variant class census mismatched")
        source_an=n["QTL_source_variant_N_and_residual_proxy"]
        require(source_an["original_strict_nonpal_SNPs"]==strict
                and source_an["positive_even_genotyped_AN_count"]==strict
                and source_an["invalid_AN_AC_count"]==0
                and source_an["AN_over_2_higher_than_catalogue_study_N"]==0,
                "Original QTL genotype AN/2 source requirements drift")
        full_per_tissue[key[0]]+=count
        strict_per_tissue[key[0]]+=strict
        excluded_classes={
            "source_INDEL_all_cis_not_in_nonpal_subset":classes.get("INDEL",0),
            "source_SNP_palindromic_or_noncanonical_or_filtered":
                classes["SNP"]-strict,
        }
        require(sum(excluded_classes.values())==count-strict,
                "Original cis variant class exclusions do not partition full cis")
        for ancestry,accession in ANCESTRIES.items():
            original=s["per_ancestry"][ancestry]
            require(original["dataset"]==accession,"Ancestry and source GWAS accession differ")
            reasons=original["counts_by_signed_GWAS_allele_match"]
            require(set(reasons).issubset(set(EXCLUSIVE_GWAS_STATES)|{
                        "AF_ABS_DELTA_GT_0_2_DESCRIPTIVE"
                    }),"Unknown GWAS reason present")
            counts={name:reasons.get(name,0) for name in EXCLUSIVE_GWAS_STATES}
            require(all(isinstance(v,int) and v>=0 for v in counts.values()),
                    "Invalid GWAS source count")
            matched=counts["MATCH_ALT_EFFECT"]+counts["MATCH_REF_EFFECT_SIGN_FLIPPED"]
            absent=counts["NO_GWAS_POSITION"]
            discordant=sum(counts[x] for x in EXCLUSIVE_GWAS_STATES
                           if x.startswith("HOLD_"))
            require(matched==original["n_nonpal_QTL_GWAS_recipient_allele_matches"]
                    and matched+absent+discordant==strict
                    and original["n_valid_QTL_ALT_AF_vs_GWAS_ALT_EAF_pairs"]<=matched,
                    "GWAS allele signed-match vs original QTL strict SNP denominator drift")
            af_flag=reasons.get("AF_ABS_DELTA_GT_0_2_DESCRIPTIVE",0)
            require(0<=af_flag<=matched,
                    "AF difference must be SUBSET of matched records, not exclusive category")
            if discordant:
                mismatches.append({
                    "QTL_dataset":key[0],
                    "gene_symbol":s["gene_symbol"],
                    "ancestry":ancestry,
                    "nonpal_discordant_original_assay_rows":discordant,
                    "not_independent_variant_or_locus_count":True,
                })
            by_ancestry[ancestry].update({
                "original_all_cis":count,
                "strict":strict,
                "matched":matched,
                "absent":absent,
                "discordant":discordant,
                "allele_alt":counts["MATCH_ALT_EFFECT"],
                "allele_ref_flipped":counts["MATCH_REF_EFFECT_SIGN_FLIPPED"],
                "af_gt_0_2_descriptive":af_flag,
                "SNP_excluded":excluded_classes["source_SNP_palindromic_or_noncanonical_or_filtered"],
                "INDEL_excluded":excluded_classes["source_INDEL_all_cis_not_in_nonpal_subset"],
            })
            nonpal_fraction=matched/strict
            all_fraction=matched/count
            results.append({
                "QTL_dataset":key[0],
                "gene_id":key[1],
                "gene_symbol":s["gene_symbol"],
                "ancestry":ancestry,
                "GWAS_accession":accession,
                "original_permutation_full_cis_tested_variant_assay_records":count,
                "original_nominal_nonpal_canonical_SNP_assay_records":strict,
                "source_cis_variants_excluded_before_GWAS_matching":excluded_classes,
                "source_original_reciprocal_liftover_allele_REF_success":strict,
                "matched_nonpal_SNP_assay_records":matched,
                "matched_original_GWAS_effect_ALTs":counts["MATCH_ALT_EFFECT"],
                "matched_original_GWAS_effect_REF_requires_beta_flip":counts[
                    "MATCH_REF_EFFECT_SIGN_FLIPPED"],
                "original_GWAS_position_absent_for_QTL_nonpal_SNP":absent,
                "source_GWAS_allele_discordance_HOLD":discordant,
                "original_allele_mismatch_categories":{
                    name:counts[name] for name in counts if name.startswith("HOLD_")
                    and counts[name]},
                "original_QTL_vs_GWAS_ALT_AF_delta_gt_0_2_DESCRIPTIVE_nonexclusive":af_flag,
                "original_SNP_AN_over_2_min":source_an["min_observed_SNP_N"],
                "original_SNP_AN_over_2_max":source_an["max_observed_SNP_N"],
                "original_SNP_genotype_N_is_not_approved_regression_N":True,
                "fraction_QTL_strict_nonpal_with_GWAS_match":round(nonpal_fraction,7),
                "fraction_QTL_ALL_original_cis_with_nonpal_GWAS_match":round(all_fraction,7),
                "excluded_palin_noncanonical_or_indel_is_not_missing_GWAS_position":True,
                "threshold_for_technical_nonpal_coverage_review_only":
                    PREDECLARED_NONPAL_COVERAGE_REVIEW_THRESHOLD,
                "at_least_80pct_strict_SNP_match_TECHNICAL_SCREEN_ONLY":
                    nonpal_fraction>=PREDECLARED_NONPAL_COVERAGE_REVIEW_THRESHOLD,
                "original_QTL_direct_sdY_gene_regression_N_approved":False,
                "complete_coloc_tested_full_cis_alleles_approved":False,
                "native_Coloc_case_control_beta_SE_requires_case_fraction_s":False,
                "formal_coloc_abf_eligible":False,
                "formal_coloc_susie_eligible":False,
                "scientific_status":"HOLD_QTL_MODEL_AND_UNMATCHED_FULL_CIS_AND_PRIOR_SENSITIVITY",
            })
    require({k:full_per_tissue[k] for k in DATASETS}=={
                k:v["all_cis"] for k,v in DATASETS.items()}
            and {k:strict_per_tissue[k] for k in DATASETS}=={
                k:v["strict"] for k,v in DATASETS.items()},
            "Original per-tissue complete cis + nonpal original source sum drift")
    totals={a:dict(x) for a,x in by_ancestry.items()}
    for ancestry in ANCESTRIES:
        t=totals[ancestry]
        require(t["original_all_cis"]==EXPECTED_TOTAL["cis"]
                and t["strict"]==EXPECTED_TOTAL["strict_nonpal"]
                and t["matched"]==EXPECTED_TOTAL[ancestry+"_signed"]
                and t["absent"]==EXPECTED_TOTAL[ancestry+"_absent"],
                "Original two-ancestry 36 group signed-allele summary count drift")
        require(t["matched"]+t["absent"]+t["discordant"]==t["strict"]
                and t["matched"]+t["absent"]+t["discordant"]+
                    t["SNP_excluded"]+t["INDEL_excluded"]==t["original_all_cis"],
                "All original cis QTL gene–assay variants must be accounted for")
        t["matched_fraction_among_strict_nonpal"] = round(t["matched"]/t["strict"],7)
        t["matched_fraction_among_original_all_cis"] = round(
            t["matched"]/t["original_all_cis"],7)
        t["fraction_excluded_before_matching_not_GWAS_absence"] = round(
            (t["SNP_excluded"]+t["INDEL_excluded"])/t["original_all_cis"],7)
        t["assay_groups_with_at_least_80pct_nonpal_GWAS_overlap"]=sum(
            x["at_least_80pct_strict_SNP_match_TECHNICAL_SCREEN_ONLY"]
            for x in results if x["ancestry"]==ancestry)
        subset=[x for x in results if x["ancestry"]==ancestry]
        t["min_per_assay_fraction_nonpal_matched"]=min(
            x["fraction_QTL_strict_nonpal_with_GWAS_match"] for x in subset)
        t["max_per_assay_fraction_nonpal_matched"]=max(
            x["fraction_QTL_strict_nonpal_with_GWAS_match"] for x in subset)
    require(totals["EAS"]["discordant"]==0
            and totals["EUR"]["discordant"]==EXPECTED_TOTAL["EUR_conflict"]
            and len(results)==72 and len(mismatches)==8
            and sum(x["nonpal_discordant_original_assay_rows"] for x in mismatches)==8,
            "Original EUR 8 ROSMAP allele discordances lost or misclassified")
    for ancestry in ANCESTRIES:
        require(source_snp["matched_nonpal_summary"][ancestry+"_match"]==
                totals[ancestry]["matched"]
                and source_snp["matched_nonpal_summary"][ancestry+"_no_GWAS_pos"]==
                    totals[ancestry]["absent"],
                "Original audited top-line matched summary differs from groups")
    return {
        "schema":"is-p3a-original-36x2-gwas-qtl-fullcis-nonpal-coverage-gate-v1",
        "technical_status":"PASS_PINNED_3_ORIGINAL_SOURCE_AUDITS_AND_72_RECONCILED_GROUPS",
        "scientific_status":"HOLD_ACTUAL_QTL_MODEL_N_SDY_COMPLETE_CIS_PRIORS",
        "qtl_genome_build":"GRCh38",
        "gwas_genome_build":"GRCh37",
        "n_original_gene_tissue_assays":36,
        "n_original_unique_genes":20,
        "n_original_gwas_loci":1,
        "n_original_gwas_ancestries":2,
        "n_assay_ancestry_rows":len(results),
        "n_permutation_vs_nominal_variant_count_exact_groups":36,
        "dataset_source_counts":DATASETS,
        "aggregate_by_ancestry":totals,
        "technical_predeclared_strict_nonpal_coverage_threshold":
            PREDECLARED_NONPAL_COVERAGE_REVIEW_THRESHOLD,
        "technical_only_highest_5_fraction_each_ancestry": {
            ancestry:rank_summary([x for x in results if x["ancestry"]==ancestry])
            for ancestry in ANCESTRIES},
        "EUR_only_eight_discordant_gene_assay_rows":mismatches,
        "notes_on_denominators":{
            "all_cis":"Publisher original 36 permuted n_variants, includes SNP and INDEL, deduplicated by gene assay",
            "nonpal":"Canonical nonpal single SNPs, excludes palindromes, INDEL, noncanonical records and other filters",
            "GWAS_absent":"GWAS target GRCh37 chrom position missing ONLY among original candidate nonpal QTLs after verified reciprocal mapping",
            "allele_AF_discrepancy":"NONEXCLUSIVE descriptive category among signed allele matched; not an additional failed SNP",
            "records_not_independent_SNPs":"One genomic SNP can be present in multiple gene×tissue assays; one chr10 locus only",
            "prior_original_QTL_GWAS_and_sdy_est_source_evidence_not_rerun":True,
        },
        "original_gwas_chr10_beta_log_OR_source_verified_prior_725707_rows":True,
        "original_binary_cc_beta_varbeta_s_required":False,
        "original_gene_specific_expression_model_N_verified":False,
        "original_genotype_N_equals_regression_N_proven":False,
        "original_transformed_expression_phenotype_sdY_directly_verified":False,
        "scientifically_full_cis_all_variant_harmonization_APPROVED":False,
        "source_independent_EAS_BBJ_replication_APPROVED":False,
        "formal_coloc_abf_eligible":0,
        "formal_coloc_abf_executed":0,
        "formal_coloc_susie_eligible":0,
        "formal_coloc_susie_executed":0,
        "no_per_SNP_effect_values_or_patient_ids_exported":True,
        "rows":results,
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--manifest",type=Path,required=True)
    p.add_argument("--reciprocal",type=Path,required=True)
    p.add_argument("--cis",type=Path,required=True)
    p.add_argument("--qtl-n",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    outcome={
        "schema":"is-p3a-original-36x2-gwas-qtl-fullcis-nonpal-coverage-gate-v1",
        "technical_status":"HOLD_BASELINE_ARTIFACT_OR_GROUP_QC",
        "scientific_status":"HOLD",
        "formal_coloc_abf_eligible":0,
        "formal_coloc_abf_executed":0
    }
    try:
        manifest=json.loads(args.manifest.read_text(encoding="utf8"))
        sources={}
        provenance={}
        for kind,path in (("reciprocal",args.reciprocal),
                          ("cis",args.cis),("qtl_n",args.qtl_n)):
            sources[kind],provenance[kind]=load_verified_artifact(path,kind)
        outcome=compute_36x2(manifest,sources["reciprocal"],
                             sources["cis"],sources["qtl_n"])
        outcome["original_source_artifact_provenance"]=provenance
    except (ValueError,OSError,UnicodeError,KeyError,TypeError,ArithmeticError,
            OverflowError) as error:
        outcome["error_class_only"]=type(error).__name__
        outcome["safe_error_reason"]=str(error)[:180]
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(outcome,indent=2,ensure_ascii=False)+"\n",
                        encoding="utf8")
    print(json.dumps({
        "technical_status":outcome["technical_status"],
        "scientific_status":outcome["scientific_status"],
        "n_original_groups":outcome.get("n_original_gene_tissue_assays"),
        "n_ancestry_group_rows":outcome.get("n_assay_ancestry_rows"),
        "ancestry_aggregates":outcome.get("aggregate_by_ancestry"),
        "n_formal_coloc":0,
        "safe_error_reason":outcome.get("safe_error_reason"),
    },ensure_ascii=False),flush=True)
    return 0 if outcome["technical_status"].startswith("PASS_") else 1


if __name__=="__main__":
    raise SystemExit(main())
