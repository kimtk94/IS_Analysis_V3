#!/usr/bin/env python3
"""Read-only audit of previously produced FASTA-anchored IS evidence.

Archived file fields are evidence *claims*, not an independent new reference
genome check. Report denominators by rows and unique genomic variants
separately. No coloc, SuSiE, MR, or causal gene conclusions.
"""
from __future__ import annotations
import collections
import csv
import gzip
import math
from pathlib import Path

ARCHIVED_RELATIVE=(
    "is_eas5_98_egene_fullcis_gwas_overlap_20261011_v3_fastaanchored/"
    "IS_EAS5_98_FDR_EGENE_FULLCIS_BBJ_EAS_MATCHED_VARIANTS.tsv.gz"
)
READY_RELATIVE=(
    "is_eas5_qtl_gwas_official_alt_harmonization_20261011_v1/"
    "IS_EAS5_GENE_TISSUE_ALLELE_HARMONIZED_READINESS.tsv"
)
OVERLAP_RELATIVE=(
    "is_eas5_qtl_gwas_snp_ref37_20261011_v2/"
    "IS_EAS5_GENE_TISSUE_GWAS_SNP_OVERLAP_QC.tsv"
)
TRUE={"true","1","yes"}
FALSE={"false","0","no","na","none",""}
MAX_ROWS=2_000_000
REQUIRED_ARCHIVE={
 "region_id","QTL_dataset","gene_id","mapped_genomic_variant_GRCh37",
 "BBJ_allele_exact_Fasta37","EAS_AIS_allele_exact_Fasta37",
 "QTL_p","BBJ_p","EAS_AIS_p","QTL_beta_GRCh37_ALT",
 "BBJ_beta_ALT","EAS_AIS_beta_ALT",
 "QTL_se",
}
REQUIRED_OVERLAP={
 "region_id","QTL_dataset","gene_id",
 "GRCh37_FASTA_VERIFIED_SNP_rows","FASTA_REF_CONFLICT_rows",
 "both_GWAS_exact_variant_rows","both_GWAS_exact_nonpalindromic_rows",
 "eligible_for_coloc","blocked_on",
}
REQUIRED_READY={
 "region_id","QTL_dataset","gene_id","possible_single_signal_ABF_followup",
 "can_now_assert_source_genome_GRCh37_REF_after_liftover",
 "can_now_assert_OFFICIAL_ALT_EFFECT_DIRECTION","QTL_source_ancestry_verified",
 "molecular_QTL_matched_LD_verified",
}

def flag(v):
    k=str(v or "").strip().lower()
    if k in TRUE:return True
    if k in FALSE:return False
    return None

def nonnegative(v):
    try:
        f=float(v)
        if not f.is_integer() or f<0:return None
        return int(f)
    except (ValueError,TypeError,OverflowError):
        return None

def valid_p(value):
    try:
        p=float(value)
        return math.isfinite(p) and 0<=p<=1
    except (ValueError,TypeError):return False

def valid_effect(value):
    try:return math.isfinite(float(value))
    except (ValueError,TypeError):return False

def valid_se(value):
    try:
        v=float(value)
        return math.isfinite(v) and v>0
    except (ValueError,TypeError):return False

def rows(path:Path,required:set[str]):
    if not path.is_file():raise FileNotFoundError(str(path))
    op=gzip.open if path.suffix==".gz" else open
    with op(path,"rt",encoding="utf-8-sig",errors="strict",newline="") as stream:
        reader=csv.DictReader(stream,delimiter="\t")
        if not reader.fieldnames or len(reader.fieldnames)!=len(set(reader.fieldnames)):
            raise ValueError("missing/duplicate headers: "+path.name)
        if not required.issubset(reader.fieldnames):
            raise ValueError("missing columns: "+str(sorted(required-set(reader.fieldnames))))
        for i,row in enumerate(reader):
            if i>=MAX_ROWS:raise ValueError("row cap exceeded")
            if None in row or None in row.values():raise ValueError("malformed TSV row")
            yield row

def audit(audits:Path):
    src=audits/ARCHIVED_RELATIVE
    ready=audits/READY_RELATIVE
    overlap=audits/OVERLAP_RELATIVE
    qc=collections.Counter()
    groups=collections.defaultdict(collections.Counter)
    unique_variants=set()
    source_joint_variants=set()
    source_joint_gene_ids=set()
    for row in rows(src,REQUIRED_ARCHIVE):
        qc["archived_rows"]+=1
        key=(row["region_id"],row["QTL_dataset"],row["gene_id"])
        groups[key]["archived_rows"]+=1
        g=groups[key]
        variant=row["mapped_genomic_variant_GRCh37"].strip()
        if variant:
            unique_variants.add(variant)
        b=flag(row["BBJ_allele_exact_Fasta37"])
        e=flag(row["EAS_AIS_allele_exact_Fasta37"])
        if b is None or e is None:
            qc["unknown_fasta_flag_encoding"]+=1
            continue
        if b:qc["bbj_archived_fasta_flag_rows"]+=1
        if e:qc["eas_archived_fasta_flag_rows"]+=1
        if b and e:
            g["both_fasta_source_flags"]+=1
            if (valid_p(row["QTL_p"]) and valid_p(row["BBJ_p"]) and
                valid_p(row["EAS_AIS_p"]) and valid_effect(row["QTL_beta_GRCh37_ALT"]) and
                valid_effect(row["BBJ_beta_ALT"]) and valid_effect(row["EAS_AIS_beta_ALT"]) and
                valid_se(row["QTL_se"])):
                g["source_flags_and_valid_qtl_se_effect_p"]+=1
            else:
                g["source_flags_missing_or_invalid_statistics"]+=1
            qc["both_archived_fasta_flag_rows"]+=1
            groups[key]["both_archived_fasta_flag_rows"]+=1
            if variant:source_joint_variants.add(variant)
            source_joint_gene_ids.add(row["gene_id"])
    qc["archived_gene_dataset_groups"]=len(groups)
    qc["archived_unique_variant_strings"]=len(unique_variants)
    qc["both_fasta_flag_unique_variant_strings"]=len(source_joint_variants)
    qc["both_fasta_flag_unique_genes"]=len(source_joint_gene_ids)
    if not qc["archived_rows"]:raise ValueError("empty archived overlap")
    overlap_counter=collections.Counter()
    overlap_by_group={}
    for row in rows(overlap,REQUIRED_OVERLAP):
        key=(row["region_id"],row["QTL_dataset"],row["gene_id"])
        if key in overlap_by_group:raise ValueError("duplicate overlap group")
        verified=nonnegative(row["GRCh37_FASTA_VERIFIED_SNP_rows"])
        conflict=nonnegative(row["FASTA_REF_CONFLICT_rows"])
        both=nonnegative(row["both_GWAS_exact_variant_rows"])
        nonpal=nonnegative(row["both_GWAS_exact_nonpalindromic_rows"])
        if any(x is None for x in (verified,conflict,both,nonpal)):
            raise ValueError("invalid nonnegative overlap count")
        if nonpal>both:raise ValueError("nonpal exceeds both overlap rows")
        overlap_counter["groups"]+=1
        overlap_counter["reported_fasta_verified_rows_sum"]+=verified
        overlap_counter["reported_fasta_conflict_rows_sum"]+=conflict
        overlap_counter["both_exact_rows_sum"]+=both
        overlap_counter["both_nonpal_rows_sum"]+=nonpal
        eligible=flag(row["eligible_for_coloc"])
        overlap_counter["eligible_flag_true" if eligible else "eligible_flag_not_true"]+=1
        overlap_by_group[key]={"fasta_verified":verified,"fasta_conflicts":conflict,
                               "both_exact":both,"both_nonpal":nonpal,
                               "blocked_on":row["blocked_on"][:150]}
    review=collections.Counter()
    flags=collections.Counter()
    source_candidates=[]
    for row in rows(ready,REQUIRED_READY):
        key=(row["region_id"],row["QTL_dataset"],row["gene_id"])
        review["groups"]+=1
        qtl_ld=flag(row["molecular_QTL_matched_LD_verified"])
        qtl_ancestry=flag(row["QTL_source_ancestry_verified"])
        abf=flag(row["possible_single_signal_ABF_followup"])
        reference=flag(row["can_now_assert_source_genome_GRCh37_REF_after_liftover"])
        direction=flag(row["can_now_assert_OFFICIAL_ALT_EFFECT_DIRECTION"])
        for name,value in (("abf_followup",abf),("reference_asserted",reference),
                           ("direction_asserted",direction),("qtl_ld_verified",qtl_ld),
                           ("qtl_ancestry_verified",qtl_ancestry)):
            flags[name+("_true" if value is True else "_false" if value is False else "_unknown")]+=1
        if abf:
            source_candidates.append({
                "region_id":key[0],"QTL_dataset":key[1],"gene_id":key[2],
                "reference_asserted":reference,"effect_direction_asserted":direction,
                "matched_qtl_ld_verified":qtl_ld,
                "ancestry_verified":qtl_ancestry,
                "archive_fasta_flag_both_rows":groups[key]["both_archived_fasta_flag_rows"],
                "overlap_reported_fasta_rows":overlap_by_group.get(key,{}).get("fasta_verified",0),
                "fullcis_archived_rows":groups[key]["archived_rows"],
                "source_flagged_stat_rows":groups[key]["source_flags_and_valid_qtl_se_effect_p"],
                "source_flagged_invalid_stat_rows":groups[key]["source_flags_missing_or_invalid_statistics"],
                "p3a_status":"EXPLORATORY_INPUT_SCREEN_ONLY",
                "still_required":["complete_region_variant_coverage","independent_ref_alt_harmonization",
                    "QTL_and_GWAS_sample_sizes","stroke_case_fraction","trait_type",
                    "BBJ_EAS_effect_standard_errors_or_valid_alternate_inputs",
                    "prior_sensitivity","single_causal_signal_assumption"],
                "scientific_status":"ABF_FOLLOWUP_FLAG_ONLY_NOT_COLOC_VALIDATION"
            })
    return {
       "schema":"is-archived-fasta-abf-audit-v1","integrity":"PASS",
       "archived":dict(qc),"overlap_report":dict(overlap_counter),
       "readiness":{**dict(review),**dict(flags),"abf_candidate_rows":len(source_candidates)},
       "abf_flagged_groups":source_candidates,
       "abf_screen":{"groups_with_archived_rows":sum(x["fullcis_archived_rows"]>0 for x in source_candidates),
          "groups_with_source_flagged_stats":sum(x["source_flagged_stat_rows"]>0 for x in source_candidates),
          "groups_without_source_flagged_stats":sum(x["source_flagged_stat_rows"]==0 for x in source_candidates),
          "screened":len(source_candidates),"coloc_abf_executed":0},
       "limitations":[
         "Archived FASTA flags and reported SNP overlap counts are not newly revalidated against a reference FASTA.",
         "Rows and unique variant strings have different denominators.",
         "A possible ABF followup flag is neither a completed coloc.abf nor SuSiE result.",
         "Single-signal coloc.abf does not inherently require an LD matrix, but complete compatible summary stats, allele and trait metadata, prior sensitivity and scientific review are prerequisites.",
       ]
    }
