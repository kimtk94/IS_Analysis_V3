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
import hashlib
import json
import re
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
# The registry is optional and server-local. Paths cannot be supplied from GitHub
# Issues or untrusted external comments.
GWAS_REGISTRY_RELATIVE="github_actions_runs/evidence/science_inputs/P3A_APPROVED_GWAS_SOURCES.tsv"
GWAS_ALLOWED_ROOT=Path("/srv/is-analysis")
GWAS_REGISTRY_COLUMNS={
    "source","accession","file_path","file_sha256","genome_build",
    "trait_type","n_cases","n_controls","reference_url",
    "independent_match_status","matched_variant_count",
    "reviewed_by","reviewed_at",
}
GWAS_SOURCES={"BBJ","EAS_AIS"}
HEX64=re.compile(r"^[a-f0-9]{64}$")
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

def audit_approved_gwas_registry(audits:Path, *, allowed_root:Path=GWAS_ALLOWED_ROOT):
    """Fail-closed external GWAS identity audit; never adopts sample sizes
    without independent matched-variant verification and actual file digest."""
    path=audits.parent/GWAS_REGISTRY_RELATIVE
    result={"schema":"is-p3a-gwas-source-identity-v1",
            "registry_status":"MISSING_APPROVED_REGISTRY",
            "accepted_sources":[],
            "blocked_sources":{},
            "required_sources":sorted(GWAS_SOURCES),
            "scientific_status":"NO_GWAS_METADATA_ADOPTED"}
    if not path.is_file():
        return result
    if path.is_symlink() or path.stat().st_size>1_000_000:
        raise ValueError("Invalid approved source registry")
    with path.open(encoding="utf-8-sig",newline="") as stream:
        reader=csv.DictReader(stream,delimiter="\t")
        if not reader.fieldnames or not GWAS_REGISTRY_COLUMNS.issubset(reader.fieldnames):
            raise ValueError("Approved GWAS registry required fields missing")
        records=list(reader)
    if len(records)>20 or any(None in r or None in r.values() for r in records):
        raise ValueError("Malformed GWAS source registry")
    seen=set()
    result["registry_status"]="REVIEW_REQUIRED"
    for record in records:
        source=record["source"]
        if source not in GWAS_SOURCES or source in seen:
            raise ValueError("Duplicate or unknown GWAS source")
        seen.add(source)
        reasons=[]
        if not record["accession"].strip() or not record["reference_url"].startswith("https://"):
            reasons.append("SOURCE_ACCESSION_PROVENANCE_MISSING")
        if record["genome_build"]!="GRCh37" or record["trait_type"]!="cc":
            reasons.append("UNVERIFIED_BUILD_OR_TRAIT_TYPE")
        try:
            cases=int(record["n_cases"])
            controls=int(record["n_controls"])
            if cases<1 or controls<1:
                raise ValueError()
        except ValueError:
            reasons.append("MISSING_VALID_CASE_CONTROL_COUNTS")
            cases=controls=0
        if record["independent_match_status"]!="VERIFIED" or not record["reviewed_by"] or not record["reviewed_at"]:
            reasons.append("INDEPENDENT_SOURCE_MATCH_NOT_APPROVED")
        try:
            match_count=int(record["matched_variant_count"])
            if match_count<20:
                reasons.append("INSUFFICIENT_INDEPENDENT_VARIANT_MATCHES")
        except ValueError:
            reasons.append("MISSING_MATCH_COUNT")
        source_path=Path(record["file_path"])
        safe=False
        try:
            if source_path.is_absolute() and not source_path.is_symlink() and source_path.is_file():
                safe=source_path.resolve().is_relative_to(allowed_root.resolve())
        except OSError:
            safe=False
        if not safe or not HEX64.fullmatch(record["file_sha256"]):
            reasons.append("MISSING_SAFE_FILE_AND_SHA256")
        else:
            h=hashlib.sha256()
            with source_path.open("rb") as handle:
                for chunk in iter(lambda:handle.read(1024*1024),b""):
                    h.update(chunk)
            if h.hexdigest()!=record["file_sha256"]:
                reasons.append("GWAS_SOURCE_FILE_DIGEST_MISMATCH")
        if reasons:
            result["blocked_sources"][source]=sorted(set(reasons))
            continue
        # A reviewed register is still not a replacement for retained per-variant
        # comparison evidence. Do NOT treat the registry alone as source proof.
        result["blocked_sources"][source]=["MATCH_EVIDENCE_NOT_REPRODUCED_BY_THIS_AUDIT"]
    for missing in sorted(GWAS_SOURCES-seen):
        result["blocked_sources"][missing]=["NO_SOURCE_RECORD"]
    return result


def source_header_inventory(path:Path):
    """Report presence of declared source fields without claiming metadata validity."""
    opener=gzip.open if path.suffix==".gz" else open
    with opener(path,"rt",encoding="utf-8-sig",newline="") as stream:
        fields=next(csv.reader([stream.readline()],delimiter="\t"))
    present=set(fields)
    expected={
        "qtl_beta_and_se":["QTL_beta_GRCh37_ALT","QTL_se"],
        "bbj_beta_and_p":["BBJ_beta_ALT","BBJ_p"],
        "eas_beta_and_p":["EAS_AIS_beta_ALT","EAS_AIS_p"],
        "bbj_effect_se":["BBJ_se","BBJ_SE","BBJ_beta_se"],
        "eas_effect_se":["EAS_AIS_se","EAS_AIS_SE","EAS_AIS_beta_se"],
        "bbj_sample_size":["BBJ_N","BBJ_n","BBJ_sample_size"],
        "eas_sample_size":["EAS_AIS_N","EAS_AIS_n","EAS_AIS_sample_size"],
        "stroke_case_fraction":["BBJ_case_fraction","EAS_AIS_case_fraction","case_fraction"],
        "genomic_variant":["mapped_genomic_variant_GRCh37"],
    }
    return {"source":"ARCHIVED_FASTA_QTL_GWAS_FULLCIS",
            "field_presence":{k:all(x in present for x in v) if k.endswith("_and_se") or k.endswith("_and_p") else any(x in present for x in v) for k,v in expected.items()},
            "field_count":len(fields),
            "qualification":"Header audit only; missing sample metadata must not be inferred."}


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
    header_audit=source_header_inventory(src)
    gwas_source=audit_approved_gwas_registry(audits)
    qc=collections.Counter()
    groups=collections.defaultdict(collections.Counter)
    unique_variants=set()
    source_joint_variants=set()
    source_joint_gene_ids=set()
    group_unique=collections.defaultdict(set)
    group_dual_unique=collections.defaultdict(set)
    for row in rows(src,REQUIRED_ARCHIVE):
        qc["archived_rows"]+=1
        key=(row["region_id"],row["QTL_dataset"],row["gene_id"])
        groups[key]["archived_rows"]+=1
        g=groups[key]
        variant=row["mapped_genomic_variant_GRCh37"].strip()
        if variant:
            unique_variants.add(variant)
            group_unique[key].add(variant)
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
            if variant:
                source_joint_variants.add(variant)
                group_dual_unique[key].add(variant)
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
            total=groups[key]["archived_rows"]
            valid=groups[key]["source_flags_and_valid_qtl_se_effect_p"]
            unique=len(group_unique[key])
            dual_unique=len(group_dual_unique[key])
            reasons=[]
            if total==0:reasons.append("NO_FULLCIS_ARCHIVE_ROWS")
            if valid==0:reasons.append("NO_DUAL_FLAGGED_VALID_SUMMARY_STATS")
            if reference is not True:reasons.append("SOURCE_REF_ASSERTION_UNVERIFIED")
            if direction is not True:reasons.append("SOURCE_EFFECT_DIRECTION_UNVERIFIED")
            # These are not inferred from nominal P-values or prior pipeline flags.
            reasons += ["INDEPENDENT_ALLELE_REFERENCE_PENDING",
                        "COMPLETE_CIS_WINDOW_COVERAGE_UNVERIFIED",
                        "QTL_GWAS_SAMPLE_SIZES_UNVERIFIED",
                        "CASE_CONTROL_FRACTION_UNVERIFIED",
                        "GWAS_EFFECT_VARIANCE_UNVERIFIED",
                        "TRAIT_TYPE_AND_PRIORS_UNVERIFIED",
                        "SINGLE_SIGNAL_ASSUMPTION_UNVERIFIED"]
            source_candidates.append({
                "region_id":key[0],"QTL_dataset":key[1],"gene_id":key[2],
                "reference_asserted":reference,"effect_direction_asserted":direction,
                "matched_qtl_ld_verified":qtl_ld,
                "ancestry_verified":qtl_ancestry,
                "archive_fasta_flag_both_rows":groups[key]["both_archived_fasta_flag_rows"],
                "overlap_reported_fasta_rows":overlap_by_group.get(key,{}).get("fasta_verified",0),
                "fullcis_archived_rows":total,
                "fullcis_unique_variant_strings":unique,
                "dual_fasta_source_unique_variant_strings":dual_unique,
                "source_flagged_stat_rows":valid,
                "source_flagged_invalid_stat_rows":groups[key]["source_flags_missing_or_invalid_statistics"],
                "source_stat_row_fraction":round(valid/total,6) if total else None,
                "source_stat_screen":"HAS_VALID_SOURCE_STATS" if valid else "NO_VALID_SOURCE_STATS",
                "coloc_abf_eligible":False,
                "p3a_status":"EXPLORATORY_INPUT_SCREEN_ONLY",
                "blocking_reasons":reasons,
                "scientific_status":"ABF_FOLLOWUP_FLAG_ONLY_NOT_COLOC_VALIDATION"
            })
    return {
       "schema":"is-archived-fasta-abf-audit-v1","integrity":"PASS",
       "p3a_source_header_audit":header_audit,
       "p3a_gwas_source_identity":gwas_source,
       "archived":dict(qc),"overlap_report":dict(overlap_counter),
       "readiness":{**dict(review),**dict(flags),"abf_candidate_rows":len(source_candidates)},
       "abf_flagged_groups":source_candidates,
       "abf_gate_summary":{
           "eligible_for_coloc_abf":0,
           "groups_with_source_stat_rows":sum(x["source_flagged_stat_rows"]>0 for x in source_candidates),
           "groups_with_no_source_stat_rows":sum(x["source_flagged_stat_rows"]==0 for x in source_candidates),
           "source_stat_rows_total":sum(x["source_flagged_stat_rows"] for x in source_candidates),
           "source_stat_unique_variants_total_over_groups":sum(x["dual_fasta_source_unique_variant_strings"] for x in source_candidates),
           "pending_independent_allele_and_sample_metadata":len(source_candidates),
       },
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
