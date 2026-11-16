#!/usr/bin/env python3
"""EAS ImmuNexUT genuine first-gene nominal p/beta to signed-Z methods pilot.

Reads actual original E-GEAD-420 ZIP64 member by <=2MiB byte-range requests,
decodes first 835 source rows + start of next gene, and outputs no patient IDs.
Signed z=sign(slope(ALT))*Phi^{-1}(1-p/2); a two-sided asymptotic assumption
which MUST NOT be described as publisher-released original Z/SE.

WASH7P is source-format sentinel, NOT prioritized ischemic-stroke locus.
NO GWAS join, eQTL effective sample N, LD R, eCAVIAR, coloc, or clinical claim.
"""
from __future__ import annotations
import argparse
import csv
from datetime import datetime,timezone
import json
import math
import os
from pathlib import Path
import re
import zlib

from scipy.stats import norm
from audit_EGEAD420_nominal_zip64 import (
    ZIP,FILELIST,IDF,MAX_TAIL,MAX_PREFIX,read_range,
    zip64_eocd,zip_directory_members,source_member_prefix,
    source_filelist,get_small,IMMUNEXUT_13_COL_HEADER,
    first_gene_full_cis_block_probe,
)


def first_gene_source_rows(tar_gzip_prefix):
    if not tar_gzip_prefix.startswith(b"\x1f\x8b"):
        raise ValueError("Actual nested gene association TAR must be gzip")
    inflated=zlib.decompressobj(16+zlib.MAX_WBITS).decompress(
        tar_gzip_prefix,2*MAX_PREFIX)
    if len(inflated)<2048:
        raise ValueError("Unusable original nominal TAR prefix")
    tar_header=inflated[:512]
    member=tar_header[:100].split(b"\0",1)[0].decode("utf-8","strict")
    if member!="CD16p_Mono_nominal.txt":
        raise ValueError("Unexpected first EAS immune-cell nominal source")
    all_text=inflated[512:].decode("utf-8","strict").split("\n")
    complete=all_text[:-1]
    if not complete or tuple(complete[0].split())!=IMMUNEXUT_13_COL_HEADER:
        raise ValueError("Source 13-column original header mismatch")
    values=[]
    first_gene=None
    seen_next=False
    num_reported=None
    for line in complete[1:]:
        if not line.strip():continue
        fields=line.split()
        if len(fields)!=13:raise ValueError("Actual source row not 13 fields")
        gene=fields[0]
        if first_gene is None:
            first_gene=gene
            num_reported=int(fields[4])
        if gene!=first_gene:
            seen_next=True
            break
        if int(fields[4])!=num_reported:
            raise ValueError("Different cis tested variant count within first gene")
        values.append(fields)
    if not seen_next or len(values)!=num_reported:
        raise ValueError("First original gene block incomplete or count differs from source")
    return values,{"source_gene":first_gene,
                   "source_declared_cis_variants":num_reported,
                   "source_next_gene_boundary_observed":seen_next}


def signed_z_from_slope_p(beta,p):
    """Do not fabricate signed Z for p=0, invalid model or 0 slope."""
    if not math.isfinite(beta) or not math.isfinite(p) or not 0<p<=1:
        raise ValueError("Original slope / p not usable for finite signed Z")
    if p<1 and beta==0:
        raise ValueError("Original nominal association has significant p with zero beta")
    z=float(norm.isf(p/2))
    if not math.isfinite(z):
        raise ValueError("Normal inverse survival not finite")
    return math.copysign(z,beta) if beta!=0 else 0.0


def convert_block(rows):
    """Convert only original rows with unambiguous sequence alleles and finite p.

    Count all source rows; do not fabricate a sequence allele for non-ACGT,
    symbolic, dash or missing markers. These are SOURCE-format exclusions,
    not failure of the original 835-row nominal association completeness.
    """
    from collections import Counter
    output=[]
    seen=set()
    reasons=Counter()
    raw_nonsig=0
    class_counts=Counter()
    for row in rows:
        if len(row)!=13:
            raise ValueError("Original row not 13 source fields")
        _,gene_symbol,chrom,tss,n_cis,rsid,ref,alt,vchr,beg,end,pval,slope=row
        if chrom!=vchr or not chrom.startswith("chr"):
            raise ValueError("Original variant chromosome inconsistent")
        try:
            beg=int(beg);end=int(end);reported=int(n_cis)
        except ValueError as exc:
            raise ValueError("Original gene cis count or coordinate not numeric") from exc
        if beg<1 or end<beg or reported!=len(rows):
            raise ValueError("Original gene cis count or genomic coordinates invalid")
        key=f"{vchr}:{beg}-{end}:{ref}>{alt}"
        if key in seen:
            raise ValueError("Duplicated original variant coordinate and allele key")
        seen.add(key)
        try:
            p=float(pval);beta=float(slope)
        except ValueError:
            reasons["non_numeric_p_or_slope"]+=1
            continue
        if not math.isfinite(p) or not 0<=p<=1 or not math.isfinite(beta):
            reasons["non_finite_or_invalid_p_slope"]+=1
            continue
        raw_nonsig+=(p>0.05)
        if not ref or not alt or ref==alt:
            reasons["missing_or_identical_REF_ALT"]+=1
            continue
        if not (set(ref)<=set("ACGT") and set(alt)<=set("ACGT")):
            reasons["noncanonical_OR_symbolic_REF_ALT"]+=1
            continue
        try:
            signed_z=signed_z_from_slope_p(beta,p)
        except ValueError:
            reasons["no_finite_signed_z_from_p_and_slope"]+=1
            continue
        if len(ref)==len(alt)==1:
            variant_class="SNV"
        elif len(ref)!=len(alt):
            variant_class="INDEL_UNNORMALIZED"
        else:
            variant_class="MNV_UNNORMALIZED"
        class_counts[variant_class]+=1
        output.append({
            "original_grch38_variant_alleles":key,
            "original_gene_symbol":gene_symbol,
            "original_nominal_p":p,
            "original_ALT_effect_slope":beta,
            "derived_signed_z_two_sided_normal_assumption":signed_z,
            "variant_class_without_ref_genome_normalization":variant_class,
            "downstream_GWAS_effect_orientation_ready":False,
        })
    if not output:
        raise ValueError("No usable original nominal rows for signed-Z source sentinel")
    excluded=sum(reasons.values())
    if len(output)+excluded!=len(rows) or len(seen)!=len(rows):
        raise ValueError("Source row conservation failure during signed-Z classification")
    return output,{
        "real_original_gene_cis_rows":len(rows),
        "source_rows_in_pilot_eligible_signed_z":len(output),
        "source_rows_excluded_from_signed_z":excluded,
        "source_exclusion_reasons":dict(reasons),
        "source_variant_classes_used":dict(class_counts),
        "nonsignificant_p_gt_0_05":raw_nonsig,
        "indel_records":class_counts["INDEL_UNNORMALIZED"],
        "positive_signed_z":sum(x["derived_signed_z_two_sided_normal_assumption"]>0 for x in output),
        "negative_signed_z":sum(x["derived_signed_z_two_sided_normal_assumption"]<0 for x in output),
        "all_original_cis_records_accounted_for":True,
        "all_original_cis_records_output":excluded==0,
        "excluded_rows_are_not_implicitly_recoded":True,
        "original_publisher_z_scores":False,
        "two_sided_normal_p_to_z_assumption":True,
        "effect_allele":"ALT",
        "has_eQTL_study_MAF":False,
        "has_eQTL_study_effective_N":False,
        "LD_correlation_matrix_available":False,
        "GWAS_harmonized":False,
        "eCAVIAR_CLPP_calculated":False
    }

def read_original_first_gene_prefix():
    manifest=source_filelist(get_small(FILELIST,6000))
    tail,info=read_range(ZIP,f"-{MAX_TAIL}",MAX_TAIL)
    cd=zip64_eocd(tail,info)
    st,stop=cd["central_directory_start"],cd["central_directory_end_inclusive"]
    if st>=info["start"] and stop<=info["end"]:
        directory=tail[st-info["start"]:st-info["start"]+cd["central_directory_byte_count"]]
    else:
        directory,_=read_range(ZIP,f"{st}-{stop}",2*MAX_PREFIX)
    members=zip_directory_members(directory,cd["members_expected"])
    file=[m for m in members if m["name"].endswith("nominal_eQTL_2.tar.gz")]
    if len(file)!=1 or file[0]["uncompressed_size"]!=manifest["nominal_eQTL_2.tar.gz"]["size_bytes"]:
        raise ValueError("Original nominal TAR publisher filelist mismatch")
    prefix,download=source_member_prefix(file[0])
    rows,counts=first_gene_source_rows(prefix)
    return rows,counts,download


def run(out):
    out.mkdir(parents=True,exist_ok=True)
    result={
        "created_utc":datetime.now(timezone.utc).isoformat(),
        "source":"Original ImmuNexUT E-GEAD-420 CD16p_Mono WASH7P source-method pilot",
        "status":"SOURCE_OR_QC_HOLD",
        "stroke_locus_prioritization":False,
        "GIGASTROKE_alignment":"NOT_RUN",
        "eCAVIAR_CLPP":"NOT_RUN",
        "whole_source_publisher_MD5":"NOT_VERIFIED",
        "eQTL_variants_are_GRCh38":True,
    }
    try:
        records,counts,prefix=read_original_first_gene_prefix()
        output,stats=convert_block(records)
        result.update({"source_first_gene_block":counts,
                       "source_range_evidence":prefix,
                       "derived_z_method_evidence":stats,
                       "status":"REAL_EAS_EQTL_GENE_SIGNED_Z_METHOD_SENTINEL_PASS_NO_LD"})
        table=out/"CD16p_Mono_WASH7P_source_derived_signed_z_SENTINEL.tsv"
        with table.open("w",encoding="utf-8",newline="") as f:
            writer=csv.DictWriter(f,fieldnames=list(output[0]),delimiter="\t")
            writer.writeheader()
            writer.writerows(output)
    except (OSError,ValueError,TimeoutError,KeyError,UnicodeDecodeError,zlib.error) as exc:
        result["source_issue"]=f"{type(exc).__name__}: {str(exc)[:220]}"
    (out/"z_method_go_hold.json").write_text(
        json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    lines=["# ImmuNexUT first-gene original eQTL p, beta -> signed Z method feasibility",
           "", "**"+result["status"]+"**","",
           "This is **WASH7P a SOURCE FORMAT SENTINEL**, not an ischemic stroke gene discovery."]
    if details:=result.get("derived_z_method_evidence"):
        lines.extend([
            f"- Original first-gene records: {details['real_original_gene_cis_rows']}",
            f"- Nonsignificant original p-values: {details['nonsignificant_p_gt_0_05']}",
            f"- Indel records retained in source summary: {details['indel_records']}",
            f"- Original source rows eligible for signed Z: {details['source_rows_in_pilot_eligible_signed_z']}",
            f"- Source rows excluded without recoding: {details['source_rows_excluded_from_signed_z']}",
            f"- Exclusion categories: {details['source_exclusion_reasons']}",
            f"- Original source classes retained: {details['source_variant_classes_used']}",
            f"- Positive signed Z: {details['positive_signed_z']}; negative signed Z: {details['negative_signed_z']}",
            "- Z signs from original ALT slope; magnitudes from two-sided normal approximation of original nominal P.",
        ])
    if issue:=result.get("source_issue"):
        lines.append("Source or method QC HOLD: "+issue)
    lines.extend(["",
        "No original per-celltype eQTL sample N, MAF, study-genotype LD R, GWAS matching or eCAVIAR CLPP.",
        "GRCh38 source alleles include indels and multi-rsID; not normalized to VCF or lifted to GIGASTROKE GRCh37.",
        "Do not interpret the signed Z transform as a validated replacement for original eQTL beta SE or a biological replication."
    ])
    summary="\n".join(lines)+"\n"
    (out/"summary.md").write_text(summary,encoding="utf-8")
    print(summary)
    if p:=os.environ.get("GITHUB_STEP_SUMMARY"):
        with open(p,"a",encoding="utf-8") as f:f.write(summary)
    if result["status"]!="REAL_EAS_EQTL_GENE_SIGNED_Z_METHOD_SENTINEL_PASS_NO_LD":
        raise SystemExit(1)


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--out-dir",default="artifacts/is_ImmuNexUT_signed_Z_sentinel")
    a=parser.parse_args()
    run(Path(a.out_dir))
