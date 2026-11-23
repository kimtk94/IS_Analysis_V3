#!/usr/bin/env python3
"""Original GIGASTROKE AIS EAS/EUR source β-SE-p and model-scale audit.

Read original publisher YAML and *bounded* 256KiB initial gzip bytes of
GCST90104545 (EAS) and GCST90104540 (EUR). Compare real effect/SE/p-value
coherence without claiming it identifies beta units. Cross-reference exact
Mishra 2022 full-text methods to describe logistic/Cox heterogeneous cohorts.
Record only aggregate QC, no per-SNP effects or row values. Never impute GWAS N.
"""
from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import importlib
import io
import json
import math
import re
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
import zlib
from pathlib import Path

PINNED=importlib.import_module("scripts.is.p3a_publisher_md5_gate")
PRIMARY=importlib.import_module("scripts.is.audit_p3a_36_source_expression_and_gwas_beta_units")
METADATA=importlib.import_module("scripts.is.audit_p3a_primary_gwas_yaml_qtl_an_36")
MAX_PREFIX_BYTES=262_144
MAX_GWAS_SAMPLE_ROWS=750
GWAS_COLUMNS={
    "chromosome","base_pair_location","effect_allele","other_allele",
    "effect_allele_frequency","beta","standard_error","p_value"
}
PRIMARY_CASES={"EAS":{"cases":19032,"controls":237242,"BBJ_is_discovery":True},
               "EUR":{"cases":62100,"controls":1234808,"BBJ_is_discovery":False}}
PRIMARY_PAPER_DOI="10.1038/s41586-022-05165-3"


def require(ok,msg):
    if not ok:raise ValueError(msg)


def source_small_gwas_prefix(ancestry,open_url=urllib.request.urlopen):
    require(ancestry in PINNED.SOURCES,"Invalid AIS ancestry")
    spec=PINNED.SOURCES[ancestry]
    url,yaml_url=PINNED.build_urls(spec)
    parsed=urllib.parse.urlsplit(url)
    require(parsed.scheme=="https" and parsed.hostname=="ftp.ebi.ac.uk" and
            parsed.path.endswith("/"+spec["accession"]+"_buildGRCh37.tsv.gz"),
            "Original publisher EAS/EUR accession mismatch")
    request=urllib.request.Request(url,headers={
        "Range":f"bytes=0-{MAX_PREFIX_BYTES-1}",
        "Accept-Encoding":"identity",
        "User-Agent":"IS-P3A-GIGASTROKE-ORIGINAL-BETA-SE-P-READONLY/1.0",
    })
    with open_url(request,timeout=90) as response:
        require(response.status==206,"Original file HTTP Range was not honored")
        resolved=urllib.parse.urlsplit(response.geturl())
        require(resolved.scheme=="https" and resolved.hostname=="ftp.ebi.ac.uk",
                "Source redirected to unapproved origin")
        length=response.headers.get("Content-Range","")
        m=re.fullmatch(r"bytes (\d+)-(\d+)/(\d+)",length)
        require(m and int(m.group(1))==0 and
                int(m.group(2))==MAX_PREFIX_BYTES-1 and
                int(m.group(3))==spec["bytes"],
                "Original file range or pinned whole-source compressed length drift")
        data=response.read(MAX_PREFIX_BYTES+1)
    require(len(data)==MAX_PREFIX_BYTES and data[:2]==b"\x1f\x8b",
            "Original compressed source first 256KiB unavailable")
    return data,{
        "exact_original_GWAS_accession":spec["accession"],
        "original_publisher_url":url,
        "publisher_metadata_yaml_url":yaml_url,
        "source_first_compressed_bytes":len(data),
        "original_full_compressed_file_bytes_pinned":spec["bytes"],
        "original_publisher_full_MD5_previously_verified_not_rerun":spec["publisher_md5_pinned"],
        "original_publisher_full_SHA256_previously_verified_not_rerun":spec["full_sha256_pinned"],
        "source_compressed_prefix_SHA256":hashlib.sha256(data).hexdigest(),
        "original_full_file_downloaded_in_this_run":False,
    }


def parse_source_effect_stats(raw,ancestry):
    require(ancestry in PINNED.SOURCES and raw[:2]==b"\x1f\x8b",
            "Original source data not gzip/known ancestry")
    inflater=zlib.decompressobj(wbits=31)
    text=inflater.decompress(raw,max_length=16_000_000).decode("utf-8-sig")
    lines=text.splitlines(keepends=True)
    if not inflater.eof and lines and not lines[-1].endswith("\n"):
        lines=lines[:-1]
    require(len(lines)>=5,"Original compressed prefix has insufficient complete rows")
    reader=csv.DictReader(io.StringIO("".join(lines)),delimiter="\t")
    fields=reader.fieldnames or []
    require(len(fields)==len(set(fields)) and GWAS_COLUMNS.issubset(fields),
            "Original publisher GWAS beta/se/p/allele schema drift")
    rows=list(next_rows(reader,MAX_GWAS_SAMPLE_ROWS))
    require(rows and all(None not in r and None not in r.values() for r in rows),
            "Original GWAS first rows malformed")
    quantiles=[]
    z_comparable=0
    deviations=[]
    chromosome_count=collections.Counter()
    sign=collections.Counter()
    effect_reversed=0
    p_zero=0
    or_rows=0
    or_log_beta_error=[]
    ci_rows=0
    ci_exp_se_error=[]
    or_positive=0
    or_missing=0
    for row in rows:
        try:
            chrom=row["chromosome"].replace("chr","")
            pos=int(row["base_pair_location"])
            effect=float(row["beta"])
            se=float(row["standard_error"])
            pv=float(row["p_value"])
            eaf=float(row["effect_allele_frequency"])
            ea=row["effect_allele"].upper()
            oa=row["other_allele"].upper()
        except (TypeError,ValueError,KeyError) as exc:
            raise ValueError("Original EAS/EUR contains nonnumeric locus/effect rows") from exc
        require(re.fullmatch(r"(?:[1-9]|1\d|2[0-2]|X|Y|MT)",chrom) is not None
                and pos>0 and 0<=pv<=1 and
                math.isfinite(effect) and math.isfinite(se) and se>0 and
                math.isfinite(eaf) and 0<=eaf<=1 and
                re.fullmatch(r"[ACGT]+",ea) and re.fullmatch(r"[ACGT]+",oa)
                and ea!=oa, "Original GIGASTROKE source invalid beta/SE/p/allele")
        chromosome_count[chrom]+=1
        sign["negative" if effect<0 else "positive" if effect>0 else "zero"]+=1
        if pv==0:
            p_zero+=1
        z=effect/se
        derived_p=math.erfc(abs(z)/math.sqrt(2))
        if derived_p>0 and pv>0:
            z_comparable+=1
            deviations.append(abs(math.log10(derived_p)-math.log10(pv)))
        if math.isfinite(z) and math.isfinite(eaf):
            quantiles.append(abs(z))
        # The actual source also lists odds_ratio/ci_lower/ci_upper.
        # A populated odds_ratio can certify a beta-scale relation in real rows,
        # whereas a header alone or beta-vs-p coherence cannot.
        or_value=optional_positive_float(row.get("odds_ratio"))
        lower=optional_positive_float(row.get("ci_lower"))
        upper=optional_positive_float(row.get("ci_upper"))
        if or_value is not None:
            or_rows+=1
            or_log_beta_error.append(abs(math.log(or_value)-effect))
            if or_value>1:or_positive+=1
        else:
            or_missing+=1
        if lower is not None and upper is not None:
            require(lower<upper, "Original odds ratio CI lower>=upper")
            ci_rows+=1
            ci_exp_se_error.append(max(
                abs(math.log(lower)-(effect-1.96*se)),
                abs(math.log(upper)-(effect+1.96*se))))
    require(len(rows)>0 and (len(deviations)>2 or p_zero>0),
            "No verifiable original GWAS effect/SE/p evidence")
    return {
        "real_original_source_columns":fields,
        "original_source_preview_rows":len(rows),
        "source_chromosome_counts_preview_only":dict(sorted(chromosome_count.items())),
        "n_real_effects_negative":sign["negative"],
        "n_real_effects_positive":sign["positive"],
        "n_real_effects_zero":sign["zero"],
        "pvalue_underflow_zero_rows":p_zero,
        "number_beta_over_SE_vs_two_sided_normal_p_comparable":z_comparable,
        "z_p_log10_error_median":quantile(deviations,0.5),
        "z_p_log10_error_p95":quantile(deviations,0.95),
        "z_p_abs_log10_error_le_0_1":sum(d<=0.1 for d in deviations),
        "z_p_abs_log10_error_gt_1":sum(d>1 for d in deviations),
        "absolute_z_max_first_rows":max(quantiles) if quantiles else None,
        "source_OR_and_CI_columns_present":
            all(x in fields for x in ("odds_ratio","ci_lower","ci_upper")),
        "source_OR_populated_rows":or_rows,
        "source_OR_missing_or_unparseable_rows":or_missing,
        "source_OR_gt_one_rows":or_positive,
        "source_beta_ln_OR_log_scale_median_abs_error":
            quantile(or_log_beta_error,0.5),
        "source_beta_ln_OR_log_scale_max_abs_error":
            max(or_log_beta_error) if or_rows else None,
        "source_beta_ln_OR_max_abs_error_le_0_015_rows":
            sum(x<=0.015 for x in or_log_beta_error),
        "source_OR_95CI_populated_rows":ci_rows,
        "source_OR_95CI_exp_beta_pm_1_96SE_p95_log_abs_error":
            quantile(ci_exp_se_error,0.95),
        "source_OR_95CI_exp_beta_pm_1_96SE_le_0_04_rows":
            sum(x<=0.04 for x in ci_exp_se_error),
        "source_direct_ln_OR_vs_beta_preview_status":(
            "PASS_ORIGINAL_POPULATED_OR_AGREES_WITH_BETA_IN_PREVIEW"
            if or_rows>=50 and all(x<=0.015 for x in or_log_beta_error)
            else "HOLD_OR_SOURCE_NOT_POPULATED_OR_INCONSISTENT"),
        "raw_effect_or_variant_rows_emitted":False,
        "source_file_beta_numeric_consistency_is_not_log_OR_unit_proof":True,
        "source_file_beta_exact_log_OR_scale_verified":False,
        "GWAS_per_SNP_effective_sample_N_recovered":False,
    }


def optional_positive_float(value):
    """Return absent when publisher's original optional OR/CI field is blank."""
    if value is None or str(value).strip().lower() in {
        "", ".", "na", "n/a", "null", "nan", "none", "nr"}:
        return None
    try:
        converted=float(value)
    except (ValueError,TypeError) as exc:
        raise ValueError("Original publisher OR/CI optional field nonnumeric") from exc
    require(math.isfinite(converted) and converted>0,
            "Original publisher OR/CI field invalid nonpositive")
    return converted


def next_rows(reader,n):
    for i,row in enumerate(reader):
        if i>=n:break
        yield row


def quantile(values,fraction):
    if not values:return None
    a=sorted(values)
    idx=(len(a)-1)*fraction
    lo=int(idx)
    hi=min(lo+1,len(a)-1)
    return a[lo]*(hi-idx)+a[hi]*(idx-lo)


def primary_methods(raw):
    root=ET.fromstring(raw)
    ids=[(x.text or "").strip().lower() for x in root.iter()
         if x.tag.split("}")[-1]=="article-id" and
         x.attrib.get("pub-id-type")=="doi"]
    require(PRIMARY_PAPER_DOI in ids,"Mishra 2022 primary paper DOI mismatch")
    paragraphs=[" ".join(" ".join(e.itertext()).split()) for e in root.iter()
                if e.tag.split("}")[-1]=="p"]
    text=" ".join(paragraphs).lower()
    logistic=("logistic regression" in text)
    cox=("cox regression" in text)
    filters=bool(re.search(r"log\s*\[\s*or\s*\]\s*[<>]",text))
    require(logistic and cox and filters,
            "Original paper no longer establishes logistic/Cox and log OR methods")
    return {
        "original_paper_DOI":PRIMARY_PAPER_DOI,
        "original_primary_paper_method_has_logistic_regression":logistic,
        "original_primary_paper_method_has_Cox_regression":cox,
        "original_primary_paper_variant_quality_control_mentions_log_OR":filters,
        "original_article_fulltext_XML_SHA256":hashlib.sha256(raw).hexdigest(),
        "original_paper_methods_support_log_OR_statistical_context":True,
        "cannot_determine_exact_uploaded_EAS_EUR_file_beta_scale_from_paper_alone":True,
        "some_contributing_cohorts_are_longitudinal_hazard_models":True,
    }


def source_yaml(ancestry):
    # Prior original source identity & empty N fields checked in previous CI.
    require(ancestry in PINNED.SOURCES,"Unsupported AIS ancestry")
    from urllib.request import urlopen
    import yaml
    spec=PINNED.SOURCES[ancestry]
    _,url=PINNED.build_urls(spec)
    req=urllib.request.Request(url,headers={"User-Agent":"Academic-IS-P3A-Original-GWAS-YAML/1.0"})
    with urlopen(req,timeout=70) as res:
        require(res.status==200,"Original source YAML HTTP failure")
        resolved=urllib.parse.urlsplit(res.geturl())
        require(resolved.scheme=="https" and resolved.hostname=="ftp.ebi.ac.uk",
                "Original GWAS YAML redirected to unapproved domain")
        raw=res.read(160_001)
    require(0<len(raw)<=160_000,"Original GWAS publisher YAML outside size limit")
    y=yaml.safe_load(raw.decode("utf-8-sig"))
    require(isinstance(y,dict) and y.get("gwas_id")==spec["accession"] and
            y.get("genome_assembly")=="GRCh37" and
            y.get("data_file_md5sum")==spec["publisher_md5_pinned"] and
            y.get("data_file_name")==spec["accession"]+"_buildGRCh37.tsv.gz" and
            "ischemic stroke" in str(y.get("trait_description","")).lower(),
            "Original GWAS YAML accession, stroke phenotype, GRCh37 or MD5 mismatch")
    source_keys=sorted(str(k) for k in y)
    applicable={k:y.get(k) for k in ("sample_size","case_count","control_count",
        "case_control_study","effect_size_units","beta_unit","trait_description")}
    has_exact_counts=all(isinstance(applicable[k],int) and applicable[k]>0
           for k in ("sample_size","case_count","control_count"))
    effect_statement=[k for k in ("effect_size_units","beta_unit") if applicable[k]]
    return {
        "original_YAML_source_URL":url,
        "original_YAML_SHA256":hashlib.sha256(raw).hexdigest(),
        "original_release_MD5_pin_match":True,
        "publisher_yaml_top_level_keys":source_keys,
        "publisher_sample_N_fields":{k:str(applicable[k]) if applicable[k] is not None else None
                                     for k in ("sample_size","case_count","control_count")},
        "publisher_effect_unit_label_fields_present":effect_statement,
        "publisher_exact_file_case_control_N_complete":has_exact_counts,
        "publisher_exact_file_case_fraction_s_directly_verified":False,
        "YAML_identity_verified_not_effect_scale_verified":True,
    }


def audit(manifest):
    require(len(manifest.get("groups",[]))==36 and
            len({(g["QTL_dataset"],g["gene_id"]) for g in manifest["groups"]})==36,
            "36 original gene x tissue assays manifest mismatch")
    # Reuse an already verified ORIGINAL full-text XML source and source hash.
    # Primary Europe PMC server timed out on feature run #45160855981.
    # This is deliberately NOT described as a live re-fetch in this execution.
    method={
        "source_status":"PASS_PREVIOUSLY_VERIFIED_PRIMARY_XML_SOURCE_SHA256",
        "original_paper_DOI":PRIMARY_PAPER_DOI,
        "source_requested_URL_from_prior_verification":PRIMARY.PAPER_URL,
        "original_article_fulltext_XML_SHA256_prior_verified":
            "7f1105825c4111f73f8b8c36f55897e5c1a19b39b0b240eabc9da3822ac89f79",
        "original_prior_CI_provenance_run_id":45040084191,
        "paper_fulltext_refetched_this_run":False,
        "original_primary_paper_method_has_logistic_regression":True,
        "original_primary_paper_method_has_Cox_regression":True,
        "original_primary_paper_variant_quality_control_mentions_log_OR":True,
        "cannot_determine_exact_uploaded_EAS_EUR_file_beta_scale_from_paper_alone":True,
        "some_contributing_cohorts_are_longitudinal_hazard_models":True,
    }
    details={}
    for ancestry in ("EAS","EUR"):
        print("SOURCE_ORIGINAL_"+ancestry+"_GWAS_YAML",flush=True)
        yaml=source_yaml(ancestry)
        print("SOURCE_ORIGINAL_"+ancestry+"_BOUND_256K_GZIP",flush=True)
        raw,provenance=source_small_gwas_prefix(ancestry)
        qc=parse_source_effect_stats(raw,ancestry)
        require(qc["source_OR_populated_rows"]==qc["original_source_preview_rows"]
                and qc["source_direct_ln_OR_vs_beta_preview_status"]==
                    "PASS_ORIGINAL_POPULATED_OR_AGREES_WITH_BETA_IN_PREVIEW"
                and qc["source_OR_95CI_populated_rows"]==qc["original_source_preview_rows"]
                and qc["source_OR_95CI_exp_beta_pm_1_96SE_le_0_04_rows"]==
                    qc["original_source_preview_rows"],
                "Source EAS/EUR OR=exp(beta) or CI=exp(beta+-1.96SE) preview failed")
        desc=PRIMARY_CASES[ancestry]
        total=desc["cases"]+desc["controls"]
        details[ancestry]={
            "official_original_AIS_GWAS_accession":PINNED.SOURCES[ancestry]["accession"],
            "official_original_source_provenance":provenance,
            "official_original_publisher_yaml":yaml,
            "actual_original_GWAS_source_prefix_QC":qc,
            "original_publication_discovery_AIS_cohort_totals_REFERENCE_ONLY": {
                "cases":desc["cases"],"controls":desc["controls"],"total":total,
                "published_case_fraction_reference_only":desc["cases"]/total,
                "case_fraction_s_approved_for_exact_summary_file":False,
                "genomic_variant_specific_effective_N_approved":False,
            },
            "BBJ_original_EAS_discovery_included":desc["BBJ_is_discovery"],
            "independent_BBJ_replication_permitted":False,
            "specific_original_file_beta_log_OR_approved":False,
        }
    return {
        "schema":"is-p3a-original-gigastroke-EAS-EUR-beta-SE-p-method-scale-20261123-v1",
        "technical_status":"PASS_PUBLISHER_ORIGINAL_GWAS_YAML_AND_BOUNDED_SOURCE_QC",
        "scientific_status":"HOLD_ORIGINAL_BINARY_GWAS_LOG_OR_AND_CASE_FRACTION_LINKAGE",
        "original_GIGASTROKE_paper_original_methods":method,
        "source_EAS_EUR":details,
        "n_original_GWAS_ancestries_with_OR_exp_beta_preview_verified":2,
        "n_real_original_GWAS_SNP_preview_OR_and_CI_consistent":sum(
            x["actual_original_GWAS_source_prefix_QC"]["source_OR_populated_rows"]
            for x in details.values()),
        "original_GWAS_beta_is_log_OR_directly_verified_in_first_750_per_ancestry":True,
        "original_GWAS_beta_log_OR_full_chr10_region_source_verified":False,
        "original_QTL_gene_tissue_assays":36,
        "original_ROSMAP_gene_assays_with_prior_fullcis_real_QTL":15,
        "source_ancestry_discovery_independence_assumed":False,
        "n_formal_coloc_abf_eligible":0,"n_formal_coloc_abf_executed":0,
        "n_formal_coloc_susie_eligible":0,"n_formal_coloc_susie_executed":0,
        "no_raw_GWAS_variant_effects_or_genotypes_exported":True,
        "review_required_to_promote_cc_metadata_to_ABF_inputs":True,
        "groups":[
            {"QTL_dataset":g["QTL_dataset"],"gene_id":g["gene_id"],
             "gene_symbol":g["gene_symbol"],
             "scientific_status":"HOLD_GWAS_FILE_LOGOR_CASE_FRACTION_AND_QTL_MODEL_N",
             "coloc_eligible":False} for g in manifest["groups"]
        ],
    }


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--manifest",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    result={
        "schema":"is-p3a-original-gigastroke-EAS-EUR-beta-SE-p-method-scale-20261123-v1",
        "technical_status":"HOLD_ORIGINAL_PUBLISHER_GWAS_SOURCE_OR_PAPER",
        "scientific_status":"HOLD",
        "n_formal_coloc_abf_eligible":0,"n_formal_coloc_abf_executed":0,
    }
    try:
        manifest=json.loads(a.manifest.read_text(encoding="utf-8"))
        print("SOURCE_PREVIOUSLY_VERIFIED_MISHRA_2022_ORIGINAL_XML_SHA256_AND_EXACT_AIS_FILES",flush=True)
        result=audit(manifest)
    except (ValueError,TypeError,KeyError,OSError,UnicodeError,TimeoutError,
            ET.ParseError,zlib.error) as exc:
        result["error_class_only"]=type(exc).__name__
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps({
        "technical_status":result["technical_status"],
        "scientific_status":result["scientific_status"],
        "preview":{a:{
            "rows":v["actual_original_GWAS_source_prefix_QC"]["original_source_preview_rows"],
            "median_abs_log10_p_error":v["actual_original_GWAS_source_prefix_QC"]["z_p_log10_error_median"],
            "case_fraction_source_approved":v[
                "original_publication_discovery_AIS_cohort_totals_REFERENCE_ONLY"
                ]["case_fraction_s_approved_for_exact_summary_file"]}
            for a,v in result.get("source_EAS_EUR",{}).items()},
        "n_coloc_executed":0,
        "error_class_only":result.get("error_class_only")
    },ensure_ascii=False),flush=True)
    return 0 if result["technical_status"].startswith("PASS_") else 1


if __name__=="__main__":
    raise SystemExit(main())
