#!/usr/bin/env python3
"""P3A primary GWAS-SSF YAML / QTL per-SNP an/2 audit on authentic 36 assays.

Publisher GWAS YAML case/control phenotype fields are observational, and can
certify those source-published declarations. They do not certify per-SNP GWAS
effective N, BBJ independence, or independent stroke replication.

Official eQTL Catalogue Columns.md says *usually* autosomal an/2 is sample
size. Count and compare it per variant/assay; never silently substitute an/2
for regression effective N or automatically authorize coloc.
"""
from __future__ import annotations

import argparse
import collections
import hashlib
import importlib
import json
import math
import time
import urllib.request
from pathlib import Path

QTL=importlib.import_module("scripts.is.audit_p3a_indexed_qtl_real_chr10")
CIS=importlib.import_module("scripts.is.audit_p3a_official_cis_union_36")
PINNED=importlib.import_module("scripts.is.p3a_publisher_md5_gate")

EXPECTED_SAMPLES={"QTD000131":387,"QTD000171":205,"QTD000434":560}
EXPECTED_LITERATURE={"EAS":{"case_count":19032,"control_count":237242},
                     "EUR":{"case_count":62100,"control_count":1234808}}
GWAS_META_KEYS=("gwas_id","trait_description","sample_ancestry_category",
                "sample_ancestry","sample_size","case_control_study",
                "case_count","control_count","genome_assembly",
                "coordinate_system","data_file_name","data_file_md5sum",
                "date_metadata_last_modified","author_notes")
MAX_YAML_BYTES=160_000


def require(ok,reason):
    if not ok:raise ValueError(reason)


def fetch_primary_gwas_metadata(ancestry,open_url=urllib.request.urlopen):
    import yaml
    meta=PINNED.SOURCES[ancestry]
    original_url,metadata_url=PINNED.build_urls(meta)
    req=urllib.request.Request(metadata_url,headers={
        "User-Agent":"Academic-IS-P3A-Primary-GWAS-Case-Control-Metadata-QC/1.0"})
    with open_url(req,timeout=65) as response:
        require(getattr(response,"status",200)==200,
                "primary GWAS YAML missing HTTP 200")
        payload=response.read(MAX_YAML_BYTES+1)
    require(0<len(payload)<=MAX_YAML_BYTES,
            "primary GWAS YAML empty/too large")
    doc=yaml.safe_load(payload.decode("utf-8-sig"))
    require(isinstance(doc,dict),"primary GWAS YAML not mapping")
    require(doc.get("gwas_id")==meta["accession"]
            and doc.get("genome_assembly")=="GRCh37"
            and doc.get("data_file_name")==
                meta["accession"]+"_buildGRCh37.tsv.gz"
            and doc.get("data_file_md5sum")==meta["publisher_md5_pinned"],
            "published GWAS YAML accession/build/pinned original MD5 changed")
    declared={k:doc.get(k) for k in GWAS_META_KEYS}
    for key in GWAS_META_KEYS:
        val=declared[key]
        if not isinstance(val,(type(None),int,float,str,bool,list)):
            declared[key]=str(val)[:180]
    count_keys=("sample_size","case_count","control_count")
    integers={}
    for key in count_keys:
        val=doc.get(key)
        if isinstance(val,int) and not isinstance(val,bool):
            integers[key]=val
        elif isinstance(val,str) and val.strip().replace(",","").isdigit():
            integers[key]=int(val.strip().replace(",",""))
        else:
            integers[key]=None
    case=integers["case_count"]
    control=integers["control_count"]
    total=integers["sample_size"]
    counts_reconcile=(case is not None and control is not None
                      and total is not None and case>0 and control>0
                      and case+control==total)
    compare=EXPECTED_LITERATURE[ancestry]
    reference_matches=(counts_reconcile
                and case==compare["case_count"]
                and control==compare["control_count"])
    actual_trait=str(doc.get("trait_description") or "")
    phenotype_mentioned=any(x in actual_trait.lower()
          for x in ("ischemic stroke","ischaemic stroke","ais"))
    return {
        "source_accession":meta["accession"],
        "original_publisher_YAML":metadata_url,
        "original_compressed_gwas_file":original_url,
        "original_publisher_YAML_SHA256":hashlib.sha256(payload).hexdigest(),
        "original_release_source_identity_verified":True,
        "publisher_fields":declared,
        "declared_sample_counts":integers,
        "case_control_sum_reconciles_primary_file":counts_reconcile,
        "followup_publication_case_control_numbers_match":reference_matches,
        "primary_phenotype_text_mentions_ischemic_stroke":phenotype_mentioned,
        "primary_GWAS_N_metadata_technically_complete":bool(counts_reconcile),
        "scientific_case_fraction_and_N_approved_for_coloc":False,
        "BBJ_independence_approved":False,
        "source_status":("PASS_PUBLISHER_CASE_CONTROL_DESCRIPTIVE"
                         if counts_reconcile else "HOLD_PUBLISHER_CASE_CONTROL_MISSING"),
    }


def new_gene_stats():
    return {"count_valid_nonpal_original_QTL":0,
       "count_an_valid_even_autosomal":0,
       "count_an_missing_or_invalid":0,
       "count_an_above_dataset_N":0,
       "count_ac_an_invalid":0,
       "count_ma_samples_invalid":0,
       "count_QTL_maf_vs_minor_ac_an_delta_gt_0_01":0,
       "min_variant_an_over_2":None,
       "max_variant_an_over_2":None,
       "variant_an_over_2_histogram":collections.Counter()}


def aggregate_an(reg,genes,tabix):
    result={g:new_gene_stats() for g in genes}
    seen={g:set() for g in genes}
    actual_datasetN=EXPECTED_SAMPLES[reg]
    for line in tabix:
        fields=line.rstrip("\r\n").split("\t")
        require(len(fields)==len(QTL.OFFICIAL_COLUMN_ORDER),
                "original source physical 19-column layout changed")
        raw=dict(zip(QTL.OFFICIAL_COLUMN_ORDER,fields))
        gene=raw["gene_id"]
        if gene not in result:continue
        parsed=importlib.import_module(
          "scripts.is.audit_p3a_36_real_gwas_qtl_reciprocal_snp").qtl_record(raw)
        if parsed is None:continue
        key=(parsed["pos38"],parsed["ref"],parsed["alt"],raw["molecular_trait_id"])
        if key in seen[gene]:continue
        seen[gene].add(key)
        stats=result[gene]
        stats["count_valid_nonpal_original_QTL"]+=1
        raw_an=raw["an"]
        raw_ac=raw["ac"]
        raw_ma=raw["ma_samples"]
        try:
            an=float(raw_an)
            if not (math.isfinite(an) and an>0 and an==int(an)
                    and int(an)%2==0):
                raise ValueError("odd or invalid original an")
            N=int(an)//2
        except (TypeError,ValueError,OverflowError):
            stats["count_an_missing_or_invalid"]+=1
            continue
        stats["count_an_valid_even_autosomal"]+=1
        stats["variant_an_over_2_histogram"][N]+=1
        stats["min_variant_an_over_2"]=min(
            stats["min_variant_an_over_2"] or N,N)
        stats["max_variant_an_over_2"]=max(
            stats["max_variant_an_over_2"] or N,N)
        if N>actual_datasetN:
            stats["count_an_above_dataset_N"]+=1
        try:
            ac=float(raw_ac)
            if not math.isfinite(ac) or not 0<=ac<=an:
                raise ValueError("ac invalid")
            observed_minor=min(ac/an,1-ac/an)
            if abs(observed_minor-parsed["maf"])>0.01:
                stats["count_QTL_maf_vs_minor_ac_an_delta_gt_0_01"]+=1
        except (TypeError,ValueError):
            stats["count_ac_an_invalid"]+=1
        try:
            ma=float(raw_ma)
            if not math.isfinite(ma) or not 0<=ma<=N:
                raise ValueError("ma_samples larger than available donors")
        except (TypeError,ValueError):
            stats["count_ma_samples_invalid"]+=1
    for r in result.values():
        counter=r.pop("variant_an_over_2_histogram")
        r["n_distinct_variant_an_over_2"]=len(counter)
        r["top_variant_an_over_2_count"]=[
          {"variant_n":n,"records":c} for n,c in counter.most_common(8)]
        r["all_valid_SNP_have_an_over_2"]=(
            r["count_valid_nonpal_original_QTL"]>0 and
            r["count_an_missing_or_invalid"]==0)
        r["has_single_constant_n_across_SNPs"]=len(counter)==1
        r["QTL_model_regression_effective_N_independently_approved"]=False
    return result


def scan_dataset(dataset,source,groupgenes,open_tabix=None):
    if open_tabix is None:
        import pysam
        open_tabix=pysam.TabixFile
    with open_tabix(source) as index:
        chrom="10" if "10" in index.contigs else "chr10" if "chr10" in index.contigs else None
        require(chrom is not None,"original source QTL chr10 missing")
        stream=index.fetch(chrom,
                CIS.EXPECTED_WINDOW["start"]-1,CIS.EXPECTED_WINDOW["end"])
        return aggregate_an(dataset,groupgenes,stream)


def check_inputs(manifest,union,official):
    groups=manifest.get("groups",[])
    require(manifest.get("schema")=="is-p3a-36-real-group-keys-from-legacy-393-v1"
            and manifest.get("scientific_status")=="HOLD"
            and len(groups)==36
            and len({(g["region_id"],g["QTL_dataset"],g["gene_id"])
                     for g in groups})==36
            and len({g["gene_id"] for g in groups})==20,
            "original archived assay groups not confirmed")
    require(union.get("schema")=="is-p3a-36-expanded-original-cis-reference-union-v1"
            and union.get("technical_status")==
                "PASS_THREE_REAL_INDEXED_REFERENCE_UNION_SURVEYS"
            and union.get("scientific_status")=="HOLD"
            and union.get("n_genuine_assay_groups")==36
            and union.get("n_approved_coloc_abf")==0
            and union.get("cis_union_reference_Ensembl105_GRCh38_1based")==
                CIS.EXPECTED_WINDOW
            and len(union.get("groups",[]))==36,
            "actual QTL reference window source baseline missing")
    require(official.get("technical_status")==
            "PASS_3_SOURCE_IDENTITIES_IN_IMMUTABLE_CATALOGUE"
            and set(official["QTL_datasets"])==set(EXPECTED_SAMPLES),
            "immutable source dataset provenance mismatch")
    old={(x["QTL_dataset"],x["gene_id"]):x for x in union["groups"]}
    require(len(old)==36
            and set(old)=={(x["QTL_dataset"],x["gene_id"]) for x in groups},
            "old group source registry differs")
    for qtl,expected in EXPECTED_SAMPLES.items():
        require(official["QTL_datasets"][qtl]["nominal_catalogue_sample_size"]==expected
                and official["QTL_datasets"][qtl]["QTL_effect_allele"]=="ALT_GRCh38",
                "individual tissue source N or ALT orientation drift")
    return old


def main():
    ap=argparse.ArgumentParser()
    for p in ("manifest","union","official","out"):
        ap.add_argument("--"+p,type=Path,required=True)
    a=ap.parse_args()
    result={"schema":"is-p3a-primary-gwas-yaml-QTL-allele-N-36-v1",
       "technical_status":"HOLD_SOURCE_INCOMPLETE","scientific_status":"HOLD",
       "formal_coloc_abf_eligible":0,"formal_coloc_abf_executed":0,
       "formal_coloc_susie_eligible":0,"formal_coloc_susie_executed":0}
    try:
        load=lambda p:json.loads(getattr(a,p).read_text(encoding="utf8"))
        manifest,union,official=(load(s) for s in ("manifest","union","official"))
        baseline=check_inputs(manifest,union,official)
        primary={}
        for ancestry in ("EAS","EUR"):
            primary[ancestry]=fetch_primary_gwas_metadata(ancestry)
        all_qtl={}
        for i,qtl in enumerate(EXPECTED_SAMPLES):
            if i:time.sleep(8)  # one region per source, no automatic retries
            qurl=official["QTL_datasets"][qtl]["source_https_probe"]
            genes={g["gene_id"] for g in manifest["groups"]
                   if g["QTL_dataset"]==qtl}
            all_qtl[qtl]=scan_dataset(qtl,qurl,genes)
        rows=[]
        for x in manifest["groups"]:
            qtl,gene=x["QTL_dataset"],x["gene_id"]
            stats=all_qtl[qtl][gene]
            expected=baseline[qtl,gene]["extended_context_valid_nonpal_SNP_assay_rows"]
            require(stats["count_valid_nonpal_original_QTL"]==expected,
                    "independently repeated original full QTL SNP count drift")
            rows.append({
                "QTL_dataset":qtl,"gene_id":gene,"gene_symbol":x["gene_symbol"],
                "original_variant_allele_N_summary":stats,
                "study_specific_QTL_model_effective_N_approved":False,
                "source_GWAS_stroke_case_fraction_approved":False,
                "coloc_abf_eligible":False,"coloc_susie_eligible":False,
                "status":"HOLD_EXPOSURE_MODEL_N_AND_OUTCOME_COHORT_PROOF",
            })
        require(len(rows)==36,"lost QTL groups")
        result.update({
           "technical_status":"PASS_PRIMARY_PUBLISHER_METADATA_AND_36_QTL_AN_SOURCE_READ",
           "n_true_assay_groups":36,"n_distinct_genes":20,
           "primary_gwas_yaml":primary,"groups":rows,
           "QTL_an_documentation":("Source eQTL Catalogue release7 tabix/Columns.md: "
               "for autosomal variants an/2 is usually sample size, not always "
               "regression effective N or study phenotype-specific N"),
           "study_genotyped_allele_N_is_not_coloc_regression_N":True,
           "publication_counts_not_automatically_approved":True,
           "all_36_Coloc_HOLD":True,
        })
    except (ValueError,TypeError,KeyError,OSError,UnicodeError,ImportError,
            TimeoutError,OverflowError) as exc:
        result["technical_status"]="HOLD_PRIMARY_METADATA_OR_QTL_SOURCE"
        result["error"]=str(exc)[:360]
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",
                     encoding="utf8")
    print(json.dumps({"technical_status":result["technical_status"],
        "primary_gwas_fields":{k:{"counts":v["declared_sample_counts"],
          "trait":v["publisher_fields"]["trait_description"],
          "identity":v["original_release_source_identity_verified"],
          "case_control_reconciled":v["case_control_sum_reconciles_primary_file"]}
          for k,v in result.get("primary_gwas_yaml",{}).items()},
        "n_assay_groups":result.get("n_true_assay_groups"),
        "formal_coloc_abf_eligible":0,
        "error":result.get("error")},ensure_ascii=False))
    return 0 if result["technical_status"].startswith("PASS_") else 1


if __name__=="__main__":
    raise SystemExit(main())
