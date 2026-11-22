#!/usr/bin/env python3
"""Primary source trail for QTL normalized phenotype SD and GIGASTROKE effect unit.

Audit exact 36 original gene×tissue source identities; query only publicly
released source metadata and original methods (no individual-level data).
Raw/processed count matrices are not interchangeable with normalized
association-test expression matrices. A genome-wide binary GWAS beta field
is not proven log(OR) by signed numerical correspondence alone.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from pathlib import Path

ZENODO_URL="https://zenodo.org/api/records/4678936"
QTLMAP_URL="https://raw.githubusercontent.com/eQTL-Catalogue/qtlmap/master/docs/usage.md"
PAPER_URL="https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9524349/fullTextXML"
ORIGINAL_PAPER_DOI="10.1038/s41586-022-05165-3"
EXPECTED_QTLS={"QTD000131":("QTS000015","GTEx","artery_aorta",387,10),
               "QTD000171":("QTS000015","GTEx","brain_cortex",205,11),
               "QTD000434":("QTS000025","ROSMAP","brain_naive",560,15)}
MAX_SIZE=11_000_000


def require(expr,reason):
    if not expr:raise ValueError(reason)


def get_source(url,allowed_hosts,cap=MAX_SIZE,fetch=urllib.request.urlopen):
    req=urllib.request.Request(url,headers={
        "User-Agent":"Academic-IS-P3A-primary-expression-and-beta-model-source-review/1.0",
        "Accept":"application/json, application/xml, text/plain; q=0.9, */*; q=0.1"})
    with fetch(req,timeout=65) as r:
        require(getattr(r,"status",200)==200,"original public source HTTP status != 200")
        actual=r.geturl()
        parsed=urllib.parse.urlsplit(actual)
        require(parsed.scheme=="https" and parsed.hostname in allowed_hosts,
                "original source unexpected redirected host")
        data=r.read(cap+1)
    require(0<len(data)<=cap,"original public metadata outside bounded byte size")
    return data,{"requested_url":url,"resolved_url":actual,
                 "raw_bytes":len(data),"SHA256":hashlib.sha256(data).hexdigest()}


def parse_zenodo(data):
    doc=json.loads(data.decode("utf-8"))
    require(isinstance(doc,dict) and
            (doc.get("id")==4678936 or
             "10.5281/zenodo.4678936" in str(doc.get("doi",""))),
            "not original 2021 eQTL Catalogue count-matrix public repository")
    items=doc.get("files",[])
    if isinstance(items,dict):
        items=[dict(info, key=key) for key,info in items.get("entries",{}).items()]
    require(isinstance(items,list),"Zenodo file list not an array")
    names=[]
    for item in items:
        if not isinstance(item,dict):continue
        name=str(item.get("key") or item.get("filename") or "")
        if not name:continue
        names.append({"name":name[:180],
                      "bytes":item.get("size"),
                      "md5":item.get("checksum"),
                      "mentions_study":([s for s in ("GTEx","ROSMAP")
                                         if s.casefold() in name.casefold()])})
    require(names,"Zenodo record contains no filename list; cannot certify inventory")
    exact=[v for v in names if re.search(r"(qcnorm|normalized|normalised|invnorm|qtlmap)",
                                          v["name"],flags=re.I)]
    relevant=[v for v in names if v["mentions_study"]]
    return {
        "original_record_id":4678936,"doi":"10.5281/zenodo.4678936",
        "source_dataset_release":"2021 original compendium, not certified Release 7 studyFiles",
        "n_public_files":len(names),
        "n_files_mentioning_study_names":len(relevant),
        "names_mentioning_GTEx_or_ROSMAP":relevant[:35],
        "n_names_suggesting_normalized_qtlmap_inputs":len(exact),
        "names_suggesting_normalized_qtlmap_inputs":exact[:35],
        "public_raw_or_processed_counts_cannot_prove_same_qtlmap_matrix":True,
        "direct_Release7_phenotype_expression_sdY_verified":False,
        "direct_Release7_sample_meta_verified":False,
        "source_status":"PASS_PUBLIC_ZENODO_RECORD_FILE_INVENTORY",
    }


def parse_qtlmap(data):
    raw=data.decode("utf-8-sig")
    # Source Markdown wraps input column names in inline code/backticks.
    # Normalise only formatting, never waive the three semantic requirements.
    text=re.sub(r"_\\*\\*([A-Za-z_][A-Za-z0-9_]*)\\*\\*_", r"\\1", raw)
    text=re.sub(r"[`*]", "", text)  # preserve identifiers like count_matrix, sample_meta
    required=(
        ("normalized_phenotype_matrix",
         r"count_matrix\s*[-:]\s*normali[sz]ed\s+phenotype\s+matrix"),
        ("source_sample_metadata",
         r"sample_meta\s*[-:]\s*specifies\s+the\s+location"),
        ("exact_sample_intersection",
         r"(?:qtlmap\s+uses|uses)\s+the\s+intersection\s+of\s+samples\s+present\s+in\s+count_matrix\s+and\s+sample_meta")
    )
    checks={name:bool(re.search(pattern,text,re.I))
            for name,pattern in required}
    require(all(checks.values()),
            "official qtlmap source formatting/semantics drift: "+
            repr(checks)+"; source snippets="+repr(
                [x.strip()[:190] for x in raw.splitlines()
                 if any(term in x for term in ("count_matrix","sample_meta"))][:8]))
    return {
      "source_status":"PASS_OFFICIAL_QTLMAP_NORMALIZED_MATRIX_AND_SAMPLE_INTERSECTION",
      "original_required_fields_verified":{
          "normalized_count_matrix":True,
          "sample_meta":True,
          "phenotype_metadata":("pheno_meta -" in text),
          "VCF":("5. vcf -" in text),
      },
      "study_expression_and_sample_intersection_required":True,
      "does_not_assert_genotyped_AN_over_2_equals_regression_N":True,
      "does_not_establish_exact_target_source_expression_sdY":True,
    }


def paper_evidence(data):
    root=ET.fromstring(data)
    authors_doi=[(x.text or "").strip().lower()
       for x in root.iter() if x.tag.split("}")[-1]=="article-id"
       and x.attrib.get("pub-id-type")=="doi"]
    require(ORIGINAL_PAPER_DOI in authors_doi,
            "original 2022 Mishra GIGASTROKE primary paper DOI not found")
    evidence=[]
    patterns=re.compile(
       r"(logistic regress|odds ratio|log[- ]odds|inverse[- ]variance|"
       r"fixed.effect|standard error|meta.analys|effect size|"
       r"effect estimate|binary trait|logarithm of)",flags=re.I)
    paragraphs=[]
    for elem in root.iter():
        if elem.tag.split("}")[-1]!="p":continue
        raw=" ".join(" ".join(elem.itertext()).split())
        if len(raw)>40:paragraphs.append(raw)
    for p in paragraphs:
        if not patterns.search(p):continue
        relevant=[(m.group(0),m.start()) for m in patterns.finditer(p)]
        score=2*sum(bool(re.search(q,p,re.I)) for q in
                     (r"logistic regress",r"log[- ]odds",r"inverse[- ]variance"))
        score+=sum(bool(re.search(q,p,re.I)) for q in
                    (r"GWAS",r"meta.analys",r"cohort",r"stroke"))
        if score<2:continue
        for hit,start in relevant[:2]:
            a=max(0,start-110);b=min(len(p),start+260)
            evidence.append({"keyword":hit,"score":score,
               "primary_paper_text_excerpt":p[a:b][:370]})
    evidence=sorted(evidence,key=lambda x:x["score"],reverse=True)[:20]
    exact_unit_statements=[x for x in evidence
        if re.search(r"(log[- ]odds|logarithm of.*odds|log.*odds ratio)",
                     x["primary_paper_text_excerpt"],re.I)]
    return {
        "source_status":"PASS_ORIGINAL_MISHRA_2022_PAPER_DOI_AND_METHOD_EXCERPTS",
        "verified_publication_DOI":ORIGINAL_PAPER_DOI,
        "original_primary_paper_method_excerpt_count":len(evidence),
        "top_original_primary_paper_method_excerpts":evidence,
        "mentions_explicit_log_odds_around_method":bool(exact_unit_statements),
        "paper_supplement_not_equals_GWAS_Catalog_file_per_SNP_beta_scale_proof":True,
        "per_exact_GWAS_release_GCST90104540_beta_log_OR_verified":False,
        "per_exact_GWAS_release_GCST90104545_beta_log_OR_verified":False,
    }


def check_inputs(manifest,original_qtl,official,cohort):
    groups=manifest.get("groups",[])
    keys={(g["QTL_dataset"],g["gene_id"]) for g in groups}
    require(manifest.get("schema")=="is-p3a-36-real-group-keys-from-legacy-393-v1"
        and manifest.get("scientific_status")=="HOLD"
        and len(groups)==36 and len(keys)==36
        and len({g["gene_id"] for g in groups})==20,
        "36 original QTL phenotype group identifiers lost")
    require(original_qtl.get("schema")==
         "is-p3a-36-primary-qtlsdY-proxy-GWAS-case-control-effect-size-path-v1"
        and original_qtl.get("technical_status")==
         "PASS_3_REAL_ORIGINAL_QTL_SD_Y_PROXY_INPUT_DIAGNOSTICS"
        and original_qtl.get("n_real_gene_tissue_assays")==36
        and original_qtl.get("n_strict_nonpal_source_gene_assay_SNP_records")==148344
        and original_qtl.get("n_formal_coloc_abf_eligible")==0,
        "original raw QTL SE and AN observed data source not certified")
    require({(x["QTL_dataset"],x["gene_id"]) for x in original_qtl["groups"]}==keys,
            "original QTL 36 source arrays not same")
    require(official.get("technical_status")==
         "PASS_3_SOURCE_IDENTITIES_IN_IMMUTABLE_CATALOGUE"
        and set(official.get("QTL_datasets",{}))==set(EXPECTED_QTLS),
        "official original catalogue QTD tissue release identity missing")
    require(cohort.get("schema")==
         "is-p3a-primary-2022-gigastroke-AIS-supplement-BBJ-cohort-36-v1"
        and cohort.get("technical_status")==
         "PASS_ORIGINAL_AIS_SUPPLEMENT_BOTH_ANCESTRY_TOTALS_AND_BBJ"
        and cohort.get("n_real_QTL_assays")==36
        and cohort.get("scientific_status")=="HOLD"
        and cohort["two_ancestries"]["EAS"]["original_total_from_source_K_AIS_cases"]==19032
        and cohort["two_ancestries"]["EUR"]["original_total_from_source_K_AIS_cases"]==62100
        and cohort["two_ancestries"]["EAS"]["BBJ_is_EAS_discovery_contributing_cohort"] is True,
        "source-original GIGASTROKE cohort and BBJ identity unverified")
    for q,(study,label,tissue,N,n_group) in EXPECTED_QTLS.items():
        d=official["QTL_datasets"][q]
        require(d["study_id"]==study and d["study_label"]==label
                and d["sample_group"]==tissue
                and d["nominal_catalogue_sample_size"]==N
                and sum(g["QTL_dataset"]==q for g in groups)==n_group,
                "original QTL Catalogue study ID/tissue/nominal N incompatible")
    return keys


def assemble(manifest,original_qtl,official,cohort,methods):
    keys=check_inputs(manifest,original_qtl,official,cohort)
    byqtl={(g["QTL_dataset"],g["gene_id"]):g for g in original_qtl["groups"]}
    rows=[]
    for g in manifest["groups"]:
        q,gene=g["QTL_dataset"],g["gene_id"]
        src=byqtl[q,gene]
        rows.append({
          "region_id":g["region_id"],"QTL_dataset":q,
          "gene_id":gene,"gene_symbol":g["gene_symbol"],
          "QTL_study_id":official["QTL_datasets"][q]["study_id"],
          "QTL_phenotype_type":"gene_expression_quantitative",
          "QTL_strict_nonpal_original_SNP_assay_rows":
              src["QTL_source_variant_N_and_residual_proxy"]["original_strict_nonpal_SNPs"],
          "QTL_variant_AN_over_2_heterogeneous":
              src["QTL_source_variant_N_and_residual_proxy"]["n_distinct_observed_SNP_N"]>1,
          "QTL_exact_release7_normalized_expression_matrix_verified":False,
          "QTL_exact_release7_sample_inclusion_metadata_verified":False,
          "QTL_direct_per_gene_expression_sdY_measured":False,
          "QTL_gene_regression_effective_N_verified":False,
          "GWAS_beta_log_OR_exact_file_source_verified":False,
          "original_case_control_GWAS_beta_varbeta_minimum_format_possible":True,
          "EAS_BBJ_can_be_independent_outcome_replication":False,
          "formal_coloc_abf_eligible":False,"formal_coloc_susie_eligible":False,
          "status":"HOLD_EXACT_EXPRESSION_INPUT_AND_GWAS_BETA_UNIT",
        })
    return {
        "schema":"is-p3a-36-primary-phenotype-matrix-and-stroke-beta-units-v1",
        "technical_status":"PASS_36_GROUP_SOURCE_METHOD_AND_PUBLIC_MATRIX_INVENTORY",
        "scientific_status":"HOLD","n_real_gene_tissue_assays":36,
        "n_distinct_genes":20,"n_distinct_stroke_loci":1,
        "n_public_zenodo_files":methods["zenodo"]["n_public_files"],
        "n_groups_with_direct_source_expression_sdY_approved":0,
        "n_groups_with_original_GWAS_beta_units_approved":0,
        "per_exact_QTL_model_N_imputation":False,
        "per_exact_GWAS_sample_N_eff_imputation":False,
        "n_formal_coloc_abf_eligible":0,"n_formal_coloc_abf_executed":0,
        "n_formal_coloc_susie_eligible":0,"n_formal_coloc_susie_executed":0,
        "official_public_source_method_evidence":methods,
        "groups":rows,
    }


def main():
    ap=argparse.ArgumentParser()
    for k in ("manifest","original-qtl","official","original-cohort","out"):
        ap.add_argument("--"+k,type=Path,required=True)
    ap.add_argument("--paper-verified-artifact",type=Path,default=None)
    args=ap.parse_args()
    result={"schema":"is-p3a-36-primary-phenotype-matrix-and-stroke-beta-units-v1",
            "technical_status":"HOLD_OFFICIAL_INPUT_SOURCE",
            "scientific_status":"HOLD","n_formal_coloc_abf_eligible":0,
            "n_formal_coloc_abf_executed":0}
    try:
        load=lambda k:json.loads(getattr(args,k.replace("-","_")).read_text(encoding="utf8"))
        manifest,qtl,official,cohort=[load(x) for x in
            ("manifest","original-qtl","official","original-cohort")]
        check_inputs(manifest,qtl,official,cohort)
        methods={}
        for key,url,hosts,cap,parser in (
            ("qtlmap",QTLMAP_URL,{"raw.githubusercontent.com"},300_000,parse_qtlmap),
            ("zenodo",ZENODO_URL,{"zenodo.org","www.zenodo.org"},1_800_000,parse_zenodo),
            ("gigastroke_paper",PAPER_URL,{"www.ebi.ac.uk"},11_000_000,paper_evidence)):
            try:
                data,source=get_source(url,hosts,cap)
                methods[key]={"source":source,**parser(data)}
            except (ValueError,KeyError,TypeError,OSError,TimeoutError,
                    UnicodeError,ET.ParseError,urllib.error.URLError) as exc:
                methods[key]={"source_status":"HOLD_ORIGINAL_PUBLIC_SOURCE",
                              "requested_url":url,
                              "error_type":type(exc).__name__,"error":str(exc)[:230]}
        paper=methods["gigastroke_paper"]
        if (not paper.get("source_status","").startswith("PASS_")
                and args.paper_verified_artifact is not None
                and args.paper_verified_artifact.is_file()):
            archived=json.loads(args.paper_verified_artifact.read_text(encoding="utf8"))
            original=archived.get("official_public_source_method_evidence",{}).get(
                "gigastroke_paper",{})
            archived_meta=original.get("source",{})
            if (original.get("source_status")==
                  "PASS_ORIGINAL_MISHRA_2022_PAPER_DOI_AND_METHOD_EXCERPTS"
                and original.get("verified_publication_DOI")==ORIGINAL_PAPER_DOI
                and archived_meta.get("requested_url")==PAPER_URL
                and archived_meta.get("SHA256")==
                    "7f1105825c4111f73f8b8c36f55897e5c1a19b39b0b240eabc9da3822ac89f79"):
                methods["gigastroke_paper"]={
                  **original,
                  "source_status":"PASS_PREVIOUSLY_VERIFIED_PRIMARY_PAPER_SHA256",
                  "fallback_original_GitHub_run_id":45040084191,
                  "fallback_original_unchanged_primary_article_SHA256":archived_meta["SHA256"],
                  "current_public_source_attempt_status":paper.get("source_status"),
                  "current_public_source_error_type":paper.get("error_type"),
                }
        if all(methods[k].get("source_status","").startswith("PASS_") for k in methods):
            result=assemble(manifest,qtl,official,cohort,methods)
        else:
            result["technical_status"]="HOLD_ONE_OR_MORE_PRIMARY_PUBLIC_METHOD_SOURCES"
            result["official_public_source_method_evidence"]=methods
            result["error"]="Original public source missing: "+",".join(
                 k for k in methods if not methods[k].get("source_status","").startswith("PASS_"))
    except (ValueError,KeyError,TypeError,OSError,UnicodeError) as exc:
        result["technical_status"]="HOLD_PREVIOUS_ORIGINAL_SOURCE_EVIDENCE"
        result["error_type"]=type(exc).__name__
        result["error"]=str(exc)[:250]
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    print(json.dumps({"status":result["technical_status"],
        "public_sources":{k:v["source_status"] for k,v in
            result.get("official_public_source_method_evidence",{}).items()},
        "n_groups":result.get("n_real_gene_tissue_assays"),
        "n_zenodo_public_files":result.get("n_public_zenodo_files"),
        "n_coloc_eligible":0,"error":result.get("error")},ensure_ascii=False))
    return 0 if result["technical_status"].startswith("PASS_") else 1


if __name__=="__main__":
    raise SystemExit(main())
