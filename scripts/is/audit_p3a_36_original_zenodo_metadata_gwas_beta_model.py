#!/usr/bin/env python3
"""P3A 36-group original Zenodo sample metadata and GIGASTROKE beta scale gate.

Audit *public* 2021 eQTL Catalogue metadata files GTExV8.tsv and ROSMAP.tsv.
Verify Zenodo 4678936 API original source SHA256 and publisher MD5 and byte
length for each file. Derive only aggregate tissue/sample counts, never export
sample/donor IDs. These source metadata are not release-7 qtlmap sample_meta
or analysis-stage normalized expression. No sdY/N_eff or log-OR approval.

Query original GWAS Catalog study record only for public effect-unit metadata
where possible; a file named "beta" is not equivalent to verified log(OR).
"""
from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import importlib
import io
import json
import re
import urllib.parse
import urllib.request
from pathlib import Path

PRIOR=importlib.import_module("scripts.is.audit_p3a_36_source_expression_and_gwas_beta_units")
RECORD=PRIOR.ZENODO_URL
RECORD_SHA256="818a5ebd7f448dd685f3f23d90d42f719a520a8fc27db233cb624bee21619682"
FILES={
    "GTExV8.tsv":{"bytes":1769117,"md5":"0e29ac0664ad10403edf7a833d0a4503",
                  "study":"GTExV8","candidate_qtd":["QTD000131","QTD000171"]},
    "ROSMAP.tsv":{"bytes":53064,"md5":"038f2fc144129d401144f3924109e1df",
                  "study":"ROSMAP","candidate_qtd":["QTD000434"]},
}
REST="https://www.ebi.ac.uk/gwas/rest/api/v2/studies/"
GWAS={"EAS":"GCST90104545","EUR":"GCST90104540"}
MAX_RECORDBYTES=150_000
MAX_FILEBYTES=2_200_000
MAX_META_ROWS=80_000
ALLOWED_RECORD_HOSTS={"zenodo.org","www.zenodo.org"}
ALLOWED_GWAS_HOSTS={"www.ebi.ac.uk"}
SELECTABLE_FIELDS={"sample_id","genotype_id","qtl_group","tissue",
                   "tissue_id","tissue_label","smtsd","smts","sample_group",
                   "condition","condition_label","subject_id","donor_id",
                   "run_id","study","run","genotype","sample",
                   "rna_seq","file_id","phenotype_id","sex","age","batch"}
UNIQUE_SAMPLE_HINTS={"sample_id","genotype_id","subject_id","donor_id","pseudo_id"}
GROUP_HINTS={"qtl_group","sample_group","tissue","tissue_label","smtsd","smts",
             "tissue_id"}
TARGET_GROUPS={
 "QTD000131":(("artery_aorta", "artery - aorta", "artery aorta"),387),
 "QTD000171":(("brain_cortex", "brain - cortex", "brain cortex"),205),
 "QTD000434":(("brain_naive","brain dlpfc","dlpfc","brain_naive_dlpfc"),560),
}


def require(cond,msg):
    if not cond:raise ValueError(msg)


def get(url,hosts,maxlen,open_url=urllib.request.urlopen):
    req=urllib.request.Request(url,headers={
        "User-Agent":"MASTER-DEGREE-P3A-primary-Zenodo-small-source-and-AIS-beta-units/1.0",
        "Accept":"application/json, text/plain; q=0.9, */*; q=0.1"})
    with open_url(req,timeout=75) as response:
        require(getattr(response,"status",200)==200,"unexpected original source HTTP status")
        resolved=response.geturl()
        loc=urllib.parse.urlsplit(resolved)
        require(loc.scheme=="https" and loc.hostname in hosts,
                "publisher response redirected to unapproved host")
        data=response.read(maxlen+1)
    require(10<len(data)<=maxlen,"original source empty or unexpectedly oversized")
    return data,{"requested_url":url,"resolved_url":resolved,
                 "read_bytes":len(data),"SHA256":hashlib.sha256(data).hexdigest()}


def require_original_record(data):
    sha=hashlib.sha256(data).hexdigest()
    require(sha==RECORD_SHA256,"official 2021 Zenodo API original source SHA256 changed")
    original=PRIOR.parse_zenodo(data)
    require(original["source_status"]=="PASS_PUBLIC_ZENODO_RECORD_FILE_INVENTORY"
            and original["n_public_files"]==32,"Zenodo original published inventory changed")
    j=json.loads(data.decode("utf-8"))
    listing=j.get("files",[])
    require(isinstance(listing,list),"original Zenodo files not an array")
    byname={x.get("key"):x for x in listing if isinstance(x,dict)}
    for name,expected in FILES.items():
        item=byname.get(name)
        require(item is not None,"original source file not in publisher Zenodo: "+name)
        require(item.get("size")==expected["bytes"]
                and str(item.get("checksum"))=="md5:"+expected["md5"],
                "published source metadata checksum or length drift: "+name)
    return byname,{"official_record_url":RECORD,
                   "official_record_original_SHA256":sha,
                   "official_record_n_public_files":original["n_public_files"]}


def official_file_url(item,name):
    """Use publisher link from pinned record, not file URLs from search snippets."""
    links=item.get("links",{})
    require(isinstance(links,dict),"missing official file URL mapping")
    options=[links.get("self"),links.get("download"),links.get("content")]
    for link in options:
        if not isinstance(link,str):continue
        parsed=urllib.parse.urlsplit(link)
        if parsed.scheme=="https" and parsed.hostname in ALLOWED_RECORD_HOSTS:
            # strict unambiguous exact original filename decoded as URL path
            if urllib.parse.unquote(parsed.path).endswith("/"+name) or (
                urllib.parse.unquote(parsed.path).endswith("/"+name+"/content")):
                return link
    raise ValueError("No publisher-linked HTTPS source download for "+name)


def get_official_file(name,item,fetch=get):
    expected=FILES[name]
    official_url=official_file_url(item,name)
    data,source=fetch(official_url,ALLOWED_RECORD_HOSTS,MAX_FILEBYTES)
    md5=hashlib.md5(data,usedforsecurity=False).hexdigest()
    require(len(data)==expected["bytes"] and md5==expected["md5"],
            "original small public sample metadata file MD5/size mismatch: "+name)
    source["original_publisher_deposited_MD5"]=expected["md5"]
    source["original_MD5_confirmed_on_actual_download"]=True
    return data,source


def norm(s):
    return re.sub(r"[^a-z0-9]+","_",str(s).strip().lower()).strip("_")


def source_qc_flag(value):
    """Classify public metadata QC boolean, not any patient data."""
    val=norm(value)
    if val in ("true","t","yes","y","1","pass","passed"):return True
    if val in ("false","f","no","n","0","fail","failed"):return False
    return None


def safe_count_source_metadata(data,name):
    """Only numerical aggregate counts; never any original sample identifiers."""
    text=data.decode("utf-8-sig")
    require(text and "\n" in text and not "\x00" in text,
            "original source sample metadata not text")
    first=text.splitlines()[0]
    tab=first.count("\t")
    comma=first.count(",")
    delim="\t" if tab>=comma and tab>0 else "," if comma>0 else None
    require(delim is not None,"original source sample table delimiter unsupported")
    reader=csv.reader(io.StringIO(text),delimiter=delim)
    original_header=next(reader,None)
    require(original_header and 2<=len(original_header)<=250,
            "original source metadata missing 2..250-column header")
    canonical=[norm(x) for x in original_header]
    require(len(canonical)==len(set(canonical)),
            "original metadata duplicate normalized column names")
    recognized=[c for c in canonical if c in SELECTABLE_FIELDS]
    locations={c:i for i,c in enumerate(canonical)}
    require(all(k in locations for k in ("pseudo_id","rna_qc_passed","genotype_qc_passed")),
            "official source metadata missing pseudo_id or exact RNA/genotype QC status")
    unique_sets={c:set() for c in UNIQUE_SAMPLE_HINTS if c in locations}
    relevant_groups={q:0 for q in FILES[name]["candidate_qtd"]}
    qc_by_qtl={q:{"rows":0,"unique_ids":set(),"both_pass_ids":set(),
                  "joint_qc":collections.Counter(),"unparsed_flags":0}
               for q in relevant_groups}
    categorical=collections.defaultdict(collections.Counter)
    known_groups=[c for c in GROUP_HINTS if c in locations]
    n=0;malformed=0;blanks=0
    for row in reader:
        if not row or (len(row)==1 and not row[0].strip()):
            blanks+=1;continue
        n+=1
        require(n<=MAX_META_ROWS,"original metadata row count exceeded bound")
        if len(row)!=len(canonical):
            malformed+=1;continue
        for c,seen in unique_sets.items():
            v=row[locations[c]].strip()
            if v:seen.add(v)
        for c in known_groups:
            v=row[locations[c]].strip()
            if v:
                key=norm(v)
                if len(categorical[c])<300 or key in categorical[c]:
                    categorical[c][key]+=1
        for q in relevant_groups:
            aliases=TARGET_GROUPS[q][0]
            for c in known_groups:
                v=norm(row[locations[c]])
                if any(norm(alias)==v for alias in aliases):
                    relevant_groups[q]+=1
                    record=qc_by_qtl[q]
                    record["rows"]+=1
                    pseudo=row[locations["pseudo_id"]].strip()
                    require(bool(pseudo),"public original metadata missing pseudonymous sample ID")
                    record["unique_ids"].add(pseudo)
                    rna=source_qc_flag(row[locations["rna_qc_passed"]])
                    geno=source_qc_flag(row[locations["genotype_qc_passed"]])
                    rlabel="PASS" if rna is True else "FAIL" if rna is False else "UNKNOWN"
                    glabel="PASS" if geno is True else "FAIL" if geno is False else "UNKNOWN"
                    record["joint_qc"]["RNA_"+rlabel+"_GENOTYPE_"+glabel]+=1
                    if rna is None or geno is None:
                        record["unparsed_flags"]+=1
                    if rna is True and geno is True:
                        record["both_pass_ids"].add(pseudo)
                    break
    require(n>5 and malformed==0,"public metadata sample rows not parsed completely")
    qc_summary={}
    for q,v in qc_by_qtl.items():
        both=v["joint_qc"]["RNA_PASS_GENOTYPE_PASS"]
        expected_n=TARGET_GROUPS[q][1]
        qc_summary[q]={
            "2021_metadata_rows_in_tissue":v["rows"],
            "2021_distinct_pseudo_ids_in_tissue_count_only":len(v["unique_ids"]),
            "2021_both_rna_genotype_QC_passed_rows":both,
            "2021_both_QC_passed_distinct_pseudo_ids":len(v["both_pass_ids"]),
            "2021_public_RNA_genotype_QC_joint_flags":dict(sorted(v["joint_qc"].items())),
            "unparsed_source_QC_flags":v["unparsed_flags"],
            "catalogue_release7_declared_tissue_N_reference_only":expected_n,
            "2021_both_QC_rows_equal_Release7_declared_N_DESCRIPTIVE":
                both==expected_n,
            "2021_QC_pass_is_not_release7_qtlmap_sample_intersection_proof":True,
        }
    return {
        "source_original_file":name,"original_deposit_study":FILES[name]["study"],
        "original_2021_metadata_row_count":n,
        "original_table_number_columns":len(canonical),
        "original_source_schema_header_names_only_no_sample_values":canonical,
        "recognized_source_column_names_no_values":recognized,
        "n_recognized_schema_columns":len(recognized),
        "recognized_sample_or_donor_identifier_columns":
            sorted(unique_sets),
        "distinct_sample_donor_identifier_count_only":
            {key:len(value) for key,value in unique_sets.items()},
        "recognized_tissue_or_qtl_group_columns":sorted(known_groups),
        "source_original_target_label_row_matches_lower_bound_only":relevant_groups,
        "public_2021_qtl_group_RNA_and_genotype_QC_flag_aggregates":qc_summary,
        "n_blank_physical_rows_ignored":blanks,
        "source_column_group_distinct_label_counts_only":
            {k:len(v) for k,v in categorical.items()},
        "source_2021_metadata_n_equal_Release7_qtlmap_intersection_verified":False,
        "sample_identifier_values_uploaded":False,
        "exact_Release7_normalised_expression_matrix_present":False,
        "direct_original_36_gene_quantitative_sdY_computed":False,
        "reusable_raw_data_sample_meta_for_release7_approved":False,
    }


def gwas_record(accession,fetch=get):
    require(accession in GWAS.values(),"unsupported GIGASTROKE publisher accession")
    url=REST+accession
    raw,source=fetch(url,ALLOWED_GWAS_HOSTS,MAX_RECORDBYTES)
    j=json.loads(raw.decode("utf-8"))
    require(isinstance(j,dict),"GWAS Catalog original study record not JSON mapping")
    # Official 2026 GWAS Catalog API v2 uses study_accession; v1 is HTTP 410.
    got=j.get("study_accession") or j.get("accession_id") or j.get("accessionId")
    require(got==accession,"GWAS Catalog returned incorrect original accession")
    descriptive={}
    for key in ("study_accession","disease_trait","initial_sample_size",
                "replication_sample_size","genotyping_technology","cohort",
                "pubmed_id","first_author","study_design_comment",
                "studyDesignComment","initialSampleSize","diseaseTrait"):
        v=j.get(key)
        if isinstance(v,(str,int,float,bool)) and v is not None:
            descriptive[key]=str(v)[:220]
    # No "OR" inference from prose, study summary or beta header.
    return {"publisher_rest_accession":accession,
        "original_public_GWAS_Catalog_study_record_SHA256":source["SHA256"],
        "original_GWAS_Catalog_official_rest_url":url,
        "published_study_descriptor_fields":descriptive,
        "effect_size_beta_exact_original_summary_file_explicit_log_OR_verified":False,
        "hazard_ratio_source_contributions_undone_or_harmonized_into_logOR_proven":False,
        "source_publisher_yml_beta_scale_identified":False,
        "outcome_case_control_coloc_beta_varbeta_eligibility_approved":False,
        "status":"PASS_ORIGINAL_GWAS_STUDY_METADATA_V2_NOT_BETA_SCALE"}


def check_source_groups(m,previous,original_official,original_gwas):
    groups=m.get("groups",[])
    k={(x["QTL_dataset"],x["gene_id"]) for x in groups}
    require(m.get("schema")=="is-p3a-36-real-group-keys-from-legacy-393-v1"
            and m.get("scientific_status")=="HOLD" and len(k)==36
            and len(groups)==36 and len({x["gene_id"] for x in groups})==20,
            "36 original QTL groups not unique and verified")
    require(previous.get("schema")==
              "is-p3a-36-primary-phenotype-matrix-and-stroke-beta-units-v1"
            and previous.get("technical_status")==
              "PASS_36_GROUP_SOURCE_METHOD_AND_PUBLIC_MATRIX_INVENTORY"
            and previous.get("n_real_gene_tissue_assays")==36
            and previous.get("n_groups_with_direct_source_expression_sdY_approved")==0
            and previous.get("n_groups_with_original_GWAS_beta_units_approved")==0,
            "previous 36 source matrix/method publication audit missing")
    require({(x["QTL_dataset"],x["gene_id"]) for x in previous["groups"]}==k,
            "36 gene-assay keys changed since previous source integrity audit")
    require(original_official.get("technical_status")==
            "PASS_3_SOURCE_IDENTITIES_IN_IMMUTABLE_CATALOGUE"
            and set(original_official.get("QTL_datasets",{}))==
               set(PRIOR.EXPECTED_QTLS),
            "original QTL Catalogue source study identities changed")
    require(original_gwas.get("schema")==
            "is-p3a-primary-2022-gigastroke-AIS-supplement-BBJ-cohort-36-v1"
            and original_gwas.get("technical_status")==
            "PASS_ORIGINAL_AIS_SUPPLEMENT_BOTH_ANCESTRY_TOTALS_AND_BBJ"
            and original_gwas.get("scientific_status")=="HOLD",
            "2022 publisher original AIS discovery cohort evidence missing")
    for q,(study,label,tissue,N,ng) in PRIOR.EXPECTED_QTLS.items():
        actual=original_official["QTL_datasets"][q]
        require(actual["study_id"]==study and actual["study_label"]==label
                and actual["sample_group"]==tissue and
                actual["nominal_catalogue_sample_size"]==N
                and sum(x["QTL_dataset"]==q for x in groups)==ng,
                "QTL release7 study/tissue reference metadata drift")
    return k


def evaluate(m,prev,official,cohort,files_by_study,gwas,provenance):
    check_source_groups(m,prev,official,cohort)
    require(set(files_by_study)==set(FILES),
            "public original small GTEx/ROSMAP metadata files unavailable")
    out=[]
    for g in m["groups"]:
        q=g["QTL_dataset"];gene=g["gene_id"]
        name="ROSMAP.tsv" if q=="QTD000434" else "GTExV8.tsv"
        source=files_by_study[name]
        out.append({
          "region_id":g["region_id"],"QTL_dataset":q,
          "gene_id":gene,"gene_symbol":g["gene_symbol"],
          "source_2021_aggregated_metadata_file":name,
          "source_2021_original_metadata_rows_total":source["original_2021_metadata_row_count"],
          "source_original_source_tabular_schema_verified":True,
          "source_original_qtl_reference_Release7_nominal_N":
              official["QTL_datasets"][q]["nominal_catalogue_sample_size"],
          "source_2021_sample_group_matches_2026_Release7_qtlmap_verified":False,
          "source_exact_release7_expression_sample_meta_linked":False,
          "quantitative_original_expression_sdY_directly_observed":False,
          "original_GIGASTROKE_EUR_or_EAS_beta_logOR_unit_approved":False,
          "can_calculate_formal_coloc_abf":False,
          "can_calculate_formal_coloc_susie":False,
          "status":"HOLD_NEEDS_RELEASE7_NORMALISED_EXPRESSION_AND_GWAS_BETA_SCALE",
        })
    require(len(out)==36 and len({x["gene_id"] for x in out})==20,
            "36 source groups missing after actual small file audit")
    return {
      "schema":"is-p3a-36-official-2021-GTEx-ROSMAP-small-sample-metadata-and-2022-gwas-beta-source-v1",
      "technical_status":"PASS_OFFICIAL_SMALL_GTEx_ROSMAP_PUBLIC_SAMPLE_METADATA_ORIGINAL_MD5",
      "scientific_status":"HOLD",
      "n_actual_gene_tissue_groups":36,"n_actual_unique_genes":20,
      "n_stroke_loci":1,
      "n_original_official_small_metadata_sources_read_and_MD5_checked":2,
      "official_Zenodo_publisher_source_provenance":provenance,
      "public_metadata_aggregate_sources":files_by_study,
      "official_GWAS_study_source_lookups":gwas,
      "direct_Release7_phenotype_sdY_approved":0,
      "exact_Release7_QTL_effective_N_approved":0,
      "exact_GIGASTROKE_beta_logOR_effect_unit_approved":0,
      "n_formal_coloc_abf_eligible":0,"n_formal_coloc_abf_executed":0,
      "n_formal_coloc_susie_eligible":0,"n_formal_coloc_susie_executed":0,
      "sample_ids_or_individual_level_data_in_artifact":False,
      "groups":out,
    }


def main():
    p=argparse.ArgumentParser()
    for name in ("manifest","prior","official","gwas-cohort","out"):
        p.add_argument("--"+name,type=Path,required=True)
    args=p.parse_args()
    result={"schema":"is-p3a-36-official-2021-GTEx-ROSMAP-small-sample-metadata-and-2022-gwas-beta-source-v1",
            "technical_status":"HOLD_ORIGINAL_PUBLIC_SMALL_METADATA",
            "scientific_status":"HOLD",
            "n_formal_coloc_abf_eligible":0,"n_formal_coloc_abf_executed":0}
    try:
        files=[json.loads(getattr(args,name.replace("-","_")).read_text(encoding="utf8"))
               for name in ("manifest","prior","official","gwas-cohort")]
        manifest,previous,official,cohort=files
        check_source_groups(manifest,previous,official,cohort)
        raw,api=get(RECORD,ALLOWED_RECORD_HOSTS,MAX_RECORDBYTES)
        items,source=require_original_record(raw)
        summaries={};prov={}
        for name in FILES:
            content,info=get_official_file(name,items[name])
            summaries[name]=safe_count_source_metadata(content,name)
            prov[name]=info
        optional={}
        for ancestry,accession in GWAS.items():
            try:
                optional[ancestry]=gwas_record(accession)
            except (OSError,ValueError,TypeError,KeyError,TimeoutError,
                    json.JSONDecodeError,UnicodeError) as err:
                optional[ancestry]={
                    "accession":accession,
                    "status":"HOLD_GWAS_STUDY_API_OR_SOURCE",
                    "error_type":type(err).__name__,
                    "error":str(err)[:180],
                    "effect_size_beta_exact_original_summary_file_explicit_log_OR_verified":False,
                }
        result=evaluate(manifest,previous,official,cohort,summaries,
                        optional,{"api":api,**source,"files":prov})
    except (OSError,ValueError,TypeError,KeyError,TimeoutError,
            json.JSONDecodeError,UnicodeError) as err:
        result["technical_status"]="HOLD_SOURCE_MD5_OR_SAMPLE_METADATA_LAYOUT"
        result["error_type"]=type(err).__name__
        result["error"]=str(err)[:400]
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    print(json.dumps({
         "technical_status":result["technical_status"],
         "n_real_groups":result.get("n_actual_gene_tissue_groups"),
         "source_metadata":{
            k:{"rows":v.get("original_2021_metadata_row_count"),
               "columns":v.get("original_table_number_columns"),
               "schema":v.get("recognized_source_column_names_no_values"),
               "target_counts":v.get("source_original_target_label_row_matches_lower_bound_only"),
               "source_QC":v.get("public_2021_qtl_group_RNA_and_genotype_QC_flag_aggregates")}
            for k,v in result.get("public_metadata_aggregate_sources",{}).items()},
         "gwas_metadata":{k:v["status"] for k,v in result.get(
             "official_GWAS_study_source_lookups",{}).items()},
         "n_coloc_approved":0,"error":result.get("error")
      },ensure_ascii=False))
    return 0 if result["technical_status"].startswith("PASS_") else 1


if __name__=="__main__":
    raise SystemExit(main())
