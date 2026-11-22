#!/usr/bin/env python3
"""Audits source GIGASTROKE 2022 Supplementary Table 1 original XLSX.

The Xlsx is a public primary supplemental source. Read OOXML directly with
standard library, not a private cohort/genotype file. Case field K is any
ischemic stroke (AIS); field P is shared discovery controls. Distinguish
EAS AIS discovery, incidental-only and PGS validation; NEVER assert
participant-level independence or automatically approve coloc sample N.

Pair this primary evidence with original GWAS Catalog publisher YAML and the
authentic existing 36 QTL gene/tissue keys in one machine-verifiable report.
"""
from __future__ import annotations

import argparse
import hashlib
import io
import json
import re
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
import zipfile
from collections import Counter
from pathlib import Path

FILENAME="41586_2022_5165_MOESM4_ESM.xlsx"
URLS=(
    "https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-022-05165-3/MediaObjects/"+FILENAME,
    "https://pmc.ncbi.nlm.nih.gov/articles/instance/9524349/bin/"+FILENAME,
)
EXPECTED={"EAS":{"cases":19032,"controls":237242},
          "EUR":{"cases":62100,"controls":1234808}}
EAS_AIS_CONTRIBUTORS={"KCPS2_ASN**":46,"BIOBANK_JAPAN*":17651,
                      "HISAYAMA_FSR":1113,"INTERSTROKE_ASN":222}
EAS_AIS_CONTROL_CONTRIBUTORS={"KCPS2_ASN**":12366,"BIOBANK_JAPAN*":159610,
                      "HISAYAMA_FSR":901,"INTERSTROKE_ASN":698}
MAX_XLSX_BYTES=28_000_000
MAX_OOXML_UNCOMPRESSED=95_000_000
NS="{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
REL="{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"
PKG="{http://schemas.openxmlformats.org/package/2006/relationships}"


def require(value,reason):
    if not value:raise ValueError(reason)


def download_workbook(fetch=urllib.request.urlopen):
    errors=[]
    for url in URLS:
        try:
            req=urllib.request.Request(url,headers={
                "User-Agent":"IS-Masters-GIGASTROKE-2022-primary-supp-source-QC/1.0",
                "Accept":"application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"})
            with fetch(req,timeout=75) as resp:
                require(getattr(resp,"status",200)==200,"XLSX publisher HTTP status")
                resolved=resp.geturl()
                host=urllib.parse.urlsplit(resolved).hostname
                require(urllib.parse.urlsplit(resolved).scheme=="https" and
                        host in {"static-content.springer.com","pmc.ncbi.nlm.nih.gov"},
                        "unexpected redirected publisher host")
                data=resp.read(MAX_XLSX_BYTES+1)
            require(5000<len(data)<=MAX_XLSX_BYTES,
                    "publisher source XLSX outside expected bounded size")
            require(zipfile.is_zipfile(io.BytesIO(data)),"publisher XLSX invalid OOXML")
            return data,{"original_official_url":url,"resolved_url":resolved,
                "publisher_spreadsheet_size_bytes":len(data),
                "sha256":hashlib.sha256(data).hexdigest(),
                "source_is_primary_2022_supplement":True,
                "2022_publisher_correction_fully_audited":False}
        except (ValueError,OSError,TimeoutError,zipfile.BadZipFile) as exc:
            errors.append(type(exc).__name__+":"+str(exc)[:150])
    raise ValueError("official GIGASTROKE Supplementary Table 1 unavailable: "+
                     "; ".join(errors))


def xlsx_table1_rows(content):
    with zipfile.ZipFile(io.BytesIO(content)) as z:
        names=set(z.namelist())
        require(all(".." not in p.split("/") for p in names),
                "unexpected path traversal in official XLSX")
        require(sum(info.file_size for info in z.infolist())<=MAX_OOXML_UNCOMPRESSED,
                "published xlsx uncompressed contents exceed bound")
        wb=ET.fromstring(z.read("xl/workbook.xml"))
        rel=ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
        rmap={r.attrib.get("Id"):r.attrib.get("Target")
              for r in rel.findall(PKG+"Relationship")}
        target=None
        for sheet in wb.findall(".//"+NS+"sheet"):
            if sheet.attrib.get("name")=="Suppl.Table1":
                target=rmap.get(sheet.attrib.get(REL+"id"))
                break
        require(target is not None,"original supplementary Table1 missing")
        target=target.lstrip("/")
        if not target.startswith("xl/"):
            target="xl/"+target
        target=target.replace("xl/xl/","xl/")
        require(target in names,"Suppl.Table1 relationship path missing")
        shared=[]
        if "xl/sharedStrings.xml" in names:
            root=ET.fromstring(z.read("xl/sharedStrings.xml"))
            shared=["".join(t.text or "" for t in si.iter(NS+"t"))
                    for si in root.findall(NS+"si")]
        tree=ET.fromstring(z.read(target))
        rows={}
        for r in tree.findall(".//"+NS+"sheetData/"+NS+"row"):
            try:rowidx=int(r.attrib["r"])
            except (TypeError,ValueError):continue
            if rowidx>300:continue
            item={}
            for cell in r.findall(NS+"c"):
                coord=cell.attrib.get("r","")
                match=re.fullmatch(r"([A-Z]+)([0-9]+)",coord)
                if not match or match.group(1) not in ("A","B","J","K","P"):continue
                typ=cell.attrib.get("t")
                if typ=="inlineStr":
                    value="".join(t.text or "" for t in cell.iter(NS+"t"))
                else:
                    node=cell.find(NS+"v")
                    value=node.text if node is not None else ""
                    if typ=="s" and value:
                        value=shared[int(value)]
                item[match.group(1)]=str(value or "").strip()
            if item:rows[rowidx]=item
        require(len(rows)>30,"Suppl.Table1 unexpectedly short")
        return rows


def cell_count(val):
    if not val:return None
    s=str(val).strip().replace(",","")
    return int(s) if re.fullmatch(r"[0-9]+(?:\.0+)?",s) and float(s).is_integer() else None


def extract_ancestry(rows,ancestry):
    start=[i for i,x in rows.items() if x.get("A","").upper()==ancestry]
    total=[i for i,x in rows.items()
           if x.get("A","").upper()==("TOTAL "+ancestry)]
    require(len(start)==1 and len(total)==1 and start[0]<total[0],
            f"missing original {ancestry} and TOTAL {ancestry} anchors")
    begin,end=start[0],total[0]
    participants=[]
    for line in range(begin,end):
        values=rows.get(line,{})
        name=values.get("B","")
        if not name:continue
        case=cell_count(values.get("K"))
        control=cell_count(values.get("P"))
        participants.append({
           "xlsx_row":line,"published_study_name":name,
           "AIS_cases_K":case,"discovery_controls_P":control,
           "other_stroke_cases_J":cell_count(values.get("J")),
           "primary_table_ancestry":ancestry,
           # In 2022 primary Table1 EAS row60 is incident-only not part
           # of global AIS meta-analysis, per archived cohort audit.
           "source_row_is_incident_only_separate_from_AIS_meta":
               ancestry=="EAS" and line==60,
        })
    require(len(participants)>3,f"no original {ancestry} study rows")
    total_cases=cell_count(rows[end].get("K"))
    total_controls=cell_count(rows[end].get("P"))
    require(total_cases is not None and total_controls is not None,
            f"original TOTAL {ancestry} AIS cases or controls unreadable")
    usable=[p for p in participants
            if not p["source_row_is_incident_only_separate_from_AIS_meta"]]
    calculated_cases=sum(p["AIS_cases_K"] or 0 for p in usable)
    calculated_ctrl=sum(p["discovery_controls_P"] or 0 for p in usable)
    exact=EXPECTED[ancestry]
    return {
        "original_ancestry_header_row":begin,
        "original_ancestry_total_row":end,
        "original_total_from_source_K_AIS_cases":total_cases,
        "original_total_from_source_P_discovery_controls":total_controls,
        "recomputed_usable_AIS_cases":calculated_cases,
        "recomputed_usable_discovery_controls":calculated_ctrl,
        "published_expected_cases":exact["cases"],
        "published_expected_controls":exact["controls"],
        "primary_total_matches_target_GWAS_case_control_counts":
            (total_cases,total_controls)==(exact["cases"],exact["controls"]),
        "source_individual_cohort_rows_reconcile_total":
            (calculated_cases,calculated_ctrl)==(total_cases,total_controls),
        "cohort_rows":participants,
        "n_source_cohort_rows_including_separate_incident":len(participants),
    }


def validate_eas(eas):
    src={p["published_study_name"]:p for p in eas["cohort_rows"]
         if not p["source_row_is_incident_only_separate_from_AIS_meta"]}
    contributors={name:p["AIS_cases_K"] for name,p in src.items()
                  if p["AIS_cases_K"] is not None and p["AIS_cases_K"]>0}
    require(contributors==EAS_AIS_CONTRIBUTORS,
            "original EAS AIS case contributor identities changed")
    for name,control in EAS_AIS_CONTROL_CONTRIBUTORS.items():
        require(src[name]["discovery_controls_P"]==control,
                "original EAS AIS contributor controls changed "+name)
    ckbb=[p for p in src.values() if p["published_study_name"].startswith("CKBB_ASN")]
    require(len(ckbb)==10 and all(p["AIS_cases_K"] is None for p in ckbb),
            "original CKBB 10 regional controls used as AIS cases")
    require(eas["source_individual_cohort_rows_reconcile_total"] and
            eas["primary_total_matches_target_GWAS_case_control_counts"],
            "original EAS AIS total not reconciled to discovery cohort rows")
    require("BIOBANK_JAPAN*" in src,
            "BBJ contributor to original EAS discovery missing")
    eas["official_BBJ_member_source_study"]="BIOBANK_JAPAN*"
    eas["original_BBJ_AIS_cases"]=src["BIOBANK_JAPAN*"]["AIS_cases_K"]
    eas["original_BBJ_shared_discovery_controls"]=src["BIOBANK_JAPAN*"]["discovery_controls_P"]
    eas["BBJ_is_EAS_discovery_contributing_cohort"]=True
    eas["BBJ_reuse_as_independent_EAS_replication_prohibited"]=True
    eas["sample_level_independence_of_other_datasets_proven"]=False
    eas["AIS_discovery_source_crosswalk_status"]="PASS_EAS_AIS_ORIGINAL_SUPP_COORT_SUMS_AND_BBJ"
    return eas


def evidence(manifest,yaml,official,content,source):
    require(manifest.get("schema")=="is-p3a-36-real-group-keys-from-legacy-393-v1"
            and len(manifest.get("groups",[]))==36
            and len({(g["region_id"],g["QTL_dataset"],g["gene_id"])
                     for g in manifest["groups"]})==36
            and len({g["gene_id"] for g in manifest["groups"]})==20
            and manifest.get("scientific_status")=="HOLD",
            "original 36 source groups invalid")
    require(yaml.get("schema")=="is-p3a-primary-gwas-yaml-QTL-allele-N-36-v1"
            and yaml.get("technical_status")==
              "PASS_PRIMARY_PUBLISHER_METADATA_AND_36_QTL_AN_SOURCE_READ"
            and yaml.get("scientific_status")=="HOLD"
            and len(yaml.get("groups",[]))==36
            and yaml.get("formal_coloc_abf_eligible")==0,
            "source original primary GWAS YAML and QTL variant AN audit invalid")
    keys={(x["QTL_dataset"],x["gene_id"]) for x in manifest["groups"]}
    require({(x["QTL_dataset"],x["gene_id"]) for x in yaml["groups"]}==keys,
            "original source QTL AN groups disagree with real manifest")
    for ancestry,accession in (("EAS","GCST90104545"),("EUR","GCST90104540")):
        x=yaml["primary_gwas_yaml"][ancestry]
        require(x["source_accession"]==accession
                and x["original_release_source_identity_verified"] is True
                and x["primary_GWAS_N_metadata_technically_complete"] is False
                and x["scientific_case_fraction_and_N_approved_for_coloc"] is False,
                "GWAS YAML missing-N scientific blocker cannot be overridden")
    rows=xlsx_table1_rows(content)
    actual={}
    for ancestry in ("EAS","EUR"):
        actual[ancestry]=extract_ancestry(rows,ancestry)
    actual["EAS"]=validate_eas(actual["EAS"])
    # EUR source Table1 cohort aggregation is checked independently; do
    # not approve it until rows and total match, not just a cited headline.
    eur=actual["EUR"]
    eur["source_cohort_exact_sum_verified"]=bool(
       eur["source_individual_cohort_rows_reconcile_total"] and
       eur["primary_total_matches_target_GWAS_case_control_counts"])
    for ancestry,a in (("EAS",actual["EAS"]),("EUR",eur)):
        a["model_sample_N_eff_approved"]=False
        a["case_fraction_approved_for_coloc"]=False
        a["study_original_sumstat_release_sample_mapping_checked"]=False
        a["individual_participant_overlap_across_cohorts_proven"]=False
    groups=[]
    nmap={(r["QTL_dataset"],r["gene_id"]):r for r in yaml["groups"]}
    for g in manifest["groups"]:
        n=nmap[(g["QTL_dataset"],g["gene_id"])]
        groups.append({
            "region_id":g["region_id"],
            "QTL_dataset":g["QTL_dataset"],
            "gene_id":g["gene_id"],
            "gene_symbol":g["gene_symbol"],
            "observed_QTL_nonpal_source_SNP_assay_rows":
                n["original_variant_allele_N_summary"]["count_valid_nonpal_original_QTL"],
            "QTL_variant_genotyped_an_over_2_technically_observed":True,
            "QTL_regression_sample_N_approved":False,
            "EAS_AIS_original_discovery_case_control_crosswalk":
                "PASS_COHORT_AGGREGATE_NOT_MODEL_N",
            "EAS_GIGASTROKE_contains_BBJ":True,
            "EUR_AIS_original_discovery_case_control_crosswalk":
                "PASS_COHORT_AGGREGATE_NOT_MODEL_N" if eur["source_cohort_exact_sum_verified"]
                else "HOLD_EUR_COHORT_ROW_RECONCILIATION",
            "QTL_nominal_full_tested_cis_approved":False,
            "BBJ_independent_EAS_replication_approved":False,
            "formal_coloc_abf_eligible":False,
            "formal_coloc_susie_eligible":False,
            "science_status":"HOLD_QTL_MODEL_N_AND_SUMSTAT_COHORT_APPROVAL",
        })
    status=("PASS_ORIGINAL_AIS_SUPPLEMENT_BOTH_ANCESTRY_TOTALS_AND_BBJ"
            if eur["source_cohort_exact_sum_verified"] else
            "PASS_ORIGINAL_EAS_AIS_BBJ_EUR_TOTAL_REVIEW_REQUIRED")
    require(len(groups)==36,"37th group inserted/missing")
    return {
        "schema":"is-p3a-primary-2022-gigastroke-AIS-supplement-BBJ-cohort-36-v1",
        "technical_status":status,"scientific_status":"HOLD",
        "official_publication":"https://doi.org/10.1038/s41586-022-05165-3",
        "official_original_supplement":source,
        "original_supplement_sheet":"Suppl.Table1",
        "original_case_definition_col":"K",
        "original_control_col":"P",
        "two_ancestries":actual,
        "n_real_QTL_assays":36,"n_distinct_genes":20,
        "n_distinct_IS_G0015_loci":1,
        "original_publisher_YAML_has_approved_N":False,
        "N_eff_or_case_fraction_imputed":False,
        "QTL_effective_N_approved":False,
        "n_coloc_abf_eligible":0,"n_coloc_abf_executed":0,
        "n_coloc_susie_eligible":0,"n_coloc_susie_executed":0,
        "groups":groups,
    }


def main():
    p=argparse.ArgumentParser()
    for field in ("manifest","qtl-gwas-yaml","out"):
        p.add_argument("--"+field,type=Path,required=True)
    args=p.parse_args()
    out={"schema":"is-p3a-primary-2022-gigastroke-AIS-supplement-BBJ-cohort-36-v1",
         "technical_status":"HOLD_PRIMARY_SOURCE_NOT_VERIFIED",
         "scientific_status":"HOLD","n_coloc_abf_eligible":0,
         "n_coloc_abf_executed":0,"n_coloc_susie_eligible":0,
         "n_coloc_susie_executed":0}
    try:
        manifest=json.loads(args.manifest.read_text(encoding="utf8"))
        yaml=json.loads(args.qtl_gwas_yaml.read_text(encoding="utf8"))
        content,provenance=download_workbook()
        out=evidence(manifest,yaml,None,content,provenance)
    except (ValueError,TypeError,KeyError,OSError,OverflowError,
            zipfile.BadZipFile,ET.ParseError,UnicodeError,TimeoutError) as exc:
        out["technical_status"]="HOLD_PRIMARY_SOURCE_OR_COHORT_ROW_FAILURE"
        out["error"]=str(exc)[:400]
        out["error_type"]=type(exc).__name__
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(out,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print(json.dumps({"technical_status":out["technical_status"],
       "ancestries":{k:{"total_cases":v.get("original_total_from_source_K_AIS_cases"),
              "total_controls":v.get("original_total_from_source_P_discovery_controls"),
              "sum_reconciles":v.get("source_individual_cohort_rows_reconcile_total")}
              for k,v in out.get("two_ancestries",{}).items()},
       "bbj_eas_cases":out.get("two_ancestries",{}).get("EAS",{}).get("original_BBJ_AIS_cases"),
       "coloc_abf_eligible":0,"error":out.get("error")},ensure_ascii=False))
    return 0 if out["technical_status"].startswith("PASS_") else 1


if __name__=="__main__":
    raise SystemExit(main())
