#!/usr/bin/env python3
"""FinnGen R13 official endpoint definition, GWAS manifest and sample-N review.

Only public R13 release metadata is read.  Endpoint *definition* is distinct
from a publicly listed R13 GWAS *phenotype*; a clinical case description and
case-control totals are distinct from SNP-specific N_eff.  Fail closed when
published source bytes/schemas cannot establish a claimed result.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
from pathlib import Path
import re
import urllib.request
from urllib.parse import urlsplit
import xml.etree.ElementTree as ET
import zipfile

from scripts.ckm_stroke.audit_gigastroke_supplement_cohort_neff import (
    EXPECTED_SHA256 as COHORT_SOURCE_SHA,
    EXPECTED_SUBTYPES, EXPECTED_CONTROLS,
)
from scripts.ckm_stroke.audit_alternative_stroke_gwas_neff_sources import (
    FINN_R13_MANIFEST, classify_finn_manifest,
)

ENDPOINT_URL="https://www.finngen.fi/sites/default/files/inline-files/FINNGEN_ENDPOINTS_DF13_Final_2025-08-14_public.xlsx"
CONTROL_URL="https://www.finngen.fi/sites/default/files/inline-files/Controls_FINNGEN_ENDPOINTS_DF13_Final_2025-08-14_public.xlsx"
OFFICIAL_CLINICAL_PAGE="https://www.finngen.fi/en/researchers/clinical-endpoints"
R13_MANIFEST_SHA="c9caafa9b98ee5ef050a705766bd1c451fc2154fb8e81b069022012a2f75ae20"
DF13_ENDPOINT_SHA256="0004624f727b96b4a56c6fac75ff182eb8bc3daf78cdf690bdad276b31fca6c7"
DF13_CONTROL_SHA256="7cbbf9072296e34f6bde44b2958c7520c7caedeade924d941e5bbc4fa89cc29b"
SOURCE_CODES=("I9_STR_EXH","I9_STR","I9_STR_EMBOLIC",
              "I9_STR_EMBOLIC_WITH_FLIMMER","I9_ICH","I9_SAH")
INSPECTED_FIELDS=("NAME","LONGNAME","DESCRIPTION","PHENOSTRING","HD_ICD_10","COD_ICD_10",
                  "HD_ICD_9","HD_ICD_8","COD_ICD_9","COD_ICD_8",
                  "INCLUDE","CONDITIONS","CONTROL_EXCLUDE","CONTROL_CONDITIONS",
                  "PRE_CONDITIONS","MAINONLY","SEX")
NS="{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
CAP_DOC=13_000_000
CAP_XML_TOTAL=115_000_000
CAP_PER_XML=38_000_000

def official_bytes(url,limit,opener=urllib.request.urlopen):
    u=urlsplit(url)
    if u.scheme!="https" or u.hostname not in ("www.finngen.fi","storage.googleapis.com"):
        raise ValueError("untrusted official FinnGen release hostname")
    req=urllib.request.Request(url,headers={"User-Agent":"CKM-Stroke-FinnGen-R13-Endpoint-Definitions/1.0"})
    with opener(req,timeout=65) as response:
        actual=urlsplit(response.geturl())
        if response.status!=200 or actual.scheme!="https" or actual.hostname not in (
            "www.finngen.fi","storage.googleapis.com"):
            raise ValueError("official FinnGen source HTTP or redirected host invalid")
        data=response.read(limit+1)
        if not data or len(data)>limit:
            raise ValueError("FinnGen source absent or exceeds byte safety cap")
    return data,hashlib.sha256(data).hexdigest()

def manifest_review(data):
    h=hashlib.sha256(data).hexdigest()
    if h!=R13_MANIFEST_SHA:
        raise ValueError("R13 actual public GWAS manifest SHA256 drift")
    try:
        raw=data.decode("utf-8-sig","strict")
    except UnicodeError as e:
        raise ValueError("R13 official manifest text encoding invalid") from e
    overview=classify_finn_manifest(raw)
    if overview["manifest_entries"]!=2754:
        raise ValueError("R13 public GWAS manifest expected 2754 source endpoints")
    reader=csv.DictReader(io.StringIO(raw),delimiter="\t")
    cols=reader.fieldnames
    if cols is None or len(cols)!=len(set(cols)):
        raise ValueError("R13 manifest header missing or duplicated")
    if "phenocode" not in cols:
        raise ValueError("R13 GWAS manifest missing phenocode")
    selected={}
    for row in reader:
        ph=row.get("phenocode")
        if ph in SOURCE_CODES:
            if ph in selected:
                raise ValueError("duplicate R13 GWAS phenocode")
            selected[ph]={k:row[k][:400] if row.get(k) is not None else None for k in cols}
    return {"official_R13_manifest_sha256":h,
            "official_manifest_columns":cols,
            "official_manifest_entries":overview["manifest_entries"],
            "GWAS_manifest_present_endpoints":sorted(selected),
            "requested_endpoint_rows":selected,
            "I9_STR_EXH_in_R13_GWAS_manifest":"I9_STR_EXH" in selected,
            "I9_STR_is_not_ischaemic_specific_by_name":True,
            "manifest_N_case_column_candidates":[x for x in cols if re.search(
                r"(?:^|_)(?:n|num)_?cases?$|(?:^|_)(?:n|num)_?controls?$|case_?count|control_?count",x,re.I)],
            "manifest_is_not_per_variant_Neff":True}

def _colnum(cellref):
    letters=re.match(r"([A-Z]{1,3})[0-9]+",cellref or "")
    if not letters:
        raise ValueError("invalid OOXML cell coordinate")
    v=0
    for ch in letters.group(1):
        v=v*26+ord(ch)-64
    return v

def _workbook_rows(payload):
    if not payload.startswith(b"PK"):
        raise ValueError("FinnGen official endpoint workbook not ZIP")
    z=zipfile.ZipFile(io.BytesIO(payload))
    names=z.namelist()
    if len(names)>220 or len(names)!=len(set(names)):
        raise ValueError("workbook member count or duplicated entries exceeds cap")
    total=0
    for info in z.infolist():
        name=info.filename
        if (name.startswith("/") or "\\" in name or ".." in name.split("/")
            or info.file_size>CAP_PER_XML):
            raise ValueError("unsafe or oversized FinnGen OOXML member")
        total+=info.file_size
    if total>CAP_XML_TOTAL:
        raise ValueError("FinnGen workbook decompression cap exceeded")
    shared=[]
    if "xl/sharedStrings.xml" in names:
        root=ET.fromstring(z.read("xl/sharedStrings.xml"))
        for si in root.findall(f"{NS}si"):
            shared.append("".join(el.text or "" for el in si.iter(f"{NS}t"))[:20000])
    sheets=[]
    for name in names:
        if not re.fullmatch(r"xl/worksheets/sheet[0-9]+\.xml",name):
            continue
        root=ET.fromstring(z.read(name))
        rows=[]
        for row in root.findall(f".//{NS}sheetData/{NS}row"):
            if len(rows)>60000:
                raise ValueError("too many source definition rows")
            values={}
            for cell in row.findall(f"{NS}c"):
                index=_colnum(cell.get("r"))
                if index>300:
                    raise ValueError("spreadsheet cell column cap exceeded")
                value=cell.find(f"{NS}v")
                raw=value.text if value is not None and value.text is not None else ""
                if cell.get("t")=="s" and raw:
                    idx=int(raw)
                    if not 0<=idx<len(shared):
                        raise ValueError("source OOXML shared string index invalid")
                    raw=shared[idx]
                elif cell.get("t")=="inlineStr":
                    inline=cell.find(f"{NS}is")
                    raw="".join(t.text or "" for t in inline.iter(f"{NS}t")) if inline is not None else ""
                if raw:
                    values[index]=raw.strip()[:1200]
            rows.append((int(row.get("r","0")),values))
        sheets.append((name,rows))
    if not sheets:
        raise ValueError("FinnGen official workbook has no worksheet")
    return sheets

def extract_endpoint_rows(payload,source_type):
    results={key:[] for key in SOURCE_CODES}
    sheets=_workbook_rows(payload)
    inventory=[]
    for name,rows in sheets:
        header=None
        header_line=None
        for rn,vals in rows[:25]:
            if "NAME" in [v.upper() for v in vals.values()]:
                header={key:val.strip().upper() for key,val in vals.items()}
                header_line=rn
                break
        inventory.append({"worksheet_xml":name,"rows":len(rows),
                          "header_found":bool(header),
                          "header_labels":list(header.values())[:36] if header else []})
        if not header:
            continue
        name_cols=[k for k,v in header.items() if v=="NAME"]
        if len(name_cols)!=1:
            raise ValueError("ambiguous NAME field in source endpoint workbook")
        name_col=name_cols[0]
        for rn,vals in rows:
            if rn<=header_line:
                continue
            phenocode=vals.get(name_col)
            if phenocode in results:
                selected={col:vals.get(k) for k,col in header.items()
                          if (col in INSPECTED_FIELDS or
                              col.startswith("HD_ICD") or col.startswith("COD_ICD")) and vals.get(k)}
                results[phenocode].append({"worksheet":name,"row":rn,"source_fields":selected,
                                           "source_columns_exist":list(header.values())[:50]})
    return {"workbook_type":source_type,
            "worksheet_inventory":inventory,
            "endpoint_definition_rows":results,
            "endpoints_with_unique_keyed_definition_rows":sorted(k for k,v in results.items() if len(v)==1),
            "endpoint_definition_rows_not_necessarily_GWAS_result_rows":True}

def validated_cohort_sensitivity(source):
    if (source.get("technical_status")!="PASS_PRIMARY_CORRECTED_COHORT_AND_TABLE3_ONLY"
        or source.get("publisher_source_sha256")!=COHORT_SOURCE_SHA
        or source.get("publisher_original_table1_eur_cohort_count")!=44
        or source.get("ready_for_ldsc") is not False
        or source.get("LDSC_h2_rg_executed") is not False):
        raise ValueError("GIGASTROKE publisher-supplement cohort evidence mismatch")
    outcomes=source.get("outcomes",[])
    if not isinstance(outcomes,list) or len(outcomes)!=4:
        raise ValueError("four phenotype supplement evidence required")
    plans=[]
    for r in outcomes:
        code=r["phenotype"]
        expected=EXPECTED_SUBTYPES.get(code)
        if expected is None or r.get("accession")!=expected["accession"]:
            raise ValueError("unrecognized published stroke outcome sample counts")
        cases=r["EUR_case_count_from_original_table1"]
        controls=r["EUR_control_count_from_original_table1"]
        if cases!=expected["cases"] or controls!=EXPECTED_CONTROLS:
            raise ValueError("published original EUR case/control table changed")
        pooled=r["illustrative_harmonic_neff"]
        detail=r["illustrative_cohortwise_neff"]
        cohort=detail["cohortwise_sum_harmonic_neff_illustrative"]
        missing=detail["cases_missing_cohort_control_count"]
        if (not math.isfinite(pooled) or not math.isfinite(cohort) or
            not 0<cohort<=pooled or type(missing) is not int or missing<0 or
            detail.get("not_per_variant_N_and_not_approved_for_LDSC") is not True):
            raise ValueError("invalid / unapproved N_eff scenario")
        plans.append({"phenotype":code,"accession":expected["accession"],
                      "published_cases":cases,"published_controls":controls,
                      "pooled_harmonic_scenario_Neff_NOT_APPROVED":pooled,
                      "known_cohort_harmonic_scenario_Neff_NOT_APPROVED":cohort,
                      "relative_scenario_difference_percent":round(100*(pooled-cohort)/pooled,5),
                      "published_cases_missing_cohort_control_number":missing,
                      "neither_value_is_variant_specific_Neff":True,
                      "not_validated_lower_upper_bounds_for_true_meta_Neff":True,
                      "primary_LDSC_allowed":False})
    if {x["phenotype"] for x in plans}!={"AIS","CE","LAA","SVD"}:
        raise ValueError("missing original four subtypes")
    return {
        "scenario_type":"PRESPECIFIED_INPUT_MAGNITUDE_COMPARISON_ONLY_NOT_LDSC",
        "GIGASTROKE_44_cohort_publisher_SHA256_verified":True,
        "cohort_and_pooled_scenarios_are_not_true_per_variant_Neff_bounds":True,
        "per_SNP_contributing_cohorts_and_meta_weights_unknown":True,
        "cases_controls_by_SNP_known":False,
        "never_approve_primary_or_sensitivity_LDSC_from_this_number_comparison":True,
        "scenarios":plans,
        "actual_LDSC_h2_intercept_rg_computed":False,
        "scientific_status":"HOLD",
    }

def review(manifest,definitions,control_definitions,cohort):
    listing=manifest_review(manifest)
    endpoint=extract_endpoint_rows(definitions,"FinnGen_DF13_case_definitions")
    controls=extract_endpoint_rows(control_definitions,"FinnGen_DF13_control_definitions")
    for k in ("I9_STR","I9_STR_EMBOLIC","I9_STR_EMBOLIC_WITH_FLIMMER"):
        if k not in listing["GWAS_manifest_present_endpoints"]:
            raise ValueError("prior R13 GWAS phenotype selection changed: "+k)
        entry=listing["requested_endpoint_rows"][k]
        for field in ("num_cases","num_controls"):
            raw=entry.get(field)
            if not isinstance(raw,str) or not raw.isdecimal() or int(raw)<=0:
                raise ValueError("official R13 endpoint-specific case/control count missing: "+k)
        if int(entry["num_cases"])+int(entry["num_controls"])>500186:
            raise ValueError("R13 phenotype cohort count exceeds official total participants")
    if listing["I9_STR_EXH_in_R13_GWAS_manifest"]:
        raise ValueError("prior no-R13-ischemic-specific source evidence changed; review primary outcome")
    srcdefs=endpoint["endpoint_definition_rows"]
    ctrl=controls["endpoint_definition_rows"]
    if any(len(srcdefs.get(k,[]))!=1 or len(ctrl.get(k,[]))!=1 for k in SOURCE_CODES):
        raise ValueError("missing unique source-keyed DF13 case/control definition")
    def code_fields(name):
        return srcdefs[name][0]["source_fields"]
    broad=code_fields("I9_STR").get("HD_ICD_10","").split("|")
    ischemic=code_fields("I9_STR_EXH").get("HD_ICD_10","").split("|")
    embolic=code_fields("I9_STR_EMBOLIC").get("HD_ICD_10","")
    af=code_fields("I9_STR_EMBOLIC_WITH_FLIMMER")
    af_ctrl=ctrl["I9_STR_EMBOLIC_WITH_FLIMMER"][0]["source_fields"]
    if (not {"I61","I63","I64"}.issubset(set(broad)) or
        "I61" in ischemic or not {"I63","I64"}.issubset(set(ischemic)) or
        embolic not in ("I631","I63.1") or
        af.get("INCLUDE")!="I9_STR_EMBOLIC" or
        af.get("CONDITIONS")!="I9_AF" or
        af_ctrl.get("CONTROL_CONDITIONS")!="!I9_AF"):
        raise ValueError("official DF13 clinical ischemic versus hemorrhagic/embolic definition drift")
    # Clinical endpoint may exist in DF13 case definition file while not
    # being listed among released GWAS; do not infer data availability.
    targets=[]
    for code in SOURCE_CODES:
        d=endpoint["endpoint_definition_rows"].get(code,[])
        cd=controls["endpoint_definition_rows"].get(code,[])
        mr=listing["requested_endpoint_rows"].get(code)
        cls=("AIS_CLINICAL_CANDIDATE" if code=="I9_STR_EXH" else
             "BROAD_STROKE_MAY_INCLUDE_HEMORRHAGE" if code=="I9_STR" else
             "EMBOLIC_CLINICAL_SUBTYPE_NOT_TOAST_EQUIVALENT" if code.startswith("I9_STR_EMBOLIC") else
             "EXCLUDE_HEMORRHAGE_FROM_AIS")
        targets.append({
            "phenocode":code,"interpretation_category":cls,
            "present_in_public_R13_GWAS_manifest":mr is not None,
            "case_definition_rows_found":len(d),
            "control_definition_rows_found":len(cd),
            "case_definition":d[0]["source_fields"] if len(d)==1 else None,
            "control_definition":cd[0]["source_fields"] if len(cd)==1 else None,
            "manifest_public_data_row":mr,
            "released_R13_GWAS_endpoint_case_count":int(mr["num_cases"]) if mr is not None else None,
            "released_R13_GWAS_endpoint_control_count":int(mr["num_controls"]) if mr is not None else None,
            "original_manifest_counts_are_not_variantwise_N_eff":True,
            "case_or_control_count_source_quality":(
                "MANIFEST_NUMBER_ONLY_IF_EXPLICIT_SOURCE_FIELDS" if mr is not None else
                "NOT_AVAILABLE_NO_R13_GWAS_MANIFEST_ENTRY"),
            "per_variant_effective_N_present":False,
            "clinical_equal_to_TOAST_AIS_or_CE_proven":False,
        })
    return {
      "schema":"ckm-FinnGen-R13-official-endpoint-definition-versus-manifest-and-Neff-scenarios-v1",
      "technical_status":"PASS_FG_R13_PUBLISHED_CASE_CONTROL_DEFINITIONS_AND_GWAS_AVAILABILITY_AUDIT",
      "manifest":listing,
      "case_definitions":{"worksheet_inventory":endpoint["worksheet_inventory"]},
      "control_definitions":{"worksheet_inventory":controls["worksheet_inventory"]},
      "endpoints":targets,
      "ischemic_stroke_I9_STR_EXH_preferred_over_I9_STR_on_clinical_label_only":True,
      "I9_STR_EXH_GWAS_data_unavailable_if_not_listed_in_R13_manifest":not listing["I9_STR_EXH_in_R13_GWAS_manifest"],
      "public_manifest_n_cases_and_n_controls_are_not_SNP_specific":True,
      "previously_reported_three_candidate_list_was_incomplete":True,
      "GIGASTROKE_meta_N_eff_input_sensitivity":validated_cohort_sensitivity(cohort),
      "scientific_status":"HOLD","ready_for_ldsc":False,
      "LDSC_h2_intercept_rg_MR_coloc_performed":False,
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--original-cohort-evidence",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    a=ap.parse_args()
    try:
        manifest,msha=official_bytes(FINN_R13_MANIFEST,6_000_000)
        df13,dsha=official_bytes(ENDPOINT_URL,CAP_DOC)
        cdf13,csha=official_bytes(CONTROL_URL,CAP_DOC)
        previous=json.loads(a.original_cohort_evidence.read_text("utf-8"))
        if dsha!=DF13_ENDPOINT_SHA256 or csha!=DF13_CONTROL_SHA256:
            raise ValueError("original official DF13 clinical definition workbook SHA256 drift")
        result=review(manifest,df13,cdf13,previous)
        result["official_sources"]={"public_GWAS_manifest":FINN_R13_MANIFEST,
           "DF13_endpoint_definitions":ENDPOINT_URL,
           "DF13_control_definitions":CONTROL_URL,
           "official_definition_website":OFFICIAL_CLINICAL_PAGE}
        result["official_full_source_SHA256"]={
            "manifest":msha,"DF13_endpoint_definitions":dsha,
            "DF13_control_definitions":csha}
        result["no_raw_person_or_genotype_data_fetched"]=True
    except Exception as exc:
        result={"technical_status":"HOLD_FG13_CASE_CONTROL_SOURCES_UNVERIFIED",
                "error_type":type(exc).__name__,"error_detail":str(exc)[:420],
                "scientific_status":"HOLD","ready_for_ldsc":False,
                "LDSC_h2_intercept_rg_MR_coloc_performed":False}
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps({
        "status":result["technical_status"],
        "manifest_endpoints":result.get("manifest",{}).get("official_manifest_entries"),
        "listed_stroke_codes":result.get("manifest",{}).get("GWAS_manifest_present_endpoints"),
        "DF13_definition_rows":{x["phenocode"]:x["case_definition_rows_found"] for x in result.get("endpoints",[])},
        "case_control_sources":{x["phenocode"]:{
           "n_cases":x["released_R13_GWAS_endpoint_case_count"],
           "n_controls":x["released_R13_GWAS_endpoint_control_count"],
           "status":x["case_or_control_count_source_quality"]}
           for x in result.get("endpoints",[])},
        "sensitivity_n_scenarios":len(result.get("GIGASTROKE_meta_N_eff_input_sensitivity",{}).get("scenarios",[])),
        "error":result.get("error_detail"),"science":result["scientific_status"]
    }))
    if result["technical_status"]!="PASS_FG_R13_PUBLISHED_CASE_CONTROL_DEFINITIONS_AND_GWAS_AVAILABILITY_AUDIT":
        raise SystemExit(1)

if __name__=="__main__":
    main()
