#!/usr/bin/env python3
"""Validate corrected primary GIGASTROKE Supplementary Tables 1 & 3 against EUR outcome manifest.

A single 756,910-byte, SHA256-pinned 2022 publisher XLSX is parsed read-only
using Python stdlib XLSX/OOXML. Publishes aggregate evidence, never cohort IDs.
N_eff scenarios are *illustrations* only; no LDSC or scientific GO is run.
"""
from __future__ import annotations
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import io
import json
import math
from pathlib import Path
import re
import urllib.request
from xml.etree import ElementTree as ET
from zipfile import ZipFile

SOURCE_URL=("https://static-content.springer.com/esm/"
            "art%3A10.1038%2Fs41586-022-05492-5/MediaObjects/"
            "41586_2022_5492_MOESM2_ESM.xlsx")
EXPECTED_SIZE=756_910
EXPECTED_SHA256="7c688f864c73ec2e9624bde9f53ac48a9880a6a08e59f5e1266690e243618ba4"
NS="{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
EXPECTED_SUBTYPES={
    "AIS":{"accession":"GCST90104540","cases":62_100,"total":1_296_908,"supp_col":"K","ldsc_col":"C"},
    "LAA":{"accession":"GCST90104542","cases":6_399,"total":1_241_207,"supp_col":"N","ldsc_col":"F"},
    "CE":{"accession":"GCST90104541","cases":10_804,"total":1_245_612,"supp_col":"L","ldsc_col":"E"},
    "SVD":{"accession":"GCST90104543","cases":6_811,"total":1_241_619,"supp_col":"M","ldsc_col":"D"},
}
EXPECTED_CONTROLS=1_234_808
MAX_ZIP_MEMBER=9_000_000
MAX_EXTRACTED_BYTES=24_000_000
SOURCE_COHORT_TABLE=("xl/worksheets/sheet2.xml","xl/worksheets/sheet4.xml")


def download_source(fetcher=urllib.request.urlopen)->bytes:
    req=urllib.request.Request(SOURCE_URL,headers={"User-Agent":"CKM-Stroke-GIGASTROKE-Cohort-Audit/1.0"})
    with fetcher(req,timeout=45) as response:
        if getattr(response,"status",200)!=200:
            raise ValueError("publisher source not HTTP 200")
        if not response.geturl().startswith("https://static-content.springer.com/"):
            raise ValueError("publisher source changed hostname")
        sz=response.headers.get("Content-Length","")
        if not sz.isdecimal() or int(sz)!=EXPECTED_SIZE:
            raise ValueError("corrected supplement content length changed")
        data=response.read(EXPECTED_SIZE+1)
    if len(data)!=EXPECTED_SIZE or hashlib.sha256(data).hexdigest()!=EXPECTED_SHA256:
        raise ValueError("corrected supplement SHA256/byte length invalid")
    return data


def _read_entry(z:ZipFile,name:str)->bytes:
    meta=z.getinfo(name)
    if meta.file_size>MAX_ZIP_MEMBER:
        raise ValueError("XLSX source XML member exceeds cap")
    payload=z.read(name)
    if len(payload)!=meta.file_size:
        raise ValueError("XLSX uncompressed member incomplete")
    return payload


def workbook_cells(blob:bytes):
    if len(blob)!=EXPECTED_SIZE or hashlib.sha256(blob).hexdigest()!=EXPECTED_SHA256:
        raise ValueError("corrected supplement pin mismatch")
    with ZipFile(io.BytesIO(blob)) as z:
        if any(x.file_size>MAX_ZIP_MEMBER for x in z.infolist()):
            raise ValueError("oversized XLSX member")
        if sum(x.file_size for x in z.infolist())>MAX_EXTRACTED_BYTES:
            raise ValueError("XLSX excessive total inflation")
        root=ET.fromstring(_read_entry(z,"xl/workbook.xml"))
        sheets=[s.get("name") for s in root.findall(f".//{NS}sheet")]
        if sheets[:4]!=["Table.of.content","Suppl.Table1","Suppl.Table2","Suppl.Table3"]:
            raise ValueError("corrected publisher supplement sheet order changed")
        shared_root=ET.fromstring(_read_entry(z,"xl/sharedStrings.xml"))
        shared=["".join(x.text or "" for x in item.iter(f"{NS}t"))
                for item in shared_root.findall(f"{NS}si")]
        tables=[]
        for path in SOURCE_COHORT_TABLE:
            table=ET.fromstring(_read_entry(z,path))
            rows={}
            for r in table.findall(f".//{NS}sheetData/{NS}row"):
                rownum=int(r.get("r","0"))
                cells={}
                for c in r.findall(f"{NS}c"):
                    match=re.fullmatch(r"([A-Z]+)[0-9]+",c.get("r",""))
                    if not match:
                        raise ValueError("invalid spreadsheet coordinate")
                    v=c.find(f"{NS}v")
                    value=v.text if v is not None else ""
                    if c.get("t")=="s" and value:
                        i=int(value)
                        if not 0<=i<len(shared):
                            raise ValueError("invalid XLSX shared string index")
                        value=shared[i]
                    elif c.get("t")=="inlineStr":
                        inline=c.find(f"{NS}is")
                        value="".join(x.text or "" for x in inline.iter(f"{NS}t")) if inline is not None else ""
                    if value is not None and value!="":
                        cells[match.group(1)]=str(value)
                rows[rownum]=cells
            tables.append(rows)
    return tuple(tables)


def _integer(x,context):
    if not re.fullmatch(r"[0-9]+",str(x)):
        raise ValueError("invalid integer "+context)
    return int(x)


def cohortwise_neff_scenario(table1:dict,case_col:str,expected_cases:int):
    """Illustration only: sum case-control harmonic N over rows with both counts.

    This does not account for per-SNP missingness or meta-analysis weights.
    Rows lacking controls (e.g. HVH2) are flagged as MISSING, not zero.
    """
    observed_cases=0
    supported_cases=0
    supported_cohorts=0
    missing_controls=[]
    total=0.0
    for rowno in range(5,49):
        row=table1.get(rowno,{})
        label=row.get("B","")
        if not label:
            continue
        value=row.get(case_col,"")
        if not str(value).isdigit():
            continue
        cases=int(value)
        observed_cases+=cases
        ctrl=row.get("P","")
        if not str(ctrl).isdigit() or int(ctrl)<=0:
            missing_controls.append({"row":rowno,"cohort":label,"cases":cases})
            continue
        controls=int(ctrl)
        if cases<=0:
            continue
        supported_cohorts+=1
        supported_cases+=cases
        total+=4*cases*controls/(cases+controls)
    if observed_cases!=expected_cases:
        raise ValueError("per-cohort reported case sum differs from original EUR total")
    return {
        "reported_case_total_from_rows":observed_cases,
        "cases_with_reported_cohort_controls":supported_cases,
        "cases_missing_cohort_control_count":observed_cases-supported_cases,
        "cohort_count_with_cases_and_reported_controls":supported_cohorts,
        "missing_cohort_control_rows":missing_controls,
        "cohortwise_sum_harmonic_neff_illustrative":round(total,3),
        "not_per_variant_N_and_not_approved_for_LDSC":True
    }


def source_evidence(blob:bytes,manifest_text:str):
    table1,table3=workbook_cells(blob)
    header=table1.get(4,{})
    if [header.get(k) for k in ("J","K","L","M","N","O","P")]!=[
        "AS","AIS","CES","SVS","LAS","ICH","Controls"
    ]:
        raise ValueError("Table1 combined incidence/prevalence column mapping changed")
    if table1.get(3,{}).get("J")!="Incident & prevalent combined":
        raise ValueError("wrong cohort sample-count table panel")
    cohorts=[table1.get(i,{}) for i in range(5,49) if table1.get(i,{}).get("B")]
    ukb=[r for r in cohorts if r.get("B")=="UKBiobank"]
    if len(cohorts)!=44 or len(ukb)!=1 or table1.get(15,{}).get("B")!="UKBiobank":
        raise ValueError("original EUR cohort composition differs")
    ukb=ukb[0]
    if (_integer(ukb.get("K"),"UKB AIS cases")!=2211 or
        _integer(ukb.get("P"),"UKB controls")!=369419 or
        [ukb.get(k) for k in ("L","M","N")]!=["NA","NA","NA"]):
        raise ValueError("original UK Biobank subtype composition changed")
    total=table1.get(49,{})
    if total.get("A")!="TOTAL EUR":
        raise ValueError("original EUR summary row missing")
    controls=_integer(total.get("P"),"EUR controls")
    if controls!=EXPECTED_CONTROLS:
        raise ValueError("original EUR control total changed")
    if (table3.get(70,{}).get("C")!="AIS" or
        table3.get(79,{}).get("C")!="AIS" or
        table3.get(71,{}).get("A")!="EUR" or
        table3.get(80,{}).get("A")!="EUR"):
        raise ValueError("original Table3 LDSC EUR schema changed")
    metadata=list(csv.DictReader(io.StringIO(manifest_text),delimiter="\t"))
    if len(metadata)!=4 or {x.get("phenotype") for x in metadata}!=set(EXPECTED_SUBTYPES):
        raise ValueError("manifest missing or duplicated phenotype")
    results=[]
    for phenotype,ref in EXPECTED_SUBTYPES.items():
        row=next(x for x in metadata if x.get("phenotype")==phenotype)
        if row.get("accession")!=ref["accession"] or row.get("build")!="GRCh37":
            raise ValueError("official publisher accession mapping mismatch "+phenotype)
        ncases=_integer(total.get(ref["supp_col"]),"EUR "+phenotype+" cases")
        if (ncases!=ref["cases"] or
            ncases+controls!=ref["total"] or
            _integer(row.get("sample_total"),"GWAS YAML total")!=ref["total"]):
            raise ValueError("official EUR case/controls not aligned with GWAS YAML "+phenotype)
        col=ref["ldsc_col"]
        intercept=float(table3[71][col])
        match=re.fullmatch(r"\s*([0-9.eE+-]+)\s*\(([0-9.eE+-]+)\)\s*",table3[80][col])
        if not match:
            raise ValueError("original Cov-LDSC h2(SE) not parseable "+phenotype)
        h2,stderr=(float(match.group(1)),float(match.group(2)))
        if not(all(map(math.isfinite,(intercept,h2,stderr))) and 0<intercept<2 and
               0<h2<1 and stderr>0):
            raise ValueError("invalid original Table3 LDSC metrics")
        arithmetic_neff=4*ncases*controls/(ncases+controls)
        cohort_scenario=cohortwise_neff_scenario(table1,ref["supp_col"],ncases)
        results.append({
            "phenotype":phenotype,"accession":ref["accession"],
            "EUR_case_count_from_original_table1":ncases,
            "EUR_control_count_from_original_table1":controls,
            "EUR_total_N_from_original_table1":ncases+controls,
            "illustrative_harmonic_neff":round(arithmetic_neff,3),
            "illustrative_cohortwise_neff":cohort_scenario,
            "illustrative_neff_not_per_snp_or_approved":True,
            "original_CovLDSC_intercept_from_table3B":intercept,
            "original_CovLDSC_observed_h2_from_table3C":h2,
            "original_CovLDSC_observed_h2_SE":stderr,
            "original_CovLDSC_h2_SE_ratio":round(h2/stderr,3),
            "this_pipeline_ldsc_executed":False
        })
    return {
        "schema":"ckm-gigastroke-corrected-primary-cohort-ldsc-gate-v1",
        "publisher_source_url":SOURCE_URL,"publisher_source_sha256":EXPECTED_SHA256,
        "corrected_supplement_size_bytes":EXPECTED_SIZE,
        "publisher_original_table1_eur_cohort_count":len(cohorts),
        "uk_biobank":{
            "row":15,"study_name":"UKBiobank",
            "AIS_case_count":2211,"UKB_reported_control_count":369419,
            "CES_SVS_LAS":"NA in UKBiobank row; absence of subtype column not cohort-independence proof",
            "incident_and_prevalent_combined":True,
            "AIS_cases_pct_of_reported_EUR_total":round(100*2211/62100,4),
            "control_pct_of_reported_EUR_total":round(100*369419/controls,4),
            "overlap_in_BMI_actual_count_known":False,
            "overlap_status":"POSSIBLE_CONFIRMED_SHARED_COHORT_BUT_INDIVIDUAL_COUNT_UNKNOWN",
        },
        "outcomes":results,
        "accession_discrepancy_warning":"Official EBI YAML CE=GCST90104541 and LAS=GCST90104542; some secondary tables reverse this.",
        "sample_overlap_estimated_count":None,
        "scientific_status":"HOLD","ready_for_ldsc":False,
        "LDSC_h2_rg_executed":False,"original_full_sumstats_downloaded":False,
        "technical_status":"PASS_PRIMARY_CORRECTED_COHORT_AND_TABLE3_ONLY"
    }


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--manifest",type=Path,default=Path("config/ckm_stroke/gigastroke_eur_official.tsv"))
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--local-source",type=Path,default=None,
                   help="Use a previously SHA-pinned ignored local publisher supplement")
    args=p.parse_args()
    result={"technical_status":"HOLD_SOURCE_OR_COHORT_QC",
            "scientific_status":"HOLD","ready_for_ldsc":False,
            "LDSC_h2_rg_executed":False}
    try:
        raw=args.local_source.read_bytes() if args.local_source else download_source()
        result=source_evidence(raw,args.manifest.read_text(encoding="utf-8-sig"))
        result["checked_at_utc"]=datetime.now(timezone.utc).isoformat()
    except Exception as exc:
        result["error_type"]=type(exc).__name__
        result["error_detail"]=str(exc)[:220]
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps({"status":result["technical_status"],
                      "uk_biobank":result.get("uk_biobank"),
                      "outcomes":result.get("outcomes"),
                      "error":result.get("error_detail"),
                      "scientific_status":result["scientific_status"]},ensure_ascii=False))
    if result["technical_status"]!="PASS_PRIMARY_CORRECTED_COHORT_AND_TABLE3_ONLY":
        raise SystemExit(1)


if __name__=="__main__":
    main()
