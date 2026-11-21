#!/usr/bin/env python3
"""First-party original-cohort *name* concordance, not individual sample overlap.

Sources:
Locke 2015 Supplementary Table 1 via EMBL-EBI BioStudies S-EPMC4382211.
GIGASTROKE corrected 2022 Supplementary Table 1 via Springer Nature.
SHA256-pinned first-party published XLSX, fail-closed read-only ZIP/OOXML parser.
A shared recruitment cohort name is NEVER evidence of shared participant identities.
"""
from __future__ import annotations
import argparse
from collections import Counter
from datetime import datetime,timezone
import hashlib
import io
import json
from pathlib import Path
import re
from urllib.parse import urlsplit
import urllib.request
from xml.etree import ElementTree as ET
from zipfile import ZipFile

SOURCE={
 "LOCKE":{
    "url":"https://ftp.ebi.ac.uk/pub/databases/biostudies/S-EPMC/211/S-EPMC4382211/Files/NIHMS668049-supplement-2.xlsx",
    "sha256":"dadd561bba67dbda9cbe33331a2fea195e638f0951ace75706149d137d2171f8",
    "bytes":1_626_884,"sheet":"xl/worksheets/sheet1.xml",
    "name":"S1-Study QC","host":"ftp.ebi.ac.uk",
 },
 "GIGASTROKE":{
    "url":"https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-022-05492-5/MediaObjects/41586_2022_5492_MOESM2_ESM.xlsx",
    "sha256":"7c688f864c73ec2e9624bde9f53ac48a9880a6a08e59f5e1266690e243618ba4",
    "bytes":756_910,"sheet":"xl/worksheets/sheet2.xml",
    "name":"Suppl.Table1","host":"static-content.springer.com",
 },
}
NS="{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
REL="{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"
MAX_INFLATED=22_000_000
MAX_XML_ENTRY=10_000_000
MAX_FILES=160

# Alias is permitted ONLY when the actual first-party workbook contains each name.
# This is a classification of cohort NAME/FAMILY reuse, NOT shared individuals.
ALIASES={
 "ARIC_EA*":("ARIC","COHORT_LABEL_SHARED","ARIC EUR label"),
 "CHS_EA*":("CHS","COHORT_LABEL_SHARED","Cardiovascular Health Study EUR label"),
 "FHS_EA*":("FRAM","COHORT_LABEL_SHARED","Framingham Heart Study / FRAM"),
 "RS1_EA*":("RS-I","COHORT_LABEL_SHARED","Rotterdam Study I"),
 "RS3_EA**":("RSIII","COHORT_LABEL_SHARED","Rotterdam Study III"),
 "WGHS_EA*":("WGHS","COHORT_LABEL_SHARED","Women's Genome Health Study"),
 "SHIP*":("SHIP","COHORT_LABEL_SHARED","Study of Health in Pomerania"),
 "HABC*":("Health ABC","COHORT_LABEL_SHARED","Health Aging and Body Composition"),
 "AGES":("AGES","COHORT_LABEL_SHARED","Age Gene/Environment Susceptibility Reykjavik"),
 "DECODE":("deCODE","COHORT_LABEL_SHARED","deCODE source cohort label"),
 "PROSPERlogistic_EA*":("PROSPER/PHASE","STUDY_FAMILY_SHARED","PROSPER/PHASE analysis subgroup uncertain"),
 "HUNT**":("HUNT 2","STUDY_FAMILY_SHARED","HUNT waves/participants not verified"),
 "RS2_EA*":("RS-I","STUDY_FAMILY_SHARED","Different Rotterdam Study waves; cannot infer subject equality"),
 "EPICCVD*":("EPIC-Obesity Study","STUDY_FAMILY_SHARED","EPIC umbrella with different subcohorts"),
 "EstBB**":("EGCUT","STUDY_FAMILY_SHARED","Estonian Biobank/EGCUT umbrella"),
 "UTARTU-EstBB**":("EGCUT-370","STUDY_FAMILY_SHARED","Tartu Estonian Biobank/EGCUT umbrella"),
}
CATEGORY={"COHORT_LABEL_SHARED","STUDY_FAMILY_SHARED"}
EXPECTED_EUR={
 "AIS":62100, "CE":10804,"SVD":6811,"LAA":6399,"controls":1234808
}


def download(name:str, opener=urllib.request.urlopen):
    ref=SOURCE[name]
    req=urllib.request.Request(ref["url"],headers={"User-Agent":"CKM-Stroke-Original-Cohort-Evidence/1.0"})
    with opener(req,timeout=50) as r:
        p=urlsplit(r.geturl())
        if getattr(r,"status",200)!=200 or p.scheme!="https" or p.hostname!=ref["host"]:
            raise ValueError("official original source redirected unexpectedly "+name)
        length=r.headers.get("Content-Length","")
        if not length.isdecimal() or int(length)!=ref["bytes"]:
            raise ValueError("official corrected supplement content length changed "+name)
        blob=r.read(ref["bytes"]+1)
    assert_publisher_file(name,blob)
    return blob


def assert_publisher_file(name:str,blob:bytes):
    ref=SOURCE[name]
    if len(blob)!=ref["bytes"] or hashlib.sha256(blob).hexdigest()!=ref["sha256"]:
        raise ValueError("publisher-hosted original XLSX SHA256 or bytes changed "+name)
    if not blob.startswith(b"PK\x03\x04"):raise ValueError("not XLSX/ZIP")
    return True


def source_rows(name:str,blob:bytes):
    assert_publisher_file(name,blob)
    ref=SOURCE[name]
    with ZipFile(io.BytesIO(blob)) as z:
        members=z.infolist()
        if len(members)>MAX_FILES or sum(x.file_size for x in members)>MAX_INFLATED:
            raise ValueError("oversized XLSX container")
        if any(x.file_size>MAX_XML_ENTRY for x in members):
            raise ValueError("oversized XLSX XML member")
        wb=ET.fromstring(z.read("xl/workbook.xml"))
        sheets=[(s.get("name"),s.get(REL)) for s in wb.iter(NS+"sheet")]
        if not any(s[0]==ref["name"] for s in sheets):
            raise ValueError("expected first-party study sheet missing")
        shared=[]
        if "xl/sharedStrings.xml" in z.namelist():
            shared=["".join(x.text or "" for x in element.iter(NS+"t"))
                    for element in ET.fromstring(z.read("xl/sharedStrings.xml")).findall(NS+"si")]
        root=ET.fromstring(z.read(ref["sheet"]))
        rows={}
        for r in root.iter(NS+"row"):
            index=int(r.get("r","0"))
            cells={}
            for c in r.findall(NS+"c"):
                match=re.fullmatch(r"([A-Z]+)\d+",c.get("r",""))
                if match is None:raise ValueError("invalid XLSX cell coordinate")
                el=c.find(NS+"v")
                value=el.text if el is not None else ""
                if c.get("t")=="s" and value:
                    si=int(value)
                    if not 0<=si<len(shared):raise ValueError("bad shared string")
                    value=shared[si]
                elif c.get("t")=="inlineStr":
                    v=c.find(NS+"is")
                    value="".join(x.text or "" for x in v.iter(NS+"t")) if v is not None else ""
                if value is not None and value!="":
                    cells[match.group(1)]=str(value).strip()
            rows[index]=cells
    return rows


def parse_n_from_cell(v:str):
    match=re.match(r"^\s*([0-9][0-9,]*)",v)
    return int(match.group(1).replace(",","")) if match else None


def extract_locke(rows:dict):
    if (rows.get(2,{}).get("A")!="Study" or
        rows.get(2,{}).get("B")!="Full name" or
        "Samples in analyses" not in rows.get(2,{}).get("H","")):
        raise ValueError("Locke Supplement 1 cohort table columns changed")
    out={}
    for idx,r in rows.items():
        if not 5<=idx<=143 or not r.get("B") or not r.get("A"):
            continue
        name=r["A"]
        if name in out:
            # Some primary source labels recur under GWAS vs Metabochip.
            # Do not silently equate participants; record additional rows explicitly.
            out[name]["additional_source_rows"].append({
                "row":idx,"sample_count_cell":r.get("H"),
                "ethnicity":r.get("D"),
                "source_section":("PRIOR_GWAS" if idx<53 else
                                  "UPDATED_GWAS" if idx<93 else
                                  "METABOCHIP" if idx<142 else "ADDITIONAL_STUDY"),
            })
            continue
        out[name]={
            "row":idx,"study":name,"full_name":r["B"],
            "sample_count_cell":r.get("H"),
            "reported_analysis_n_leading_numeric":parse_n_from_cell(r.get("H","")),
            "ethnicity":r.get("D"),
            "source_section":(
                "PRIOR_GWAS" if idx<53 else
                "UPDATED_GWAS" if idx<93 else
                "METABOCHIP" if idx<142 else "ADDITIONAL_STUDY"
            ),
            "additional_source_rows":[],
        }
    if len(out)<105 or len(out)>145 or out.get("ARIC",{}).get("reported_analysis_n_leading_numeric")!=8108:
        raise ValueError("Locke EUR cohort table profile missing")
    if (out.get("CHS",{}).get("reported_analysis_n_leading_numeric")!=3228 or
        out.get("FRAM",{}).get("reported_analysis_n_leading_numeric")!=8094 or
        out.get("AGES",{}).get("reported_analysis_n_leading_numeric")!=3207):
        raise ValueError("Locke source anchor rows changed")
    if any(out[k]["ethnicity"]!="white European" for k in ("ARIC","CHS","FRAM","AGES")):
        raise ValueError("Locke European ancestry anchors changed")
    return out


def extract_gig(rows:dict):
    head=rows.get(4,{})
    if [head.get(x) for x in ("K","L","M","N","P")]!=["AIS","CES","SVS","LAS","Controls"]:
        raise ValueError("GIGASTROKE corrected source phenotype column order changed")
    total=rows.get(49,{})
    if total.get("A")!="TOTAL EUR":
        raise ValueError("GIGASTROKE corrected EUR total missing")
    got={k:parse_n_from_cell(total.get(col,"")) for k,col in
         (("AIS","K"),("CE","L"),("SVD","M"),("LAA","N"),("controls","P"))}
    if got!=EXPECTED_EUR:
        raise ValueError("GIGASTROKE corrected subtype totals changed")
    result={}
    for n in range(5,49):
        cells=rows.get(n,{})
        name=cells.get("B")
        if not name or name in result:
            raise ValueError("missing/duplicated GIGASTROKE study name in EUR")
        result[name]={
            "row":n,"study":name,
            "AIS_cases":parse_n_from_cell(cells.get("K","")),
            "CE_cases":parse_n_from_cell(cells.get("L","")),
            "SVD_cases":parse_n_from_cell(cells.get("M","")),
            "LAA_cases":parse_n_from_cell(cells.get("N","")),
            "controls":parse_n_from_cell(cells.get("P","")),
            "design":cells.get("T",""),
        }
    if len(result)!=44:
        raise ValueError("expected exactly 44 source EUR cohorts")
    return result


def compare(locke:dict,gig:dict):
    joined=[]
    for name in sorted(ALIASES):
        lookup,tier,reason=ALIASES[name]
        if tier not in CATEGORY or name not in gig or lookup not in locke:
            raise ValueError("unsupported or stale alias: "+name)
        l=locke[lookup]; g=gig[name]
        joined.append({
            "GIGASTROKE_row":g["row"],
            "GIGASTROKE_source_label":name,
            "Locke_2015_row":l["row"],
            "Locke_2015_source_label":lookup,
            "Locke_2015_full_cohort_name":l["full_name"],
            "source_label_classification":tier,
            "interpretation":reason,
            "Locke_reported_analysis_N_cell":l["sample_count_cell"],
            "Locke_source_section":l["source_section"],
            "Locke_ethnicity":l["ethnicity"],
            "GIGASTROKE_reported_AIS_cases":g["AIS_cases"],
            "GIGASTROKE_reported_LAA_cases":g["LAA_cases"],
            "GIGASTROKE_reported_CE_cases":g["CE_cases"],
            "GIGASTROKE_reported_SVD_cases":g["SVD_cases"],
            "GIGASTROKE_reported_controls":g["controls"],
            "same_individuals_verified":False,
            "individual_overlap_count":None,
        })
    matched={x["GIGASTROKE_source_label"] for x in joined}
    unmatched=sorted(set(gig)-matched)
    counts=Counter(x["source_label_classification"] for x in joined)
    if ("UKBiobank" not in unmatched or
        set(gig)!=matched|set(unmatched) or
        counts["COHORT_LABEL_SHARED"]<10):
        raise ValueError("shared-cohort names changed or overlap validation needs revision")
    return {
        "schema":"ckm-locke2015-gigastroke-eur-firstparty-cohortname-audit-v1",
        "technical_status":"PASS_FIRSTPARTY_TWO_SUPPLEMENT_COHORT_NAME_RECONCILIATION",
        "source_provenance":{
            name:{"url":SOURCE[name]["url"],"sha256":SOURCE[name]["sha256"],
                  "bytes":SOURCE[name]["bytes"],"sheet":SOURCE[name]["name"]}
            for name in SOURCE
        },
        "Locke_source_study_labels":len(locke),
        "GIGASTROKE_EUR_cohorts":len(gig),
        "direct_cohort_label_matches":counts["COHORT_LABEL_SHARED"],
        "study_family_only_matches":counts["STUDY_FAMILY_SHARED"],
        "unmatched_GIGASTROKE_study_labels":len(unmatched),
        "matched_study_rows":joined,
        "unmatched_GIGASTROKE_labels":unmatched,
        "UK_Biobank_in_Locke2015_source_labels":False,
        "UK_Biobank_in_GIGASTROKE_source_labels":True,
        "exact_individual_overlap_count_identified":False,
        "LDScore_overlap_intercept_assumption_approved":False,
        "ready_for_ldsc":False,"scientific_status":"HOLD",
        "LDSC_h2_rg_MR_run":False,
        "limitations":[
            "Matching cohorts by name does NOT establish any shared individuals, visits or genotype files",
            "Labels for linked cohorts (e.g. HUNT waves, Estonian Biobank, EPIC umbrella) are possible families only",
            "Some GIGASTROKE labels have NA subtype counts; not a zero individual-overlap proof",
            "Locke Table 1 includes GWAS, Metabochip and additional-study rows; not per-SNP participating cohort indicators",
            "No sample-identifiers or cohort-specific SNP effects available from these two published tables",
            "No individual overlap fraction may be computed by summing or comparing these cohort sample counts"
        ],
    }


def audit(*,local_locke=None,local_gig=None,opener=urllib.request.urlopen):
    out={"technical_status":"HOLD_FIRSTPARTY_COHORT_SOURCE_OR_ALIAS_QC",
         "scientific_status":"HOLD","ready_for_ldsc":False,"LDSC_h2_rg_MR_run":False}
    try:
        locke=(Path(local_locke).read_bytes() if local_locke else download("LOCKE",opener))
        gig=(Path(local_gig).read_bytes() if local_gig else download("GIGASTROKE",opener))
        a=extract_locke(source_rows("LOCKE",locke))
        b=extract_gig(source_rows("GIGASTROKE",gig))
        out=compare(a,b)
        out["checked_at_utc"]=datetime.now(timezone.utc).isoformat()
    except Exception as e:
        out["error_type"]=type(e).__name__
        out["error_detail"]=str(e)[:240]
    return out


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--local-locke",type=Path,default=None)
    parser.add_argument("--local-gig",type=Path,default=None)
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    out=audit(local_locke=args.local_locke,local_gig=args.local_gig)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(out,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps({
        "status":out["technical_status"],
        "matched_source_labels":out.get("direct_cohort_label_matches"),
        "family_labels":out.get("study_family_only_matches"),
        "unmatched":out.get("unmatched_GIGASTROKE_study_labels"),
        "error":out.get("error_detail"),"scientific_status":out["scientific_status"]
    }))
    if out["technical_status"]!="PASS_FIRSTPARTY_TWO_SUPPLEMENT_COHORT_NAME_RECONCILIATION":
        raise SystemExit(1)


if __name__=="__main__":main()
