#!/usr/bin/env python3
"""Official FinnGen older-release AIS manifest and R13 access/provenance triage.

This is a public metadata-only CI gate. It does not download GWAS sumstats, use
protected credentials, infer per-SNP N_eff, run LDSC, or claim that separate
consortium names establish nonoverlapping participants.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import re
import urllib.error
import urllib.request
from urllib.parse import urlsplit
from scripts.ckm_stroke.audit_finngen_r13_endpoint_definitions_neff import (
    FINN_R13_MANIFEST, R13_MANIFEST_SHA, manifest_review,
)
from scripts.ckm_stroke.audit_gigastroke_supplement_cohort_neff import (
    download_source, workbook_cells, EXPECTED_SHA256 as SUPPLEMENT_SHA
)

MANIFESTS={
    "R9":"https://storage.googleapis.com/finngen-public-data-r9/summary_stats/R9_manifest.tsv",
    "R10":"https://storage.googleapis.com/finngen-public-data-r10/summary_stats/R10_manifest.tsv",
    "R13":FINN_R13_MANIFEST,
}
RELEVANT=("I9_STR_EXH","I9_STR","I9_STR_EMBOLIC","I9_STR_EMBOLIC_WITH_FLIMMER")
PUBLIC_GCS="storage.googleapis.com"
LIMIT=5_000_000
COHORT_KEYWORDS=re.compile(r"finn|finrisk|fing|helsin|turku|suomi|sisu|kuopio|nordic",re.I)

def get_text(url,opener=urllib.request.urlopen):
    initial=urlsplit(url)
    if initial.scheme!="https" or initial.hostname!=PUBLIC_GCS:
        raise ValueError("only official FinnGen public GCS host allowed")
    request=urllib.request.Request(url,headers={"User-Agent":"CKM-Stroke-FinnGen-AIS-R9-R10-Metadata/1.0"})
    with opener(request,timeout=40) as response:
        actual=urlsplit(response.geturl())
        if getattr(response,"status",None)!=200 or actual.scheme!="https" or actual.hostname!=PUBLIC_GCS:
            raise ValueError("not official FinnGen public GCS manifest")
        data=response.read(LIMIT+1)
        if not data or len(data)>LIMIT:
            raise ValueError("public manifest empty or over maximum size")
    return data,hashlib.sha256(data).hexdigest()

def parse_release(release,content,sha):
    if release not in MANIFESTS:
        raise ValueError("unsupported FinnGen public release")
    if release=="R13" and sha!=R13_MANIFEST_SHA:
        raise ValueError("R13 previously verified public manifest SHA256 changed")
    lines=content.decode("utf-8-sig","strict").splitlines()
    if len(lines)<2 or len(lines)>10000:
        raise ValueError("FinnGen release manifest missing or excessive entries")
    reader=csv.DictReader(lines,delimiter="\t")
    columns=reader.fieldnames
    if columns is None or len(set(columns))!=len(columns) or "phenocode" not in columns:
        raise ValueError("official release manifest phenotype code column missing")
    # R9/R10 schema might differ from R13. Inspect instead of assuming.
    selected={}
    count=0
    for row in reader:
        count+=1
        code=row.get("phenocode")
        if not code or len(code)>130:
            raise ValueError("invalid public manifest endpoint code")
        if code in RELEVANT:
            if code in selected:
                raise ValueError("duplicated candidate endpoint in release manifest")
            ncase=row.get("num_cases") or row.get("n_cases")
            nctrl=row.get("num_controls") or row.get("n_controls")
            for x in (ncase,nctrl):
                if x is not None and not str(x).isdigit():
                    raise ValueError("noninteger public endpoint sample metadata")
            url=row.get("path_https")
            if url:
                u=urlsplit(url)
                if u.scheme!="https" or u.hostname!=PUBLIC_GCS:
                    raise ValueError("endpoint file URL points outside official public GCS")
            selected[code]={
                "phenocode":code,"phenotype":row.get("phenotype") or row.get("phenostring"),
                "release_cases":int(ncase) if ncase is not None else None,
                "release_controls":int(nctrl) if nctrl is not None else None,
                "exact_path_https":url,
                "file_link_is_not_evidence_of_access_or_per_variant_Neff":True,
            }
    if release=="R13":
        if count!=2754 or "I9_STR_EXH" in selected or not all(code in selected for code in RELEVANT[1:]):
            raise ValueError("previous pinned R13 phenotype availability changed")
    return {
        "release":release,"manifest_url":MANIFESTS[release],
        "manifest_sha256":sha,"publisher_independent_expected_SHA256_known":release=="R13",
        "total_manifest_endpoints":count,"manifest_columns":columns,
        "requested_endpoints":{k:{"present":k in selected,**(selected.get(k) or {})} for k in RELEVANT},
        "source_manifest_verified":True,
        "per_SNP_Neff_in_manifest":False,
        "older_release_does_not_establish_independent_cohort_from_R13":True,
    }

def head_check(url,opener=urllib.request.urlopen):
    u=urlsplit(url)
    if u.scheme!="https" or u.hostname!=PUBLIC_GCS:
        raise ValueError("not original public manifest-linked FinnGen GCS URL")
    req=urllib.request.Request(url,method="HEAD",headers={
        "User-Agent":"CKM-Stroke-FinnGen-HEAD-Only-No-Download/1.0"})
    try:
        with opener(req,timeout=25) as res:
            actual=urlsplit(res.geturl())
            if actual.scheme!="https" or actual.hostname!=PUBLIC_GCS:
                raise ValueError("unexpected redirected GWAS resource host")
            status=res.status
            length=res.headers.get("Content-Length")
            return {"HTTP_status":status,
                "declared_resource_bytes":int(length) if length and length.isdecimal() else None,
                "request_method":"HEAD_ONLY_NO_SUMSTATS_DOWNLOAD",
                "remote_file_header_or_variant_Neff_verified":False,
                "access_form_still_required_per_FinnGen_official_instructions":True,
                "clinical_data_access_approval_inferred":False}
    except urllib.error.HTTPError as exc:
        if exc.code not in (401,403,404,405,429):
            raise
        return {"HTTP_status":exc.code,"declared_resource_bytes":None,
            "request_method":"HEAD_ONLY_NO_SUMSTATS_DOWNLOAD",
            "remote_file_header_or_variant_Neff_verified":False,
            "access_form_still_required_per_FinnGen_official_instructions":True,
            "clinical_data_access_approval_inferred":False}

def cohort_review(blob):
    # The reusable upstream parser verifies the original 756910-byte publisher
    # supplement SHA256, workbook sheet names and source sheet cells.
    cohort,_=workbook_cells(blob)
    rows=[]
    for i in range(5,49):
        a=cohort.get(i,{})
        label=a.get("B","")
        if label:
            rows.append({"source_row":i,"name":label})
    if len(rows)!=44 or not any(x["name"]=="UKBiobank" for x in rows):
        raise ValueError("not original 44 EUR cohort source rows")
    hits=[x for x in rows if COHORT_KEYWORDS.search(x["name"])]
    return {
       "corrected_GIGASTROKE_supplement_SHA256":SUPPLEMENT_SHA,
       "EUR_cohort_rows":44,
       "explicit_Nordic_Finnish_named_cohort_candidates":hits,
       "explicit_label_FinnGen":any(x["name"].lower()=="finngen" for x in rows),
       "absence_of_FinnGen_label_proves_zero_overlap":False,
       "participant_ID_overlap_quantified":False,
       "independent_FinnGen_replication_cohort_claim_approved":False,
       "UKBiobank_cohort_present":True,
       "cohort_names_not_individual_participant_IDs":True,
    }

def audit(manifest_sources,cohort,head_checker=head_check):
    if set(manifest_sources)!=set(MANIFESTS):
        raise ValueError("all pinned release manifest keys required")
    reports={}
    for release in ("R9","R10","R13"):
        data,sha=manifest_sources[release]
        reports[release]=parse_release(release,data,sha)
    if not reports["R13"]["requested_endpoints"]["I9_STR_EMBOLIC"]["present"]:
        raise ValueError("R13 official GWAS stroke case no longer published")
    r13=reports["R13"]
    access={}
    for code in ("I9_STR_EMBOLIC","I9_STR_EMBOLIC_WITH_FLIMMER"):
        url=r13["requested_endpoints"][code].get("exact_path_https")
        if not url:
            raise ValueError("R13 manifest GWAS path_https missing")
        access[code]=head_checker(url)
    return {"schema":"ckm-FinnGen-R9-R10-R13-published-AIS-availability-and-overlap-v1",
        "technical_status":"PASS_PUBLIC_RELEASE_MANIFEST_AND_GIGASTROKE_COHORT_PROVENANCE",
        "release_manifest_evidence":reports,
        "FinnGen_R13_source_GWAS_HEAD_access_only":access,
        "GIGASTROKE_EUR_original_cohort_labels":cohort,
        "clinical_I9_STR_EXH_is_not_exact_TOAST_AIS":True,
        "I9_STR_includes_ICH_not_valid_as_AIS":True,
        "R9_R10_same_FinnGen_project_independent_replication_between_releases":False,
        "MEGASTROKE_GIGASTROKE_participant_overlap_resolved":False,
        "all_FinnGen_public_downloads_subject_to_official_access_form":True,
        "primary_LDSC_SNP_specific_Neff_resolved":False,
        "primary_LDSC_h2_rg_MR_coloc_executed":False,
        "scientific_status":"HOLD","ready_for_ldsc":False}

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True,type=Path)
    a=p.parse_args()
    output=a.out
    try:
        records={}
        for k,u in MANIFESTS.items():
            records[k]=get_text(u)
        supplement=download_source()
        result=audit(records,cohort_review(supplement))
    except Exception as exc:
        result={"technical_status":"HOLD_OFFICIAL_RELEASE_AVAILABILITY_OR_COHORT_PROVENANCE",
            "error_type":type(exc).__name__,"error_detail":str(exc)[:380],
            "scientific_status":"HOLD","ready_for_ldsc":False,
            "primary_LDSC_SNP_specific_Neff_resolved":False}
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"status":result["technical_status"],
        "r9_ischaemic_exists":result.get("release_manifest_evidence",{}).get("R9",{}).get("requested_endpoints",{}).get("I9_STR_EXH",{}).get("present"),
        "r10_ischaemic_exists":result.get("release_manifest_evidence",{}).get("R10",{}).get("requested_endpoints",{}).get("I9_STR_EXH",{}).get("present"),
        "r13_gwas_head":{x:y["HTTP_status"] for x,y in result.get("FinnGen_R13_source_GWAS_HEAD_access_only",{}).items()},
        "cohort_label_candidates":result.get("GIGASTROKE_EUR_original_cohort_labels",{}).get("explicit_Nordic_Finnish_named_cohort_candidates"),
        "error":result.get("error_detail"),"science":result["scientific_status"]}))
    if result["technical_status"]!="PASS_PUBLIC_RELEASE_MANIFEST_AND_GIGASTROKE_COHORT_PROVENANCE":
        raise SystemExit(1)
if __name__=="__main__":
    main()
