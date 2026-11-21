#!/usr/bin/env python3
"""Check original TPMI and CKB PheWeb public *phenotype pages*, not GWAS data.

Tested endpoints are not assumed valid TOAST AIS just because the URL resolves.
A PGS Catalog scoring file is NOT an original GWAS summary association table.
Disjointness of GIGASTROKE vs TPMI participants remains unverified.
"""
from __future__ import annotations
import argparse
import hashlib
from html.parser import HTMLParser
import json
from pathlib import Path
import re
from urllib.parse import urlsplit
from urllib.request import Request,urlopen
from urllib.error import HTTPError,URLError

TPMI="pheweb.ibms.sinica.edu.tw"
CKB="pheweb.ckbiobank.org"
PGS="www.pgscatalog.org"
URLS={
 "tpmi_about":"https://pheweb.ibms.sinica.edu.tw/about",
 "tpmi_433_21":"https://pheweb.ibms.sinica.edu.tw/pheno/433.21",
 "tpmi_433_2":"https://pheweb.ibms.sinica.edu.tw/pheno/433.2",
 "tpmi_433_3":"https://pheweb.ibms.sinica.edu.tw/pheno-filter/433.3",
 "ckb_i63":"https://pheweb.ckbiobank.org/pheno/i63",
 "pgs_019250":"https://www.pgscatalog.org/score/PGS019250/",
 "pgs_019234":"https://www.pgscatalog.org/score/PGS019234/",
}
MAX_BYTES=1_500_000
MAX_EXPECTED_CASES=1_000_000
EXPECTED_PHEWEB_BUILD="GRCh38"
CLINICAL={
 "tpmi_433_21":{"label":"Cerebral artery occlusion, with cerebral infarction","category":"POTENTIAL_ISCHEMIC_INFARCTION_NOT_IMAGING_ADJUDICATED"},
 "tpmi_433_2":{"label":"Occlusion of cerebral arteries","category":"NON_SPECIFIC_OCCLUSION_NOT_PROVEN_INFARCTION"},
 "tpmi_433_3":{"label":"Cerebral ischemia","category":"BROAD_ISCHEMIA_NOT_PROVEN_ACUTE_INFARCTION"},
 "ckb_i63":{"label":"i63","category":"ICD10_ISCHEMIC_STROKE_BUT_GIGASTROKE_CKB_COHORT_OVERLAP_RISK"},
}
CASES=re.compile(r"(?<!\w)([0-9][0-9,]*(?:\.0)?)\s+cases?\s*,?\s*([0-9][0-9,]*(?:\.0)?)\s+controls?",re.I)
H1=re.compile(r"<h1[^>]*>(.*?)</h1>",re.I|re.S)
TAGS=re.compile(r"<[^>]+>")
class HtmlText(HTMLParser):
    def __init__(self):
        super().__init__()
        self.parts=[]
    def handle_data(self,d):
        self.parts.append(d)
    def text(self):
        return " ".join(" ".join(self.parts).split())

def safe_fetch(url,fetcher=urlopen):
    target=urlsplit(url)
    if target.scheme!="https" or target.hostname not in (TPMI,CKB,PGS):
        raise ValueError("not a supported official GWAS source hostname")
    req=Request(url,headers={"User-Agent":"CKM-Stroke-EAS-Phenotype-Provenance/1.0",
                              "Accept":"text/html"})
    try:
        with fetcher(req,timeout=40) as result:
            actual=urlsplit(result.geturl())
            if result.status!=200 or actual.scheme!="https" or actual.hostname!=target.hostname:
                raise ValueError("source HTTP status or host redirect unexpected")
            data=result.read(MAX_BYTES+1)
            if len(data)>MAX_BYTES or not data:
                raise ValueError("source page empty or over size cap")
            try:
                page=data.decode("utf-8-sig","strict")
            except UnicodeError:
                raise ValueError("source not valid UTF8")
            return {"status":"HTTP_200","source_sha256":hashlib.sha256(data).hexdigest(),"text":page}
    except HTTPError as exc:
        if exc.code in (403,404,429,503):
            return {"status":"HTTP_"+str(exc.code),"source_sha256":None,"text":None}
        raise
    except URLError as exc:
        return {"status":"NETWORK_ERROR","source_sha256":None,"text":None,
                "error_type":type(exc).__name__}

class DownloadAnchors(HTMLParser):
    """Metadata-only discovery of public anchor labels; never open large files."""
    def __init__(self):
        super().__init__()
        self.entries=[]
        self.current=None
    def handle_starttag(self,tag,attrs):
        if tag=="a":
            href=dict(attrs).get("href")
            self.current={"href":href,"label":""} if href else None
    def handle_data(self,data):
        if self.current is not None:
            self.current["label"]+=data
    def handle_endtag(self,tag):
        if tag=="a" and self.current is not None:
            self.entries.append(self.current)
            self.current=None

def advertised_download_anchors(key,html):
    if key not in URLS:
        raise ValueError("unknown official phenotype URL")
    from urllib.parse import urljoin
    host=urlsplit(URLS[key]).hostname
    p=DownloadAnchors()
    p.feed(html)
    selected=[]
    for item in p.entries:
        if "summary statistics" not in item["label"].lower():
            continue
        url=urljoin(URLS[key],item["href"])
        parts=urlsplit(url)
        if parts.scheme!="https" or parts.hostname!=host:
            raise ValueError("public source download anchor points off official host")
        selected.append({"label":" ".join(item["label"].split())[:90],"url":url,
                         "sumstats_response_or_header_has_NOT_been_downloaded":True})
    if len(selected)>6:
        raise ValueError("unexpected excess summary statistics download anchors")
    return selected


def clean_html(source):
    p=HtmlText()
    p.feed(source)
    return p.text()

def parse_gwas_page(key,web):
    if key not in CLINICAL:
        raise ValueError("invalid GWAS phenotype reference")
    if web["status"]!="HTTP_200":
        return {"phenotype":key,"status":web["status"],
                "publicly_indexed_endpoint_with_source_verified":False,
                "exact_phenotype_definition_approved":False,
                "n_cases":None,"n_controls":None,
                "not_per_variant_Neff":True}
    plain=clean_html(web["text"])
    match=CASES.search(plain)
    pattern=CLINICAL[key]
    label_found=pattern["label"].lower() in plain.lower()
    if not label_found or not match:
        return {"phenotype":key,"status":"HOLD_SOURCE_PAGE_PARSING",
                "publicly_indexed_endpoint_with_source_verified":False,
                "exact_phenotype_definition_approved":False,
                "n_cases":None,"n_controls":None,
                "not_per_variant_Neff":True,"source_sha256":web["source_sha256"]}
    def asint(x):
        n=int(float(x.replace(",","")))
        if n<=0 or n>MAX_EXPECTED_CASES:
            raise ValueError("public phenotype N exceeds source guard")
        return n
    cases,controls=(asint(match.group(1)),asint(match.group(2)))
    download_anchors=advertised_download_anchors(key,web["text"])
    return {
        "phenotype":key,"clinical_category":pattern["category"],
        "status":"PASS_ORIGINAL_SOURCE_PHENOTYPE_COUNT",
        "publicly_indexed_endpoint_with_source_verified":True,
        "reported_label":pattern["label"],"n_cases":cases,
        "n_controls":controls,"not_per_variant_Neff":True,
        "not_proven_same_original_2026_TPMI_early_onset_GWAS":True,
        "not_proven_TOAST_or_imaging_adjudicated":True,
        "source_sha256":web["source_sha256"],
        "on_page_advertised_summary_stats_links_metadata_only":download_anchors,
        "raw_GWAS_full_summary_statistics_downloaded":False,
        "sumstats_variant_schema_or_SNP_Neff_verified":False,
        "exact_phenotype_definition_approved":False,
    }

def parse_about(web):
    if web["status"]!="HTTP_200":
        raise ValueError("official TPMI portal methods/build cannot be verified")
    page=clean_html(web["text"])
    if not("GRCh38" in page and "GWAS" in page and "TPMI" in page):
        raise ValueError("TPMI original PheWeb GRCh38 GWAS provenance missing")
    return {"genome_build":EXPECTED_PHEWEB_BUILD,"source_sha256":web["source_sha256"],
            "GWAS_and_PGS_resources_described":True,
            "original_2026_TPMI_paper_files_equal_this_portal_unverified":True}

def parse_pgs(key,web):
    if key not in ("pgs_019250","pgs_019234"):
        raise ValueError("unsupported PGS Catalog source")
    if web["status"]!="HTTP_200":
        return {"status":web["status"],"is_original_GWAS_summary_statistics":False}
    text=clean_html(web["text"])
    pgsid="PGS019250" if key=="pgs_019250" else "PGS019234"
    if pgsid not in text or "Polygenic Score" not in text or "GRCh38" not in text:
        raise ValueError("PGS Catalog identifier or source build unverified")
    return {"status":"PASS_PUBLIC_PGS_METADATA_ONLY",
            "PGS_ID":pgsid,"is_original_GWAS_summary_statistics":False,
            "not_validated_per_SNP_effect_estimates_of_outcome":True,
            "source_sha256":web["source_sha256"]}

def head_only_download_preflight(anchors,opener=urlopen):
    """Metadata-only HEAD, at most 2 links. No GWAS file content is read."""
    if len(anchors)>2:
        raise ValueError("too many TPMI original candidate HEAD links")
    out=[]
    for item in anchors:
        url=item["url"].strip()
        p=urlsplit(url)
        if (p.scheme!="https" or p.hostname!=TPMI
            or not p.path.startswith("/download/")
            or p.query or p.fragment):
            raise ValueError("TPMI source download path not an approved origin")
        req=Request(url,method="HEAD",
                    headers={"User-Agent":"CKM-Stroke-Original-TPMI-HEAD-Only/1.0"})
        try:
            with opener(req,timeout=25) as result:
                actual=urlsplit(result.geturl())
                is_official=actual.scheme=="https" and actual.hostname==TPMI
                length=result.headers.get("Content-Length","")
                out.append({
                  "label":item["label"],"original_link":url,
                  "http_status":result.status,
                  "final_redirect_hostname_official":is_official,
                  "content_length":int(length) if length.isdecimal() else None,
                  "request_method":"HEAD_ONLY",
                  "download_body_bytes_read":0,
                  "source_GWAS_variant_schema_verified":False,
                  "license_and_source_use_requirements_still_apply":True
                })
        except HTTPError as exc:
            if exc.code not in (401,403,404,405,429,503):
                raise
            out.append({
              "label":item["label"],"original_link":url,
              "http_status":exc.code,"request_method":"HEAD_ONLY",
              "download_body_bytes_read":0,
              "source_GWAS_variant_schema_verified":False,
              "license_and_source_use_requirements_still_apply":True
            })
        except URLError as exc:
            out.append({
              "label":item["label"],"original_link":url,
              "http_status":None,"network_error_type":type(exc).__name__,
              "request_method":"HEAD_ONLY","download_body_bytes_read":0,
              "source_GWAS_variant_schema_verified":False,
              "license_and_source_use_requirements_still_apply":True
            })
    return out


def evidence(source):
    about=parse_about(source["tpmi_about"])
    targets={x:parse_gwas_page(x,source[x]) for x in CLINICAL}
    pgs={x:parse_pgs(x,source[x]) for x in ("pgs_019250","pgs_019234")}
    return {
      "schema":"ckm-EAS-AIS-source-phenotype-overlap-and-access-v1",
      "technical_status":"PASS_OFFICIAL_PHEWEB_SOURCE_PROVENANCE_WITH_PER_PHENOTYPE_GATES",
      "original_sources":URLS,"TPMI_original_PheWeb":about,
      "original_public_endpoint_checks":targets,"PGS_Catalog_not_outcome_sumstats":pgs,
      "external_TPMI_2026_study":{
         "publication":"Neurology 2026 doi:10.1212/WNL.0000000000218126",
         "source_cohort":"Taiwan Precision Medicine Initiative",
         "published_2026_study_reported_cases":5546,
         "published_2026_study_reported_controls":143017,
         "reported_study_phenotype":"early-onset ischemic stroke analysis",
         "equivalence_to_2025_TPMI_PheWeb_phenotypes_verified":False,
         "SNP_specific_Neff_obtained":False,
         "original_2026_study_GWAS_full_sumstats_downloaded":False,
      },
      "study_overlap_risk":{
         "BBJ_prior_GIGASTROKE_or_MEGASTROKE_contribution_known":True,
         "CKB_prior_GIGASTROKE_contribution_risk":True,
         "TPMI_no_GIGASTROKE_participant_ID_overlap_proven":False,
         "2025_PheWeb_and_2026_TPMI_are_same_recruitment_program":True,
         "GIGASTROKE_EUR_vs_TPMI_EAS_populations_need_separate_LD_and_phenotype_review":True,
      },
      "research_decision":{
         "top_target":"TPMI PheWeb source-backed phenotype with documented infarction definition, not generic cerebral ischemia",
         "primary_EAS_AIS_replication_GO":False,
         "TPMI_PheWeb_source_candidate_to_study":"SOURCE_FOUND_PHENOTYPE_NEEDS_VALIDATION",
         "independent_replication_claim":False,
         "new_h2_rg_MR_coloc_executed":False,
         "original_GIGASTROKE_SNP_specific_Neff_approved":False,
      },
      "scientific_status":"HOLD","ready_for_ldsc":False
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",required=True,type=Path)
    args=p.parse_args()
    try:
        source={key:safe_fetch(url) for key,url in URLS.items()}
        report=evidence(source)
        selected=report["original_public_endpoint_checks"]["tpmi_433_21"]
        links=selected.get("on_page_advertised_summary_stats_links_metadata_only",[])
        report["TPMI_433_21_original_GWAS_access_HEAD_only"]=head_only_download_preflight(links)
        report["TPMI_full_GWAS_downloaded"]=False
        report["original_outcome_beta_se_p_and_SNP_Neff_schema_verified"]=False
    except Exception as exc:
        report={"technical_status":"HOLD_ORIGINAL_EAS_SOURCE_UNVERIFIED",
                "error_type":type(exc).__name__,"error_detail":str(exc)[:250],
                "scientific_status":"HOLD","ready_for_ldsc":False}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n","utf-8")
    print(json.dumps({"status":report["technical_status"],
                      "phenotype_evidence":{k:{"status":v["status"],"n_cases":v["n_cases"],"n_controls":v["n_controls"]} for k,v in report.get("original_public_endpoint_checks",{}).items()},
                      "data_source_only_no_LDSC":report.get("ready_for_ldsc") is False,
                      "tpmi_433_21_HEAD_results":[{"http_status":x.get("http_status"),"label":x["label"],"size":x.get("content_length")} for x in report.get("TPMI_433_21_original_GWAS_access_HEAD_only",[])],
                      "error":report.get("error_detail")}))
    if report["technical_status"]!="PASS_OFFICIAL_PHEWEB_SOURCE_PROVENANCE_WITH_PER_PHENOTYPE_GATES":
        raise SystemExit(1)
if __name__=="__main__":
    main()
