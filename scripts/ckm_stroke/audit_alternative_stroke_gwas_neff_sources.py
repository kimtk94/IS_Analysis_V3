#!/usr/bin/env python3
"""Official GWAS alternative-source schema and SNP-wise N_eff eligibility gate.

Fetches only public consortium information/schema, R13 manifest where accessible,
and a pinned prior GIGASTROKE source run. Does NOT download protected summary
statistics or replace phenotype-level N with per-variant case/control N.
No LDSC, MR, beta/SE manipulation, or scientific GO.
"""
from __future__ import annotations
import argparse
from html.parser import HTMLParser
import hashlib
import json
from pathlib import Path
import re
import urllib.request
from urllib.parse import urlsplit

MEGA_URL="https://www.megastroke.org/"
FINN_SCHEMA_URL="https://finngen.gitbook.io/documentation/data-description.md"
FINN_ACCESS_URL="https://finngen.gitbook.io/documentation/data-download.md"
FINN_R13_MANIFEST="https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_manifest.tsv"
MAX_DOC_BYTES=280_000
MAX_MANIFEST_BYTES=6_000_000
GWAS_HEADERS=("chromosome","base_pair_location","effect_allele_frequency","beta",
              "standard_error","p_value","odds_ratio","ci_lower","ci_upper","effect_allele","other_allele")
MEGA_COLS=("MarkerName","Allele1","Allele2","Freq1","Effect","StdErr","P-value")
FINN_COLS=("#chrom","pos","ref","alt","rsids","nearest_genes","pval","mlogp","beta",
           "sebeta","af_alt","af_alt_cases","af_alt_controls")
N_FIELDS=("n","n_eff","n_cases","n_controls","num_cases","num_controls",
          "ncase","ncontrol","samplesize","sample_size","effective_n")
STROKE_PAT=re.compile(r"(?:stroke|cerebr|infarct|intracrani|I9_STR|I9_CER|I9_SAH|I9_ICH|I9_AIS)",re.I)
REPO="kimtk94/IS_Analysis_V3"
ORIGINAL_RUN=44860586028
ORIGINAL_SHA="7c1ca53a57f3e44b1cdb8e2cf19753b943cf45ba"
REFERENCES={
    "GIGASTROKE":"https://www.ebi.ac.uk/gwas/docs/summary-statistics-format/",
    "MEGASTROKE":MEGA_URL,
    "FinnGen_R13_schema":FINN_SCHEMA_URL,
    "FinnGen_R13_access":FINN_ACCESS_URL,
}

class TextOnly(HTMLParser):
    def __init__(self):
        super().__init__()
        self.parts=[]
    def handle_data(self,data):
        self.parts.append(data)
    def value(self):
        return " ".join(" ".join(self.parts).split())

def fetch(url,limit,fetcher=urllib.request.urlopen,hosts=None):
    parsed=urlsplit(url)
    if parsed.scheme!="https" or (hosts and parsed.hostname not in hosts):
        raise ValueError("non-approved source host")
    req=urllib.request.Request(url,headers={
        "User-Agent":"CKM-Stroke-Official-Neff-Source-Triage/1.0",
        "Accept":"text/plain, text/markdown, text/html, application/json; q=0.8"})
    with fetcher(req,timeout=35) as response:
        actual=urlsplit(response.geturl())
        if (response.status!=200 or actual.scheme!="https"
            or (hosts and actual.hostname not in hosts)):
            raise ValueError("official source status/redirect mismatch")
        data=response.read(limit+1)
        if len(data)>limit:
            raise ValueError("official bounded documentation/manifest exceeded limit")
        # The legacy consortium HTML page may use Windows-1252 smart quotes.
        # Preserve SHA of original bytes; fail closed on malformed official
        # Markdown/manifest rather than silently replacing N/schema characters.
        if parsed.hostname in ("www.megastroke.org","megastroke.org"):
            try:
                plain=data.decode("utf-8-sig","strict")
            except UnicodeDecodeError:
                plain=data.decode("cp1252","strict")
        else:
            plain=data.decode("utf-8-sig","strict")
        return plain,hashlib.sha256(data).hexdigest()

def check_finn_schema(doc):
    columns=[x for x in FINN_COLS if re.search(r"(?<![A-Za-z0-9_])"+re.escape(x)+r"(?![A-Za-z0-9_])",doc)]
    if len(columns)!=len(FINN_COLS) or "GRCh38" not in doc or "regenie" not in doc.lower():
        raise ValueError("FinnGen official summary schema or genome build changed")
    # The first table only: don't include fine-mapping field names or prose in N test.
    first=doc.split("## Fine-mapping results",1)[0]
    if not first or any(re.search(r"\|\s*\*{0,2}\x60?"+re.escape(x)+r"\x60?\*{0,2}\s*\|",first,re.I) for x in ("n_eff","n_cases","n_controls","num_cases","num_controls","sample_size")):
        raise ValueError("unexpected FinnGen per-variant sample N column in first GWAS table")
    return {"verified_summary_columns":list(FINN_COLS),
            "variant_N_or_Neff_column_present":False,
            "reference_genome":"GRCh38",
            "published_effect_allele":"ALT",
            "variant_specific_cases_controls_in_source_columns":False,
            "not_allele_frequency_based_sample_size_inference":True}

def check_mega_description(document):
    text=document
    if "<html" in document[:500].lower() or "<!doctype" in document[:500].lower():
        parser=TextOnly()
        parser.feed(document)
        text=parser.value()
    needed=["MEGASTROKE","Europeans","40,585","406,111","n_cases > 50%"]
    if any(item.lower() not in text.lower() for item in needed):
        raise ValueError("official MEGASTROKE EUR GWAS description changed")
    for col in MEGA_COLS:
        if not re.search(re.escape(col)+r"\s*:",text,re.I):
            raise ValueError("MEGASTROKE documented association column missing: "+col)
    descriptions=re.findall(r"(MarkerName|Allele1|Allele2|Freq1|Effect|StdErr|P-value)\s*:",text,re.I)
    if len(set(x.lower() for x in descriptions))!=7:
        raise ValueError("MEGASTROKE official output field partition drift")
    return {
        "documented_summary_columns":list(MEGA_COLS),
        "public_per_variant_Neff_column_documented":False,
        "variant_level_n_cases_used_for_filter_only":True,
        "EUR_cases":40585,"EUR_controls":406111,
        "stroke_pheno_set":["AS","AIS","LAS","CES","SVS"],
        "EUR_only_files_supported":True,
        "meta_analysis_is_EUR_independent_of_GIGASTROKE_participants_verified":False,
    }

def check_finn_access(document):
    if not ("Release 13" in document and "finngen_R13_manifest.tsv" in document
            and "online form" in document.lower()):
        raise ValueError("FinnGen R13 release/access provenance differs from official docs")
    return {
        "release":"R13",
        "manifest_source_url":FINN_R13_MANIFEST,
        "official_download_requires_request_or_access_step":True,
        "user_access_credentials_available":False,
        "can_assert_full_sumstats_file_header_checked":False,
    }

def classify_finn_manifest(body):
    lines=body.splitlines()
    if not lines or len(lines)>40000 or len(lines[0])>8192:
        raise ValueError("FinnGen R13 manifest header/line count invalid")
    header=[x.strip().lower().lstrip("\ufeff") for x in lines[0].split("\t")]
    if not any(x in header for x in ("phenocode","phenotype","trait","endpoint","phenotype_code")):
        raise ValueError("FinnGen manifest no phenotype code column; cannot safely classify")
    if len(set(header))!=len(header):
        raise ValueError("FinnGen manifest duplicated column headers")
    codeidx=next(i for i,h in enumerate(header)
                 if h in ("phenocode","phenotype","trait","endpoint","phenotype_code"))
    candidates=[]
    matched_metadata=[]
    description_idx=next((i for i,h in enumerate(header)
                          if h in ("phenostring","description","phenotype_description","phenotype_name","trait_name")),None)
    for index,line in enumerate(lines[1:],start=1):
        parts=line.split("\t")
        if len(parts)!=len(header):
            raise ValueError("FinnGen manifest row width changed")
        ph=parts[codeidx].strip()
        if not ph or len(ph)>160:
            raise ValueError("invalid FinnGen phenotype code")
        if STROKE_PAT.search(ph):
            candidates.append(ph)
            matched_metadata.append({
                "phenocode":ph,
                "manifest_label":parts[description_idx][:160].strip() if description_idx is not None else None,
                "clinical_case_definition_verified":False,
            })
    if len(candidates)>1000:
        raise ValueError("candidate endpoints exceed source review cap")
    return {"manifest_retrieved":True,"manifest_entries":len(lines)-1,
            "likely_stroke_endpoint_codes":sorted(set(candidates))[:80],
            "likely_stroke_endpoint_labels":sorted(matched_metadata,key=lambda d:d["phenocode"])[:80],
            "likely_stroke_endpoint_count":len(set(candidates)),
            "ischemic_stroke_candidate_codes_require_definition_review":[x for x in
                ("I9_STR","I9_STR_EMBOLIC","I9_STR_EMBOLIC_WITH_FLIMMER") if x in candidates],
            "hemorrhagic_endpoint_codes_EXCLUDED_from_ischemic_stroke":[x for x in
                ("I9_ICH","I9_SAH","I9_STR_SAH","I9_SAHANEUR") if x in candidates],
            "no_verified_LAA_or_SVD_specific_endpoint_claim":True,
            "endpoint_code_match_is_not_validated_stroke_case_definition":True,
            "per_variant_N_eff_verified_from_manifest":False,
            "per_variant_full_sumstats_header_read":False}

def review_sources(mega,finn_schema,finn_access,manifest=None,manifest_error=None):
    m=check_mega_description(mega)
    f=check_finn_schema(finn_schema)
    access=check_finn_access(finn_access)
    stat={"manifest_retrieved":False,
          "manifest_status":"NOT_RETRIEVED_OR_ACCESS_RESTRICTED",
          "likely_stroke_endpoint_codes":[],
          "per_variant_N_eff_verified_from_manifest":False,
          "per_variant_full_sumstats_header_read":False}
    if manifest is not None:
        stat=classify_finn_manifest(manifest)
        stat["manifest_status"]="SOURCE_MANIFEST_PARSED_NOT_SUMSTATS"
    if manifest_error:
        stat["manifest_unavailable_reason"]=manifest_error[:150]
    return {
      "schema":"ckm-stroke-alternative-gwas-per-variant-Neff-official-availability-v1",
      "technical_status":"PASS_OFFICIAL_CONSORTIUM_SCHEMA_NEFF_GAP_REVIEWED",
      "official_sources":{"MEGASTROKE":MEGA_URL,
                          "FinnGen_R13_schema":FINN_SCHEMA_URL,
                          "FinnGen_R13_access":FINN_ACCESS_URL,
                          "FinnGen_R13_manifest":FINN_R13_MANIFEST},
      "GIGASTROKE": {
        "source":"publisher-MD5-verified GCST90104540/41/42/43 EUR whole original files",
        "source_schema_columns":list(GWAS_HEADERS),
        "variant_specific_effective_N":False,
        "original_source_integrity_run":ORIGINAL_RUN,
        "reference_build":"GRCh37",
        "primary_scientific_gate":"HOLD_SNP_SPECIFIC_META_N_EFFECTIVE"},
      "MEGASTROKE_2018":{
        **m,
        "per_variant_N_eff_available_from_documented_schema":False,
        "possible_role":"older EUR subtype comparison, not automatically sample-independent",
        "genome_build_for_files_not_verified_by_this_schema_probe":True,
        "scientific_approval":"HOLD_NEFF_AND_PARTICIPANT_OVERLAP",
      },
      "FinnGen_R13":{
        **f,**access,"manifest":stat,
        "per_variant_N_eff_available_from_documented_schema":False,
        "case_count_endpoint_metadata_is_NOT_SNP_specific_N":True,
        "possible_role":"Finnish phenotype-defined replication candidate, subject to source phenotype and sample overlap validation",
        "requires_GRCh38_to_GRCh37_harmonization_before_cross_GWAS_merge":True,
        "scientific_approval":"HOLD_ENDPOINT_SELECTION_N_VALIDATION_AND_ACCESS",
      },
      "eligibility_rules":{
        "requires_source_documented_per_variant_case_and_control_N_or_validated_per_SNP_effective_N":True,
        "cohortwise_illustrative_N_not_original_effective_N":True,
        "constant_N_sensitivity_requires_prior_written_justification":True,
        "never_infer_per_variant_N_from_EAF_or_beta_SE_without_calibration":True,
        "source_provenance_and_overlapping_cohorts_must_be_checked":True,
        "independent_replication_is_not_a_substitute_for_Neff_of_primary_GIGASTROKE":True,
      },
      "findings":{
        "verified_alternative_outcomes_satisfying_primary_LDSC_Neff_criterion":0,
        "candidate_review_priority":["FinnGen R13 for replication once phenotype metadata and access confirmed",
                                     "MEGASTROKE EUR for subtype direction comparison, accounting for sample overlap",
                                     "GIGASTROKE original for primary effects pending an approved N_eff policy"],
        "alternative_per_variant_N_eff_problem_solved":False,
      },
      "LDSC_h2_intercept_rg_MR_coloc_executed":False,
      "scientific_status":"HOLD","ready_for_ldsc":False
    }

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",required=True,type=Path)
    args=ap.parse_args()
    try:
        a,ah=fetch(MEGA_URL,MAX_DOC_BYTES,hosts={"www.megastroke.org","megastroke.org"})
        b,bh=fetch(FINN_SCHEMA_URL,MAX_DOC_BYTES,hosts={"finngen.gitbook.io"})
        c,ch=fetch(FINN_ACCESS_URL,MAX_DOC_BYTES,hosts={"finngen.gitbook.io"})
        manifest=None
        manifestsha=None
        manifesterror=None
        try:
            manifest,manifestsha=fetch(FINN_R13_MANIFEST,MAX_MANIFEST_BYTES,
                                      hosts={"storage.googleapis.com"})
        except Exception as exc:
            manifesterror=type(exc).__name__+": "+str(exc)[:120]
        result=review_sources(a,b,c,manifest=manifest,manifest_error=manifesterror)
        result["source_document_sha256"]={"MEGASTROKE":ah,"FinnGen_R13_schema":bh,
                                           "FinnGen_R13_access":ch,
                                           "FinnGen_R13_manifest":manifestsha}
    except Exception as exc:
        result={"technical_status":"HOLD_OFFICIAL_ALTERNATIVE_SOURCE_SCHEMA_NOT_VERIFIED",
                "error_type":type(exc).__name__,"error_detail":str(exc)[:240],
                "scientific_status":"HOLD","ready_for_ldsc":False,
                "LDSC_h2_intercept_rg_MR_coloc_executed":False}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"status":result["technical_status"],
                      "manifest":result.get("FinnGen_R13",{}).get("manifest",{}).get("manifest_status"),
                      "likely_stroke_endpoint_count":result.get("FinnGen_R13",{}).get("manifest",{}).get("likely_stroke_endpoint_count"),
                      "alternative_Neff_solved":result.get("findings",{}).get("alternative_per_variant_N_eff_problem_solved"),
                      "scientific_status":result["scientific_status"],
                      "error":result.get("error_detail")},ensure_ascii=False))
    if result["technical_status"]!="PASS_OFFICIAL_CONSORTIUM_SCHEMA_NEFF_GAP_REVIEWED":
        raise SystemExit(1)
if __name__=="__main__":
    main()
