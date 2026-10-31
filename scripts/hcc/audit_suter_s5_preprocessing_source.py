#!/usr/bin/env python3
"""Verify Suter 2022 S5 preprocessing PDF, preserving source provenance.

No patient/feature data, PDFs or verbatim extracted prose in CI artifacts.
Transport failure is HOLD, not 'preprocessing not performed'.
"""
import argparse
import datetime as dt
import hashlib
import io
import json
import pathlib
import re
import shutil
import subprocess
import tempfile
import urllib.error
import urllib.request

DOI="10.1371/journal.pcbi.1009767.s013"
SOURCES=(
    "https://journals.plos.org/ploscompbiol/article/file?id="+DOI+"&type=supplementary",
    "https://pmc.ncbi.nlm.nih.gov/articles/PMC9481159/bin/pcbi.1009767.s013.pdf",
    "https://doi.org/"+DOI,
)
MAX_BYTES=3_000_000
PATTERNS={
  "rna_or_transcriptome":r"\b(?:RNA|transcriptom\w*)\b",
  "protein_or_proteome":r"\b(?:protein|proteom\w*)\b",
  "imputation":r"\bimput\w*\b",
  "knn":r"\b(?:k[ -]?nearest|k[ -]?NN|impute\.knn)\b",
  "batch_or_combat":r"\b(?:batch|ComBat|removeBatchEffect)\b",
  "normalization":r"\bnormali[sz]\w*\b",
  "logarithm":r"\blog\s*2\b|logarithm",
  "normal_reference":r"\b(?:healthy|normal livers?|healthy liver)\b",
  "transformation":r"\b(?:center(?:ing|ed)?|standardiz\w*|scale|z.score)\b",
  "missingness":r"\bmissing\b",
}
def summarize_text(s):
    if len(s)<150 or not re.search(r"data\s+pre[ -]?processing|pre[ -]?processing\s+of|appendix",s,re.I):
        raise ValueError("PDF text is not identifiable as an appendix")
    return {k:len(re.findall(p,s,flags=re.I)) for k,p in PATTERNS.items()}

def fetch_pdf(url):
    req=urllib.request.Request(url,headers={"Accept":"application/pdf",
        "User-Agent":"HCC-method-provenance-audit/1.0"})
    with urllib.request.urlopen(req,timeout=35) as resp:
        b=resp.read(MAX_BYTES+1)
    if len(b)>MAX_BYTES or len(b)<500 or not b.startswith(b"%PDF-"):
        raise ValueError("not a bounded PDF")
    return b

def extract_pdf_text(b):
    if shutil.which("pdftotext"):
        with tempfile.TemporaryDirectory(prefix="hcc-s5-") as d:
            p=pathlib.Path(d)/"s5.pdf";p.write_bytes(b)
            run=subprocess.run(["pdftotext","-layout",str(p),"-"],
                capture_output=True,text=True,timeout=25,check=True)
            return run.stdout,"pdftotext"
    # GitHub-hosted runners may lack poppler. Pin pypdf in the workflow.
    from pypdf import PdfReader
    reader=PdfReader(io.BytesIO(b),strict=True)
    if len(reader.pages)==0 or len(reader.pages)>20:
        raise ValueError("PDF page count out of range")
    return "\n".join((p.extract_text() or "") for p in reader.pages),"pypdf_6.1.3"


def layer_method_keyword_sections(text):
    """Extract *mentions within headings*, not assert application to HCCfull.

    Unlike globally counting 'impute' or 'batch', separate each assay section
    and keep provenance in the official final-S5 PDF. Derived flags are hints
    that still require a researcher's source-method interpretation.
    """
    flat=re.sub(r"\\s+"," ",text)
    headings=[
      ("transcriptome",r"for\\s+transcriptome\\s+data[^.]{0,70}"),
      ("proteome",r"for\\s+proteome\\s+data[^.]{0,70}"),
      ("phosphoproteome",r"for\\s+phosphoproteome\\s+data[^.]{0,70}"),
      ("cna",r"(?:the\\s+)?CNA\\s+data[^.]{0,80}"),
    ]
    found=[]
    for name,regex in headings:
        m=re.search(regex,flat,flags=re.I)
        if m:found.append((m.start(),name))
    found.sort()
    result={}
    for i,(start,name) in enumerate(found):
        end=found[i+1][0] if i+1<len(found) else len(flat)
        segment=flat[start:end]
        result[name]={
          "section_heading_detected":True,
          "segment_character_count":len(segment),
          "upper_quartile_mentioned":bool(re.search(r"upper[- ]quartile",segment,re.I)),
          "normalization_to_1000_mentioned":bool(re.search(r"\\b1000\\b",segment,re.I)),
          "log2_mentioned":bool(re.search(r"log\\s*2",segment,re.I)),
          "median_subtraction_mentioned":bool(re.search(r"median.{0,30}subtract|subtract.{0,30}median",segment,re.I)),
          "half_sample_detection_filter_mentioned":bool(re.search(r"50\\s*%|half.{0,15}sample",segment,re.I)),
          "imputation_mentioned":bool(re.search(r"imput",segment,re.I)),
          "clustering_only_qualifier_mentioned":bool(re.search(
              r"for\\s+clustering\\s+only|only\\s+for\\s+clustering",segment,re.I)),
          "unimputed_differential_expression_mentioned":bool(re.search(
              r"unimputed.{0,60}differential|differential.{0,100}unimputed",segment,re.I)),
          "batch_correction_mentioned":bool(re.search(r"batch.{0,30}correct",segment,re.I)),
          "edgeR_mentioned":bool(re.search(r"\\bedgeR\\b",segment,re.I)),
        }
    return result


def audit():
    report={
        "schema":"hcc-stage25-suter-s5-source-audit-v1",
        "checked_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
        "expected_supplement_doi":DOI,
        "source_paper_doi":"10.1371/journal.pcbi.1009767",
        "original_ng_doi":"10.1038/s41467-022-29960-8",
        "status":"HOLD_UNVERIFIED_S5_SOURCE",
        "verified_source_url":None,
        "source_pdf_sha256":None,
        "source_pdf_bytes":None,
        "pdftotext_available":shutil.which("pdftotext") is not None,
        "pdf_extractor":None,
        "method_keyword_counts":{},
        "layer_method_section_indicators":{},
        "layer_method_indicators_are_text_mentions_only":True,
        "annotation":"Counts identify text mentions only, not the layer-specific transformation, order, parameters, or validation.",
        "s5_all_processing_parameters_manually_reviewed":False,
        "raw_ng_counts_reanalysis_performed":False,
        "hcc_processed_rds_identical_to_ng_raw_data":False,
        "original_ng_gene_set_identically_reproduced":False,
        "automatic_V03_pass":False,
        "pdf_or_full_verbatim_text_exported":False,
        "failures_redacted":[],
    }
    for idx,url in enumerate(SOURCES):
        try:
            b=fetch_pdf(url)
            extracted,extractor=extract_pdf_text(b)
            words=summarize_text(extracted)
            layers=layer_method_keyword_sections(extracted)
            report.update({
                "status":"SOURCE_PDF_TEXT_VERIFIED_NEEDS_LAYER_METHOD_REVIEW",
                "verified_source_url":url,
                "source_pdf_sha256":hashlib.sha256(b).hexdigest(),
                "source_pdf_bytes":len(b),
                "method_keyword_counts":words,
                "layer_method_section_indicators":layers,
                "pdf_extractor":extractor,
            })
            return report
        except (OSError,ValueError,TimeoutError,subprocess.SubprocessError,
                urllib.error.URLError) as exc:
            report["failures_redacted"].append(
                "SOURCE_"+str(idx+1)+"_"+type(exc).__name__)
    return report

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",required=True)
    a=ap.parse_args()
    r=audit()
    out=pathlib.Path(a.out);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(r,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(r,indent=2))

if __name__=="__main__":main()
