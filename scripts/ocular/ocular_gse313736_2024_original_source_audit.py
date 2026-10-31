#!/usr/bin/env python3
"""Audit public 2024 author reports preceding GEO GSE313736.

Extract only compact evidence from three PUBLIC institutional PDFs and the GEO
Series. No genotype data, no author 80-member set invented, no full-text
redistribution. The final state is a source status, not independent validation.
"""
import argparse,hashlib,json,re,time,urllib.parse,urllib.request
from pathlib import Path
import fitz
SOURCES={
 "AVR_2024_annual_report":{
   "url":"https://australianvisionresearch.org/wp-content/uploads/2025/08/2024_AVR-Annual-Report.pdf",
   "limit":28_000_000,
   "terms":["Trabecular Meshwork Cell-Specific Mapping of Intraocular Pressure",
            "63 upregulated","17 downregulated","164,378","6,178","74"]},
 "ANZGS_2024_official_program":{
   "url":"https://anzgsconference.com/wp-content/uploads/2024/02/ANZGS-2024-Program_FINAL.pdf",
   "limit":12_000_000,
   "terms":["FRIDAY 16 FEBRUARY 2024","Expression Quantitative Trait Mapping",
            "Thomas Guinan","17:00"]},
 "UTAS_2024_research_showcase":{
   "url":"https://mooc-media.utas.edu.au/files/7/5/171877818175/CoHM%20Research%20Showcase%202024%20-%20Presentation%20Abstracts%20Booklet.pdf",
   "limit":15_000_000,
   "terms":["Thomas Guinan","Expression Quantitative Trait Loci mapping",
            "Fifty nine samples","demultiplexed","nine single cell datasets"]}
}
GEO="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE313736"
def download(url,cap):
    last=""
    for i in range(3):
        try:
            with urllib.request.urlopen(urllib.request.Request(
                url,headers={"User-Agent":"ocular-primary-source-institutional-archive-audit/1.0"}),
                timeout=44) as resp:
                buf=resp.read(cap+1)
                final=resp.url
            if len(buf)>cap or len(buf)<300:raise ValueError("SIZE_OUT_OF_BOUNDS")
            return buf,final
        except (OSError,ValueError) as e:
            last=type(e).__name__+": "+str(e)[:160]
            if i<2:time.sleep(2)
    raise RuntimeError(last)
def nearby(text,term,width=130):
    m=re.search(re.escape(term),text,re.I)
    return text[max(0,m.start()-width):min(len(text),m.end()+width)].strip() if m else None
def source_audit(key,cfg):
    record={"source":key,"original_url":cfg["url"],
      "status":"HOLD_SOURCE_NOT_RETRIEVED","contains_complete_author_80_gene_list":False}
    try:
        payload,final=download(cfg["url"],cfg["limit"])
        doc=fitz.open(stream=payload,filetype="pdf")
        if len(doc)>160:raise ValueError("UNEXPECTED_PDF_PAGE_COUNT")
        record.update({"status":"GO_OFFICIAL_PUBLIC_PDF_PROVENANCE",
                       "final_url":final,"sha256":hashlib.sha256(payload).hexdigest(),
                       "pdf_bytes":len(payload),"pages":len(doc)})
        evidence=[]
        for term in cfg["terms"]:
            found=[]
            for i,page in enumerate(doc):
                body=re.sub(r"\s+"," ",page.get_text())
                snippet=nearby(body,term)
                if snippet:found.append({"pdf_page_number":i+1,"short_excerpt":snippet[:360]})
            evidence.append({"term":term,"evidence":found[:2]})
        record["verified_short_excerpts"]=evidence
        record["all_required_terms_found"]=all(x["evidence"] for x in evidence)
        if not record["all_required_terms_found"]:
            record["status"]="REVIEW_PDF_TEXT_EXTRACTION_INCOMPLETE"
        entire=" ".join(re.sub(r"\s+"," ",p.get_text()) for p in doc)
        if key=="AVR_2024_annual_report":
            record["author_74_316799_already_reported_2024"]=bool(
                re.search(r"\b74\b",entire) and re.search(r"316,799",entire))
            record["author_63up_17down_already_reported_2024"]=bool(
                re.search(r"63\s+upregulated",entire,re.I) and
                re.search(r"17\s+downregulated",entire,re.I))
            record["threshold_author_down_direction_report_ambiguous"]=True
        elif key=="UTAS_2024_research_showcase":
            record["original_single_cell_genotype_QC_59"]=bool(
                re.search(r"Fifty\s+nine\s+samples\s+passed",entire,re.I))
            record["original_eQTL_work_at_2024_presented_not_source_80_table"]=True
        else:
            first_page=re.sub(r"\\s+"," ",doc[0].get_text())
            scheduled=all(re.search(pattern,first_page,re.I) for pattern in
                [r"16\\s+FEBRUARY\\s+2024",r"17:00",r"17:15",
                 r"Expression Quantitative Trait Mapping",r"Thomas Guinan"])
            record["ANZGS_Feb16_2024_1700_1715_talk_scheduled"]=bool(scheduled)
            record["conference_program_contains_original_80_DEG_IDs"]=False
            record["conference_schedule_proves_80_DEG_reported_in_February"]=False
            if not scheduled:
                record["status"]="REVIEW_CONFERENCE_TIME_OR_SPEAKER_NOT_VERIFIED"
        doc.close()
    except Exception as e:
        record.update({"status":"HOLD_OFFICIAL_PDF_DOWNLOAD_OR_TEXT",
            "error":type(e).__name__+": "+str(e)[:220]})
    return record
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",type=Path,required=True)
    a=ap.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    records=[source_audit(k,v) for k,v in SOURCES.items()]
    r={x["source"]:x for x in records}
    all_three=all(x["status"]=="GO_OFFICIAL_PUBLIC_PDF_PROVENANCE" for x in records)
    summary={
      "status":"GO_THREE_AUTHOR_PRE_GEO_SOURCE_REPORTS_VERIFIED" if all_three else
        "HOLD_ONE_OR_MORE_PRE_GEO_REPORTS_NOT_VERIFIED",
      "public_GEO_2026_series":GEO,
      "public_institutional_sources":records,
      "GEO_74_donors_and_316799_singlets_already_2024":r["AVR_2024_annual_report"].get(
          "author_74_316799_already_reported_2024",False),
      "author_63up_17down_already_2024":r["AVR_2024_annual_report"].get(
          "author_63up_17down_already_reported_2024",False),
      "ANZGS_Feb16_2024_official_program_confirmed":r["ANZGS_2024_official_program"].get(
          "ANZGS_Feb16_2024_1700_1715_talk_scheduled",False),
      "ANZGS_schedule_does_not_prove_80_DEG_results_presented":True,
      "UTAS_genotype_scRNA_QC_59_in_2024":r["UTAS_2024_research_showcase"].get(
          "original_single_cell_genotype_QC_59",False),
      "our_GSE313736_BeamA_57_paired_donors_is_different_pipeline_filter":True,
      "published_59_identical_to_our_57_donors_not_asserted":True,
      "author_full_80_gene_ID_list_available_in_sources":False,
      "original_decrease_threshold_mathematically_defined_by_source":False,
      "no_degnames_or_batch_independence_inferred":True,
      "action":"Obtain primary author differential gene list with gene ID, log2FC, FDR, contrast specification, normalization and original unit-of-replication."}
    (a.out/"GSE313736_2024_ORIGINAL_AUTHOR_PRE_GEO_AUDIT.json").write_text(
      json.dumps(summary,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print("GSE313736_PRE_GEO",summary["status"],"all_three",all_three,
      "original_63_17",summary["author_63up_17down_already_2024"],
      "original_qc_59",summary["UTAS_genotype_scRNA_QC_59_in_2024"],flush=True)
    if not all_three:raise SystemExit("AT_LEAST_ONE_OFFICIAL_PDF_UNVERIFIED")
if __name__=="__main__":main()
