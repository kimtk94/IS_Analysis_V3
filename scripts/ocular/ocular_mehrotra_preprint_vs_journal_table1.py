#!/usr/bin/env python3
"""Does original Mehrotra preprint Table 1 resolve published TM155/TM151 6h discrepancy?

Primary sources only: official EuropePMC public fullTextXML.
A missing preprint Table 1 is HOLD, never a fabricated correction.
"""
import argparse,hashlib,json,re,time,urllib.request,xml.etree.ElementTree as ET
from pathlib import Path
SOURCES={
 "journal":"PMC12882097",
 "preprint":"PMC12776393"}
def text(x):return re.sub(r"\s+"," "," ".join(x.itertext())).strip()
def fetch(pmc):
    url="https://www.ebi.ac.uk/europepmc/webservices/rest/"+pmc+"/fullTextXML"
    last=None
    for i in range(2):
        try:
            with urllib.request.urlopen(urllib.request.Request(
                url,headers={"User-Agent":"ocular-published-vs-preprint-experiment-table-audit/1.0"}),
                timeout=40) as f:payload=f.read(8_000_001)
            if len(payload)>8_000_000:raise ValueError("OVER_LIMIT")
            if b"<article" not in payload:raise ValueError("NOT_ARTICLE_XML")
            return payload,url
        except (OSError,ValueError) as e:
            last=type(e).__name__+": "+str(e)[:160]
            if i==0:time.sleep(2)
    raise RuntimeError(last)
def extract_table1(payload,version):
    root=ET.fromstring(payload)
    results=[]
    for wrapper in root.findall(".//table-wrap"):
        cap=text(wrapper.find("./caption")) if wrapper.find("./caption") is not None else ""
        if "demographic data of donors" in cap.lower() or "demographic and" in cap.lower():
            results.append(wrapper)
    detection="CAPTION_MATCH"
    if not results:
        # Publisher and preprint can label the same table differently. Identify
        # candidate using actual primary-source TM151/TM155 + exposure-duration
        # rows; don't infer a missing table from caption changes.
        for wrapper in root.findall(".//table-wrap"):
            table=wrapper.find(".//table")
            if table is None:continue
            source_rows=[]
            for tr in table.findall(".//tr"):
                cells=[text(x) for x in tr if x.tag in ("th","td")]
                if len(cells)>1:source_rows.append(cells)
            whole=" ".join(" ".join(row) for row in source_rows)
            has_strains=bool(re.search(r"TM\s*151",whole,re.I)) and bool(re.search(r"TM\s*155",whole,re.I))
            has_time=bool(re.search(r"\b(2\s*d|6\s*h|1\s*h)\b",whole,re.I))
            if has_strains and has_time:results.append(wrapper)
        detection="SOURCE_STRAIN_AND_EXPOSURE_ROWS"
    if not results:
        return {"status":"HOLD_NOT_FOUND_PAPER_TABLE1",
                "table_captions":[text(w.find("./caption"))[:160] if w.find("./caption") is not None else ""
                                  for w in root.findall(".//table-wrap")][:12]}
    if len(results)!=1:raise ValueError("MULTIPLE_SOURCE_TABLES")
    table=results[0].find(".//table")
    if table is None:raise ValueError("NO_TABLE_TAG")
    strains={}
    for tr in table.findall(".//tr"):
        cells=[text(x) for x in tr if x.tag in ("th","td")]
        if not cells:continue
        m=re.search(r"\b(TM|SC)\s*(\d{2,3})\b",cells[0],re.I)
        if not m:continue
        name=m.group(1).upper()+m.group(2)
        if name in strains:raise ValueError("DUPLICATE_SOURCE_STRAIN")
        raw=cells[-1]
        units=re.findall(r"\b(?:2\s*(?:days?|d)|48\s*(?:hours?|h)|1\s*(?:hours?|h)|6\s*(?:hours?|h))\b",raw,re.I)
        times=sorted(set("2Days" if (v.startswith("2d") or v.startswith("48")) else
                         "1Hr" if v.startswith("1") else "6Hrs"
                         for v in (re.sub(r"\s+","",t.lower()) for t in units)))
        strains[name]={"timepoints":times,"raw_exposure_text":raw[:120]}
    return {"status":"GO_PAPER_STRAIN_TABLE1_READ","detection":detection,"strains":strains,
            "TM_count":sum(x.startswith("TM") for x in strains),
            "SCE_count":sum(x.startswith("SC") for x in strains)}
def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=Path,required=True)
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    pages={}
    for key,pmcid in SOURCES.items():
        d={"source_kind":key,"PMC":pmcid}
        try:
            raw,url=fetch(pmcid)
            d.update({"fullTextXML_url":url,"sha256":hashlib.sha256(raw).hexdigest()})
            d.update(extract_table1(raw,key))
        except Exception as e:
            d.update({"status":"HOLD_SOURCE_XML_UNAVAILABLE",
                      "error":type(e).__name__+": "+str(e)[:150]})
        pages[key]=d
    j=pages["journal"];pre=pages["preprint"]
    compared=[]
    if j.get("status")=="GO_PAPER_STRAIN_TABLE1_READ" and pre.get("status")=="GO_PAPER_STRAIN_TABLE1_READ":
        for id in ("TM151","TM155"):
            compared.append({"strain":id,"preprint":pre["strains"].get(id),
                             "journal":j["strains"].get(id),
                             "changed":pre["strains"].get(id,{}).get("timepoints")!=j["strains"].get(id,{}).get("timepoints"),
                             "raw_duration_phrase_changed_only":pre["strains"].get(id,{}).get("raw_exposure_text")!=j["strains"].get(id,{}).get("raw_exposure_text")})
    result={"status":"GO_PREPRINT_JOURNAL_TABLE1_COMPARISON" if compared else
                        "HOLD_PREPRINT_TABLE1_NOT_COMPARABLE",
            "sources":pages,
            "review_genes":["TM151","TM155"],
            "strain_timepoint_comparison":compared,
            "source_s003_6h_TM151_only":True,
            "biological_duration_labels_compared_after_h_hour_days_d_normalization":True,
            "sample_specific_SRA_timepoint_evidence":False,
            "cannot_infer_journal_table_typos_without_author_confirmation":True,
            "no_private_or_raw_FASTQ_data":True}
    (a.out/"MEHROTRA_PREPRINT_JOURNAL_TABLE1_PROVENANCE.json").write_text(
        json.dumps(result,indent=2,ensure_ascii=False)+"\n")
    print("PAPER_VERSION_TABLE1",result["status"],"changes",compared,flush=True)
if __name__=="__main__":main()
