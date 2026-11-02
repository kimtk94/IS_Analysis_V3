#!/usr/bin/env python3
"""Official PMC article supplementary-source inventory; bounded samples/headers only.

Fetch 4 exact filenames declared by PMC12882097 fullTextXML, and scan DOCX
table row patterns and XLSX sheet metadata for explicit sample/time links and
sample-level raw counts. No unpublished counts or full text redistributed.
"""
import argparse,collections,hashlib,io,json,re,time,urllib.request,zipfile
import xml.etree.ElementTree as ET
from pathlib import Path

PMC="PMC12882097";NUM="12882097"
FILES=["iovs-67-2-13_s001.docx","iovs-67-2-13_s002.docx",
       "iovs-67-2-13_s003.xlsx","iovs-67-2-13_s004.docx"]
BASE="https://pmc.ncbi.nlm.nih.gov/articles"
DEFAULT_LIMIT=48_000_000
TIME=re.compile(r"(?i)(?<![a-z0-9])(?:1|6|48|2)\s*(?:hr|hrs|hour|hours|h|day|days|d)\b")
SAMPLES=re.compile(r"(?i)(?:\b(?:h?TM|SC)[_-]?[0-9]{2,3}(?:[._-][0-9]+)?\b)")
XNS="{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
def get(url,limit=DEFAULT_LIMIT):
    last=None
    for retry in range(2):
        try:
            req=urllib.request.Request(url,headers={
                "User-Agent":"academic-ocular-PMC-supplement-header-integrity/1.0",
                "Accept":"application/vnd.openxmlformats-officedocument.spreadsheetml.sheet,application/zip,*/*"})
            with urllib.request.urlopen(req,timeout=45) as f:
                body=f.read(limit+1);final=f.url;ctype=f.headers.get("Content-Type","")
            if len(body)>limit:raise ValueError("SUPPLEMENT_EXCEEDS_BOUNDED_BYTES")
            if not zipfile.is_zipfile(io.BytesIO(body)):raise ValueError("NOT_OFFICE_ZIP_DOCUMENT")
            with zipfile.ZipFile(io.BytesIO(body)) as z:
                if z.testzip():raise ValueError("SUPPLEMENT_OFFICE_ZIP_CRC_FAILURE")
            return body,final,ctype
        except (OSError,ValueError,TimeoutError) as e:
            last=type(e).__name__+": "+str(e)[:140]
            if retry==0:time.sleep(2)
    raise RuntimeError(last)
def text(node):
    return re.sub(r"\s+"," "," ".join(node.itertext())).strip()
def relevant(x,max_chars=260):
    s=re.sub(r"\s+"," ",x).strip()
    return s[:max_chars]
def docx(data):
    rows=[];headings=[]
    with zipfile.ZipFile(io.BytesIO(data)) as z:
        names=z.namelist()
        if "word/document.xml" not in names:raise ValueError("NO_DOCX_DOCUMENT_XML")
        xml=z.read("word/document.xml")
        if len(xml)>23_000_000:raise ValueError("DOCX_XML_OVER_BOUND")
    root=ET.fromstring(xml)
    W="{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
    ptexts=[]
    for p in root.findall(".//"+W+"p"):
        s="".join((x.text or "") for x in p.findall(".//"+W+"t"))
        if s.strip():ptexts.append(s.strip())
    for tbl in root.findall(".//"+W+"tbl"):
        trs=tbl.findall("./"+W+"tr")
        for row in trs[:300]:
            cells=[]
            for tc in row.findall("./"+W+"tc"):
                pieces=[x.text or "" for x in tc.findall(".//"+W+"t")]
                cells.append(relevant("".join(pieces),180))
            joined=" | ".join(cells)
            if SAMPLES.search(joined) and (TIME.search(joined) or "control" in joined.lower() or "dex" in joined.lower()):
                rows.append({"source":"DOCX_TABLE","cell_count":len(cells),"evidence_row":relevant(joined,300),
                             "sample_tokens":list(dict.fromkeys(SAMPLES.findall(joined)))[:7],
                             "explicit_time_terms":list(dict.fromkeys(TIME.findall(joined)))[:5]})
    examples=[relevant(x,250) for x in ptexts
              if SAMPLES.search(x) and TIME.search(x)]
    return {"paragraph_count":len(ptexts),
        "table_sample_time_or_treatment_rows":rows[:25],
        "paragraph_sample_time_examples":examples[:12],
        "has_sample_level_time_crosswalk":bool(rows or examples),
        "raw_count_sample_matrix_declared":False}
def xlsx(data):
    # Only workbook metadata and first rows, not entire (potentially giant)
    # copyrighted DE expression workbook. We DO NOT publish original cells
    # except limited header labels needed for provenance.
    with zipfile.ZipFile(io.BytesIO(data)) as z:
        files=set(z.namelist())
        if "xl/workbook.xml" not in files:raise ValueError("XLSX_MISSING_WORKBOOK")
        root=ET.fromstring(z.read("xl/workbook.xml"))
        sheets=[{"name":s.get("name","")[:90],"sheetId":s.get("sheetId","")}
                for s in root.findall(".//"+XNS+"sheet")]
        shared=[]
        if "xl/sharedStrings.xml" in files:
            b=z.read("xl/sharedStrings.xml")
            if len(b)>30_000_000:raise ValueError("SHARED_STRINGS_OVER_BOUND")
            sr=ET.fromstring(b)
            for si in sr.findall("./"+XNS+"si"):
                shared.append("".join(t.text or "" for t in si.findall(".//"+XNS+"t")))
        samples=[];sheetstats=[]
        ws=sorted((x for x in files if re.fullmatch(r"xl/worksheets/sheet\d+\.xml",x)),
                  key=lambda x:int(re.search(r"\d+",x).group()))
        for i,name in enumerate(ws):
            headers=[];candidate=[]
            count=0
            with z.open(name) as fh:
                for event,el in ET.iterparse(fh,events=("end",)):
                    if el.tag!=XNS+"row":continue
                    count+=1
                    cells=[]
                    for c in el.findall("./"+XNS+"c"):
                        v=c.findtext(XNS+"v")
                        typ=c.get("t","")
                        if typ=="inlineStr":
                            v="".join(t.text or "" for t in c.findall(".//"+XNS+"t"))
                        elif typ=="s" and v and v.isdigit() and int(v)<len(shared):
                            v=shared[int(v)]
                        if v is not None:
                            cells.append(relevant(v,75))
                    if count<=12:
                        headers.append(cells[:25])
                    joined=" | ".join(cells[:45])
                    if count<=250 and (SAMPLES.search(joined) or (TIME.search(joined) and "dex" in joined.lower())):
                        candidate.append({"row_index":count,
                          "preview":relevant(joined,270),"has_sample_token":bool(SAMPLES.search(joined)),
                          "has_time_and_treatment":bool(TIME.search(joined) and "dex" in joined.lower())})
                    el.clear()
                    if count>=300:break
            sheetstats.append({"sheet_xml":name,
                "workbook_sheet_name_by_ordinal":sheets[i]["name"] if i<len(sheets) else "",
                "sampled_rows":count,"header_rows":headers[:4],
                "sample_or_design_rows":candidate[:20]})
        return {"sheets":sheets,
          "sheet_first_300_rows_summary":sheetstats,
          "sample_level_count_matrix_confirmed":False,
          "note":"Only first 300 worksheet rows audited. Full sample-level raw counts must be positively identified; numerical DE result workbook does not prove such availability."}
def main():
    p=argparse.ArgumentParser();p.add_argument("--out",required=True,type=Path)
    args=p.parse_args();args.out.mkdir(parents=True,exist_ok=True)
    results=[]
    archive={}
    archive_summary={"official_europepmc_url":"https://www.ebi.ac.uk/europepmc/webservices/rest/"+PMC+"/supplementaryFiles",
                     "status":"HOLD_OFFICIAL_SUPPLEMENT_ARCHIVE_UNAVAILABLE"}
    try:
        body,source,ctype=get(archive_summary["official_europepmc_url"],limit=90_000_000)
        with zipfile.ZipFile(io.BytesIO(body)) as z:
            if z.testzip():raise ValueError("ARCHIVE_CRC_FAIL")
            for filename in FILES:
                options=[n for n in z.namelist() if n.rsplit("/",1)[-1]==filename]
                if len(options)==1:
                    raw=z.read(options[0])
                    if not zipfile.is_zipfile(io.BytesIO(raw)):
                        raise ValueError("NOT_OFFICE_DOCUMENT_"+filename)
                    archive[filename]=raw
        archive_summary.update({"status":"PASS_OFFICIAL_ARCHIVE_CRC",
            "sha256":hashlib.sha256(body).hexdigest(),
            "bytes":len(body),"final_url":source,
            "found_named_supplements":sorted(archive)})
    except Exception as e:
        archive_summary["error"]=type(e).__name__+": "+str(e)[:200]
    for name in FILES:
        x={"filename":name,
           "source_declared_in_publisher_PMC_fullTextXML":True,
           "exposed_only_short_table_excerpts_not_full_supplement":True}
        urls=[
            BASE+"/instance/"+NUM+"/bin/"+name,
            BASE+"/"+PMC+"/bin/"+name,
        ]
        data=archive.get(name)
        if data is not None:
            x.update({"official_url":archive_summary.get("final_url"),
                      "content_type":"application/zip; archived "+name,
                      "sha256":hashlib.sha256(data).hexdigest(),
                      "compressed_bytes":len(data),"crc_pass":True})
        else:
            for url in urls:
                try:
                    data,final,ctype=get(url)
                    x.update({"official_url":final,"content_type":ctype,
                              "sha256":hashlib.sha256(data).hexdigest(),
                              "compressed_bytes":len(data),"crc_pass":True})
                    break
                except Exception as e:
                    x.setdefault("attempts",[]).append({"url":url,"error":str(e)[:150]})
        if data:
            try:
                x["source_status"]="PASS_OFFICIAL_SUPPLEMENT_DOWNLOADED_AND_CRC"
                x["audit"]=xlsx(data) if name.endswith("xlsx") else docx(data)
            except Exception as e:
                x["source_status"]="HOLD_SUPPLEMENT_PARSE"
                x["parse_error"]=type(e).__name__+": "+str(e)[:200]
        else:x["source_status"]="HOLD_SUPPLEMENT_NOT_RETRIEVABLE"
        results.append(x)
        print("SUPPLEMENT",name,x["source_status"],
              "bytes",x.get("compressed_bytes"),
              "sheets",len(x.get("audit",{}).get("sheets",[])),flush=True)
    summary={"status":"PASS_SUPPLEMENT_INVENTORY_WITH_LIMITATIONS" if any(
        r["source_status"]=="PASS_OFFICIAL_SUPPLEMENT_DOWNLOADED_AND_CRC" for r in results
        ) else "HOLD_ALL_FOUR_SUPPLEMENTS_NOT_ACCESSIBLE",
      "study":"Mehrotra 2026 IOVS 67(2):13","PMC":PMC,
      "official_supplement_archive_audit":archive_summary,
      "official_supplement_names":FILES,
      "retrieved_crc_passed_n":sum(r["source_status"]=="PASS_OFFICIAL_SUPPLEMENT_DOWNLOADED_AND_CRC" for r in results),
      "direct_processed_raw_counts_confirmed":False,
      "one_to_one_48h_TM_run_strain_timepoint_crosswalk_confirmed":False,
      "raw_FASTQ_downloaded":0,
      "data_source_software_version_pinned_via_checksum":True,
      "source_files_redacted_to_metadata_only":True,
      "manual_action_if_no_sample_sheet":"Search exact official public sample sheet or request source authors' run/strain/timepoint/QC crosswalk"}
    (args.out/"MEHROTRA_4_SUPPLEMENT_OFFICIAL_ACCESS_QC.json").write_text(
        json.dumps(results,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    (args.out/"MEHROTRA_PROCESSED_COUNTS_AND_TIMEPOINT_GATE.json").write_text(
        json.dumps(summary,indent=2,ensure_ascii=False)+"\n",encoding="utf8")
    print("SUPPLEMENT_GATE",json.dumps(summary,ensure_ascii=False),flush=True)
if __name__=="__main__":main()
