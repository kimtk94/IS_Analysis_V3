#!/usr/bin/env python3
"""Stage26 Ng 2022 publisher xlsx Figure2/Source Data gene-set provenance probe.

Only aggregate source file SHA / sheet row counts and ENSG-like evidence.
No patient, gene, protein IDs, row-level values or Excel payloads in artifact.
Counts resembling 5,481 are NOT claimed to reproduce the published universe.
"""
import argparse
import datetime as dt
import hashlib
import io
import json
import pathlib
import re
import urllib.error
import urllib.request
import xml.etree.ElementTree as ET
import zipfile

BASE="https://media.springernature.com/original/springer-static/esm/art%3A10.1038%2Fs41467-022-29960-8/MediaObjects/"
PMC="https://pmc.ncbi.nlm.nih.gov/articles/PMC9068765/bin/"
# Figure 2 has downloadable Source Data XLSX; Ng Supplementary Data 4 is
# specifically referenced in main text for gene-correlation pathways.
FILES={
    "publisher_Source_Data":"41467_2022_29960_MOESM19_ESM.xlsx",
    "publisher_Supp_Data_4":"41467_2022_29960_MOESM7_ESM.xlsx",
    "publisher_Supp_Data_2":"41467_2022_29960_MOESM5_ESM.xlsx",
}
NS="{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
MAX_FILE=55*1024*1024
MAX_SHEET_XML=190*1024*1024
ENSG=re.compile(r"ENSG[0-9]{11}(?:\.[0-9]+)?")
def bounded_download(url):
    req=urllib.request.Request(url,headers={
        "User-Agent":"HCC-stage26-public-paper-supplement-integrity/1.0",
        "Accept":"application/octet-stream"})
    with urllib.request.urlopen(req,timeout=65) as resp:
        payload=resp.read(MAX_FILE+1)
    if len(payload)>MAX_FILE or len(payload)<500:
        raise ValueError("source XLSX is empty or over bounded size")
    if not zipfile.is_zipfile(io.BytesIO(payload)):
        raise ValueError("download not XLSX ZIP")
    return payload
def workbook_profile(payload, private_candidate_out=None):
    sheets=[]
    with zipfile.ZipFile(io.BytesIO(payload)) as z:
        members=set(z.namelist())
        names=sorted(p for p in members
                     if re.fullmatch(r"xl/worksheets/sheet\d+\.xml",p))
        if not names:raise ValueError("no actual XLSX worksheets")
        sheet_names={}
        if "xl/workbook.xml" in members:
            root=ET.fromstring(z.read("xl/workbook.xml"))
            for node in root.iter(NS+"sheet"):
                try:
                    sheet_index=int(node.get("sheetId","0"))
                except ValueError:continue
                name=node.get("name","")
                if 0<sheet_index<=100 and len(name)<=100:
                    sheet_names[sheet_index]=name
        shared=[]
        figure2_ids=None
        figure2_candidate_sheet=None
        cna_15272_corrob=False
        if "xl/sharedStrings.xml" in members:
            if z.getinfo("xl/sharedStrings.xml").file_size>MAX_SHEET_XML:
                raise ValueError("shared strings entry too big")
            with z.open("xl/sharedStrings.xml") as handle:
                for _,el in ET.iterparse(handle,events=("end",)):
                    if el.tag==NS+"si":
                        shared.append("".join(t.text or "" for t in el.iter(NS+"t")))
                        el.clear()
        # Unstructured scan rather than making up a gene identifier column.
        for entry in names[:40]:
            if z.getinfo(entry).file_size>MAX_SHEET_XML:
                sheets.append({"sheet_index":len(sheets)+1,
                   "status":"SHEET_XML_EXCEEDS_BOUND"})
                continue
            ids=set()
            ensg_rows=0
            observed_rows=0
            nonempty_rows=0
            with z.open(entry) as handle:
                for _,row in ET.iterparse(handle,events=("end",)):
                    if row.tag!=NS+"row":continue
                    observed_rows+=1
                    if observed_rows>300000:
                        raise ValueError("oversized sheet row count")
                    row_genes=set()
                    nonempty=False
                    for cell in row.findall(NS+"c"):
                        typ=cell.get("t")
                        if typ=="inlineStr":
                            value="".join(t.text or "" for t in cell.iter(NS+"t"))
                        else:
                            vn=cell.find(NS+"v")
                            v=vn.text if vn is not None else ""
                            if typ=="s" and v is not None and v.isdecimal():
                                n=int(v)
                                value=shared[n] if n<len(shared) else ""
                            else:value=v or ""
                        if value:nonempty=True
                        for g in ENSG.findall(value):
                            row_genes.add(g.split(".")[0])
                    if row_genes:
                        ids.update(row_genes);ensg_rows+=1
                    if nonempty:nonempty_rows+=1
                    row.clear()
            sheet_index=len(sheets)+1
            if (private_candidate_out is not None and sheet_index==28 and
                observed_rows==5482 and len(ids)==5481 and ensg_rows==5481):
                figure2_ids=sorted(ids)
                figure2_candidate_sheet=sheet_index
            if (private_candidate_out is not None and sheet_index==27 and
                observed_rows==15273 and len(ids)==15272 and ensg_rows==15272):
                cna_15272_corrob=True
            sheets.append({
              "sheet_index":sheet_index,
              "workbook_sheet_label":sheet_names.get(sheet_index),
             
              "status":"SCAN_COMPLETED",
              "observed_rows_n":observed_rows,
              "nonempty_rows_n":nonempty_rows,
              "rows_with_ENSG_like_identifier_n":ensg_rows,
              "unique_ENSG_like_identifiers_n":len(ids),
              "candidate_exact_5481_count_only":len(ids)==5481,
              "candidate_exact_5490_count_only":len(ids)==5490,
              "ENSG_union_digest_not_gene_membership_confirmation":
                  hashlib.sha256("\n".join(sorted(ids)).encode()).hexdigest()
                  if ids else None,
            })
    if private_candidate_out is not None:
        if figure2_ids is None or not cna_15272_corrob:
            raise ValueError("the original 5481 IDs cannot be anchored to both figure panels")
        # The output is STRICTLY runner temporary and is NEVER an artifact.
        destination=pathlib.Path(private_candidate_out)
        destination.parent.mkdir(parents=True,exist_ok=True)
        destination.write_text("\n".join(figure2_ids)+"\n",encoding="utf-8")
        destination.chmod(0o600)
    return sheets
def probe(private_genes_out=None):
    result={"schema":"hcc-stage26-ng-paper-gene-source-probe-v1",
       "checked_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
       "paper_doi":"10.1038/s41467-022-29960-8",
       "source_5481_correlation_benchmark":5481,
       "source_5490_comparable_gene_benchmark":5490,
       "sources":{},
       "exact_published_5481_gene_membership_verified":False,
       "original_supplement_provenance_is_not_original_2022_mapping_reproduction":True,
       "patient_gene_and_protein_identifiers_exported":False,
       "raw_original_XLSX_in_artifact":False,
    }
    for role,filename in FILES.items():
        record={"filename":filename,"status":"UNVERIFIED","tested_official_urls_n":0}
        for url in (BASE+filename,PMC+filename):
            record["tested_official_urls_n"]+=1
            try:
                content=bounded_download(url)
                sheets=workbook_profile(content,
                    private_candidate_out=private_genes_out if role=="publisher_Source_Data" else None)
                record.update({"status":"PUBLIC_XLSX_CONTENT_AUDITED",
                   "source_url":url,
                   "sha256":hashlib.sha256(content).hexdigest(),
                   "size_bytes":len(content),
                   "sheet_count_inspected":len(sheets),
                   "sheets":sheets})
                break
            except (OSError,ValueError,zipfile.BadZipFile,ET.ParseError,
                    urllib.error.HTTPError,urllib.error.URLError) as e:
                record["last_error_class"]=type(e).__name__
        result["sources"][role]=record
    return result
def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--out",required=True)
    ap.add_argument("--private-gene-ids-out",default=None)
    a=ap.parse_args()
    p=pathlib.Path(a.out);p.parent.mkdir(parents=True,exist_ok=True)
    report=probe(a.private_gene_ids_out)
    p.write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))
if __name__=="__main__":main()
