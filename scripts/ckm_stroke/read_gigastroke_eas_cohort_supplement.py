#!/usr/bin/env python3
"""Inspect original GIGASTROKE Supplementary Tables (public .xlsx).

Parse XLSX OOXML directly to inventory study-specific table 1; extract only
small text excerpts and counts from aggregate public supplementary tables.
Do not infer GWAS EAS AIS participant overlap from geography or PGS usage.
"""
import argparse,io,json,re,urllib.request,urllib.parse,zipfile,xml.etree.ElementTree as ET
from pathlib import Path

FILENAME="41586_2022_5165_MOESM4_ESM.xlsx"
URLS=[
 "https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-022-05165-3/MediaObjects/"+FILENAME,
 "https://pmc.ncbi.nlm.nih.gov/articles/instance/9524349/bin/"+FILENAME,
 "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9524349/supplementaryFiles",
]
NS="{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
REL="{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"
PKG="{http://schemas.openxmlformats.org/package/2006/relationships}"
def fetch(url):
 req=urllib.request.Request(url,headers={"User-Agent":"IS_AnalysisV3-GIGASTROKE-study-supplement/1.0","Accept":"application/octet-stream"})
 with urllib.request.urlopen(req,timeout=70) as r:
  if r.status!=200:raise ValueError("source http != 200")
  final=urllib.parse.urlsplit(r.geturl())
  if final.scheme!="https" or not final.hostname:raise ValueError("source redirected outside https")
  sz=r.headers.get("Content-Length")
  if sz and int(sz)>28000000:raise ValueError("source too large")
  buf=r.read(28000001)
 if len(buf)>28000000:raise ValueError("source exceeded hard limit")
 return buf,r.geturl()
def getstr(cell,shared):
 t=cell.get("t")
 if t=="inlineStr":
  return "".join(node.text or "" for node in cell.iter(NS+"t"))
 v=cell.find(NS+"v")
 if v is None:return ""
 raw=v.text or ""
 if t=="s":
  try:return shared[int(raw)]
  except Exception:return ""
 return raw
def xml_book(b):
 with zipfile.ZipFile(io.BytesIO(b)) as z:
  names=z.namelist()
  if any(".." in n.split("/") for n in names):raise ValueError("zip path invalid")
  if sum(x.file_size for x in z.infolist())>95000000:raise ValueError("uncompressed size too large")
  if "xl/workbook.xml" not in names:raise ValueError("not an XLSX workbook")
  shared=[]
  if "xl/sharedStrings.xml" in names:
   root=ET.fromstring(z.read("xl/sharedStrings.xml"))
   shared=["".join(x.text or "" for x in si.iter(NS+"t")) for si in root.findall(NS+"si")]
  wb=ET.fromstring(z.read("xl/workbook.xml"))
  rel=ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
  relmap={r.attrib.get("Id"):r.attrib.get("Target") for r in rel.findall(PKG+"Relationship")}
  report=[]
  for sh in wb.findall(".//"+NS+"sheet"):
   name=sh.attrib.get("name","")
   relid=sh.attrib.get(REL+"id")
   target=relmap.get(relid,"")
   target=target.lstrip("/")
   if not target.startswith("xl/"):target="xl/"+target
   target=target.replace("xl/xl/","xl/")
   if target not in names:continue
   sheet=ET.fromstring(z.read(target))
   rows=[]
   for r in sheet.findall(".//"+NS+"sheetData/"+NS+"row"):
    cells=[(c.get("r"),getstr(c,shared)) for c in r.findall(NS+"c")]
    cells=[(i,(s or "")[:280]) for i,s in cells if (s or "").strip()]
    if cells:rows.append(cells)
   keywords=("east asian","eas","japan","biobank","taiwan","china kadoorie","korea","stroke","cohort","discovery","replication")
   flagged=[]
   for row in rows:
    vals=" | ".join(v for k,v in row)
    if any(t in vals.lower() for t in keywords):
     flagged.append({"range":row[0][0], "cells":row[:15]})
   if name=="Suppl.Table1":
    # Supplementary Table 1 carries separate incident / combined / follow-up cohorts.
    # Return actual source row positions rather than classify any cohort by name.
    report_eas=[]
    for row in rows:
     m=re.search(r"([0-9]+)",row[0][0])
     if m and 45<=int(m.group(1))<=76:
      report_eas.append(row)
   else:report_eas=[]
   report.append({"sheet":name,"nonempty_rows":len(rows),
                  "first_rows":[x[:15] for x in rows[:6]],
                  "keyword_rows_count":len(flagged),
                  "keyword_rows":flagged[:85], "supp_table1_rows_45_to_76":report_eas})
  return report
def main():
 p=argparse.ArgumentParser();p.add_argument("--out",type=Path,required=True);a=p.parse_args()
 r={"schema":"GIGASTROKE-supplement-2022-primary-xlsx-cohort-audit-v1",
    "source_file":FILENAME,"scientific_status":"HOLD",
    "sample_overlap_with_TPMI_confirmed":False}
 for url in URLS:
  try:
   blob,actual=fetch(url)
   if not zipfile.is_zipfile(io.BytesIO(blob)):
    r.setdefault("attempts",[]).append({"url":url,"status":"NOT_XLSX","bytes":len(blob)})
    continue
   report=xml_book(blob)
   if not report:raise ValueError("no sheets in supplement")
   r.update({"source_status":"PASS_ORIGINAL_XLSX_PARSED",
            "source_url":actual,"file_bytes":len(blob),
            "sheets":report,"eas_discovery_cohort_list_inferred":False,
            "eas_ais_discovery_phenotype_proven":False})
   break
  except Exception as e:r.setdefault("attempts",[]).append({"url":url,"error":str(e)[:200]})
 else:r["source_status"]="HOLD_SOURCE_NOT_PARSED"
 a.out.parent.mkdir(parents=True,exist_ok=True)
 a.out.write_text(json.dumps(r,indent=2,ensure_ascii=False)+"\n")
 print(json.dumps({"status":r.get("source_status"),
    "url":r.get("source_url"),
    "sheets":[{"sheet":s["sheet"],"nonempty_rows":s["nonempty_rows"],
               "first_rows":s["first_rows"][:3],"keyword_rows":s["keyword_rows"][:15],
               "EAS_segment":s.get("supp_table1_rows_45_to_76",[])} for s in r.get("sheets",[])[:12]],
    "attempts":r.get("attempts"),"science":"HOLD"},ensure_ascii=False))
if __name__=="__main__":main()
