#!/usr/bin/env python3
"""Bounded audit of Ng 2022 *public supplementary XLSX* matrix layouts.

Only sheet dimensions, raw-file SHA256, safe header-family tags and numerical
geometry are persisted. No patient/biopsy/sample IDs, raw data cells, full XLSX
or public signed URLs are written to CI artifacts. This is an RNA quantification
access screen, NOT proof of absence or presence of a matched expression matrix.
"""
import argparse
import datetime as dt
import hashlib
import io
import json
import pathlib
import re
import urllib.error
import zipfile
import xml.etree.ElementTree as ET

from probe_quant_matrix_manifest import SPRINGER, read_bounded
from probe_ng_gao_public import _cell_value, _column_index

NS="{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
MAX_ARCHIVE=14*1024*1024
MAX_XML=24*1024*1024
MAX_SHEETS=30
# Supplement 1 is patient metadata audited by a different purpose-built
# tool. Keep sensitive ID tables out of this expression-file audit.
FILES={f"Supplementary_Data_{i}":f"41467_2022_29960_MOESM{i+3}_ESM.xlsx"
       for i in range(2,15)}
FILES["Source_Data"]="41467_2022_29960_MOESM19_ESM.xlsx"
TAGS={
 "gene_feature":r"\bgene(?:s| id| symbol| name)?\b|\bensembl\b|\bentrez\b|\btranscript\b",
 "protein_feature":r"\bproteins?\b|\baccession\b|\buniprot\b|\bpeptide\b|\bphospho",
 "sample_metadata":r"\bpatient\b|\bbiopsy\b|\bsample\b|\baliquot\b",
 "expression_measure":r"\bcount(?:s)?\b|\btpm\b|\bfpkm\b|\brpkm\b|\brsem\b|\bexpression\b",
 "summary_stat":r"\blog2?(?:fc| fold)\b|\bp[- ]?value\b|\bpadj\b|\bq.value\b|\bfdr\b|\bfold.?change\b",
 "pathway_enrichment":r"\bpathway\b|\bgsea\b|\benrichment\b|\bgo term\b",
}
def header_categories(strings):
    return sorted({tag for value in strings for tag,pat in TAGS.items()
                   if re.search(pat,str(value).lower())})

def parse_worksheets(raw):
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        pathlist=sorted((x for x in z.namelist()
                         if re.fullmatch(r"xl/worksheets/sheet\d+\.xml",x)),
                        key=lambda x:int(re.search(r"sheet(\d+)",x).group(1)))
        shared=[]
        if "xl/sharedStrings.xml" in z.namelist() and z.getinfo("xl/sharedStrings.xml").file_size<=MAX_XML:
            strs=ET.fromstring(z.read("xl/sharedStrings.xml"))
            shared=["".join(t.text or "" for t in item.iter(NS+"t")) for item in strs.findall(NS+"si")]
        results=[]
        for path in pathlist[:MAX_SHEETS]:
            f=z.getinfo(path)
            result={"sheet_ordinal":len(results)+1,
                    "xml_uncompressed_bytes":f.file_size,
                    "status":"UNVERIFIED",
                    "matrix_quantitative_same_patient_join":False}
            if f.file_size>MAX_XML:
                result["status"]="XML_TOO_LARGE_TO_INSPECT"
                results.append(result)
                continue
            root=ET.fromstring(z.read(path))
            dim=root.find(NS+"dimension")
            sheet=(root.find(NS+"sheetData"))
            rows=list(sheet) if sheet is not None else []
            bounds=(dim.attrib.get("ref") or "") if dim is not None else ""
            # Dimensions are only used for geometry; never disclose raw labels.
            match=re.search(r":([A-Z]+)(\d+)$",bounds)
            maxcol=_column_index(match.group(1))+1 if match else None
            maxrow=int(match.group(2)) if match else None
            head=[]
            for row in rows[:2]:
                row_values=[]
                for cell in row.findall(NS+"c")[:50]:
                    label=_cell_value(cell,shared)
                    if label:row_values.append(label)
                head.extend(row_values)
            result.update({
                "status":"PUBLIC_WORKSHEET_GEOMETRY_INSPECTED",
                "excel_row_count":maxrow,
                "excel_column_count":maxcol,
                "xml_row_count":len(rows),
                "header_family_tags_allowlist":header_categories(head),
                "wide_row_hint_at_least_20_columns":bool(maxcol is not None and maxcol>=20),
                "candidate_gene_by_sample_matrix_geometry":bool(
                    maxrow is not None and maxrow>=500 and maxcol is not None and maxcol>=20),
                "header_semantics_verified":False,
                "sample_identifier_columns_verified":False
            })
            results.append(result)
        return results,len(pathlist)

def inspect_file(data):
    if not zipfile.is_zipfile(io.BytesIO(data)):
        raise ValueError("source not valid XLSX")
    sheets,total=parse_worksheets(data)
    return {
      "status":"PUBLIC_SUPPLEMENT_XLSX_LAYOUT_SCREENED",
      "bytes":len(data),
      "sha256":hashlib.sha256(data).hexdigest(),
      "sheets_total":total,
      "sheets_audited":len(sheets),
      "sheet_summaries":sheets,
      "no_patient_identifiers_or_quant_values_saved":True,
      "same_patient_rna_protein_quantitative_join_verified":False,
    }
def audit():
    records={}
    for key,name in FILES.items():
        rec={"supplement":key,"public_file":name,"status":"UNVERIFIED"}
        try:
            data=read_bounded(SPRINGER+name,limit=MAX_ARCHIVE,timeout=18)
            rec.update(inspect_file(data))
        except (urllib.error.HTTPError,urllib.error.URLError,
                OSError,ValueError,ET.ParseError,zipfile.BadZipFile) as err:
            rec["status"]="NOT_SCREENED"
            rec["error_type"]=type(err).__name__
        records[key]=rec
        print(key,rec["status"],"sheets",rec.get("sheets_total"),flush=True)
    fetched=[k for k,r in records.items() if r["status"]=="PUBLIC_SUPPLEMENT_XLSX_LAYOUT_SCREENED"]
    geometry=[
        {"supplement":k,"sheet_ordinal":s["sheet_ordinal"],
         "row_count":s.get("excel_row_count"),"column_count":s.get("excel_column_count")}
        for k,r in records.items() for s in r.get("sheet_summaries",[])
        if s.get("candidate_gene_by_sample_matrix_geometry") is True
    ]
    return {
      "source_paper":"Ng et al. 2022 Nature Communications doi:10.1038/s41467-022-29960-8",
      "scope":"Public supplementary spreadsheets only; no EGA or original per-patient raw data",
      "supplements_targeted":len(FILES),
      "supplements_downloaded_and_layout_screened":len(fetched),
      "supplements_not_screened":sorted(set(FILES)-set(fetched)),
      "wide_matrix_geometry_hints":geometry,
      "absence_of_public_rna_quant_matrix_proven":False,
      "matched_rna_protein_patient_quant_matrix_verified":False,
      "EGA_rna_access":"RESTRICTED",
      "records":records
    }
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",default="results/hcc/feasibility/ng_supplement_matrix_depth.json")
    args=p.parse_args()
    report={"checked_utc":dt.datetime.now(dt.timezone.utc).isoformat(),
            "audit":audit()}
    out=pathlib.Path(args.out);out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    # A single inaccessible supplement is a per-file UNVERIFIED state, not
    # a false claim of absence. But no file at all = execution failure.
    if report["audit"]["supplements_downloaded_and_layout_screened"]==0:
        raise SystemExit(1)
if __name__=="__main__":main()
