#!/usr/bin/env python3
"""Aggregate-only forensic audit of Ng 2022 Figure2 and SuppFig13 source XLSX.

Identifies XLSX sheets by workbook relationships, and numeric column profiles
WITHOUT uploading gene IDs, sheet cell values, human identifiers or workbooks.
A numeric median close to paper 0.287 does NOT establish gene-set identity.
"""
import argparse
import datetime as dt
import hashlib
import io
import json
import pathlib
import re
import statistics
import urllib.request
import xml.etree.ElementTree as ET
import zipfile

URL=("https://media.springernature.com/original/springer-static/esm/"
     "art%3A10.1038%2Fs41467-022-29960-8/MediaObjects/"
     "41467_2022_29960_MOESM19_ESM.xlsx")
SHA="16376438aa16a09a13db417a4ae1f4cc16b87ec2c47630aeae651004501be420"
N="{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
REL="{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"
SHEETS={"Fig2a","Fig2b","SuppFig13f","SuppFig13g"}
NUMERIC_LIMIT=12
MAX_BYTES=8000000
def header_type(text):
    if re.fullmatch(r"ENSG[0-9]{11}(?:\.[0-9]+)?",text):
        return "ENS_IDENTIFIER_REDACTED"
    if re.search(r"ensembl",text,re.I):return "ENSEMBL_HEADER"
    if re.search(r"gene|symbol",text,re.I):return "GENE_HEADER"
    if re.search(r"spearman|correl|rho",text,re.I):return "CORRELATION_HEADER"
    if re.search(r"p[-_ ]?value|adj[._-]?p|fdr|q[-_ ]?value",text,re.I):
        return "SIGNIFICANCE_HEADER"
    if re.search(r"cna|copy.number",text,re.I):return "CNA_HEADER"
    if re.search(r"rna|transcript",text,re.I):return "RNA_HEADER"
    if re.search(r"protein|proteom",text,re.I):return "PROTEIN_HEADER"
    if text.strip()=="" :return "EMPTY"
    return "UNCLASSIFIED_REDACTED"
def fetch():
    req=urllib.request.Request(URL,headers={"User-Agent":"HCC-Stage27-SourceFigure2-Sheet-Audit/1.0"})
    with urllib.request.urlopen(req,timeout=55) as resp:
        b=resp.read(MAX_BYTES+1)
    if len(b)>MAX_BYTES or hashlib.sha256(b).hexdigest()!=SHA:
        raise ValueError("original published Source Data SHA mismatch")
    return b
def describe(data):
    out={"schema":"hcc-stage27-ng-2022-figure2-xlsx-structure-v1",
         "official_source_url":URL,"source_sha256":SHA,
         "source_size_bytes":len(data),"sheet_evidence":{},
         "source_fig2b_5481_Ensembl_set_identity_confirmed":False,
         "source_2022_same_gene_feature_quantitation_reproduced":False,
         "source_patient_gene_or_cell_values_exported":False,
         "note":"Only metadata types and aggregate numeric distribution are exported."}
    with zipfile.ZipFile(io.BytesIO(data)) as z:
        shared=[]
        if "xl/sharedStrings.xml" in z.namelist():
            with z.open("xl/sharedStrings.xml") as f:
                for _,e in ET.iterparse(f,events=("end",)):
                    if e.tag==N+"si":
                        shared.append("".join(t.text or "" for t in e.iter(N+"t")))
                        e.clear()
        w=ET.fromstring(z.read("xl/workbook.xml"))
        rel=ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
        relmap={}
        for e in rel:
            m=re.search(r"(?:^|/)sheet([0-9]+)\.xml$",e.get("Target",""))
            if m:relmap[e.get("Id")]=int(m.group(1))
        mapping={}
        for sh in w.iter(N+"sheet"):
            name=sh.get("name","")
            if name in SHEETS and sh.get(REL) in relmap:
                mapping[name]=relmap[sh.get(REL)]
        if set(mapping)!=SHEETS or mapping["Fig2b"]==mapping["SuppFig13g"]:
            raise ValueError("workbook original Figure/SuppFig mapping cannot be proven")
        for name,sheetid in mapping.items():
            path=f"xl/worksheets/sheet{sheetid}.xml"
            if z.getinfo(path).file_size>75_000_000:
                raise ValueError("oversized sheet")
            header=[]
            total=0; data_rows=0; ensg_rows=0
            numerics=[[] for _ in range(NUMERIC_LIMIT)]
            data_rows_with_numeric=0
            with z.open(path) as f:
                for _,e in ET.iterparse(f,events=("end",)):
                    if e.tag!=N+"row":continue
                    total+=1
                    if total>20000:raise ValueError("unexpected sheet rows")
                    values={}
                    for cell in e.findall(N+"c"):
                        key=re.match(r"[A-Z]+",cell.get("r",""))
                        if not key:continue
                        idx=0
                        for ch in key.group():idx=idx*26+ord(ch)-64
                        col=idx-1
                        if col>=NUMERIC_LIMIT or col<0:continue
                        kind=cell.get("t")
                        if kind=="inlineStr":
                            val="".join(x.text or "" for x in cell.iter(N+"t"))
                        else:
                            v=cell.find(N+"v");v=v.text if v is not None else ""
                            if kind=="s" and v and v.isdecimal():
                                j=int(v);val=shared[j] if j<len(shared) else ""
                            else:val=v or ""
                        values[col]=val
                    if total==1:
                        header=[header_type(values.get(i,"")) for i in range(NUMERIC_LIMIT)]
                    else:
                        data_rows+=1
                        if any(re.fullmatch(r"ENSG[0-9]{11}(?:\.[0-9]+)?",x) for x in values.values()):
                            ensg_rows+=1
                        had_numeric=False
                        for i,v in values.items():
                            try:
                                x=float(v)
                                if not (x==x and abs(x)<1e20):continue
                                numerics[i].append(x);had_numeric=True
                            except ValueError:pass
                        if had_numeric:data_rows_with_numeric+=1
                    e.clear()
            columns=[]
            for i,a in enumerate(numerics):
                if len(a)==0:continue
                # Values are summarized, not exported.
                columns.append({"column_index_one_based":i+1,
                    "header_category":header[i] if i<len(header) else "UNCLASSIFIED_REDACTED",
                    "numeric_cells_n":len(a),
                    "fraction_between_minus1_plus1":round(sum(-1<=x<=1 for x in a)/len(a),6),
                    "numeric_median":round(statistics.median(a),6),
                    "potential_correlation_distribution_not_verified":
                        len(a)>100 and sum(-1<=x<=1 for x in a)/len(a)>.95})
            out["sheet_evidence"][name]={
               "worksheet_physical_index":sheetid,
               "total_xml_row_n":total,"data_rows_n":data_rows,
               "rows_with_Ensembl_strings_first12_columns":ensg_rows,
               "data_rows_with_at_least_one_numeric_cell_n":data_rows_with_numeric,
               "first12_header_categories":header,
               "numeric_column_aggregate":columns,
               "sheet_rows_cannot_be_assumed_to_match_gene_universe":True}
    return out
def main():
    p=argparse.ArgumentParser();p.add_argument("--out",required=True)
    a=p.parse_args()
    out=describe(fetch());out["checked_utc"]=dt.datetime.now(dt.timezone.utc).isoformat()
    f=pathlib.Path(a.out);f.parent.mkdir(parents=True,exist_ok=True)
    f.write_text(json.dumps(out,indent=2)+"\n",encoding="utf8")
    print(json.dumps(out,indent=2))
if __name__=="__main__":main()
