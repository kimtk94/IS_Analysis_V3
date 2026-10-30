#!/usr/bin/env python3
"""Auditable numeric patient-key join between Xing 2023 proteomics and clinical S1.

Processes publicly released xlsx files within GitHub Actions memory. Never
prints or commits patient keys, patient attributes, phenotype columns,
protein expression values or raw downloaded spreadsheet files.
This is NOT an RNA/protein join and cannot promote V02.
"""
import argparse
import collections
import datetime as dt
import io
import json
import pathlib
import re
import zipfile
import xml.etree.ElementTree as ET
import probe_pmc2026_cloud_proteomics as pmc

def discover_media(meta,target):
    original=pmc.TARGET
    try:
        pmc.TARGET=target
        return pmc.parse_article_metadata(meta)["candidate_urls_in_memory"]
    finally:
        pmc.TARGET=original

def get_file(meta,target,maxbytes=pmc.MAX_XLSX):
    urls=discover_media(meta,target)
    if len(urls)!=1:raise ValueError("expected one official media URL")
    raw,code,ctype=pmc.bounded(urls[0],maxbytes,timeout=20)
    if not zipfile.is_zipfile(io.BytesIO(raw)):
        raise ValueError("not an XLSX workbook")
    return raw

def sheet_rows(z,ordinal,limit=170):
    shared=pmc.read_shared_strings(z)
    target="xl/worksheets/sheet"+str(ordinal)+".xml"
    rows=[]
    with z.open(target) as stream:
        for _,row in ET.iterparse(stream,events=("end",)):
            if row.tag!=pmc.NS+"row":continue
            rownum=int(row.get("r") or "0")
            cells={}
            for c in row.findall(pmc.NS+"c"):
                col=pmc._excel_column_number(c.get("r",""))
                if col<1 or col>307:continue
                v=c.find(pmc.NS+"v")
                t=c.get("t")
                if t=="inlineStr":
                    val="".join(x.text or "" for x in c.iter(pmc.NS+"t"))
                elif v is None:continue
                elif t=="s":
                    idx=int(v.text)
                    val=shared[idx] if 0<=idx<len(shared) else ""
                else:val=v.text or ""
                cells[col]=val.strip() if isinstance(val,str) else str(val)
            rows.append((rownum,cells))
            row.clear()
            if len(rows)>=limit:break
    return rows

def canon_id(raw):
    s=str(raw or "").strip()
    # Avoid fuzzy stripping prefixes: only pure integer codes are considered.
    if re.fullmatch(r"\d{1,8}(?:\.0+)?",s):
        return str(int(float(s)))
    return None

def match_clinical_ids(protein_headers,clinical_rows):
    # Source RNA, protein values and actual case keys never serialized.
    pairs=collections.defaultdict(set)
    for x in protein_headers:
        c=pmc.classify_patient_label(x)
        if c:
            rawkey,kind,_=c
            norm=canon_id(rawkey)
            if norm:pairs[norm].add(kind)
    paired_keys={key for key,value in pairs.items() if value=={"T","N"}}
    fields=collections.defaultdict(set)
    records_per_field=collections.Counter()
    # Detailed clinical workbook convention: rows 1-3 title/header; 4+ data.
    for rownum,cells in clinical_rows:
        if rownum<=3:continue
        for col,raw in cells.items():
            key=canon_id(raw)
            if key:
                records_per_field[col]+=1
                fields[col].add(key)
    rankings=sorted(
        ((len(keys & paired_keys),col,len(keys),records_per_field[col]) for col,keys in fields.items()),
        reverse=True)
    highest=rankings[0] if rankings else (0,None,0,0)
    runnerup=rankings[1][0] if len(rankings)>1 else 0
    return {
        "proteome_header_paired_numeric_case_keys":len(paired_keys),
        "clinic_rows_at_or_after_4":sum(row>=4 for row,_ in clinical_rows),
        "clinical_columns_with_numeric_code_candidates":len(fields),
        "best_clinical_numeric_column_index":highest[1],
        "best_clinical_column_match_count":highest[0],
        "best_clinical_column_unique_numeric_codes":highest[2],
        "best_clinical_column_numeric_cells":highest[3],
        "next_best_clinical_column_match_count":runnerup,
        "all_protein_pairs_match_clinical_candidate":len(paired_keys)>0 and highest[0]==len(paired_keys),
        "patient_ids_or_phenotypes_logged":False,
        "clinical_patient_column_annotated_and_confirmed":False,
        "protein_rna_quantitative_join":False,
        "status":"IN_MEMORY_ID_COLUMN_CANDIDATE_COMPARISON",
        "warning":"A matching numeric ID column requires independent clinical table field validation; no RNA expression matrix or survival analysis audited."
    }

def audit():
    out={"provenance":"https://pmc-oa-opendata.s3.amazonaws.com/README.txt",
         "clinical_source":"mmc2.xlsx","proteome_source":"mmc3.xlsx",
         "status":"UNVERIFIED","V02":"BLOCKED"}
    try:
        metadata,status,_=pmc.bounded(pmc.BASE+"metadata/"+pmc.PMCID+".1.json",pmc.MAX_META)
        p=get_file(metadata,"mmc3.xlsx")
        c=get_file(metadata,"mmc2.xlsx")
        with zipfile.ZipFile(io.BytesIO(p)) as zp,zipfile.ZipFile(io.BytesIO(c)) as zc:
            header_rows=sheet_rows(zp,2,limit=2)
            if len(header_rows)<2:raise ValueError("protein header rows missing")
            protein_headers=[value for col,value in sorted(header_rows[1][1].items()) if col>=4]
            clin=sheet_rows(zc,3,limit=170)
            result=match_clinical_ids(protein_headers,clin)
        out["status"]="PUBLIC_PROTEOME_VS_CLINICAL_ID_AUDITED"
        out["id_overlap"]=result
    except Exception as exc:
        out["error_type"]=type(exc).__name__
    return out

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",default="results/hcc/feasibility/xing2023_proteome_clinical_id_overlap.json")
    args=p.parse_args()
    report={"checked_utc":dt.datetime.now(dt.timezone.utc).isoformat(),"audit":audit()}
    path=pathlib.Path(args.out);path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n")
    print(json.dumps(report,indent=2,ensure_ascii=False))
if __name__=="__main__":main()
