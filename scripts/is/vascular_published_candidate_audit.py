#!/usr/bin/env python3
"""Inventory exact published VSMC candidate mentions, not a new QTL/coloc."""
import argparse, hashlib, json, re
from pathlib import Path
from openpyxl import load_workbook
GENES={"COL4A1","COL4A2","SH3PXD2A","NEURL1","CALHM2"}
RS={"RS9515201","RS4771674","RS4918058","RS55983834","RS12415501"}
TARGETS=GENES|RS
SHEETS={3,5,6,7,8,13}
def tokens(v):
    if not isinstance(v,str):return set()
    return set(re.findall(r"(?<![A-Z0-9_])[A-Z][A-Z0-9_]*(?![A-Z0-9_])",v.upper()))&TARGETS
def sheet_number(name):
    m=re.fullmatch(r"(?:supplementary|suppl\.?|supp\.?)\s*(?:data|table)?[\s._-]*(\d+)",name.strip(),re.I)
    if not m:m=re.fullmatch(r"(?:data|table|sheet)[\s._-]*(\d+)",name.strip(),re.I)
    return int(m[1]) if m else None
def scan(path):
    with path.open("rb") as f:sha=hashlib.file_digest(f,"sha256").hexdigest()
    wb=load_workbook(path,read_only=True,data_only=True,keep_links=False)
    out={"status":"PUBLISHED_SUPPLEMENT_CANDIDATE_INVENTORY","doi":"10.1038/s41467-026-69273-8","sha256":sha,"bytes":path.stat().st_size,"workbook_sheets":wb.sheetnames,"requested_data_numbers":sorted(SHEETS),"sheets":[],"limitations":["Selected published tables, not full regional QTL data, no new colocalization or MR.","Exact token mentions require interpretation of column semantics; absence of a token in a selected table is not absence of prior literature or biological effect.","No causal gene, subtype-specific or ancestry-specific effect is inferred."]}
    try:
        for k in sorted(SHEETS):
            names=[n for n in wb.sheetnames if sheet_number(n)==k]
            if len(names)!=1:
                out["sheets"].append({"data_number":k,"status":"HOLD_SHEET_IDENTITY","matches":names});continue
            s=wb[names[0]]
            if (s.max_column or 0)>300 or (s.max_row or 0)>250000:
                out["sheets"].append({"data_number":k,"sheet":s.title,"status":"HOLD_SIZE_BOUND"});continue
            counts={t:0 for t in sorted(TARGETS)}; examples=[];n=0;header=[]
            for i,row in enumerate(s.iter_rows(values_only=True),1):
                if i>250000:raise ValueError("Row bound exceeded")
                n=i
                if i<=5:header.append({"row":i,"cells":[str(x or "")[:100] for x in row[:40]]})
                matched=[{"column":j,"value":str(v)[:160],"tokens":sorted(tokens(v))} for j,v in enumerate(row,1) if tokens(v)]
                if not matched:continue
                unique={t for v in matched for t in v["tokens"]}
                for t in unique:counts[t]+=1
                if len(examples)<60:
                    examples.append({"row":i,"matched":matched,"context_first_25":[str(v)[:120] if v is not None else "" for v in row[:25]]})
            out["sheets"].append({"data_number":k,"sheet":s.title,"status":"EXACT_TOKEN_SCAN_COMPLETE","rows_scanned":n,"first_five_rows":header,"exact_token_row_counts":counts,"examples_first_60_matching_rows":examples})
    finally:wb.close()
    if any(s["status"]!="EXACT_TOKEN_SCAN_COMPLETE" for s in out["sheets"]):out["status"]="PARTIAL_PUBLISHED_SUPPLEMENT_INVENTORY"
    return out
def main():
    p=argparse.ArgumentParser();p.add_argument("--input",type=Path);p.add_argument("--out",type=Path);p.add_argument("--selftest",action="store_true");a=p.parse_args()
    if a.selftest:
        assert tokens("COL4A2; rs9515201 / COL4A20")=={"COL4A2","RS9515201"}
        assert tokens("NEURL10")==set()
        assert sheet_number("Supplementary Data 7")==7
        assert sheet_number("Supplementary Data 13")==13
        print("Focused exact token and sheet-name tests passed; no biological inference.");return
    result=scan(a.input);a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(result,indent=2,allow_nan=False)+"\n")
    print(json.dumps({**result,"sheets":[{k:v for k,v in s.items() if k!="examples_first_60_matching_rows"} for s in result["sheets"]]},indent=2))
if __name__=="__main__":main()
