#!/usr/bin/env python3
"""Exact published rsID lookup for the entire pre-existing 24-marker panel.
A matching token is not necessarily an association, nor is absence novelty.
Do not re-select markers or assert LD independence/causal genes.
"""
import argparse,collections,hashlib,json,re
from pathlib import Path
from openpyxl import load_workbook

PIN="7c688f864c73ec2e9624bde9f53ac48a9880a6a08e59f5e1266690e243618ba4"
EXPECTED={"rs880315","rs12037987","rs13026829","rs13394970","rs58447744","rs7649237","rs2129981","rs7800053","rs480","rs4744059","rs2258457","rs2767060","rs60401382","rs1973765","rs11045241","rs671","rs7989823","rs7203740","rs1054280","rs34864700","rs76892135","rs12610411","rs6021246","rs5748761"}
def require(ok,msg):
    if not ok:raise ValueError(msg)
def run(manifest,workbook):
    m=json.loads(manifest.read_text())
    require(m["verified_count"]==24 and len(m["rows"])==24,"Complete anchor panel required")
    require(m["historical_blob"]=="b34b4ae84ea4e5b0c72e3982d530a048cc9be5ca","Wrong historical table")
    rows=m["rows"];keys={v["rsid"]:k for k,v in rows.items()}
    require(set(keys)==EXPECTED and all(v["status"]=="ANCHOR_VERIFIED" for v in rows.values()),"Frozen panel identity")
    require(workbook.stat().st_size==756910 and hashlib.sha256(workbook.read_bytes()).hexdigest()==PIN,"Publisher supplement mismatch")
    wb=load_workbook(workbook,read_only=True,data_only=True,keep_links=False)
    matches=collections.defaultdict(list);sheets={}
    try:
        for s in wb.worksheets:
            require((s.max_row or 0)<=100000 and (s.max_column or 0)<=400,"Sheet bounds")
            n=0;header=[];nmatch=0
            for i,row in enumerate(s.iter_rows(values_only=True),1):
                n=i
                vals=["" if v is None else str(v) for v in row]
                if i<=8:header.append(vals[:80])
                hit=collections.defaultdict(list)
                for j,v in enumerate(vals,1):
                    for rs in set(re.findall(r"(?<![A-Za-z0-9])rs[0-9]+(?![A-Za-z0-9])",v,re.I)):
                        rs=rs.lower()
                        if rs in keys:hit[rs].append(j)
                for rs,cols in hit.items():
                    matches[rs].append({"sheet":s.title,"row":i,"matched_columns":cols,"values":vals[:80]});nmatch+=1
            sheets[s.title]={"rows":n,"exact_rs_token_matches":nmatch,"header_first_8_rows":header}
    finally:wb.close()
    result={"status":"EXACT_PUBLISHER_TOKEN_AUDIT_NOT_NOVELTY","anchor_run":45534891070,"panel_count":24,"workbook_sha256":PIN,"sheet_count":len(sheets),"matched_markers":sum(bool(matches[rs]) for rs in EXPECTED),"markers":[{"rsid":rs,"region":keys[rs],"GRCh38_chr":rows[keys[rs]]["chr"],"GRCh38_pos":rows[keys[rs]]["pos"],"historical_ref":rows[keys[rs]]["ref"],"historical_alt":rows[keys[rs]]["alt"],"exact_matches":matches[rs]} for rs in sorted(EXPECTED)],"sheets":sheets,"limits":["Exact rsID occurrence in this one workbook is not evidence of phenotype association. Inspect actual row/header, study and model; variant not found does not establish novelty.","No LD proxy inference, gene mapping, independent replication, causal claim or new significance test. All 24 retained irrespective of observed P."]}
    return result
def main():
    p=argparse.ArgumentParser();p.add_argument("--manifest",type=Path,required=True);p.add_argument("--workbook",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    x=run(a.manifest,a.workbook);a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(x,indent=2,allow_nan=False)+"\n")
    print(json.dumps({"status":x["status"],"panel_count":x["panel_count"],"sheet_count":x["sheet_count"],"matched_markers":x["matched_markers"],"markers":[{"rsid":v["rsid"],"chr":v["GRCh38_chr"],"pos":v["GRCh38_pos"],"hits_by_sheet":dict(collections.Counter(z["sheet"] for z in v["exact_matches"])),"selected_matches":v["exact_matches"][:10] if v["rsid"] in {"rs2129981","rs7203740"} else []} for v in x["markers"]],"limits":x["limits"]},indent=2,allow_nan=False))
if __name__=="__main__":main()
