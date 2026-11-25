#!/usr/bin/env python3
"""Audit exact GIGASTROKE original discovery cohort totals for AIS/CES/LAS/SVS.
Cohort aggregates do not certify per-SNP N, equal ascertainment or independence.
"""
import argparse,hashlib,json,math,re
from pathlib import Path
from openpyxl import load_workbook
from openpyxl.utils import get_column_letter
PIN="7c688f864c73ec2e9624bde9f53ac48a9880a6a08e59f5e1266690e243618ba4"
TOTAL_N={"EUR":{"AIS":1296908,"CES":1245612,"LAS":1241207,"SVS":1241619},"EAS":{"AIS":256274,"CES":238168,"LAS":238977,"SVS":242774}}
def require(x,m):
    if not x:raise ValueError(m)
def norm(x):return re.sub(r"\s+"," ",str(x or "").strip()).upper()
def n(x):
    if isinstance(x,(int,float)) and not isinstance(x,bool) and math.isfinite(x) and x>=0 and int(x)==x:return int(x)
    if isinstance(x,str) and re.fullmatch(r"\d+(?:,\d{3})*",x.strip()):return int(x.replace(",",""))
    if norm(x) in {"","NA","N/A","-","–","—","ND"}:return None
    raise ValueError("Non-count source cell: "+str(x)[:60])
def audit(path):
    require(path.stat().st_size==756910 and hashlib.sha256(path.read_bytes()).hexdigest()==PIN,"Pinned workbook changed")
    w=load_workbook(path,read_only=True,data_only=True,keep_links=False)
    try:
        tab=w["Suppl.Table1"]
        rows={i:list(row)[:30] for i,row in enumerate(tab.iter_rows(min_row=1,max_row=86,max_col=30,values_only=True),1)}
    finally:w.close()
    headers={str(i):[str(x) if x is not None else "" for x in rows[i]] for i in range(1,5)}
    out={"status":"HOLD_HEADER_OR_COHORT_SCHEMA","workbook_sha256":PIN,"sheet":"Suppl.Table1","header_first_four_rows":headers,"ancestries":{},"limits":["Only original primary discovery block through Total ALL discovery; no PGS or replication sections.","Control totals are shared published discovery controls, not necessarily subtype-specific controls or unique persons. Do not derive a missing per-SNP N or infer non-overlap.","Match against eight previously verified summary-statistic nominal totals; even exact reconciliation is not certification of log-odds scale, phenotype ascertainment or comparable regression models."]}
    try:
        require(norm(rows[86][0])=="TOTAL ALL DISCOVERY","Discovery cutoff changed")
        # Actual blinded header readback run 45539396437: duplicated incident-only
        # and incident+prevalent columns. Lock the second primary combined block.
        require(norm(rows[3][9])=="INCIDENT & PREVALENT COMBINED" and
                [norm(v) for v in rows[4][9:16]]==["AS","AIS","CES","SVS","LAS","ICH","CONTROLS"],"Exact combined-discovery header changed")
        cols={"AIS":10,"CES":11,"SVS":12,"LAS":13}
        out["trait_columns"]={k:get_column_letter(j+1) for k,j in cols.items()}
        for anc in ("EUR","EAS"):
            starts=[i for i in range(5,86) if norm(rows[i][0])==anc]
            totals=[i for i in range(5,86) if norm(rows[i][0])=="TOTAL "+anc]
            require(len(starts)==len(totals)==1 and starts[0]<totals[0],"Ancestry block changed")
            a,b=starts[0],totals[0]
            controls=n(rows[b][15]);require(controls=={"EUR":1234808,"EAS":237242}[anc],"Control total drift")
            cohort=[]
            for i in range(a,b):
                name=str(rows[i][1] or "").strip()
                if not name:continue
                record={"row":i,"study":name,"incident_only":anc=="EAS" and i==60,"published_controls":n(rows[i][15]),"cases":{k:n(rows[i][j]) for k,j in cols.items()}}
                cohort.append(record)
            usable=[v for v in cohort if not v["incident_only"]]
            x={"total_row":b,"cohort_count":len(usable),"published_controls":controls,"control_row_sum_reconciles":sum(v["published_controls"] or 0 for v in usable)==controls,"traits":{}}
            for k,j in cols.items():
                published=n(rows[b][j]);require(published is not None,"Missing published "+k+" total")
                contributing=[{"study":v["study"],"row":v["row"],"cases":v["cases"][k]} for v in usable if (v["cases"][k] or 0)>0]
                calc=sum(v["cases"] for v in contributing)
                bbj=sum(v["cases"] for v in contributing if v["study"]=="BIOBANK_JAPAN*")
                x["traits"][k]={"published_cases":published,"sum_cases":calc,"sum_reconciles":calc==published,"nominal_gwas_total":TOTAL_N[anc][k],"published_cases_plus_shared_controls_matches_nominal":published+controls==TOTAL_N[anc][k],"case_contributing_cohorts":len(contributing),"case_contributors":contributing,"BBJ_cases":bbj,"BBJ_case_fraction":bbj/published if published else None}
            out["ancestries"][anc]=x
        require(all(x["control_row_sum_reconciles"] and all(y["sum_reconciles"] and y["published_cases_plus_shared_controls_matches_nominal"] for y in x["traits"].values()) for x in out["ancestries"].values()),"One or more original aggregate totals unresolved")
        out["status"]="PUBLISHED_DISCOVERY_TOTALS_RECONCILED_NOT_MODEL_N"
    except (ValueError,KeyError,IndexError,TypeError) as e:out["reason"]=str(e)
    return out
def main():
    p=argparse.ArgumentParser();p.add_argument("--workbook",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    x=audit(a.workbook);a.out.parent.mkdir(parents=True,exist_ok=True);a.out.write_text(json.dumps(x,indent=2,allow_nan=False)+"\n")
    print(json.dumps(x,indent=2,allow_nan=False))
if __name__=="__main__":main()
