#!/usr/bin/env python3
"""Audit published GIGASTROKE MR/TWAS baselines and estimate conservative subtype contrasts.

This uses no private input, no new GWAS fitting, and makes no independence claim.
"""
import argparse
import hashlib
import itertools
import json
import math
import re
import urllib.request
from pathlib import Path

SOURCE_URL = ("https://static-content.springer.com/esm/"
              "art%3A10.1038%2Fs41586-022-05492-5/MediaObjects/"
              "41586_2022_5492_MOESM2_ESM.xlsx")
SOURCE_SHA256 = "7c688f864c73ec2e9624bde9f53ac48a9880a6a08e59f5e1266690e243618ba4"
SOURCE_SIZE = 756910
SUBTYPES = ("LAS", "CES", "SVS")
# Published Figure 2, Neurology 2020, DOI 10.1212/WNL.0000000000009814.
# SBP per 10 mmHg, rounded marginal odds ratio and 95% CI, not original per-SNP effects.
SBP_FIG2 = {"LAS": (1.68, 1.54, 1.84), "CES": (1.24, 1.16, 1.34),
            "SVS": (1.47, 1.36, 1.58)}
FAMILY_N = 6  # two exposures x three planned pairwise contrasts, including untested LDL-C

def normal_p(z):
    return math.erfc(abs(z)/math.sqrt(2))

def contrast(beta1, se1, beta2, se2, family=FAMILY_N):
    if not all(math.isfinite(x) for x in (beta1,se1,beta2,se2)) or min(se1,se2)<=0 or family<1:
        raise ValueError("Invalid effect or uncertainty")
    d=beta1-beta2
    # For an unknown correlation rho in [-1,1], the maximal contrast SE is se1+se2.
    # Thus p at rho=-1 bounds the marginal Wald p from above; no assumed shared-control rho.
    max_se=se1+se2
    pmax=normal_p(d/max_se)
    grid=[]
    for rho in (-1,-.5,0,.25,.5,.75,1):
        v=max(0,se1**2+se2**2-2*rho*se1*se2)
        p=None if v==0 else normal_p(d/math.sqrt(v))
        grid.append({"rho":rho,"se_difference":math.sqrt(v),"p":p})
    return {"difference_log_odds":d,"ratio_of_odds_ratios":math.exp(d),
            "max_possible_se":max_se,"worst_case_p":pmax,
            "worst_case_bonferroni":min(1,family*pmax),
            "passes_all_covariances":pmax<.05/family,"rho_sensitivity":grid}

def sbp_baseline():
    values={k:(math.log(v[0]),(math.log(v[2])-math.log(v[1]))/(2*1.95996398454))
            for k,v in SBP_FIG2.items()}
    rows=[]
    for a,b in itertools.combinations(SUBTYPES,2):
        rows.append({"contrast":a+"-"+b,**contrast(*values[a],*values[b])})
    return {"source_doi":"10.1212/WNL.0000000000009814","figure":"2",
            "exposure":"SBP per 10 mmHg","family_n":FAMILY_N,
            "published_rounded_OR_CI":SBP_FIG2,"contrasts":rows,
            "interpretation":"Exploratory secondary analysis of rounded published MR estimates. Marginal normal approximation and SE inferred from CI. Worst-case covariance bound is conservative if the marginal estimators are valid and comparably scaled. No covariance estimate, SNP-level refit, or new independence claim. Three LDL-C contrasts remain untested in this six-comparison family."}

def download():
    req=urllib.request.Request(SOURCE_URL,headers={"User-Agent":"publication-contrast/1.0","Accept-Encoding":"identity"})
    with urllib.request.urlopen(req,timeout=60) as resp:
        if resp.status!=200:raise ValueError("Provider did not return complete file")
        data=resp.read(SOURCE_SIZE+1)
    if len(data)!=SOURCE_SIZE or hashlib.sha256(data).hexdigest()!=SOURCE_SHA256:
        raise ValueError("Corrected publisher workbook byte/hash mismatch")
    return data

def cell_text(value):
    return "" if value is None else str(value).strip()

def sheet_extract(wb,name,pattern,header_limit):
    if name not in wb.sheetnames:
        raise ValueError("Expected published sheet not found: "+name)
    ws=wb[name]
    headers=[]
    matched=[]
    for i,r in enumerate(ws.iter_rows(values_only=True),1):
        if i>100000:raise ValueError("Sheet row bound")
        cells=[{"column":j+1,"value":cell_text(v)} for j,v in enumerate(r) if v is not None]
        if i<=header_limit:headers.append({"row":i,"cells":cells})
        if any(pattern.search(x["value"]) for x in cells):
            if len(matched)>=400:raise ValueError("Unexpected number of matched published rows")
            matched.append({"row":i,"cells":cells})
    return {"sheet":name,"row_count":i,"headers":headers,"matched_rows":matched}

def main():
    arg=argparse.ArgumentParser(description=__doc__)
    arg.add_argument("--out",type=Path,required=True)
    args=arg.parse_args()
    args.out.mkdir(parents=True,exist_ok=True)
    base=sbp_baseline()
    (args.out/"published_sbp_contrasts.json").write_text(json.dumps(base,indent=2,allow_nan=False)+"\n")
    print("PUBLISHED_SBP",json.dumps(base,allow_nan=False),flush=True)
    try:
        import io
        from openpyxl import load_workbook
        data=download()
        wb=load_workbook(io.BytesIO(data),read_only=True,data_only=True)
        mr_pattern=re.compile(r"\b(?:LDL(?:[-_ ]?C|[-_ ]?cholesterol)?|SBP|systolic(?: blood pressure)?)\b",re.I)
        gene_pattern=re.compile(r"\b(?:FGF5|CALHM2|ALDH2|NEURL1|INA|C4orf22)\b",re.I)
        tables={"mr":sheet_extract(wb,"Suppl.Table24",mr_pattern,12),
                "twas":sheet_extract(wb,"Suppl.Table27",gene_pattern,8)}
        wb.close()
        obj={"status":"PUBLISHER_TABLE_ROWS_READ","source":SOURCE_URL,
             "source_sha256":SOURCE_SHA256,"source_size":SOURCE_SIZE,**tables}
        # Read exact published cells, not a guessed coefficient schema. A subsequent
        # direct contrast requires validated exposure unit, beta, SE, ancestry and method.
        (args.out/"gigastroke_published_extract.json").write_text(json.dumps(obj,indent=2)+"\n")
        for k,v in tables.items():
            print("PUBLISHED_"+k.upper()+"_HEADER",json.dumps(v["headers"]),flush=True)
            print("PUBLISHED_"+k.upper()+"_MATCHES",json.dumps(v["matched_rows"][:35]),flush=True)
        print("TABLE_COUNTS",json.dumps({k:{"rows":v["row_count"],"matches":len(v["matched_rows"])} for k,v in tables.items()}))
    except (OSError,ValueError,KeyError,ImportError) as e:
        obj={"status":"HOLD_PUBLISHED_SOURCE_OR_SCHEMA","error":type(e).__name__+": "+str(e)[:200]}
        (args.out/"gigastroke_published_extract.json").write_text(json.dumps(obj,indent=2)+"\n")
        print("PUBLISHED_SOURCE_HOLD",json.dumps(obj),flush=True)

if __name__=="__main__":
    main()
