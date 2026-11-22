#!/usr/bin/env python3
"""Optional actual-source replay of legacy BBJ P01-P15 ABF sufficient statistics.

Run locally on already-authorized Google Drive source files; never upload SNP rows.
GitHub CI on its own DOES NOT have these files, and source-level audit is NOT a
formal coloc/QTL-N/allele approval. All formulas mirror coloc R/claudia.R
(sdY.est, approx.bf.estimates and combine.abf), independently implemented.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import importlib
import json
import math
from pathlib import Path

MODEL=importlib.import_module("scripts.is.replay_legacy_15_abf_sufficient")
COLUMNS=["snp","gwas_beta","gwas_se","gwas_maf","eqtl_beta","eqtl_se","eqtl_maf","qtl_n"]
FILES=[
 "FGF5_stats.tsv","P02_C4orf22_Cerebellar_Hemisphere_stats.tsv",
 "P03_C4orf22_Cerebellum_stats.tsv","P04_INA_stats.tsv",
 "P05_NEURL1_stats.tsv","P06_CALHM2_stats.tsv",
 "P07_RPL6_stats.tsv","P08_TRAFD1_stats.tsv",
 "COL4A2_stats.tsv","P10_COL4A2_Amygdala_stats.tsv",
 "ALDH2_stats.tsv","P12_ALDH2_Aorta_stats.tsv",
 "P13_ALDH2_Cortex_stats.tsv","SH3PXD2A_stats.tsv",
 "P15_COL4A1_stats.tsv"
]

def logbf(beta,se,prior_sd):
    variance=se*se
    W=prior_sd*prior_sd
    r=W/(W+variance)
    return 0.5*(math.log1p(-r)+r*(beta/se)**2)

def one_source(path, expected):
    _,gene,_,expected_rows,expected_bytes,expected_sha,original_n,expected_pal,expected_gap,_=expected
    MODEL.check(path.is_file() and not path.is_symlink() and path.stat().st_size<=1_000_000,
                "source input missing/not regular/symlink/oversized: "+gene)
    raw=path.read_bytes()
    MODEL.check(len(raw)==expected_bytes and hashlib.sha256(raw).hexdigest()==expected_sha,
                "original Drive source SHA256 mismatch: "+gene)
    parsed=list(csv.DictReader(raw.decode("utf-8-sig").splitlines(),delimiter="\t"))
    MODEL.check(len(parsed)==expected_rows and len(parsed)>0,
                "source cohort SNP input row count mismatch: "+gene)
    seen=set();effects=[];pal=0;gap=0
    for rec in parsed:
        MODEL.check(list(rec)==COLUMNS and None not in rec and None not in rec.values(),
                    "source fields/row formatting invalid: "+gene)
        snp=rec["snp"]
        tokens=snp.split(":")
        MODEL.check(len(tokens)==4 and tokens[1].isdigit()
                    and tokens[2] in "ACGT" and tokens[3] in "ACGT"
                    and len(tokens[2])==1 and len(tokens[3])==1 and tokens[2]!=tokens[3],
                    "SNP position or allele invalid: "+gene)
        MODEL.check(snp not in seen, "duplicate original SNP: "+gene)
        seen.add(snp)
        try:
            numbers=[float(rec[c]) for c in COLUMNS[1:]]
        except ValueError as err:
            raise ValueError("invalid numeric input: "+gene) from err
        MODEL.check(all(math.isfinite(z) for z in numbers),"nonfinite source numeric")
        bg,sg,mg,bq,sq,mq,N=numbers
        MODEL.check(sg>0 and sq>0 and 0<mg<=0.5 and 0<mq<=0.5
                    and N==original_n, "invalid source SE/MAF/historical scalar N")
        if frozenset([tokens[2],tokens[3]]) in (frozenset(["A","T"]),frozenset(["C","G"])):
            pal+=1
        if abs(mg-mq)>0.2:
            gap+=1
        effects.append((bg,sg,bq,sq,mq))
    MODEL.check(pal==expected_pal and gap==expected_gap,
                "source allele/QC aggregate drift: "+gene)

    oneover=[1/(e[3]**2) for e in effects]
    nscaled=[2*original_n*e[4]*(1-e[4]) for e in effects]
    slope=math.fsum(o*v for o,v in zip(oneover,nscaled))/math.fsum(o*o for o in oneover)
    MODEL.check(slope>0,"estimated QTL phenotype variance invalid")
    gwas=[logbf(bg,sg,0.2) for bg,sg,_,_,_ in effects]
    suff=[]
    for factor in MODEL.N_FACTORS:
        sdY=math.sqrt(factor*slope)
        eqtl=[logbf(bq,sq,0.15*sdY) for _,_,bq,sq,_ in effects]
        suff.append([MODEL.logsumexp(gwas),MODEL.logsumexp(eqtl),
                     MODEL.logsumexp([a+b for a,b in zip(gwas,eqtl)])])
    original=expected[-1]
    worst=max(abs(a-b) for t,u in zip(suff,original) for a,b in zip(t,u))
    MODEL.check(worst<=2e-6,
                "actual original source re-read differs from GitHub source-sufficient log sums: "+gene)
    return {"gene":gene,"panel":expected[0],"sha256_verified":True,
            "source_SNP_count":len(effects),"max_abs_logsum_difference":worst,
            "scientific_status":"EXPLORATORY_HOLD"}

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--source-root",type=Path,required=True)
    ap.add_argument("--sufficient",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    args=ap.parse_args()
    try:
        spec=json.loads(args.sufficient.read_text(encoding="utf-8"))
        MODEL.reconcile(spec,0,1)
        root=args.source_root.resolve(strict=True)
        rows=[one_source(root / file,panel) for file,panel in zip(FILES,spec["panels"])]
        status={"schema":"is-15-archive-direct-source-reread-v1",
                "technical_status":"PASS_ORIGINAL_SOURCE_15_FILES_SHA_AND_LOGSUM",
                "groups":15,"original_source_SNP_rows":sum(x["source_SNP_count"] for x in rows),
                "original_source_files_read_this_execution":True,
                "scientific_status":"EXPLORATORY_HOLD",
                "QTL_N_source_authority_independently_verified":False,
                "allele_and_donor_LD_approved":False,
                "source_input_checks":rows}
        MODEL.check(status["original_source_SNP_rows"]==29661,"full source universe drift")
    except (OSError,KeyError,TypeError,ValueError) as error:
        status={"schema":"is-15-archive-direct-source-reread-v1",
                "technical_status":"HOLD_MISSING_OR_MISMATCHED_SOURCE",
                "scientific_status":"EXPLORATORY_HOLD","error":str(error)[:240]}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(status,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({k:v for k,v in status.items() if k in
                      ("technical_status","groups","original_source_SNP_rows","error")}))
    return 0 if status["technical_status"].startswith("PASS_") else 1

if __name__=="__main__":
    raise SystemExit(main())
