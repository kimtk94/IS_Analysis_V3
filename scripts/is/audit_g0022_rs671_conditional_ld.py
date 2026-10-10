#!/usr/bin/env python3
"""Diagnostic single-variant conditioning on ALDH2 rs671 for G0022 AIS.

Conditional z=(z_j-r_jk*z_k)/sqrt(1-r_jk**2). Exploratory, external
1KG EAS n_ref=504, potential LD mismatch, no cohort-level validation.
Not a substitute for GCTA-COJO, conditional GWAS, or multi-causal SuSiE.
"""
import argparse,csv,json,math
from pathlib import Path
import numpy as np
BASE=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1")
RS671="12:112241766:G:A"
CS={
 "12:112241766:G:A":"ALDH2 rs671",
 "12:112168009:G:A":"ACAD10 rs11066015",
 "12:112468206:C:T":"NAA25-overlap",
 "12:112736118:A:G":"HECTD4-overlap"}
def p_from_z(z):return math.erfc(abs(z)/math.sqrt(2))
def residualized_z(z,ldcol,lead):
    z=np.asarray(z,dtype=float)
    ldcol=np.asarray(ldcol,dtype=float)
    if lead<0 or lead>=len(z) or len(ldcol)!=len(z):
        raise ValueError("Bad conditional dimensions")
    if max(abs(ldcol))>1.0000001 or not np.all(np.isfinite(ldcol)):
        raise ValueError("Invalid signed LD")
    if abs(ldcol[lead]-1)>1e-6:
        raise ValueError("Lead must be exactly the reference LD diagonal")
    residual=(z-ldcol*z[lead])/np.sqrt(np.maximum(1e-12,1-ldcol**2))
    # A duplicated reference genotype variant has |r|=1; skip for defensibility.
    residual[np.abs(ldcol)>=1-1e-8]=np.nan
    return residual
def build(root):
    source=json.loads((root/"FULL_LOCUS_INPUT_QC.json").read_text())
    model=json.loads((root/"G0022_FULL_AIS_SUSIE_SUMMARY.json").read_text())
    if source["status"]!="FULL_G0022_SIGNED_LD_INPUT_QC_PASS" or source["genotype_n"]!=504:
        raise ValueError("Original verified LD source required")
    with (root/"variants.tsv").open() as f:
        records=list(csv.DictReader(f,delimiter="\t"))
    n=len(records)
    if n!=source["genotype_snp_count"] or model["n_snps"]!=n or n>6500:
        raise ValueError("Locus dimension/source mismatch")
    positions={r["variant_id"]:i for i,r in enumerate(records)}
    if not set(CS).issubset(positions):
        raise ValueError("Original 4 CS variants missing")
    z=np.array([float(x["z"]) for x in records],dtype=np.float64)
    if not np.isfinite(z).all():raise ValueError("Nonfinite z")
    ld=np.memmap(root/"ld.rowmajor.f64",dtype="<f8",mode="r",shape=(n,n))
    lead=positions[RS671]
    signed_ld=ld[:,lead]
    cond=residualized_z(z,signed_ld,lead)
    counts=sum(bool(abs(v)>=5.45131) for v in cond if math.isfinite(v))
    clumped=[];report=[]
    for var,annotation in CS.items():
        i=positions[var];lr=float(signed_ld[i])
        r2=lr*lr
        clumped.append({
            "variant_id":var,"annotation":annotation,
            "study_z":float(z[i]),"unconditioned_p":p_from_z(z[i]),
            "signed_r_to_rs671":lr,"r2_to_rs671":r2,
            "conditional_z_given_rs671":float(cond[i]) if math.isfinite(cond[i]) else "",
            "conditional_p_approx":p_from_z(cond[i]) if math.isfinite(cond[i]) else "",
            "interpretation":"EAS_1000G_LD_PROXIMITY_ONLY_NOT_CAUSAL"})
    for i,r in enumerate(records):
        if math.isfinite(cond[i]):
            report.append({"variant_id":r["variant_id"],"pos":r["pos"],"study_z":z[i],
                "conditional_z":float(cond[i]),"conditional_p_approx":p_from_z(cond[i]),
                "r2_to_rs671":float(signed_ld[i]**2),
                "status":"APPROX_EXTERNAL_LD_SINGLE_LEAD_CONDITIONING"})
    report.sort(key=lambda x:abs(x["conditional_z"]),reverse=True)
    with (root/"G0022_AIS_RS671_CONDITIONAL_ALL_SNPS.tsv").open("w",newline="") as f:
        w=csv.DictWriter(f,delimiter="\t",fieldnames=list(report[0]));w.writeheader();w.writerows(report)
    with (root/"G0022_AIS_FOUR_CS_RS671_LD.tsv").open("w",newline="") as f:
        w=csv.DictWriter(f,delimiter="\t",fieldnames=list(clumped[0]));w.writeheader();w.writerows(clumped)
    details={
        "rs671_study_gwas_z":float(z[lead]),"lead_variant":RS671,
        "genotyped_reference":source["LD_reference"],"reference_n":source["genotype_n"],
        "gwas_snp_count":n,
        "conditional_snp_count":len(report),
        "conditional_gws_p_lt_5e8":sum(x["conditional_p_approx"]<5e-8 for x in report),
        "maximum_absolute_conditional_z":max(abs(x["conditional_z"]) for x in report),
        "minimum_conditional_p":min(x["conditional_p_approx"] for x in report),
        "top_conditional_variant":report[0]["variant_id"],
        "credible_set_variant_R2":{x["variant_id"]:round(x["r2_to_rs671"],6) for x in clumped},
        "status":"ONE_LEAD_CONDITIONAL_PROXY_ANALYSIS_NOT_VALIDATED",
        "caveat":"Formula assumes single causal lead/known population LD; 1KG EAS n=504 underpowered for high-LD residual analysis; not independent conditional GWAS."}
    (root/"G0022_AIS_RS671_CONDITIONAL_DIAGNOSTIC_SUMMARY.json").write_text(json.dumps(details,indent=2))
    print(json.dumps(details,indent=2))
    return details
if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--root",type=Path,default=BASE)
    build(p.parse_args().root)
