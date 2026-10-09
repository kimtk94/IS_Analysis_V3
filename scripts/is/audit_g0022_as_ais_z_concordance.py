#!/usr/bin/env python3
"""EAS GIGASTROKE AS-vs-AIS G0022 Z/EAF concordance audit.

These are overlapping cohorts/phenotypes; not replication and not a flip caller.
"""
import csv,gzip,json,argparse
from pathlib import Path
import numpy as np

DEFAULT=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2")
CENTERS={
 "AS_clump1":111629389, "AS_clump2":112930475,
 "AIS_clump1":112241766,"AIS_clump2":113031474,
 "AIS_clump3":110675363
}
AS="GCST90104544"
AIS="GCST90104545"
def compare(root,radius=250000):
    by={}
    with gzip.open(root/"IS_GWS_REFERENCE_MATCHED_VARIANTS.tsv.gz","rt") as f:
        for row in csv.DictReader(f,delimiter="\t"):
            if row["group_id"]!="IS_XDATA_G0022" or row["dataset"] not in (AS,AIS):
                continue
            if row["qc_status"] not in ("MATCH_ALT_EFFECT","MATCH_REF_EFFECT"):
                continue
            try:
                z=float(row["alt_effect_beta"])/float(row["se"])
                eaf=float(row["alt_effect_eaf"])
            except (ValueError,ZeroDivisionError):continue
            if not (np.isfinite(z) and np.isfinite(eaf)):continue
            k=(row["dataset"],row["variant_id"])
            if k not in by or abs(z)>abs(by[k][0]):
                by[k]=(z,eaf,int(row["pos"]),row["p"])
    overlap=[]
    for (study,vid),item in by.items():
        if study!=AS or (AIS,vid) not in by:continue
        alt=by[AIS,vid]
        overlap.append((vid,item[2],item[0],alt[0],item[1],alt[1],item[3],alt[3]))
    summaries=[];outliers=[]
    for name,center in CENTERS.items():
        arr=[x for x in overlap if abs(x[1]-center)<=radius]
        if len(arr)<10:raise ValueError(f"Insufficient matched AS/AIS variants: {name}")
        z1=np.asarray([x[2] for x in arr],dtype=float)
        z2=np.asarray([x[3] for x in arr],dtype=float)
        ef1=np.asarray([x[4] for x in arr],dtype=float)
        ef2=np.asarray([x[5] for x in arr],dtype=float)
        strong=(abs(z1)>=3) & (abs(z2)>=3)
        sig=(abs(z1)>=3)|(abs(z2)>=3)
        c=float(np.corrcoef(z1,z2)[0,1])
        slope=float(z1@z2/(z2@z2)) if z2@z2>0 else None
        big=sorted(arr,key=lambda x:abs(x[2]-x[3]),reverse=True)
        for rank,x in enumerate(big[:20],1):
            outliers.append({
              "window":name,"variant_id":x[0],"pos":x[1],"as_z":round(x[2],7),
              "ais_z":round(x[3],7),"abs_z_difference":round(abs(x[2]-x[3]),7),
              "as_alt_eaf":round(x[4],7),"ais_alt_eaf":round(x[5],7),
              "abs_eaf_diff":round(abs(x[4]-x[5]),7),
              "as_p":x[6],"ais_p":x[7],"rank_by_abs_z_diff":rank})
        summaries.append({
          "window":name,"center":center,"radius_bp":radius,"n_variant_pairs":len(arr),
          "as_ais_z_pearson":round(c,6),
          "as_z_regression_slope_on_ais_z":round(slope,6) if slope is not None else "",
          "strong_both_abs_z_ge_3":int(strong.sum()),
          "discordant_sign_strong_both":int(((z1[strong]*z2[strong])<0).sum()),
          "discordant_sign_either_strong":int(((z1[sig]*z2[sig])<0).sum()),
          "mean_abs_z_delta":round(float(np.abs(z1-z2).mean()),6),
          "max_abs_z_delta":round(float(np.abs(z1-z2).max()),6),
          "median_abs_eaf_delta":round(float(np.median(np.abs(ef1-ef2))),6),
          "max_abs_eaf_delta":round(float(np.abs(ef1-ef2).max()),6),
          "interpretation":"CROSS_PHENOTYPE_DESCRIPTIVE_NOT_INDEPENDENT_REPLICATION"})
    for name,rows in [("G0022_AS_VS_AIS_Z_CONCORDANCE.tsv",summaries),
                      ("G0022_AS_VS_AIS_LARGEST_Z_DIFFS.tsv",outliers)]:
        with (root/name).open("w",newline="") as h:
            w=csv.DictWriter(h,delimiter="\t",fieldnames=list(rows[0]))
            w.writeheader();w.writerows(rows)
    summary={"windows":len(summaries),"overlap_total_unique_variants":len(overlap),
         "any_strong_both_sign_discordance":sum(x["discordant_sign_strong_both"] for x in summaries),
         "min_window_pearson_r":min(x["as_ais_z_pearson"] for x in summaries),
         "warning":"AS and AIS share samples and phenotypes differ. Low correlation not proof of allele flip or error."}
    (root/"G0022_AS_VS_AIS_Z_CONCORDANCE_SUMMARY.json").write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2))
    for row in summaries:print(row)
    return summary
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--root",type=Path,default=DEFAULT)
    p.add_argument("--radius-bp",type=int,default=250000)
    args=p.parse_args()
    compare(args.root,args.radius_bp)
