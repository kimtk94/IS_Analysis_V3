#!/usr/bin/env python3
"""Integrate rs9515201 published-pair recovery into the original, fixed 40 tests.
Only original held rs9515201 rows may be filled. Other rows and tests stay fixed.
"""
import argparse, json, math, re
from pathlib import Path
from fixed_subtype_indexed import PANEL, STUDIES, contrast, require, write
RS="rs9515201"
TRAITS={"AIS":"Ischemic stroke","CES":"Cardioembolic stroke","LAS":"Large artery stroke","SVS":"Small vessel stroke"}
def metadata(x,study):
    m=x["metadata"];a,t=study.split("_")
    require(m["genome_assembly"]=="GRCh38" and m["is_harmonised"] is True,"Assembly/source")
    require(m["trait_description"]==[TRAITS[t]],"Phenotype identity")
    require({v for s in m["samples"] for v in s["sample_ancestry_category"]}=={"European" if a=="EUR" else "East Asian"},"Ancestry identity")
def validate_pair(x,study,original):
    require(x["study"]==study and x["accession"]==STUDIES[study],"Recovery study identity")
    row=x["row"]
    if row["status"]!="MATCHED_HARMONISED_ROW":return {"status":"HOLD_RECOVERY","reason":row["status"]}
    metadata(x,study)
    r=x["reference"]
    require(r["prespecified_pair"]=="A/C" and r["chr"]=="13" and r["pos"]==110388451 and r["genome_reference"]=="A","Prior pair identity")
    require(row["effect_allele"]=="A" and row["other_allele"]=="C","A versus C")
    for key in ("source","metadata_sha256","index_sha256","header_range_sha256"):
        require(x[key]==original[key],"Original source identity changed: "+key)
    b,s,p,eaf=[float(row[k]) for k in ("beta_A","se","p","A_eaf")]
    require(all(math.isfinite(v) for v in (b,s,p,eaf)) and s>0 and 0<=p<=1 and 0<eaf<1,"Invalid recovered numbers")
    require(row["hm_code"] in (10,11,12,13),"Unresolved harmonisation")
    return {"status":"MATCHED_HARMONISED_ROW","effect_allele":"A","other_allele":"C","beta_alt":b,"se":s,"p":p,"alt_eaf":eaf,"orientation":"A relative to C; A is genome reference, not ALT","recovery_source":x}
def integrate(x,recovery):
    require(x["family"]==40 and set(x["studies"])==set(STUDIES),"Original family/study set")
    d=x["studies"]
    for study,v in d.items():
        require(v["accession"]==STUDIES[study] and set(v["rows"])==set(PANEL),"Original study identity")
        metadata(v,study)
        require(v["rows"][RS]["status"]=="HOLD_REFERENCE","Refusing overwrite of prior valid row")
        v["rows"][RS]=validate_pair(recovery[study],study,v) if study in recovery else {"status":"HOLD_RECOVERY_PENDING"}
        v["matched"]=sum(r["status"]=="MATCHED_HARMONISED_ROW" for r in v["rows"].values())
    pairs=[("EAS_"+t,"EUR_"+t) for t in ("AIS","CES","LAS","SVS")]+[(a+"_SVS",a+"_"+t) for a in ("EUR","EAS") for t in ("CES","LAS")]
    require(len(x["contrasts"])==40 and {(r["rsid"],r["first"],r["second"]) for r in x["contrasts"]}=={(rs,a,b) for rs in PANEL for a,b in pairs},"Original test family changed")
    rows=[]
    for r in x["contrasts"]:
        if r["rsid"]!=RS:
            rows.append(r);continue
        a,b=r["first"],r["second"];u=d[a]["rows"][RS];v=d[b]["rows"][RS]
        z={"rsid":RS,"first":a,"second":b,"status":"HOLD_UNMATCHED"}
        if u["status"]==v["status"]=="MATCHED_HARMONISED_ROW":z.update(contrast(u,v))
        rows.append(z)
    valid=[r for r in rows if r["status"]=="CONDITIONAL_REPORTED_BETA_CONTRAST"]
    return {"status":"RECOVERY_INTEGRATED_EXPLORATORY_ONLY","family":40,"matched_cells":sum(v["matched"] for v in d.values()),"matched_contrasts":len(valid),"recovered_rs9515201_cells":sum(v["rows"][RS]["status"]=="MATCHED_HARMONISED_ROW" for v in d.values()),"significant_bonferroni_zero_covariance":[r for r in valid if r["bonferroni_40_zero"]<.05],"significant_bonferroni_max_variance":[r for r in valid if r["bonferroni_40_max"]<.05],"studies":d,"contrasts":rows,"limits":["Recovery only for independently published A/C pair; other two multi-allelic SNPs remain held. Original 40-test multiplicity denominator retained.","Reported beta comparisons conditional on scale/model compatibility; unknown covariance bound is not adjustment for ancestry LD, winner's curse or confounding.","No new MR/colocalization or independent cohort validation. Indexed ranges do not verify full publisher file digests."]}
def main():
    p=argparse.ArgumentParser();p.add_argument("--original",type=Path,required=True);p.add_argument("--recovery",type=Path,required=True);p.add_argument("--out",type=Path,required=True);a=p.parse_args()
    x=json.loads(a.original.read_text());r={}
    for k in STUDIES:
        path=a.recovery/k/(k+".json")
        if path.is_file():r[k]=json.loads(path.read_text())
    require(len(r)==len(STUDIES),"All eight recovery artifacts required; do not label partial results final")
    write(a.out,integrate(x,r))
if __name__=="__main__":main()
