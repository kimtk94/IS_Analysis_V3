#!/usr/bin/env python3
"""Cross-phenotype G0022 summary of four CS variants: descriptive, not replication."""
import csv,gzip,json,math,argparse
from pathlib import Path
SRC=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/IS_GWS_REFERENCE_MATCHED_VARIANTS.tsv.gz")
OUT=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1")
CS={"12:112241766:G:A","12:112168009:G:A","12:112468206:C:T","12:112736118:A:G"}
STUDIES={"BBJ":"Japanese_IS","GCST90104544":"EAS_AS","GCST90104545":"EAS_AIS","GCST90104546":"EAS_CES","GCST90104547":"EAS_LAS","GCST90104548":"EAS_SVS"}
def load(path):
    source={};tot={}
    with gzip.open(path,"rt") as f:
        for r in csv.DictReader(f,delimiter="\t"):
            if r["group_id"]!="IS_XDATA_G0022":continue
            ds=r["dataset"];tot[ds]=tot.get(ds,0)+1
            if ds not in STUDIES or r["variant_id"] not in CS:continue
            key=(ds,r["variant_id"])
            if key in source:raise ValueError("Duplicate study-SNP")
            source[key]=r
    if len(source)!=24 or set(tot)!=set(STUDIES):raise ValueError("Incomplete data")
    return source,tot
def run(src,out):
    source,n=load(src);data=[]
    for ds,trait in STUDIES.items():
        for var in sorted(CS):
            r=source[(ds,var)]
            chrom,pos,ref,alt=var.split(":")
            if [r[k] for k in ("chr","pos","ref","alt") ]!=[chrom,pos,ref,alt]:raise ValueError("Wrong allele locus")
            if r["build"]!="GRCh37" or r["ancestry"]!=("Japanese" if ds=="BBJ" else "EAS"):raise ValueError("Ancestry/build mismatch")
            if r["qc_status"] not in ("MATCH_ALT_EFFECT","MATCH_REF_EFFECT"):raise ValueError("Unharmonized")
            beta=float(r["beta"]);b=float(r["alt_effect_beta"]);se=float(r["se"]);p=float(r["p"]);af=float(r["alt_effect_eaf"])
            sign=1 if r["effect_allele"]==alt else -1 if r["effect_allele"]==ref else None
            source_eaf=float(r["eaf"])
            expected_alt_eaf=source_eaf if sign==1 else 1-source_eaf
            if (sign is None or abs(sign*beta-b)>1e-9 or se<=0 or
                not(0<af<1 and 0<source_eaf<1 and 0<=p<=1) or
                abs(expected_alt_eaf-af)>1e-8):
                raise ValueError("Allele beta/EAF contract")
            data.append({"study":ds,"trait":trait,"variant":var,"ref":ref,"effect_allele":"ALT",
                "original_effect_allele":r["effect_allele"],"ALT_beta":b,"se":se,
                "ALT_OR":math.exp(b),"OR_95CI_lower":math.exp(b-1.96*se),
                "OR_95CI_upper":math.exp(b+1.96*se),"ALT_EAF":af,"p":p,
                "GWS":int(p<5e-8),"nominal_significant":int(p<.05),
                "effect_direction":"protective" if b<0 else "risk",
                "scientific_gate":"DESCRIPTIVE_ONLY_OVERLAPPING_STUDIES_NOT_INDEPENDENT_REPLICATION"})
    out.mkdir(parents=True,exist_ok=True)
    for filename,rows in [("G0022_FOUR_CS_SIX_GWAS_EFFECTS.tsv",data),
         ("G0022_RS671_SIX_GWAS_SUBTYPE_EFFECTS.tsv",[x for x in data if x["variant"]=="12:112241766:G:A"])]:
        with (out/filename).open("w",newline="") as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t");w.writeheader();w.writerows(rows)
    rs=[x for x in data if x["variant"]=="12:112241766:G:A"]
    summary={"rows":len(data),"rs671_phenotypes":len(rs),"all_24_study_variant_effects_protective":all(x["ALT_beta"]<0 for x in data),
      "rs671_AIS_OR":next(x["ALT_OR"] for x in rs if x["trait"]=="EAS_AIS"),
      "rs671_SVS_p":next(x["p"] for x in rs if x["trait"]=="EAS_SVS"),
      "rs671_CES_p":next(x["p"] for x in rs if x["trait"]=="EAS_CES"),
      "rs671_LAS_p":next(x["p"] for x in rs if x["trait"]=="EAS_LAS"),
      "source_record_counts":n,"phenotype_independence_verified":False,
      "sex_stratified_summary_availability":"NOT_AVAILABLE",
      "formal_effect_difference_or_meta_analysis":"NOT_VALID_WITHOUT_COHORT_COVARIANCE",
      "causal_stroke_gene":False}
    (out/"G0022_RS671_STROKE_SUBTYPE_SUMMARY.json").write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2))
    for r in rs:print(r["trait"],r["ALT_beta"],r["se"],r["ALT_OR"],r["p"])
    return summary
if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--source",type=Path,default=SRC);p.add_argument("--out",type=Path,default=OUT);a=p.parse_args()
    run(a.source,a.out)
