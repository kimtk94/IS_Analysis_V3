#!/usr/bin/env python3
"""Offline exact-reference SNP lookup in authorized TPMI PheWeb 433.21 GWAS.

Input schema verified from previously saved 2MB download probe:
chrom,pos,ref,alt,rsids,pval,beta,sebeta,af,case_af,control_af.
Requires the full legitimately obtained original GRCh38 gzip file.
Never scrapes, retries or bypasses the TPMI website.
"""
import argparse,csv,gzip,json,math
from pathlib import Path
from collections import defaultdict,Counter

REQUIRED={"chrom","pos","ref","alt","rsids","pval","beta","sebeta","af"}
def number(x):
    try:
        v=float(x)
        return v if math.isfinite(v) else None
    except (TypeError,ValueError):return None

def extract(gwas_path,targets):
    by_pos=defaultdict(list)
    for r in targets:by_pos[(r["GRCh38"].split(":")[0],int(r["GRCh38"].split(":")[1]))].append(r)
    observed=defaultdict(list)
    at_other=defaultdict(list)
    rows_scanned=0
    opener=gzip.open if str(gwas_path).endswith(".gz") else open
    with opener(gwas_path,"rt",newline="") as f:
        reader=csv.DictReader(f,delimiter="\t")
        if not REQUIRED.issubset(reader.fieldnames or []):
            raise ValueError("TPMI PheWeb 433.21 expected GRCh38 source schema missing: "+
                repr(sorted(REQUIRED-set(reader.fieldnames or []))))
        for r in reader:
            rows_scanned+=1
            try:key=(r["chrom"].removeprefix("chr"),int(r["pos"]))
            except (ValueError,TypeError):continue
            if key not in by_pos:continue
            for target in by_pos[key]:
                if f"{key[0]}:{key[1]}:{r['ref'].upper()}:{r['alt'].upper()}"==target["GRCh38"]:
                    observed[target["GRCh38"]].append(r)
                else:at_other[target["GRCh38"]].append(r)
    results=[]
    for target in targets:
        vid=target["GRCh38"]
        matched=observed.get(vid,[])
        other=at_other.get(vid,[])
        if len(matched)>1: raise ValueError("Duplicate exact target: "+vid)
        row={"gene":target["gene"],"rsid":target["rsid"],
             "GRCh37":target["GRCh37"],"GRCh38":vid,
             "status":"NOT_REPORTED","n_same_position_other_alleles":len(other)}
        if len(matched)==1:
            m=matched[0]
            beta,se,p,af=(number(m[x]) for x in ("beta","sebeta","pval","af"))
            if beta is None or se is None or se<=0 or p is None or p<0 or p>1 or af is None or af<0 or af>1:
                raise ValueError("Nonvalid TPMI association row at "+vid)
            reported=set(str(m["rsids"]).replace(";",",").split(","))
            row.update({"status":("EXACT_GRCH38_REF_ALT_MATCH" if target["rsid"] in reported else "EXACT_ALLELES_RSID_REVIEW"),
                        "source_reported_rsids":m["rsids"],
                        "target_rsid_in_reported":target["rsid"] in reported,
                        "tpmi_beta_ALT":beta,"tpmi_se":se,"tpmi_p":p,
                        "tpmi_p_numeric_underflow":p==0,
                        "tpmi_alt_af":af,
                        "tpmi_maf":min(af,1-af),
                        "tpmi_case_af":m.get("case_af",""),
                        "tpmi_control_af":m.get("control_af","")})
        elif len(other)>0:
            row["status"]="SAME_POSITION_DIFFERENT_ALLELE_PAIR_REVIEW"
        results.append(row)
    return results,rows_scanned

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--summary-gwas",type=Path,required=True)
    p.add_argument("--crossbuild-map",type=Path,required=True)
    p.add_argument("--out-dir",type=Path,required=True)
    p.add_argument("--phenotype",choices=["433.21"],required=True)
    p.add_argument("--build",choices=["GRCh38"],required=True)
    p.add_argument("--effect-alt-confirmed",action="store_true",help="Explicitly attest TPMI source beta is for ALT allele")
    a=p.parse_args()
    if not a.effect_alt_confirmed:
        p.error("TPMI source effect orientation must be confirmed (--effect-alt-confirmed) before ALT beta reporting")
    if not a.summary_gwas.is_file():raise FileNotFoundError(a.summary_gwas)
    with a.crossbuild_map.open(newline="") as f:
        targets=list(csv.DictReader(f,delimiter="\t"))
    if len(targets)!=11 or any(x["reference_status"]!="CHAIN_UCSC_GRCH38_REF_ALT_PASS" for x in targets):
        raise ValueError("11 audited GRCh38 reference positions are required")
    if len({x["GRCh38"] for x in targets})!=11:raise ValueError("Non-unique GRCh38 IDs")
    if a.out_dir.exists() and any(a.out_dir.iterdir()):raise FileExistsError(a.out_dir)
    results,n=extract(a.summary_gwas,targets)
    a.out_dir.mkdir(parents=True,exist_ok=True)
    dest=a.out_dir/"TPMI_43321_ALCOHOL_11_SNP_ASSOCIATIONS.tsv"
    fields=list(dict.fromkeys(k for r in results for k in r))
    with dest.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter="\t")
        w.writeheader();w.writerows(results)
    summary={"status":"DESCRIPTIVE_TPMI_STROKE_GWAS_LOOKUP_ONLY",
        "phenotype":"433.21_cerebral_artery_occlusion_with_infarction",
        "build":"GRCh38","cases_from_PheWeb":9249,"controls_from_PheWeb":304660,
        "n_variants_scanned":n,"n_targets":len(results),
        "counts":dict(Counter(r["status"] for r in results)),
        "source_file":str(a.summary_gwas),
        "beta_effect_orientation":"ALT_EXPLICITLY_USER_ATTESTED_SOURCE_CONVENTION",
        "sampling_overlap_with_BBJ":"COHORT_ORIGIN_REVIEW_PENDING",
        "independent_EAS_replication_verified":False,
        "causal_inference":"NOT_PERFORMED",
        "ALDH2_finemap":"BLOCKED","ADH1B_finemap":"EXPLORATORY"}
    (a.out_dir/"TPMI_43321_ASSOCIATION_MANIFEST.json").write_text(json.dumps(summary,indent=2)+"\n")
    print("TPMI_43321_ASSOCIATION_LOOKUP_COMPLETE",json.dumps(summary["counts"]))
    for r in results:print("SNP",r["rsid"],r["status"],r.get("tpmi_p","NA"))

if __name__=="__main__":main()
