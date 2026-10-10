#!/usr/bin/env python3
"""Japanese alcohol non-ALDH2 positional sentinel discovery and EAS AIS allele QC.

Genomewide alcohol p<5e-8, n>=100k, EAF 1..99%; 500kb coarse bins and
1Mb positional separation (NOT LD-clumped independent instruments).
EAS AIS matched on exact GRCh37 chromosome-position and both alleles;
A/T and G/C palindrome excluded (even if EAF reported). No MR attempted.
"""
import argparse,csv,gzip,json,math
from collections import Counter
from pathlib import Path
ALCOHOL=Path("/srv/is-analysis/data/is/gwas_exposure/japanese_alcohol_2024/1_Alcohol_intake_Unstratified.tsv.gz")
AIS=Path("/srv/is-analysis/data/is/processed/gigastroke/eas/GCST90104545_AIS_GRCh37.canonical.tsv.gz")
ROOT=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1")
HEAD=("SNP","CHR","POS","EA","NEA","EAF","BETA","SE","P","HetP","N")
CHR={str(x) for x in range(1,23)}
def pair(a,b):return tuple(sorted((a.upper(),b.upper())))
def palindrome(a,b):return pair(a,b) in (("A","T"),("C","G"))
def candidate(r):
    chr_=r["CHR"].replace("chr","")
    if chr_ not in CHR:return None
    pos=int(r["POS"])
    # 10Mb conservative regional exclusion around 12q24:
    if chr_=="12" and 107000000<=pos<=117000000:return None
    ref_alt=r["SNP"].split("_")
    if len(ref_alt)!=4 or ref_alt[0].removeprefix("chr")!=chr_ or int(ref_alt[1])!=pos:
        return None
    ref,alt=ref_alt[2].upper(),ref_alt[3].upper()
    ea,nea=r["EA"].upper(),r["NEA"].upper()
    if any(len(x)!=1 or x not in "ACGT" for x in (ref,alt,ea,nea)):return None
    if pair(ref,alt)!=pair(ea,nea):return None
    if palindrome(ea,nea):return None
    b=float(r["BETA"]);se=float(r["SE"]);eaf=float(r["EAF"]);p=float(r["P"]);n=int(float(r["N"]))
    if not all(map(math.isfinite,(b,se,eaf,p))) or se<=0:return None
    if n<100000 or not(.01<=eaf<=.99) or not(0<=p<5e-8):return None
    sign=1 if ea==alt else -1
    het=r["HetP"]
    return {"chr":chr_,"pos":pos,"reference_allele":ref,"alternative_allele":alt,
            "variant_grch37":f"{chr_}:{pos}:{ref}:{alt}",
            "alcohol_original_effect_allele":ea,"alcohol_original_other_allele":nea,
            "alcohol_beta_ALT":sign*b,"alcohol_se":se,"alcohol_p":p,
            "alcohol_beta_over_se_abs":abs(b/se),"alcohol_ALT_eaf":eaf if sign==1 else 1-eaf,
            "alcohol_n":n,"alcohol_cohort_heterogeneity_p":het,
            "alcohol_p_zero_underflow":int(p==0)}
def select_bins(rows,limit=60):
    bins={};qc=Counter()
    for r in rows:
        qc["total_source_records"]+=1
        x=candidate(r)
        if x is None:continue
        qc["quality_gws_pass_records"]+=1
        key=(x["chr"],x["pos"]//500000)
        if key not in bins or x["alcohol_beta_over_se_abs"]>bins[key]["alcohol_beta_over_se_abs"]:
            bins[key]=x
    ranked=sorted(bins.values(),key=lambda x:-x["alcohol_beta_over_se_abs"])
    keep=[]
    for x in ranked:
        if any(x["chr"]==y["chr"] and abs(x["pos"]-y["pos"])<1000000 for y in keep):
            continue
        x["positional_spacing_to_other_selected"]="ONE_MEGABASE_NOT_LD_CLUMPING"
        x["independent_genetic_IV_certified"]=False
        x["pleiotropy_exclusion_restriction_certified"]=False
        keep.append(x)
        if len(keep)>=limit:break
    qc["subselected_positional_leads"]=len(keep)
    qc["gws_bin_count"]=len(bins)
    return keep,qc
def match_ais(candidates,ais):
    by={(r["chr"],str(r["pos"])):r for r in candidates}
    rows={};qc=Counter()
    for r in ais:
        qc["ais_total_source_rows"]+=1
        key=(r["chr"],r["pos"])
        if key not in by:continue
        x=by[key]
        if key in rows:raise ValueError("Duplicated EAS AIS locus position")
        alleles=(r["effect_allele"].upper(),r["other_allele"].upper())
        if pair(*alleles)!=pair(x["reference_allele"],x["alternative_allele"]):
            qc["allele_pair_discordant"]+=1
            rows[key]=dict(x,ais_qc="ALLELE_PAIR_DISCORDANT",ais_beta_ALT="",ais_p="")
            continue
        b=float(r["beta"]);se=float(r["se"]);p=float(r["p"]);eaf=float(r["eaf"])
        if not all(map(math.isfinite,(b,se,p,eaf))) or se<=0 or not(0<eaf<1 and 0<=p<=1):
            rows[key]=dict(x,ais_qc="INVALID_AIS_STATS",ais_beta_ALT="",ais_p="")
            continue
        sign=1 if r["effect_allele"].upper()==x["alternative_allele"] else -1
        rows[key]=dict(x,ais_qc="ALLELE_HARMONIZED",
            ais_beta_ALT=sign*b,ais_se=se,ais_p=p,
            ais_ALT_eaf=eaf if sign==1 else 1-eaf,
            ais_effect_direction_with_alcohol=("SAME" if sign*b*x["alcohol_beta_ALT"]>0 else "OPPOSITE"),
            ais_nominal_p_below_0_05=int(p<.05),
            between_study_cohort_independence_verified=False)
        qc["ais_exact_allele_match"]+=1
    for k,x in by.items():
        if k not in rows:
            rows[k]=dict(x,ais_qc="NOT_IN_EAS_AIS_SOURCE",ais_beta_ALT="",ais_p="")
    qc["no_ais_counterpart"]=sum(r["ais_qc"]=="NOT_IN_EAS_AIS_SOURCE" for r in rows.values())
    return sorted(rows.values(),key=lambda x:-x["alcohol_beta_over_se_abs"]),qc
def run(a,b,out,limit=60):
    with gzip.open(a,"rt") as h:
        r=csv.DictReader(h,delimiter="\t")
        if tuple(r.fieldnames)!=HEAD:raise ValueError("Alcohol summary schema mismatch")
        alcohol,qc=select_bins(r,limit)
    with gzip.open(b,"rt") as h:
        r=csv.DictReader(h,delimiter="\t")
        if r.fieldnames!=['dataset','phenotype','ancestry','build','chr','pos','effect_allele','other_allele','beta','se','p','eaf','or','variant_id']:
            raise ValueError("EAS stroke GWAS schema mismatch")
        joined,qc2=match_ais(alcohol,r)
    out.mkdir(exist_ok=True,parents=True)
    p=out/"IS_RS671_NON_ALDH2_ALCOHOL_POSITIONAL_GWS_AIS_QC.tsv"
    fieldnames=list(dict.fromkeys(k for row in joined for k in row.keys()))
    with p.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fieldnames,delimiter="\t",extrasaction="ignore")
        w.writeheader();w.writerows(joined)
    verdict={
      "source_alcohol":"Koyanagi_2024_Japanese_zenodo10038152_MD5_VERIFIED",
      "source_stroke":"GIGASTROKE_2022_EAS_AIS_GCST90104545",
      **dict(qc),**dict(qc2),
      "selected_loci_are_independent_IVs":False,
      "LD_clumping_performed":False,
      "source_case_control_overlap_ruled_out":False,
      "genetic_instrument_exclusion_restriction_tested":False,
      "MR_estimated":False,
      "scientific_gate":"EXPLORATORY_SPACED_SENTINELS_NOT_VALID_IV_SET"}
    (out/"IS_RS671_NON_ALDH2_ALCOHOL_POSITIONAL_GWS_AIS_QC_SUMMARY.json").write_text(json.dumps(verdict,indent=2))
    print(json.dumps(verdict,indent=2))
    for r in joined[:20]:
        print("CANDIDATE",r["variant_grch37"],"ALC_Z",round(r["alcohol_beta_over_se_abs"],1),
           "AIS",r.get("ais_beta_ALT"),r.get("ais_p"),r["ais_qc"])
    return verdict
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--alcohol",type=Path,default=ALCOHOL)
    p.add_argument("--ais",type=Path,default=AIS)
    p.add_argument("--out",type=Path,default=ROOT)
    p.add_argument("--limit",type=int,default=60)
    opt=p.parse_args();run(opt.alcohol,opt.ais,opt.out,opt.limit)
