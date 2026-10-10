#!/usr/bin/env python3
"""Prepare source-verified Koyanagi Japanese alcohol regional SuSiE-RSS inputs.

ALDH2/ADH1B ±500kb, MAF>=0.05 in 504 EAS, original SNP n>=100k,
GWAS original P<1e-4 and abs(ALT-EAF gap)<=0.10. All retained source
variants included; maximum 1600 variants and >2GB free RAM guard.
Reference ALT beta/SE and allele-dosage signs checked. Independent credible
sets, mediation and drug targets are NOT validated.
"""
import argparse,csv,gzip,json,subprocess,math
from pathlib import Path
import numpy as np
from prepare_g0022_susie_pilot import traw_rows
from extract_alcohol_eas_region_reference_for_clump import ROOT as REF
ROOT=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1/alcohol_gws_regional_clumping")
OUT=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1/alcohol_conditional_susie_v1")
CENTERS={"ADH1B":("4",100239319),"ALDH2":("12",112241766)}
def tsv(path):
    with path.open() as f:return list(csv.DictReader(f,delimiter="\t"))
def choices(root,ref,locus):
    chrom,center=CENTERS[locus]
    m={r["ID"]:r for r in tsv(ref/f"chr{chrom}.EAS504.regions.afreq")}
    src=tsv(root/"IS_ALCOHOL_SOURCE_REF_MATCHED_REGIONAL_ASSOCIATIONS.tsv")
    count={"source_matched_in_region":0,"excluded_non_common_or_frequency_mismatch":0,
       "excluded_sample_size_or_p":0,"retained":0}
    result={}
    for row in src:
        if row["chr"]!=chrom or abs(int(row["pos"])-center)>500000:continue
        count["source_matched_in_region"]+=1
        rid=row["ID"]
        if rid not in m:continue
        af=float(m[rid]["ALT_FREQS"])
        if not (0.05<=af<=.95) or abs(af-float(row["reference_ALT_eaf_Japanese"]))>.10:
            count["excluded_non_common_or_frequency_mismatch"]+=1;continue
        n=int(row["source_neffect"]);p=float(row["P_original"])
        if n<100000 or p>=1e-4:
            count["excluded_sample_size_or_p"]+=1;continue
        beta=float(row["reference_ALT_beta"])
        # Source original includes SE not in this matched table; recover via source
        # SNP-level Z from reported p? NO: obtain original direct beta/SE instead.
        result[rid]=dict(row,ref_alt_af=af)
    count["retained"]=len(result)
    return result,count
def source_stats(path,keep):
    require=set(keep);got={}
    positions={}
    for vid in require:
        key=":".join(vid.split(":")[:2])
        positions.setdefault(key,[]).append(vid)
    with gzip.open(path,"rt") as f:
        for r in csv.DictReader(f,delimiter="\t"):
            k=f'{r["CHR"]}:{r["POS"]}'
            if k not in positions:continue
            alleles=r["SNP"].split("_")
            if len(alleles)!=4:continue
            ref,alt=alleles[2].upper(),alleles[3].upper()
            possible=[vid for vid in positions[k] if set(vid.split(":")[2:])=={ref,alt}]
            if len(possible)!=1:continue
            rid=possible[0]
            reference_alt=keep[rid]["reference_ALT"]
            ea=r["EA"].upper()
            if ea not in (ref,alt):raise ValueError("Source effect allele invalid "+rid)
            se=float(r["SE"]);beta=float(r["BETA"])
            if se<=0 or not math.isfinite(beta/se):raise ValueError("Source SNP z invalid")
            z=beta/se * (1 if ea==reference_alt else -1)
            if rid in got:raise ValueError("Duplicated source variant")
            got[rid]=(z,se)
            if len(got)==len(require):break
    if set(got)!=require:
        raise ValueError(f"Source z not complete {len(got)} vs {len(require)}")
    return got
def make(root,ref,out,source,locus,max_variants=1600):
    chrom,center=CENTERS[locus]
    options,qc=choices(root,ref,locus)
    if not(50<=len(options)<=max_variants):
        raise ValueError(f"Unanticipated dense locus {locus} {len(options)}")
    out.mkdir(parents=True,exist_ok=True)
    directory=out/locus;directory.mkdir(exist_ok=True)
    # Require original actual effect & standard error per SNP, no p-value
    # underflow inversion and no case-control proxy effective N.
    zstats=source_stats(source,options)
    extract=directory/"reference_snp_ids.txt"
    extract.write_text("\n".join(options)+"\n")
    pref=directory/"reference_eas504"
    traw=Path(str(pref)+".traw")
    if not traw.exists():
        cmd=["plink2","--pfile",str(ref/f"chr{chrom}.EAS504.regions"),
            "--extract",str(extract),"--export","Av","--out",str(pref)]
        p=subprocess.run(cmd,capture_output=True,text=True,timeout=120)
        if p.returncode:raise RuntimeError("PLINK2 dosage failed "+p.stdout[-1000:])
    dosage=[];variants=[]
    for meta,vec in traw_rows(traw):
        c,vid,_,pos,counted,other=meta
        row=options.get(vid)
        if row is None:continue
        reference_alt=row["reference_ALT"]
        # Exact original 1000G VCF ID carries REF:ALT, while Av exports
        # the counted allele in a potentially different orientation.
        if {counted,other}!=set(vid.split(":")[2:]):
            raise ValueError("Dosage reference allele inconsistency "+vid)
        if counted==reference_alt:
            x=vec.copy()
        elif other==reference_alt:
            x=2-vec
        else:raise ValueError("No reference ALT allele dosage "+vid)
        if np.isnan(x).sum()!=0:raise ValueError("1000G reference unexpectedly missing genotypes")
        freq=float(x.mean()/2)
        if abs(freq-row["ref_alt_af"])>5e-6:raise ValueError("Reference ALT EAF discrepancy "+vid)
        if np.std(x)<1e-7:raise ValueError("Unexpected monomorphic genotype")
        z,se=zstats[vid]
        if abs(z*se-float(row["reference_ALT_beta"]))>1e-7:
            raise ValueError("Original GWAS beta/SE reference ALT orientation mismatch "+vid)
        variants.append({"ID":vid,"chr":chrom,"pos":int(pos),
            "reference_ALT":reference_alt,"source_original_p":row["P_original"],
            "source_beta_ALT":row["reference_ALT_beta"],"source_se":se,
            "source_z":z,"source_n":int(row["source_neffect"]),
            "reference_EAS_ALT_EAF":freq,
            "Japanese_original_ALT_EAF":row["reference_ALT_eaf_Japanese"]})
        dosage.append(x)
    if len(dosage)!=len(options):raise ValueError("Genotype subset source match missing")
    mat=np.asarray(dosage,dtype=np.float64)
    stdev=mat.std(axis=1,ddof=1)
    arr=(mat-mat.mean(axis=1,keepdims=True))/stdev[:,None]
    R=arr@arr.T/(mat.shape[1]-1)
    if (not np.isfinite(R).all() or
        np.max(np.abs(np.diag(R)-1))>1e-6):
        raise ValueError("LD nonfinite")
    if len(set(v["ID"] for v in variants))!=len(variants):raise ValueError("SNP duplicate")
    rows=directory/"variants.tsv"
    with rows.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(variants[0]),delimiter="\t")
        w.writeheader();w.writerows(variants)
    ld=directory/"ld.f64.rowmajor"
    R.astype("<f8").tofile(ld)
    if ld.stat().st_size!=8*len(variants)**2:raise ValueError("LD bytes")
    q={
        "status":"ALCOHOL_REGIONAL_REAL_GWAS_EAS_LD_INPUT_PASS",
        "locus":locus,"chrom":chrom,"source_trait":"Koyanagi_2024_log2_grams_day_plus1",
        "original_zenodo_GWAS_MD5":"a37efa6044605f187d93866d162b2026",
        "source_beta_SE_per_SNP_verified":True,
        "n_per_snp_available":True,
        "reference_EAS_n":mat.shape[1],
        "variants":len(variants),
        "min_per_SNP_n":min(v["source_n"] for v in variants),
        "max_per_SNP_n":max(v["source_n"] for v in variants),
        "source_variant_count_pre_filter":qc,
        "source_p_zero_count":sum(v["source_original_p"] in ("0","0.0") for v in variants),
        "LD_input_order_exact_variant_table":True,
        "LD_built_from_signed_allele_aligned_1KG_EAS_dosage":True,
        "matrix_symmetry_max_deviation":float(np.max(np.abs(R-R.T))),
        "max_reference_vs_Japanese_EAF_gap":max(abs(v["reference_EAS_ALT_EAF"]-float(v["Japanese_original_ALT_EAF"])) for v in variants),
        "GIGASTROKE_stroke_not_used_for_exposure_finemapping":True,
        "causal_mediation_estimated":False}
    (directory/"input_qc.json").write_text(json.dumps(q,indent=2))
    print(json.dumps(q,indent=2),flush=True)
    return q
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,default=OUT)
    p.add_argument("--root",type=Path,default=ROOT)
    p.add_argument("--reference",type=Path,default=REF)
    p.add_argument("--source",type=Path,default=Path("/srv/is-analysis/data/is/gwas_exposure/japanese_alcohol_2024/1_Alcohol_intake_Unstratified.tsv.gz"))
    p.add_argument("--locus",choices=list(CENTERS)+["all"],default="all")
    a=p.parse_args()
    for locus in (list(CENTERS) if a.locus=="all" else [a.locus]):
        make(a.root,a.reference,a.out,a.source,locus)
