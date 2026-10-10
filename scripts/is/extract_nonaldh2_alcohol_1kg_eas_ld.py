#!/usr/bin/env python3
"""Exact 1KG phase3 EAS GRCh37 seven-variant LD audit for alcohol/AIS rs671.

Remote indexed 1000G EBI release/20130502 v5b genotype VCF; 504 EAS
sample list already established, not a full genome download.
Hard requirements: exact genomic allele REF/ALT, same sample order,
biallelic SNPs, no dosage missingness. Signed ALT dosage r; multi-chrom
pairwise correlations exploratory/finite-sample, not formal LD-clumping.
No causal MR/pleiotropy conclusion inferred from r2.
"""
import argparse,csv,json,math,subprocess
from pathlib import Path
import numpy as np
ROOT=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/expanded_reference_v2/g0022_full_locus_ais_v1")
OUT=Path("/srv/is-analysis/data/is/ld_reference/nonaldh2_alcohol_eas_20261010")
SAMPLES=Path("/srv/is-analysis/data/is/ld_reference/1kg_eas_grch37/EAS.samples.txt")
URLBASE="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"
REGIONS={
  "2":[("2",27730940,"T","C")],
  "4":[("4",39413780,"A","G"),("4",100239319,"T","C")],
  "9":[("9",38395928,"T","C"),("9",75461066,"T","C")],
  "12":[("12",106750302,"A","G"),("12",112241766,"G","A")]}
def url(chrom):
    if chrom not in REGIONS:raise ValueError("Unsupported chromosome")
    return f"{URLBASE}/ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
def manifest():
    return {f"{ch}:{pos}:{ref}:{alt}":(ch,pos,ref,alt)
      for group in REGIONS.values() for ch,pos,ref,alt in group}
def get_plan(chrom,out,samples):
    targets=REGIONS[chrom]
    if len(samples)!=504 or len(set(samples))!=504:
        raise ValueError("EAS ancestry sample count/uniqueness changed")
    region=",".join(f"{c}:{pos}-{pos}" for c,pos,ref,alt in targets)
    return {"chromosome":chrom,"region":region,"source_url":url(chrom),
      "EAS_n":len(samples),"output":str(out/f"chr{chrom}.7snp_eas.grch37.vcf.gz"),
      "expected_variants":[f"{c}:{pos}:{ref}:{alt}" for c,pos,ref,alt in targets]}
def extract_chrom(chrom,out,sample_file,execute=False):
    s=[x.strip() for x in sample_file.read_text().splitlines() if x.strip()]
    plan=get_plan(chrom,out,s)
    if not execute:return plan
    out.mkdir(parents=True,exist_ok=True)
    target=Path(plan["output"])
    if not target.exists():
        tmp=Path(str(target)+".part.gz")
        if tmp.exists():raise FileExistsError("Partial tmp file present; manual inspection")
        cmd=["bcftools","view","--threads","1","-r",plan["region"],
            "-S",str(sample_file),"-m2","-M2","-v","snps","-Oz","-o",str(tmp),plan["source_url"]]
        p=subprocess.run(cmd,capture_output=True,text=True,timeout=280)
        if p.returncode:
            raise RuntimeError("Indexed official source failed "+p.stderr[-1000:])
        if not tmp.exists() or tmp.stat().st_size<1000:
            raise RuntimeError("Reference extract empty/undersized")
        tmp.replace(target)
    if not Path(str(target)+".tbi").exists():
        subprocess.run(["tabix","-p","vcf",str(target)],check=True,timeout=90)
    recorded_samples=subprocess.run(["bcftools","query","-l",str(target)],
             capture_output=True,text=True,timeout=30)
    if recorded_samples.returncode or recorded_samples.stdout.splitlines()!=s:
        raise ValueError("Sample order between EAS panel and original VCF changed")
    plan["EAS_sample_order_exact_verified"]=True
    p=subprocess.run(["bcftools","query","-f","%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n",str(target)],
                     capture_output=True,text=True,timeout=50)
    if p.returncode:raise RuntimeError(p.stderr)
    gt={}
    wanted_by_pos={(c,str(pos)):(ref,alt) for c,pos,ref,alt in REGIONS[chrom]}
    observed_alleles={}
    for raw in p.stdout.splitlines():
        t=raw.split("\t")
        if len(t)!=4+504:raise ValueError("Reference VCF selected sample n not 504")
        chrom0,pos,ref,alt=t[:4]
        key=(chrom0,pos)
        if key not in wanted_by_pos:continue
        want_ref,want_alt=wanted_by_pos[key]
        if (ref,alt)==(want_ref,want_alt):
            swap=False
        elif (ref,alt)==(want_alt,want_ref):
            swap=True
        else:
            raise ValueError("Reference allele pair mismatch for "+chrom0+":"+pos)
        vid=f"{chrom0}:{pos}:{want_ref}:{want_alt}"
        if vid in gt:raise ValueError("Duplicated genomic SNP position")
        dosage=[]
        for x in t[4:]:
            if x in ("./.",".|.","0|.","1|.",".|0",".|1","0/.","1/.","./0","./1"):
                dosage.append(float("nan"));continue
            gtparts=x.replace("|","/").split("/")
            if len(gtparts)!=2 or any(q not in ("0","1") for q in gtparts):
                raise ValueError("Unexpected genotype (phase/ploidy/multiallelic): "+x)
            reference_alt_dosage=sum(int(q) for q in gtparts)
            dosage.append(2-reference_alt_dosage if swap else reference_alt_dosage)
        gt[vid]=np.asarray(dosage,dtype=float)
        observed_alleles[vid]={"1kg_REF":ref,"1kg_ALT":alt,
            "GWAS_REF":want_ref,"GWAS_ALT":want_alt,
            "genotype_dosage_2_minus_for_swapped_REF_ALT":swap}
    for v in plan["expected_variants"]:
        if v not in gt:raise ValueError("Requested allele pair not found in 1KG EAS: "+v)
    plan["ref_alt_harmonization_status"]=observed_alleles
    plan["swapped_REF_ALT_source_count"]=sum(int(x["genotype_dosage_2_minus_for_swapped_REF_ALT"]) for x in observed_alleles.values())
    plan["source_genotypes_exact_verified"]=len(gt)
    plan["missing_genotype_per_variant"]={k:int(np.isnan(v).sum()) for k,v in gt.items()}
    data=out/f"chr{chrom}.7snp_eas.dosages.tsv"
    with data.open("w",newline="") as f:
        w=csv.writer(f,delimiter="\t")
        w.writerow(["EAS_sample_ID"]+plan["expected_variants"])
        for i,sample in enumerate(s):
            w.writerow([sample]+[gt[x][i] for x in plan["expected_variants"]])
    (out/f"chr{chrom}.source_qc.json").write_text(json.dumps(plan,indent=2))
    print("SOURCE_QC",json.dumps(plan,indent=2),flush=True)
    return plan
def build_ld(out,results):
    expect=list(manifest())
    if set(results)!=set(REGIONS):raise ValueError("Missing ancestry reference chromosome")
    arrays={};s0=None
    for chrom in REGIONS:
        src=out/f"chr{chrom}.7snp_eas.dosages.tsv"
        with src.open() as f:
            r=csv.DictReader(f,delimiter="\t")
            if set(r.fieldnames or [])!={"EAS_sample_ID"}|set(x for x in expect if x.startswith(chrom+":")):
                raise ValueError("Reference dosage table schema mismatch")
            rows=list(r)
        sample=[x["EAS_sample_ID"] for x in rows]
        if len(sample)!=504 or (s0 is not None and s0!=sample):
            raise ValueError("EAS sample order mismatch or <504")
        s0=sample
        for v in expect:
            if v.startswith(chrom+":"):
                arrays[v]=np.array([float(x[v]) for x in rows])
    pairs=[];byid={}
    for v in expect:
        x=arrays[v];present=np.isfinite(x)
        if present.sum()<450:raise ValueError("Too much missing reference dosage "+v)
        af=float(np.nanmean(x)/2)
        if not (.005<=af<=.995):raise ValueError("Monomorphic or very rare reference genotype "+v)
        byid[v]={"n_nonmissing":int(present.sum()),"1KG_EAS_ALT_allele_freq":af,
                 "n_missing":int((~present).sum())}
    for i,a in enumerate(expect):
        for b in expect[i+1:]:
            x=arrays[a];y=arrays[b];present=np.isfinite(x)&np.isfinite(y)
            xx=x[present];yy=y[present]
            if len(xx)<450 or np.std(xx)==0 or np.std(yy)==0:raise ValueError("Insufficient genotype variance")
            r=float(np.corrcoef(xx,yy)[0,1]);r2=r*r
            pairs.append({"SNP_a":a,"SNP_b":b,"chrom_a":a.split(":")[0],"chrom_b":b.split(":")[0],
               "same_chromosome":int(a.split(":")[0]==b.split(":")[0]),
               "distance_bp":abs(int(a.split(":")[1])-int(b.split(":")[1])) if a.split(":")[0]==b.split(":")[0] else "",
               "signed_ALT_dosage_r":r,"ALT_dosage_r2":r2,
               "EAS_n_pair_complete":len(xx),
               "exact_1KG_EAS_phased_genotypes":True,
               "reference_genotype_LD_provenance":"1000G_phase3_EAS_GRCh37_504",
               "MR_independence_pleiotropy_and_cohort_overlap_validated":False})
    p=out/"IS_RS671_ALCOHOL_7SNP_REAL_EAS_LD_PAIRS.tsv"
    with p.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(pairs[0]),delimiter="\t")
        w.writeheader();w.writerows(pairs)
    summary={"reference":"1KG_phase3_20130502_v5b_GRCh37_EAS",
      "reference_samples":504,
      "exact_verified_variants":len(arrays),
      "all_positional_ALCOHOL_sentinels_available":True,
      "SNPs_checked":list(arrays),
      "ALT_EAS_reference_af":byid,
      "pairwise_r2_results":len(pairs),
      "max_pair_r2":max(x["ALT_dosage_r2"] for x in pairs),
      "same_chrom_pair_r2":{f'{x["SNP_a"]}|{x["SNP_b"]}':x["ALT_dosage_r2"] for x in pairs if x["same_chromosome"]},
      "ld_ref_validated_genotype_n":504,
      "positional_IVs_pleiotropy_pass":False,
      "unbiased_MR_estimated":False,
      "status":"REAL_EAS_LD_PAIRWISE_CHECK_COMPLETE_STILL_NOT_IV_VALIDATION"}
    (out/"IS_RS671_ALCOHOL_7SNP_REAL_EAS_LD_SUMMARY.json").write_text(json.dumps(summary,indent=2))
    print("LD_SUMMARY",json.dumps(summary,indent=2))
    return summary
def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,default=OUT)
    p.add_argument("--samples",type=Path,default=SAMPLES)
    p.add_argument("--chrom",choices=list(REGIONS)+["all"],default="all")
    p.add_argument("--execute",action="store_true")
    p.add_argument("--calculate",action="store_true")
    a=p.parse_args();targets=list(REGIONS) if a.chrom=="all" else [a.chrom]
    if a.calculate:
        if not a.execute:raise ValueError("--calculate needs --execute")
        if a.chrom!="all":raise ValueError("LD table requires all 4 chromosomes")
    result={}
    for chrom in targets:
        result[chrom]=extract_chrom(chrom,a.out,a.samples,a.execute)
        if not a.execute:print(json.dumps(result[chrom],indent=2))
    if a.calculate:build_ld(a.out,result)
if __name__=="__main__":main()
