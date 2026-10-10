#!/usr/bin/env python3
"""Acquire 7 alcohol GWAS loci with authentic 1000G EAS phase3 reference.

Indexed remote extraction; 504 original EAS individuals, GRCh37, 500kb either
side of each source-GWS sentinel. Default plan-only. The raw data and PLINK
files are stored OUTSIDE Git; no existing canonical result is overwritten.
"""
import argparse,csv,json,subprocess
from pathlib import Path
ROOT=Path("/srv/is-analysis/data/is/ld_reference/alcohol_region_clump_eas_20261010")
SAMPLES=Path("/srv/is-analysis/data/is/ld_reference/1kg_eas_grch37/EAS.samples.txt")
BASE="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"
SENTINELS={"2":[27730940],"4":[39413780,100239319],
           "9":[38395928,75461066],"12":[106750302,112241766]}
def source_url(chrom):
    if chrom not in SENTINELS:raise ValueError("Chromosome outside requested source")
    return BASE+f"/ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
def build_plan(chrom,sample_file,out):
    sample_ids=[x.strip() for x in sample_file.read_text().splitlines() if x.strip()]
    if len(sample_ids)!=504 or len(set(sample_ids))!=504:
        raise ValueError("1KG EAS sample IDs must be 504 unique IDs")
    intervals=[(pos-500000,pos+500000) for pos in SENTINELS[chrom]]
    r=",".join(f"{chrom}:{l}-{u}" for l,u in intervals)
    return {"chrom":chrom,"region_intervals_500kb":r,"EAS_samples":504,
        "reference_build":"GRCh37_hg19",
        "source":source_url(chrom),"vcf":str(out/f"chr{chrom}.EAS504.regions.vcf.gz"),
        "plink_prefix":str(out/f"chr{chrom}.EAS504.regions"),
        "sentinels":SENTINELS[chrom],"sampling":"1000G_Phase3_EAS_complete_panel"}
def execute(chrom,out,sample_file):
    plan=build_plan(chrom,sample_file,out)
    out.mkdir(parents=True,exist_ok=True)
    vcf=Path(plan["vcf"])
    if not vcf.exists():
        tmp=Path(str(vcf)+".part.gz")
        if tmp.exists():raise RuntimeError("Partial existing reference extract needs manual QC")
        cmd=["bcftools","view","--threads","2","-r",plan["region_intervals_500kb"],
             "-S",str(sample_file),"-m2","-M2","-v","snps","-Oz","-o",str(tmp),
             plan["source"]]
        p=subprocess.run(cmd,capture_output=True,text=True,timeout=900,cwd=str(out))
        if p.returncode:raise RuntimeError("Official indexed source failed: "+p.stderr[-1800:])
        if not tmp.exists() or tmp.stat().st_size<1000:
            raise RuntimeError("Reference source too small or missing")
        tmp.replace(vcf)
    p=subprocess.run(["bcftools","query","-l",str(vcf)],capture_output=True,text=True,
                     timeout=60,cwd=str(out))
    if p.returncode or p.stdout.splitlines()!=sample_file.read_text().splitlines():
        raise ValueError("VCF EAS exact sample count or order different from official panel")
    idx=Path(str(vcf)+".tbi")
    if not idx.exists():
        subprocess.run(["tabix","-p","vcf",str(vcf)],check=True,timeout=120,cwd=str(out))
    count=subprocess.run(["bcftools","index","-n",str(vcf)],
         capture_output=True,text=True,timeout=60,cwd=str(out))
    if count.returncode:raise ValueError("Bad regional reference VCF index")
    n=int(count.stdout.strip())
    if n<100:raise ValueError("Reference SNP coverage anomalously small")
    pref=Path(plan["plink_prefix"])
    if not Path(str(pref)+".pgen").exists():
        p=subprocess.run(["plink2","--vcf",str(vcf),"--double-id",
                "--set-all-var-ids","@:#:$r:$a","--make-pgen","--out",str(pref)],
            capture_output=True,text=True,timeout=360,cwd=str(out))
        if p.returncode:raise RuntimeError("PLINK2 conversion: "+p.stdout[-1400:]+p.stderr[-1000:])
    if not all(Path(str(pref)+ext).is_file() for ext in (".pgen",".pvar",".psam")):
        raise RuntimeError("Regional EAS PLINK triplet incomplete")
    with Path(str(pref)+".psam").open() as f:
        actual=sum(1 for x in f if x.strip() and not x.startswith("#"))
    if actual!=504:raise ValueError("PLINK genotype panel missing EAS samples")
    plan.update({"source_total_genotyped_SNVs":n,
                 "plink_sample_count":actual,"validated":True})
    (out/f"chr{chrom}.region_source_audit.json").write_text(json.dumps(plan,indent=2))
    print(json.dumps(plan,indent=2),flush=True)
    return plan
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,default=ROOT)
    p.add_argument("--samples",type=Path,default=SAMPLES)
    p.add_argument("--chrom",choices=[*SENTINELS,"all"],default="all")
    p.add_argument("--execute",action="store_true")
    a=p.parse_args()
    targets=list(SENTINELS) if a.chrom=="all" else [a.chrom]
    for c in targets:
        if a.execute:execute(c,a.out,a.samples)
        else:print(json.dumps(build_plan(c,a.samples,a.out),indent=2))
