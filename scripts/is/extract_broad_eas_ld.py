#!/usr/bin/env python3
"""Non-destructive 1000 Genomes EAS regional LD extraction for IS GWS discoveries.

Default: plan only. --execute performs indexed remote extraction of one group.
Requires bcftools, tabix, plink2 and locally vetted 504 EAS sample IDs.
No canonical results are overwritten.
"""
import argparse,csv,json,subprocess
from pathlib import Path

ROOT=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1")
SAMPLES=Path("/srv/is-analysis/data/is/ld_reference/1kg_eas_grch37/EAS.samples.txt")
OUT=Path("/srv/is-analysis/data/is/ld_reference/broad_v1/1kg_eas")
BASE="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502"
def reference_url(chrom):
    if chrom not in [str(x) for x in range(1,23)]:
        raise ValueError("Only autosomes 1-22 supported")
    return f"{BASE}/ALL.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
def run(args):
    with (args.root/"IS_EXPANDED_REGION_EXECUTION_MANIFEST.tsv").open() as f:
        region=[r for r in csv.DictReader(f,delimiter="\t") if r["group_id"]==args.group]
    if len(region)!=1:raise ValueError("Unknown group or ambiguous input")
    r=region[0]
    if r["has_gws"]!="1" and not args.allow_suggestive:
        raise ValueError("By default only GWS groups may execute")
    if not args.samples.exists():raise FileNotFoundError(args.samples)
    n=sum(bool(x.strip()) for x in args.samples.open())
    if n!=504:raise ValueError(f"Expected 504 EAS sample IDs; found {n}")
    lead=r["lead_variant"]
    chrom=r["chr"];url=reference_url(chrom)
    start=int(r["window_start"]);end=int(r["window_end"])
    if end<start:raise ValueError("Invalid interval")
    stem=args.out/args.group
    output=Path(str(stem)+".1KG_EAS.GRCh37.vcf.gz")
    prefix=Path(str(stem)+".stable")
    pgen=Path(str(prefix)+".pgen")
    plan={"group_id":args.group,"chr":chrom,"region":f"{chrom}:{start}-{end}",
          "lead_variant":lead,"samples":n,"reference_url":url,
          "vcf":str(output),"pgen":str(pgen),
          "status":"PLAN_ONLY"}
    if not args.execute:
        print(json.dumps(plan,indent=2));return
    args.out.mkdir(parents=True,exist_ok=True)
    if output.exists():
        chk=subprocess.run(["bcftools","view","-h",str(output)],capture_output=True,text=True)
        if chk.returncode:
            raise ValueError("Existing VCF invalid: manual review required, no overwrite")
    else:
        tmp=Path(str(output)+".part.gz")
        if tmp.exists():raise FileExistsError("Partial file exists; manual inspection required: "+str(tmp))
        cmd=["bcftools","view","--threads","2",
             "-r",f"{chrom}:{start}-{end}","-S",str(args.samples),
             "-m2","-M2","-v","snps","-Oz","-o",str(tmp),url]
        proc=subprocess.run(cmd,capture_output=True,text=True,timeout=args.timeout,cwd=str(args.out))
        if proc.returncode:
            raise RuntimeError("remote indexed VCF extraction failed: "+proc.stderr[-1000:])
        if not tmp.exists() or tmp.stat().st_size<1000:
            raise RuntimeError("Empty/undersized indexed VCF; no promotion")
        subprocess.run(["bcftools","view","-h",str(tmp)],check=True,capture_output=True)
        tmp.replace(output)
    idx=Path(str(output)+".tbi")
    if not idx.exists():
        subprocess.run(["tabix","-p","vcf",str(output)],check=True,timeout=120)
    nrecords=subprocess.run(["bcftools","index","-n",str(output)],capture_output=True,text=True)
    if nrecords.returncode:raise ValueError("VCF index invalid: "+nrecords.stderr)
    count=int(nrecords.stdout.strip())
    if count<1:raise ValueError("No variants after QC: "+args.group)
    if not pgen.exists():
        proc=subprocess.run(["plink2","--vcf",str(output),
            "--double-id","--set-all-var-ids","@:#:$r:$a",
            "--make-pgen","--out",str(prefix)],capture_output=True,text=True,timeout=300)
        if proc.returncode:
            raise RuntimeError("PLINK conversion failed: "+proc.stdout[-1800:]+proc.stderr[-1800:])
    if not Path(str(prefix)+".pvar").exists() or not Path(str(prefix)+".psam").exists():
        raise RuntimeError("PGEN triplet incomplete")
    with Path(str(prefix)+".psam").open() as handle:
        loaded_samples=sum(1 for line in handle if line.strip() and not line.startswith("#"))
    if loaded_samples!=n:
        raise RuntimeError(f"PGEN sample mismatch: expected {n}, loaded {loaded_samples}")
    lead_info=lead.split(":")
    records=[]
    with Path(str(prefix)+".pvar").open() as f:
        for row in f:
            if row.startswith("#"):continue
            t=row.strip().split("\t")
            if len(t)<5:continue
            if t[0]==lead_info[0] and t[1]==lead_info[1]:
                records.append((t[3].upper(),t[4].upper()))
    a,b=lead_info[2].upper(),lead_info[3].upper()
    exact=any((x,y)==(a,b) for x,y in records)
    swapped=any((x,y)==(b,a) for x,y in records)
    plan.update({"status":"EXTRACTED_AND_CONVERTED","n_variants":count,
        "lead_allele_status":"EXACT" if exact else ("SWAPPED" if swapped else "LEAD_MISSING_OR_ALLELE_MISMATCH"),
        "fine_mapping_ready":"NO_PENDING_VARIANT_HARMONIZATION_AND_LD_QC"})
    Path(str(prefix)+".REFERENCE_QC.json").write_text(json.dumps(plan,indent=2))
    print(json.dumps(plan,indent=2))
if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--root",type=Path,default=ROOT)
    ap.add_argument("--samples",type=Path,default=SAMPLES)
    ap.add_argument("--out",type=Path,default=OUT)
    ap.add_argument("--group",required=True)
    ap.add_argument("--execute",action="store_true")
    ap.add_argument("--allow-suggestive",action="store_true")
    ap.add_argument("--timeout",type=int,default=900)
    run(ap.parse_args())
