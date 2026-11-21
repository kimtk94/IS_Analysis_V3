#!/usr/bin/env python3
"""Full normalized CKDGen GWAS variant-level QC, bounded-memory / disk sort.

Input: validated, standardized autosomal SNP summary statistics.
Outputs: JSON evidence only. POS / unordered allele-pair duplicates are
flagged; allele harmonization and cohort sample overlap are NOT proved.
"""
import argparse
import csv
import datetime as dt
import gzip
import json
import math
import os
from pathlib import Path
import subprocess
import tempfile

HEADER={"CHR","POS","EA","NEA","BETA","SE","P"}
PAL={frozenset(("A","T")),frozenset(("C","G"))}
AUTOSOMES={str(x) for x in range(1,23)}


def variant_qc(source, work_dir=None, stride=101, sort_memory="256M"):
    if stride<1:
        raise ValueError("stride must be positive")
    source=Path(source)
    if not source.is_file():
        raise FileNotFoundError(source)
    with tempfile.TemporaryDirectory(prefix="ckm_variant_qc_",dir=work_dir) as wd:
        wd=Path(wd)
        keys=wd/"variant-keys.tsv"
        sorted_keys=wd/"keys.sorted.tsv"
        counts={c:0 for c in sorted(AUTOSOMES,key=int)}
        pal=0
        n=0
        samples=0
        mismatches=0
        inconsistency_examples=[]
        with gzip.open(source,"rt",encoding="utf-8",newline="") if str(source).endswith(".gz") else open(
                source,"rt",encoding="utf-8",newline="") as inp, keys.open("w",encoding="utf-8") as out:
            reader=csv.DictReader(inp,delimiter="\t")
            if not reader.fieldnames or not HEADER.issubset(reader.fieldnames):
                raise ValueError("normalized GWAS missing columns")
            for row in reader:
                n+=1
                chrom=row["CHR"]
                pos=int(row["POS"])
                ea,nea=row["EA"].upper(),row["NEA"].upper()
                beta,se,p=map(float,(row["BETA"],row["SE"],row["P"]))
                if (chrom not in AUTOSOMES or pos<=0 or len(ea)!=1 or len(nea)!=1
                    or ea not in "ACGT" or nea not in "ACGT" or ea==nea
                    or not all(map(math.isfinite,(beta,se,p)))
                    or se<=0 or not 0<p<=1):
                    raise ValueError("invalid standardized variant at row "+str(n))
                counts[chrom]+=1
                if frozenset((ea,nea)) in PAL:
                    pal+=1
                a,b=sorted((ea,nea))
                out.write(f"{chrom}\t{pos:012d}\t{a}\t{b}\t{ea}\t{nea}\n")
                if n%stride==0:
                    samples+=1
                    z=abs(beta/se)
                    expected=math.erfc(z/math.sqrt(2))
                    if expected>=1e-150 and p>=1e-150:
                        difference=abs(math.log10(p)-math.log10(expected))
                        if difference>1.0:
                            mismatches+=1
                            if len(inconsistency_examples)<8:
                                inconsistency_examples.append({
                                    "variant":f"{chrom}:{pos}:{ea}:{nea}",
                                    "reported_p":p,"beta":beta,"se":se,
                                    "wald_p":expected,"abs_log10_p_diff":round(difference,3)})
        if n==0:
            raise ValueError("empty normalized GWAS")
        with sorted_keys.open("w",encoding="utf-8") as out:
            proc=subprocess.run(["sort","-T",str(wd),"-S",sort_memory,str(keys)],
                                stdout=out,stderr=subprocess.PIPE,text=True,
                                env={**os.environ,"LC_ALL":"C"},check=False,timeout=1800)
        if proc.returncode!=0:
            raise RuntimeError("variant external sort failed: "+proc.stderr[-400:])
        duplicate_rows=0
        duplicate_groups=0
        opposing_orientation_groups=0
        duplicate_examples=[]
        prev_key=None
        group_size=0
        group_orientations=set()
        def add_group(key,size,orientations):
            nonlocal duplicate_rows, duplicate_groups, opposing_orientation_groups
            if key is None or size<=1:
                return
            duplicate_groups+=1
            duplicate_rows+=size-1
            if len(orientations)>1:
                opposing_orientation_groups+=1
            if len(duplicate_examples)<10:
                duplicate_examples.append({"chr":key[0],"pos":int(key[1]),
                                           "alleles":key[2]+"/"+key[3],
                                           "rows":size,"orientations":sorted(orientations)})
        with sorted_keys.open(encoding="utf-8") as inp:
            for line in inp:
                fields=line.rstrip("\n").split("\t")
                if len(fields)!=6:
                    raise ValueError("variant sort key corrupted")
                key=tuple(fields[:4])
                orientation="/".join(fields[4:])
                if prev_key!=key:
                    add_group(prev_key,group_size,group_orientations)
                    group_size=0
                    group_orientations=set()
                    prev_key=key
                group_size+=1
                group_orientations.add(orientation)
            add_group(prev_key,group_size,group_orientations)
    status="HOLD" if duplicate_rows or mismatches else "PASS"
    return {
        "status":status,"scientific_status":"HOLD",
        "rows":n,"distinct_variant_keys":n-duplicate_rows,
        "duplicate_excess_rows":duplicate_rows,
        "duplicate_key_groups":duplicate_groups,
        "duplicate_opposing_effect_allele_groups":opposing_orientation_groups,
        "duplicate_examples":duplicate_examples,
        "palindromic_rows":pal,"palindromic_fraction":round(pal/n,6),
        "chromosome_counts":counts,
        "p_consistency_sample_stride":stride,
        "p_consistency_sampled":samples,
        "p_consistency_abs_log10_threshold":1.0,
        "p_consistency_mismatches":mismatches,
        "p_consistency_examples":inconsistency_examples,
        "warning":"P vs beta/SE is an approximate diagnostic, not an error proof; "
                  "study-specific phenotype scale, sample overlap, build, genomic control "
                  "and external allele harmonization remain unverified."
    }


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--input",required=True)
    parser.add_argument("--out",required=True)
    parser.add_argument("--sample-stride",type=int,default=101)
    args=parser.parse_args()
    try:
        result=variant_qc(args.input,stride=args.sample_stride)
    except (OSError,ValueError,RuntimeError,subprocess.TimeoutExpired) as exc:
        result={"status":"FAIL","scientific_status":"NOT_RUN","error":str(exc)}
    result["timestamp_utc"]=dt.datetime.now(dt.timezone.utc).isoformat()
    result["source_commit"]=os.getenv("GITHUB_SHA","")
    out=Path(args.out)
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps(result,ensure_ascii=False),flush=True)
    return 1 if result["status"]=="FAIL" else 0


if __name__=="__main__":
    raise SystemExit(main())
