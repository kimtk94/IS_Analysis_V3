#!/usr/bin/env python3
"""Non-destructive normalization of verified GIGASTROKE EUR/multi-ancestry GWAS.

The source lacks reference alleles; NEVER reinterpret effect/other allele as REF/ALT.
Stores allele-pair keys, explicit effect orientation, and QC ledger.
Only processes an existing MD5-matched raw file. No internet access in this tool.
"""
import argparse
import csv
import gzip
import hashlib
import json
import math
from collections import Counter
from pathlib import Path

MANIFEST=Path("/srv/is-analysis/results/is/stage5_functional/broad_discovery_v1/GIGASTROKE_SOURCE_VERIFIED_V3.tsv")
RAW=Path("/srv/is-analysis/data/is/reference/gigastroke/additional")
OUT=Path("/srv/is-analysis/data/is/processed/gigastroke/broad_v1")
COLS=["accession","ancestry","phenotype","build","chr","pos","variant_pair_id",
      "effect_allele","other_allele","beta","se","p","eaf","or","reported_n",
      "status","ref_alt_status"]
SOURCE_FIELDS={"chromosome","base_pair_location","effect_allele_frequency",
    "beta","standard_error","p_value","effect_allele","other_allele"}
def md5file(p):
    h=hashlib.md5()
    with p.open("rb") as f:
        for block in iter(lambda:f.read(2**20),b""):h.update(block)
    return h.hexdigest()
def clean(x):
    return str(x or "").upper().strip()
def fnum(s):
    try:return float(s)
    except (TypeError,ValueError):return float("nan")
def normalize(manifest,raw,out,accession):
    with manifest.open() as f:
        item=[r for r in csv.DictReader(f,delimiter="\t") if r["accession"]==accession]
    if len(item)!=1:raise ValueError("Unrecognized accession")
    study=item[0]
    path=raw/f"{accession}_buildGRCh37.tsv.gz"
    if not path.is_file():raise FileNotFoundError(path)
    observed=md5file(path)
    if observed!=study["source_md5"]:
        raise ValueError("Source MD5 mismatch: no output generated")
    if study["genome_build"]!="GRCh37":raise ValueError("Genome build mismatch")
    out.mkdir(parents=True,exist_ok=True)
    target=out/f"{accession}_{study['phenotype']}_{study['ancestry_verified']}_GRCh37.allele_pairs.tsv.gz"
    if target.exists():raise FileExistsError("Refusing to overwrite canonical-like output")
    temp=Path(str(target)+".part")
    summary=Counter()
    try:
        with gzip.open(path,"rt",newline="") as source,gzip.open(temp,"wt",newline="") as dst:
            rows=csv.DictReader(source,delimiter="\t")
            if not SOURCE_FIELDS.issubset(rows.fieldnames or []):
                raise ValueError("Unrecognized GWAS source schema")
            writer=csv.DictWriter(dst,fieldnames=COLS,delimiter="\t")
            writer.writeheader()
            for row in rows:
                summary["total"]+=1
                chrom=clean(row["chromosome"]).removeprefix("CHR")
                a=clean(row["effect_allele"]);b=clean(row["other_allele"])
                pos=fnum(row["base_pair_location"])
                beta=fnum(row["beta"]);se=fnum(row["standard_error"])
                p=fnum(row["p_value"]);eaf=fnum(row["effect_allele_frequency"])
                if (chrom not in {str(i) for i in range(1,23)}
                   or not math.isfinite(pos) or pos<1 or pos!=int(pos)
                   or not a or not b or a==b or any(c not in "ACGT" for c in a+b)
                   or not all(math.isfinite(x) for x in (beta,se,p,eaf))
                   or se<=0 or not 0<=p<=1 or not 0<eaf<1):
                    summary["invalid"]+=1;continue
                status="P_VALUE_UNDERFLOW_ZERO" if p==0 else "QC_VALID"
                if p==0:summary["p_underflow_zero"]+=1
                pair=":".join([chrom,str(int(pos))]+sorted((a,b)))
                data={"accession":accession,"ancestry":study["ancestry_verified"],
                    "phenotype":study["phenotype"],"build":"GRCh37","chr":chrom,
                    "pos":int(pos),"variant_pair_id":pair,"effect_allele":a,
                    "other_allele":b,"beta":row["beta"],"se":row["standard_error"],
                    "p":row["p_value"],"eaf":row["effect_allele_frequency"],
                    "or":row.get("odds_ratio",""),"reported_n":study["source_total_n"],
                    "status":status,"ref_alt_status":"UNKNOWN_REQUIRES_GENOME_REFERENCE"}
                writer.writerow(data)
                summary["valid"]+=1
                if p<=5e-8:summary["gws_variants"]+=1
                if p<=1e-6:summary["suggestive_or_better"]+=1
        temp.replace(target)
    except BaseException:
        if temp.exists():temp.unlink()
        raise
    qc={"accession":accession,"source_md5_verified":observed,
        "normalized_file":str(target),"source_total_n":int(study["source_total_n"]),
        "ancestry":study["ancestry_verified"],"phenotype":study["phenotype"],
        "counts":dict(summary),
        "limitation":"Allele-pair key is not REF/ALT. Effect allele must be harmonized against reference before meta-analysis."}
    Path(str(target)+".qc.json").write_text(json.dumps(qc,indent=2))
    return qc
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--manifest",type=Path,default=MANIFEST)
    p.add_argument("--raw",type=Path,default=RAW)
    p.add_argument("--out",type=Path,default=OUT)
    p.add_argument("--accession",required=True)
    args=p.parse_args()
    print(json.dumps(normalize(args.manifest,args.raw,args.out,args.accession),indent=2))
