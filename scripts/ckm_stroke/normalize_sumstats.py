#!/usr/bin/env python3
"""Stream CKDGen-2019 or GLGC-2021 GWAS into an autosomal SNP-only TSV.

No genome-build conversion, strand harmonization, or source provenance inferred.
Excluded variants are counted; any loss forces structural HOLD pending review.
"""
import argparse
import csv
import datetime
import gzip
import hashlib
import json
import math
import os
from pathlib import Path

MAP = {
    "ckdgen2019": {"CHR":"Chr","POS":"Pos_b37","EA":"Allele1","NEA":"Allele2",
                   "BETA":"Effect","SE":"StdErr","P":"P-value"},
    "glgc2021": {"CHR":"CHROM","POS":"POS_b37","EA":"ALT","NEA":"REF",
                 "BETA":"EFFECT_SIZE","SE":"SE","P":"pvalue"}
}
OUT_HEADER=("CHR","POS","EA","NEA","BETA","SE","P")


def open_text(path, mode="rt"):
    return gzip.open(path, mode, encoding="utf-8-sig", newline="") if str(path).endswith(".gz") else open(
        path, mode, encoding="utf-8-sig", newline="")


def sha256(path):
    h=hashlib.sha256()
    with open(path,"rb") as handle:
        while True:
            block=handle.read(1024*1024)
            if not block: break
            h.update(block)
    return h.hexdigest()


def normalize(src, dest, style, expected_sha256=""):
    src=Path(src)
    dest=Path(dest)
    if style not in MAP: raise ValueError("unsupported GWAS source")
    if src.resolve()==dest.resolve(): raise ValueError("input and output must differ")
    if not src.is_file(): raise FileNotFoundError(src)
    original_digest=sha256(src)
    if expected_sha256 and original_digest.lower()!=expected_sha256.lower():
        raise ValueError("source SHA256 mismatch")
    columns=MAP[style]
    dest.parent.mkdir(parents=True,exist_ok=True)
    tmp=dest.with_name(dest.name+".partial")
    count=0
    good=0
    excluded={"non_autosomal":0,"non_snp":0,"missing_or_invalid":0,"p_underflow":0}
    try:
        if dest.suffix!=".gz": raise ValueError("output must end .gz")
        with open_text(src) as input_fh, gzip.open(tmp,"wt",newline="",encoding="utf-8") as output_fh:
            reader=csv.DictReader(input_fh,delimiter="\t")
            if not reader.fieldnames or not set(columns.values()).issubset(reader.fieldnames):
                raise ValueError("source header missing: "+str(set(columns.values())-set(reader.fieldnames or [])))
            writer=csv.DictWriter(output_fh,fieldnames=OUT_HEADER,delimiter="\t")
            writer.writeheader()
            for row in reader:
                count+=1
                try:
                    rec={k:row[v] for k,v in columns.items()}
                    chrom=rec["CHR"].removeprefix("chr")
                    if chrom not in {str(i) for i in range(1,23)}:
                        excluded["non_autosomal"]+=1
                        continue
                    pos=int(rec["POS"])
                    ea=rec["EA"].upper()
                    nea=rec["NEA"].upper()
                    if len(ea)!=1 or len(nea)!=1 or ea not in "ACGT" or nea not in "ACGT" or ea==nea:
                        excluded["non_snp"]+=1
                        continue
                    beta,se,p=map(float,(rec["BETA"],rec["SE"],rec["P"]))
                    if p==0:
                        excluded["p_underflow"]+=1
                        continue
                    if pos<=0 or not all(map(math.isfinite,(beta,se,p))) or se<=0 or not 0<p<=1:
                        raise ValueError("position/beta/se/p")
                    writer.writerow({"CHR":chrom,"POS":pos,"EA":ea,"NEA":nea,
                                     "BETA":beta,"SE":se,"P":p})
                    good+=1
                except (ValueError,TypeError,AttributeError,KeyError):
                    excluded["missing_or_invalid"]+=1
        if good==0: raise ValueError("no usable autosomal SNPs")
        tmp.replace(dest)
    finally:
        if tmp.exists(): tmp.unlink()
    loss=sum(excluded.values())
    return {"status":"HOLD" if loss else "PASS",
            "scientific_status":"HOLD", "source":style,
            "input_file":str(src),"input_sha256":original_digest,
            "output_file":str(dest),"output_sha256":sha256(dest),
            "input_rows":count,"output_rows":good,"excluded":excluded,
            "declared_genome_build":"GRCh37",
            "note":"Source build declared by consortium documentation; no liftover, strand harmonization or independent hash verification."}


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--source",required=True,choices=sorted(MAP))
    ap.add_argument("--input",required=True)
    ap.add_argument("--output",required=True)
    ap.add_argument("--report",required=True)
    ap.add_argument("--expected-sha256",default="")
    opts=ap.parse_args()
    try:
        result=normalize(opts.input,opts.output,opts.source,opts.expected_sha256)
    except (OSError,ValueError) as exc:
        result={"status":"FAIL","scientific_status":"NOT_RUN","error":str(exc)}
    result["time_utc"]=datetime.datetime.now(datetime.timezone.utc).isoformat()
    result["source_commit"]=os.getenv("GITHUB_SHA","")
    dest=Path(opts.report)
    dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(result))
    return int(result["status"]=="FAIL")

if __name__=="__main__":
    raise SystemExit(main())
