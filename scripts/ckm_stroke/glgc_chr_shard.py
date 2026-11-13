#!/usr/bin/env python3
"""Stream one full GLGC EUR chromosome through SNP-only QC.

Uses only official HTTPS BGZF+TBI sources, requiring a source-observed
header. Writes normalized variants to ephemeral scratch. Only JSON
metadata is uploaded. No full GWAS hash, LDSC or MR inferred.
"""
import argparse
import csv
import datetime as dt
import gzip
import json
import math
import os
from pathlib import Path
import tempfile

from scripts.ckm_stroke.glgc_pysam_probe import fetch_official_header
from scripts.ckm_stroke.glgc_region_probe import BASE, FILE_NAMES, REQUIRED
from scripts.ckm_stroke.variant_qc import variant_qc

OUT_HEADER=["CHR","POS","EA","NEA","BETA","SE","P"]
AUTOSOMES={str(i) for i in range(1,23)}
BASE_DNA={"A","C","G","T"}


def chromosome_qc(trait, chromosome, factory=None, header_fetcher=None,
                  max_rows=5000000, work_dir=None, qc_runner=variant_qc):
    if trait not in FILE_NAMES or str(chromosome) not in AUTOSOMES:
        raise ValueError("invalid trait or autosome")
    if max_rows<100: raise ValueError("row cap must be >=100")
    chromosome=str(chromosome)
    url=BASE+FILE_NAMES[trait]
    source_header=(header_fetcher or fetch_official_header)(url)
    if source_header.get("status")!="PASS":
        return {"status":"HOLD","scientific_status":"NOT_RUN",
                "reason":"official header not verified", "header_evidence":source_header,
                "trait":trait,"chromosome":chromosome}
    header=source_header["header"].split("\t")
    if not REQUIRED.issubset(set(header)) or len(set(header))!=len(header):
        return {"status":"FAIL","scientific_status":"NOT_RUN",
                "reason":"invalid or duplicate source columns", "header":header}
    if factory is None:
        try:
            import pysam
        except ImportError as exc:
            return {"status":"TOOL_UNAVAILABLE","scientific_status":"NOT_RUN",
                    "reason":str(exc)}
        factory=pysam.TabixFile
    indices={k:header.index(k) for k in REQUIRED}
    counts={"original_rows":0,"retained_autosomal_snp":0,
            "non_snp":0,"malformed_columns":0,"wrong_chromosome":0,
            "invalid_coordinate":0,"invalid_alleles":0,"invalid_effect_or_se":0,
            "p_underflow_zero":0,"invalid_p":0,
            "invalid_allele_frequency":0,"palindromic_ambiguous_af_0_42_to_0_58":0}
    examples={}
    def reason(name,position,values):
        counts[name]+=1
        if name not in examples:
            examples[name]={"source_row":position,"tokens_head":values[:8]}
    try:
        with tempfile.TemporaryDirectory(prefix="ckm_glgc_chr_",dir=work_dir) as tmp:
            tmp=Path(tmp)
            normalized=tmp/"normalized.tsv.gz"
            capped=False
            with gzip.open(normalized,"wt",encoding="utf-8",newline="") as out:
                writer=csv.DictWriter(out,fieldnames=OUT_HEADER,delimiter="\t")
                writer.writeheader()
                with factory(url) as source:
                    for raw in source.fetch(chromosome):
                        counts["original_rows"]+=1
                        n=counts["original_rows"]
                        if n>max_rows:
                            capped=True
                            break
                        parts=raw.rstrip("\n").split("\t")
                        if len(parts)!=len(header):
                            reason("malformed_columns",n,parts)
                            continue
                        try:
                            chrom=parts[indices["CHROM"]].removeprefix("chr")
                            pos=int(parts[indices["POS_b37"]])
                        except (ValueError,TypeError):
                            reason("invalid_coordinate",n,parts)
                            continue
                        if chrom!=chromosome:
                            reason("wrong_chromosome",n,parts)
                            continue
                        if pos<1:
                            reason("invalid_coordinate",n,parts)
                            continue
                        ea,nea=parts[indices["ALT"]].upper(),parts[indices["REF"]].upper()
                        if len(ea)!=1 or len(nea)!=1:
                            reason("non_snp",n,parts)
                            continue
                        if ea not in BASE_DNA or nea not in BASE_DNA or ea==nea:
                            reason("invalid_alleles",n,parts)
                            continue
                        try:
                            beta=float(parts[indices["EFFECT_SIZE"]])
                            se=float(parts[indices["SE"]])
                            p=float(parts[indices["pvalue"]])
                        except (ValueError,TypeError):
                            reason("invalid_effect_or_se",n,parts)
                            continue
                        if not all(map(math.isfinite,(beta,se,p))) or se<=0:
                            reason("invalid_effect_or_se",n,parts)
                            continue
                        if p==0:
                            reason("p_underflow_zero",n,parts)
                            continue
                        if not 0<p<=1:
                            reason("invalid_p",n,parts)
                            continue
                        try:
                            af=float(parts[indices["POOLED_ALT_AF"]])
                        except (ValueError,TypeError):
                            reason("invalid_allele_frequency",n,parts)
                            continue
                        if not math.isfinite(af) or not 0<=af<=1:
                            reason("invalid_allele_frequency",n,parts)
                            continue
                        if frozenset((ea,nea)) in (frozenset(("A","T")),frozenset(("C","G"))):
                            if 0.42<=af<=0.58:
                                counts["palindromic_ambiguous_af_0_42_to_0_58"]+=1
                        writer.writerow({"CHR":chrom,"POS":pos,"EA":ea,"NEA":nea,
                                         "BETA":beta,"SE":se,"P":p})
                        counts["retained_autosomal_snp"]+=1
            if not counts["retained_autosomal_snp"]:
                return {"status":"HOLD","scientific_status":"NOT_RUN",
                        "reason":"no valid SNPs in chromosome shard",
                        "counts":counts,"examples":examples,"chromosome":chromosome,
                        "trait":trait}
            variant=qc_runner(normalized,work_dir=tmp)
            if capped:
                status="HOLD"
            elif variant.get("status")!="PASS":
                status="HOLD"
            elif any(counts[k] for k in ("malformed_columns","wrong_chromosome","invalid_coordinate",
                                         "invalid_alleles","invalid_effect_or_se","invalid_p",
                                         "invalid_allele_frequency")):
                status="HOLD"
            else:
                status="PASS"
            return {
                "status":status,"scientific_status":"HOLD",
                "source":"GLGC_Graham2021","trait":trait,"ancestry":"EUR",
                "chromosome":chromosome,"source_url":url,
                "scope":"ONE_FULL_CHROMOSOME_ONLY_NOT_WHOLE_GWAS",
                "header_evidence":source_header,
                "counts":counts,"examples":examples,"row_cap_reached":capped,
                "variant_qc":variant,
                "full_file_source_sha256_verified":False,
                "analysis_ready_data_persisted":False,
                "note":"Only chr shard and SNP subset; indels and P=0 underflow explicitly counted. "
                       "No whole-GWAS checksum, cross-GWAS harmonization, LDSC or MR performed."
            }
    except (OSError,RuntimeError,ValueError,KeyError,TimeoutError) as exc:
        return {"status":"SOURCE_UNAVAILABLE","scientific_status":"NOT_RUN",
                "reason":str(exc)[:320],"source_url":url,"chromosome":chromosome,"trait":trait}


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--trait",required=True,choices=sorted(FILE_NAMES))
    ap.add_argument("--chromosome",default="22")
    ap.add_argument("--out",required=True)
    ap.add_argument("--max-rows",type=int,default=5000000)
    args=ap.parse_args()
    result=chromosome_qc(args.trait,args.chromosome,max_rows=args.max_rows)
    result["timestamp_utc"]=dt.datetime.now(dt.timezone.utc).isoformat()
    result["source_commit"]=os.getenv("GITHUB_SHA","")
    dest=Path(args.out)
    dest.parent.mkdir(parents=True,exist_ok=True)
    dest.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(result,ensure_ascii=False),flush=True)
    if os.getenv("GITHUB_STEP_SUMMARY"):
        with open(os.environ["GITHUB_STEP_SUMMARY"],"a",encoding="utf-8") as fh:
            fh.write("\n### GLGC chr"+str(args.chromosome)+" "+args.trait+"\n")
            fh.write("- Technical: "+result.get("status","unknown")+"\n")
            fh.write("- Scientific: "+result.get("scientific_status","NOT_RUN")+"\n")
            fh.write("- Original rows: "+str(result.get("counts",{}).get("original_rows",0))+"\n")
            fh.write("- Retained SNPs: "+str(result.get("counts",{}).get("retained_autosomal_snp",0))+"\n")
    return int(result.get("status")!="PASS")


if __name__=="__main__":
    raise SystemExit(main())
