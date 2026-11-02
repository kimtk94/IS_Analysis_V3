#!/usr/bin/env python3
"""GTEx v8 historical GENCODE v26 gene TSS ±1Mb cis-window audit.

This is an independent *annotation-level* geometry screen, not proof the
original imported GTEx FastQTL used identical collapsed gene annotation,
reference FASTA, tested SNP denominator, or exact QTL sample inclusion.

Downloads a bounded official GENCODE v26 GTF in ephemeral GitHub CI. Never
uses requester-pays GCS, downloads any whole 3.4GB QTL file or calls SSH.
"""
import argparse
import gzip
import hashlib
import json
import re
import tempfile
import urllib.request
from pathlib import Path

URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_26/gencode.v26.annotation.gtf.gz"
EXPECTED={
 "FGF5":("ENSG00000138675","4",79000000,81500000),
 "C4orf22":("ENSG00000197826","4",79000000,81500000),
 "CALHM2":("ENSG00000138172","10",101500000,106500000),
 "NEURL":("ENSG00000107954","10",101500000,106500000),
}
MAX_COMPRESSED_BYTES=85_000_000
CIS_RADIUS=1_000_000
ATTR_GENE=re.compile(r'(?:^|;\s*)gene_id\s+"([^"]+)"')

def gene_records(gtf):
    found={}
    allowed={r[0] for r in EXPECTED.values()}
    with gzip.open(gtf,"rt",encoding="utf-8",newline="") as handle:
        for line in handle:
            if line.startswith("#"):continue
            bits=line.rstrip("\n").split("\t")
            if len(bits)!=9 or bits[2]!="gene":continue
            ids=ATTR_GENE.search(bits[8])
            if not ids:continue
            gene=ids.group(1).split(".")[0]
            if gene not in allowed:continue
            chrom=bits[0].removeprefix("chr")
            start,end=int(bits[3]),int(bits[4])
            strand=bits[6]
            if strand not in ("+","-") or start<1 or end<start:
                raise ValueError("Unexpected GENCODE gene feature or strand")
            if gene in found:
                raise ValueError("Duplicate GENCODE v26 gene feature")
            tss=start if strand=="+" else end
            found[gene]={"gene_id":gene,"chromosome_GRCh38":chrom,
                         "start":start,"end":end,"strand":strand,"TSS_GRCh38":tss}
    if set(found)!=allowed:
        raise ValueError("Official annotation lacks one of four stable genes")
    return found

def compare_gene_cis_boundaries(genes):
    result={}
    for name,(gene,chrom,left,right) in EXPECTED.items():
        item=genes[gene]
        start=item["TSS_GRCh38"]-CIS_RADIUS
        end=item["TSS_GRCh38"]+CIS_RADIUS
        same_chr=chrom==item["chromosome_GRCh38"]
        inside=same_chr and left<=start and right>=end
        result[name]={
            **item,"previous_screen_region":f"{chrom}:{left}-{right}",
            "GENCODE_v26_gene_TSS_plusminus1Mb_start":start,
            "GENCODE_v26_gene_TSS_plusminus1Mb_end":end,
            "both_ends_inside_predeclared_interval":inside,
            "left_margin_bp":start-left,
            "right_margin_bp":right-end,
            "gene_model_is_confirmed_original_GTEx_collapsed_model":False,
            "final_original_GTEx_cis_bed_boundaries_approved":False,
        }
    return result

def inspect(source_file):
    source_file=Path(source_file)
    if not source_file.is_file() or source_file.stat().st_size>MAX_COMPRESSED_BYTES:
        raise ValueError("GENCODE source file unavailable or exceeds bounded download")
    import zlib
    h=hashlib.sha256()
    count=0
    with open(source_file,"rb") as f:
        while True:
            part=f.read(1024*1024)
            if not part:break
            count+=len(part)
            if count>MAX_COMPRESSED_BYTES:raise ValueError("Compressed input cap")
            h.update(part)
    found=gene_records(source_file)
    comp=compare_gene_cis_boundaries(found)
    all_inside=all(x["both_ends_inside_predeclared_interval"] for x in comp.values())
    return {
        "schema":"IS_GTExV8_GENCode26_gene_TSS_screen_v1",
        "official_gencode_url":URL,
        "compressed_bytes":count,"compressed_SHA256":h.hexdigest(),
        "gene_model_source":"GENCODE release 26 gene GTF",
        "original_GTExv8_tss_1Mb_method_documented":True,
        "gene_boundary_screen":comp,
        "four_gene_TSS_windows_covered":all_inside,
        "geometry_gate":(
            "GENCODE26_TSS_1MB_PROXY_INSIDE_SCREENS_NOT_ORIGINAL_GENE_MODEL_CERTIFIED"
            if all_inside else "HOLD_PREDECLARED_INTERVAL_TRUNCATES_GENECODE26_PROXY_CIS"),
        "GTEx_original_collapsed_gene_model_independently_verified":False,
        "source_specific_GTEx_original_n_variants_independently_verified":False,
        "whole_GTEx_original_3_4GB_source_checksum_verified":False,
        "QTD2023_release6_bytes_identified":False,
        "coloc_performed":False,"scientific_status":"HOLD"
    }

def download_and_inspect():
    request=urllib.request.Request(URL,headers={"User-Agent":"IS-GTExV8-GENCODE26-TSS-Source-Audit/1.0"})
    with tempfile.TemporaryDirectory(prefix="gencode26_tss_") as temp:
        dest=Path(temp)/"gencode.v26.annotation.gtf.gz"
        total=0
        with urllib.request.urlopen(request,timeout=50) as source, dest.open("wb") as out:
            if source.url!=URL:
                raise ValueError("Unexpected GENCODE source redirect")
            while True:
                part=source.read(1024*1024)
                if not part:break
                total+=len(part)
                if total>MAX_COMPRESSED_BYTES:
                    raise ValueError("Compressed source exceeds hard bound")
                out.write(part)
        return inspect(dest)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--local-gtf",type=Path,default=None)
    args=p.parse_args()
    try:
        report=inspect(args.local_gtf) if args.local_gtf is not None else download_and_inspect()
    except (OSError,ValueError,EOFError,UnicodeError) as exc:
        report={"schema":"IS_GTExV8_GENCode26_gene_TSS_screen_v1",
                "geometry_gate":"HOLD_OFFICIAL_ANNOTATION_SOURCE_UNAVAILABLE_OR_MALFORMED",
                "failure_class":type(exc).__name__,"failure_reason_bounded":str(exc)[-180:],
                "scientific_status":"HOLD","coloc_performed":False}
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,sort_keys=True,indent=2)+"\n")
    print("IS_GTEX_TSS_GEOMETRY="+json.dumps(report,sort_keys=True,separators=(",",":")))

if __name__=="__main__":
    main()
