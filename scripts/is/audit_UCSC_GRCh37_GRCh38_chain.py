#!/usr/bin/env python3
"""Original UCSC hg19<->hg38 chain MD5, technical coordinate pilot.

Read two ~200KiB official chains + <10KiB publisher checksum manifests.
Check first ImmuNexUT WASH7P chr1 GRCh38 deletion interval coordinate
geometry; no source REF/ALT sequence, left-normalization or GWAS matches.
Do not label a variant harmonized from chain mappings alone.
"""
from __future__ import annotations
import argparse
from datetime import datetime,timezone
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import tempfile
from urllib.request import Request,urlopen
import urllib.error

from pyliftover import LiftOver

ROOT="https://hgdownload.soe.ucsc.edu/goldenPath/"
CHAIN_IDS={
    "GRCh37_to_GRCh38":("hg19","hg19ToHg38.over.chain.gz"),
    "GRCh38_to_GRCh37":("hg38","hg38ToHg19.over.chain.gz"),
}
MAX_CHAIN_BYTES=5_000_000
MAX_MANIFEST_BYTES=60_000


def get_bounded(url,cap):
    if not url.startswith(ROOT):
        raise ValueError("Source not official UCSC hgdownload")
    req=Request(url,headers={"User-Agent":"IS-analysis-V3-chain-provenance/1.0"})
    with urlopen(req,timeout=65) as resp:
        if resp.status!=200:
            raise ValueError("No public publisher original HTTP 200")
        data=resp.read(cap+1)
    if not data or len(data)>cap:
        raise ValueError("Source exceeded strict metadata cap")
    return data


def checksum_from_publisher_manifest(manifest,chainname):
    matches=[]
    for row in manifest.decode("utf-8","strict").splitlines():
        row=row.strip()
        m=re.fullmatch(r"([0-9a-fA-F]{32})\s+\*?(\S+)",row)
        if not m:
            continue
        if m.group(2)==chainname or m.group(2).endswith("/"+chainname):
            matches.append(m.group(1).lower())
    if len(matches)!=1:
        raise ValueError("UCSC chain fingerprint missing/ambiguous in original manifest")
    return matches[0]


def unique_same_strand_mapping(lifter,chrom,position1):
    if not (position1>0 and re.fullmatch(r"chr(?:[1-9]|1[0-9]|2[0-2])",chrom)):
        raise ValueError("Only positive autosomal 1-based genomic coordinates")
    possibilities=lifter.convert_coordinate(chrom,position1-1)
    if len(possibilities)!=1:
        return {"source_chr":chrom,"source_pos1":position1,
                "status":"UNMAPPED_OR_MULTIMAPPED","alternative_mappings":len(possibilities)}
    mapped=possibilities[0]
    mapped_chr,mapped0,orientation,score=mapped
    if orientation not in ("+","-") or mapped0<0:
        raise ValueError("Invalid chain coordinate mapping")
    if not re.fullmatch(r"chr(?:[1-9]|1[0-9]|2[0-2])",mapped_chr):
        return {"source_chr":chrom,"source_pos1":position1,
                "status":"UNSUPPORTED_CHROMOSOME_TARGET"}
    return {"source_chr":chrom,"source_pos1":position1,
            "mapped_chr":mapped_chr,"mapped_pos1":int(mapped0)+1,
            "orientation":orientation,
            "status":"UNIQUE_CHAIN_COORDINATE_MAPPING_NOT_ALLELE_VALIDATED"}


def map_source_interval(lifter,chrom,start1,end1):
    if end1<start1:
        raise ValueError("Original interval end before start")
    start=unique_same_strand_mapping(lifter,chrom,start1)
    end=unique_same_strand_mapping(lifter,chrom,end1)
    same=(
        start["status"]=="UNIQUE_CHAIN_COORDINATE_MAPPING_NOT_ALLELE_VALIDATED" and
        end["status"]=="UNIQUE_CHAIN_COORDINATE_MAPPING_NOT_ALLELE_VALIDATED" and
        start.get("mapped_chr")==end.get("mapped_chr") and
        start.get("orientation")==end.get("orientation") and
        abs(end["mapped_pos1"]-start["mapped_pos1"])==end1-start1
    )
    return {
        "original_source_chr":chrom,"source_interval_start1":start1,
        "source_interval_end1":end1,
        "coordinate_geometry_consistent":bool(same),
        "first_endpoint":start,"last_endpoint":end,
        "REF_ALT_base_compatibility":"NOT_TESTED",
        "indel_left_normalized":"NOT_TESTED",
        "GWAS_effect_allele_harmonized":False,
        "can_directly_run_stroke_coloc":False,
    }


def run(out):
    out.mkdir(parents=True,exist_ok=True)
    status={
        "checked_utc":datetime.now(timezone.utc).isoformat(),
        "source":"Official UCSC hg19 hg38 liftOver chains",
        "status":"SOURCE_OR_QC_HOLD",
        "allele_orientation_and_reference_validation":"NOT_RUN",
        "eQTL_to_GWAS_reconciled":"NOT_RUN",
        "original_ImmuNexUT_nominal_build":"GRCh38",
        "original_GIGASTROKE_GWAS_build":"GRCh37",
        "converted_variants_are_coloc_ready":False,
    }
    try:
        with tempfile.TemporaryDirectory(prefix="is_orig_chain_") as tmp:
            maps={}
            for direction,(build,filename) in CHAIN_IDS.items():
                parent=ROOT+build+"/liftOver/"
                manifest=get_bounded(parent+"md5sum.txt",MAX_MANIFEST_BYTES)
                wanted=checksum_from_publisher_manifest(manifest,filename)
                data=get_bounded(parent+filename,MAX_CHAIN_BYTES)
                md5=hashlib.md5(data).hexdigest()
                if md5!=wanted or not data.startswith(b"\x1f\x8b"):
                    raise ValueError("Official UCSC chain publisher MD5 or gzip proof failed")
                path=Path(tmp)/(direction+".over.chain")
                with gzip.open(__import__("io").BytesIO(data),"rb") as inp:
                    content=inp.read(24_000_000+1)
                    if len(content)>24_000_000:
                        raise ValueError("Uncompressed chain exceeds strict expected cap")
                path.write_bytes(content)
                lifter=LiftOver(str(path))
                mappings=[
                    unique_same_strand_mapping(lifter,"chr1",1000000),
                    unique_same_strand_mapping(lifter,"chr2",1000000),
                    unique_same_strand_mapping(lifter,"chr5",50000000),
                ]
                entry={
                    "publisher_MD5_chain_pass":True,
                    "original_compressed_bytes":len(data),
                    "original_compressed_md5":md5,
                    "direction":direction,
                    "three_technical_coordinate_sentinals":mappings,
                    "reference_bases_not_validated":True,
                }
                if direction=="GRCh38_to_GRCh37":
                    # First bona fide eQTL source record is an indel, coordinate
                    # conversion is NOT sequence normalization or allele check.
                    entry["WASH7P_original_indel_interval_first_source_record"]=map_source_interval(
                        lifter,"chr1",54712,54716)
                maps[direction]=entry
            status["chains"]=maps
            if not all(m["publisher_MD5_chain_pass"] for m in maps.values()):
                raise ValueError("Chain MD5 not verified")
            status["status"]="OFFICIAL_UCSC_TWO_WAY_CHAIN_MD5_AND_COORDINATE_PILOT_PASS"
    except (OSError,ValueError,urllib.error.URLError,EOFError) as e:
        status["issue"]=f"{type(e).__name__}: {str(e)[:240]}"
    (out/"original_ucsc_two_way_chains.json").write_text(
        json.dumps(status,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    lines=[
        "# Official hg19/GRCh37 ↔ hg38/GRCh38 chain source technical pilot","",
        "**"+status["status"]+"**","",
        "Both original UCSC chains independently matched their publisher MD5.",
        "Unique coordinate mapping is not REF/ALT harmonization and is not indel normalization.",
        "Original source WASH7P interval is a technical sentinel, not an IS signal."
    ]
    if status.get("issue"):
        lines.append("Source/QC issue: "+status["issue"])
    if status.get("chains"):
        interval=status["chains"]["GRCh38_to_GRCh37"]["WASH7P_original_indel_interval_first_source_record"]
        lines.extend([
            f"- chr1:54712–54716 source EAS eQTL indel geometric mapping consistent: {interval['coordinate_geometry_consistent']}",
            f"- REF/ALT verified: {interval['REF_ALT_base_compatibility']}",
            "- GWAS locus actual matching: NOT_RUN",
        ])
    summary="\n".join(lines)+"\n"
    (out/"summary.md").write_text(summary,encoding="utf-8")
    print(summary)
    if target:=os.getenv("GITHUB_STEP_SUMMARY"):
        with open(target,"a",encoding="utf-8") as f:f.write(summary)
    if status["status"]!="OFFICIAL_UCSC_TWO_WAY_CHAIN_MD5_AND_COORDINATE_PILOT_PASS":
        raise SystemExit(1)


if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--out-dir",default="artifacts/is_original_ucsc_liftover")
    a=p.parse_args()
    run(Path(a.out_dir))
