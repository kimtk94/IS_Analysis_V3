#!/usr/bin/env python3
"""Original GIGASTROKE EUR/EAS vs publisher-MD5 1000G Phase3 LD BIM pilot.

Read-only public sources; in GitHub Actions downloads ~290/323MB PLINK TGZ
plus ~260/283MB GIGASTROKE GWAS gzip to temporary local disk per ancestry.
Pilot uses chr1 and chr2 to bound RAM, streams the entire original GWAS and
verifies official publisher MD5 for both remote compressed archives.

A *coordinate+unordered allele-set* match is not harmonized effect direction,
an LD calculation, sample-size model, fine-mapping, or coloc. Palindromic
and complement cases are explicitly classified as unresolved.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import re
import tarfile
import tempfile
from urllib.parse import urlparse
from urllib.request import Request, urlopen
import urllib.error

ZENODO_API="https://zenodo.org/api/records/7768714"
ZENODO_FILE_ROOT="https://zenodo.org/records/7768714/files/"
GWAS_FTP=("https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
          "GCST90104001-GCST90105000/")
SOURCES={
    "EUR":{
        "gwas_accession":"GCST90104540",
        "gwas_publisher_md5":"1bfe4ae8a5042fb24cb0a562b05b0d2f",
        "gwas_publisher_size":283482604,
        "ld_archive":"1000G_Phase3_plinkfiles.tgz",
        "ld_publisher_md5":"a7773ab485827b533cb300c76356d76b",
    },
    "EAS":{
        "gwas_accession":"GCST90104545",
        "gwas_publisher_md5":"904e12e33ec3c2eeefe4dc93bec4995c",
        "gwas_publisher_size":260689208,
        "ld_archive":"1000G_Phase3_EAS_plinkfiles.tgz",
        "ld_publisher_md5":"abf52fba6416622ea0757ca9aca51a87",
    },
}
EXPECTED_SOURCE_MODEL={
    "gwas_genome_build":"GRCh37",
    "1000g_reference_build":"GRCh37_stated_in_S_LDSC_source_not_independently_liftover_checked",
}
VALID_BASES=frozenset("ACGT")
COMPLEMENT=dict(zip("ACGT","TGCA"))
LIMIT_BYTES=650*1024*1024


def get_json(url,limit=3_000_000):
    req=Request(url,headers={"User-Agent":"IS-Analysis-V3-SourceLD-Overlap/1.0",
                             "Accept":"application/json"})
    with urlopen(req,timeout=55) as r:
        content=r.read(limit+1)
    if len(content)>limit:
        raise ValueError("Remote metadata exceeds byte limit")
    return json.loads(content)


def download(url,target,expected_md5,expected_size=None):
    if urlparse(url).scheme!="https":
        raise ValueError("Non-HTTPS source")
    md5=hashlib.md5()
    sha=hashlib.sha256()
    total=0
    req=Request(url,headers={"User-Agent":"IS-Analysis-V3-SourceLD-Overlap/1.0"})
    with urlopen(req,timeout=120) as response:
        if response.status!=200:
            raise ValueError("Cannot download full remote source HTTP 200")
        declared=response.headers.get("Content-Length","")
        if declared.isdigit() and int(declared)>LIMIT_BYTES:
            raise ValueError("Source too large")
        with open(target,"wb") as out:
            while True:
                chunk=response.read(4*1024*1024)
                if not chunk:
                    break
                total+=len(chunk)
                if total>LIMIT_BYTES:
                    raise ValueError("Source exceeded bounded download")
                out.write(chunk)
                md5.update(chunk)
                sha.update(chunk)
    if expected_size is not None and total!=expected_size:
        raise ValueError(f"Source size mismatch {total} vs {expected_size}")
    if declared.isdigit() and total!=int(declared):
        raise ValueError("HTTP transfer size mismatch")
    if md5.hexdigest().lower()!=expected_md5.lower():
        raise ValueError("Official publisher MD5 mismatch: "+target.name)
    return {"downloaded_bytes":total,
            "publisher_MD5_verified":True,
            "MD5":md5.hexdigest(),"SHA256_local":sha.hexdigest()}


def metadata_ld_source(ancestry):
    source=SOURCES[ancestry]
    archive=source["ld_archive"]
    base={
        "source":"S-LDSC 1000G Phase3 LD PLINK EUR/EAS",
        "source_record":ZENODO_API,
        "download_file":archive,
        "record_id":7768714,
    }
    meta=get_json(ZENODO_API)
    if meta.get("id")!=7768714:
        raise ValueError("Unexpected Zenodo record metadata")
    files=[f for f in meta.get("files",[]) if f.get("key")==archive]
    if len(files)!=1:
        raise ValueError("Expected precisely one LD archive in Zenodo metadata")
    match=files[0]
    checksum=match.get("checksum","")
    if checksum.lower()!="md5:"+source["ld_publisher_md5"]:
        raise ValueError("Publisher LD checksum drift")
    length=match.get("size")
    if not isinstance(length,int) or not (100_000_000<length<LIMIT_BYTES):
        raise ValueError("Unexpected reference archive size")
    uri=(match.get("links") or {}).get("self")
    if not uri or urlparse(uri).scheme!="https":
        uri=ZENODO_FILE_ROOT+archive
    # Include stable source link and explicit publisher metadata, not a guessed genome build.
    return {**base,"download_url":uri,
            "expected_publisher_md5":source["ld_publisher_md5"],
            "published_bytes":length}


def canonical_chrom(text):
    s=text.strip().removeprefix("chr")
    return int(s) if s.isdigit() and 1<=int(s)<=22 else None


def pair(a,b):
    a=a.upper().strip()
    b=b.upper().strip()
    if len(a)!=1 or len(b)!=1 or a not in VALID_BASES or b not in VALID_BASES or a==b:
        return None
    return frozenset((a,b))


def complement_pair(p):
    return frozenset(COMPLEMENT[a] for a in p)


def reference_bim_coordinates(archive_path,selected_chromosomes):
    """Load only selected chr positional alleles; retain no genotypes or sample IDs."""
    ref=defaultdict(set)
    tally=Counter()
    bim_names=[]
    with tarfile.open(archive_path,"r:gz") as tar:
        for member in tar:
            if not member.isfile() or not member.name.lower().endswith(".bim"):
                continue
            bim_names.append(member.name)
            fileobj=tar.extractfile(member)
            if fileobj is None:
                raise ValueError("Failed to read .bim archive member")
            for raw in fileobj:
                fields=raw.decode("utf-8",errors="strict").split()
                if len(fields)!=6:
                    raise ValueError("PLINK BIM record without 6 fields")
                chrom=canonical_chrom(fields[0])
                if chrom not in selected_chromosomes:
                    continue
                tally["selected_BIM_rows"]+=1
                try:
                    position=int(fields[3])
                except ValueError as exc:
                    raise ValueError("Invalid numeric PLINK coordinate") from exc
                alleles=pair(fields[4],fields[5])
                if position<1 or alleles is None:
                    tally["excluded_invalid_biallelic_BIM_rows"]+=1
                    continue
                ref[(chrom,position)].add(alleles)
                tally["selected_valid_BIM_rows"]+=1
    if not bim_names:
        raise ValueError("LD reference archive has no PLINK BIM member")
    if not ref:
        raise ValueError("LD reference has no selected chromosome SNP records")
    tally["unique_selected_chr_pos"]=len(ref)
    tally["selected_chr_pos_with_multiple_allele_pairs"]=sum(len(x)>1 for x in ref.values())
    return ref,dict(tally),sorted(bim_names)


def classify_variant(chrom,pos,effect,other,reference):
    alleles=pair(effect,other)
    if alleles is None:
        return "UNSUPPORTED_NON_SNP_OR_ALLELE"
    sets=reference.get((chrom,pos))
    if not sets:
        return "COORDINATE_NOT_IN_LD_BIM"
    if len(sets)!=1:
        return "MULTI_ALLELIC_REFERENCE_POSITION_REVIEW"
    ref_pair=next(iter(sets))
    if alleles==ref_pair:
        if alleles in (frozenset(("A","T")),frozenset(("C","G"))):
            return "PALINDROMIC_AMBIGUOUS_NO_EAF_ORIENTATION"
        return "UNORDERED_ALLELE_SET_EXACT"
    if complement_pair(alleles)==ref_pair:
        return "STRAND_COMPLEMENT_UNRESOLVED"
    return "COORDINATE_PRESENT_ALLELES_DISCORDANT"


def classify_gwas(full_gwas,ref,chromosomes):
    counts=Counter()
    per_chr=defaultdict(Counter)
    with gzip.open(full_gwas,"rt",encoding="utf-8",newline="") as handle:
        header=handle.readline().rstrip("\n").split("\t")
        required=("chromosome","base_pair_location","effect_allele","other_allele")
        if any(f not in header for f in required):
            raise ValueError("GIGASTROKE original GWAS fields missing")
        ix={key:header.index(key) for key in required}
        for raw in handle:
            if not raw.strip():
                continue
            counts["all_original_GWAS_variant_rows"]+=1
            fields=raw.rstrip("\n").split("\t")
            chrom=canonical_chrom(fields[ix["chromosome"]])
            if chrom not in chromosomes:
                continue
            try:
                pos=int(fields[ix["base_pair_location"]])
            except ValueError as exc:
                raise ValueError("Invalid original GWAS coordinate") from exc
            if pos<1:
                raise ValueError("Nonpositive original GWAS position")
            classification=classify_variant(chrom,pos,fields[ix["effect_allele"]],
                                            fields[ix["other_allele"]],ref)
            counts["selected_chr_GWAS_variant_rows"]+=1
            counts[classification]+=1
            per_chr[chrom]["selected_GWAS_rows"]+=1
            per_chr[chrom][classification]+=1
    if not counts["selected_chr_GWAS_variant_rows"]:
        raise ValueError("No original GWAS rows in selected chromosome window")
    counts["coordinate_level_overlap"] = (
        counts["UNORDERED_ALLELE_SET_EXACT"]+counts["PALINDROMIC_AMBIGUOUS_NO_EAF_ORIENTATION"]+
        counts["STRAND_COMPLEMENT_UNRESOLVED"]+
        counts["COORDINATE_PRESENT_ALLELES_DISCORDANT"]+
        counts["MULTI_ALLELIC_REFERENCE_POSITION_REVIEW"]
    )
    counts["allele_set_match_excluding_ambiguous"] = counts["UNORDERED_ALLELE_SET_EXACT"]
    counts["selected_chr_coordinate_overlap_pct"] = round(
        counts["coordinate_level_overlap"]/counts["selected_chr_GWAS_variant_rows"]*100,4)
    counts["selected_chr_unordered_allele_match_pct"] = round(
        counts["UNORDERED_ALLELE_SET_EXACT"]/counts["selected_chr_GWAS_variant_rows"]*100,4)
    return dict(counts),{str(k):dict(v) for k,v in sorted(per_chr.items())}


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--ancestry",choices=SOURCES,required=True)
    parser.add_argument("--chromosomes",default="1,2")
    parser.add_argument("--out-dir",default="artifacts/is_gwas_ld_bim_pilot")
    args=parser.parse_args()
    chromosomes={int(x) for x in args.chromosomes.split(",")}
    if not chromosomes or not chromosomes<=set(range(1,23)):
        raise ValueError("Invalid autosome selection")
    out=Path(args.out_dir)
    out.mkdir(parents=True,exist_ok=True)
    result={
        "audit_utc":datetime.now(timezone.utc).isoformat(),
        "ancestry":args.ancestry,
        "chromosome_pilot":sorted(chromosomes),
        "source_type":"REAL_FULL_ORIGINAL_GWAS_AND_1000G_PHASE3_REFERENCE",
        "status":"SOURCE_UNAVAILABLE",
        "GWAS_to_LD_effect_alleles_harmonized":False,
        "LD_matrix_computed":False,
        "sample_size_model_ready":False,
        "MAGMA_ready":False,
        "LDSC_SEG_ready":False,
        "coloc_ready":False,
        "not_an_independent_biological_validation":True,
        "build_assumption":EXPECTED_SOURCE_MODEL,
    }
    try:
        source=SOURCES[args.ancestry]
        ref_info=metadata_ld_source(args.ancestry)
        result["LD_publisher_metadata"]=ref_info
        accession=source["gwas_accession"]
        gwas=GWAS_FTP+accession+"/"+accession+"_buildGRCh37.tsv.gz"
        with tempfile.TemporaryDirectory(prefix="is_ld_overlap_") as tmp:
            root=Path(tmp)
            ld_path=root/"LD_reference.tgz"
            gwas_path=root/"GWAS_original.tsv.gz"
            result["LD_archive"]=download(
                ref_info["download_url"],ld_path,
                ref_info["expected_publisher_md5"],ref_info["published_bytes"])
            reference,ref_tally,bim_members=reference_bim_coordinates(ld_path,chromosomes)
            result["BIM_selected_chromosomes"]=ref_tally
            result["BIM_member_names"]=bim_members[:30]
            # Reduce peak disk use; no need to retain LD archive while reading GWAS.
            ld_path.unlink()
            result["GWAS_original"]=download(
                gwas,gwas_path,source["gwas_publisher_md5"],source["gwas_publisher_size"])
            counts,by_chr=classify_gwas(gwas_path,reference,chromosomes)
            result["GWAS_vs_reference"]=counts
            result["by_chromosome"]=by_chr
        result["status"]="ORIGINAL_MD5_AND_COORDINATE_PILOT_PASS"
    except (urllib.error.URLError,OSError,TimeoutError,ValueError,
            tarfile.TarError,UnicodeDecodeError,KeyError) as exc:
        result["status"]="SOURCE_OR_QC_HOLD"
        result["error_type"]=type(exc).__name__
        result["error_short"]=str(exc)[:250]
    (out/"ld_overlap_pilot.json").write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    lines=["# Original EUR/EAS IS GWAS vs ancestry-matched 1000G Phase3 LD reference",
           "",f"**{args.ancestry}: {result['status']}** | selected chromosomes {sorted(chromosomes)}",
           "",
           "Official-source compressed MD5 verification and chromosome-position/unordered SNP alleles only.",
           "**NOT** aligned allele effects, statistical LD, MAGMA, LDSC-SEG or colocalization.",
           ""]
    for k,v in result.get("GWAS_vs_reference",{}).items():
        lines.append(f"- {k}: {v}")
    if result.get("error_short"):
        lines.append(f"- source/QC issue: {result['error_short']}")
    lines += ["","Benchmark only selected chr GWAS rows against selected chr LD BIM positions.",
              "Do not generalize chr1+2 matching to the entire genome.",
              "Palindromic and complement SNPs require EAF and effect-direction harmonization.",
              "Per-variant sample size unavailable in original GIGASTROKE standardized GWAS."]
    summary="\n".join(lines)+"\n"
    (out/"summary.md").write_text(summary,encoding="utf-8")
    print(summary)
    if os.environ.get("GITHUB_STEP_SUMMARY"):
        with open(os.environ["GITHUB_STEP_SUMMARY"],"a",encoding="utf-8") as handle:
            handle.write(summary)
    return 0 if result["status"]=="ORIGINAL_MD5_AND_COORDINATE_PILOT_PASS" else 1


if __name__=="__main__":
    raise SystemExit(main())
