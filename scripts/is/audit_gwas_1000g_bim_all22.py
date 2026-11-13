#!/usr/bin/env python3
"""All 22 autosomes of original GIGASTROKE vs ancestry-matched 1000G Phase3 BIM.

Streaming/chromosome-sharded; does NOT construct an LD matrix or align betas.
Both original compressed source MD5s checked against publisher before parsing.
No genotypes/phenotypes saved to GitHub artifacts. Temporary files are only
(position, allele1, allele2) for BIM, and (position, effect, other) for GWAS.

QC: all 22 autosomes present, whole GWAS gz CRC, exact reconciliation of
chr1/2 counts with the previous independently successful CI run #43021975612.
"""
from __future__ import annotations
import argparse
from collections import Counter, defaultdict
from contextlib import ExitStack
from datetime import datetime, timezone
import gzip
import json
import hashlib
import subprocess
import time
from urllib.parse import urlparse
import os
from pathlib import Path
import tarfile
import tempfile
import urllib.error

from audit_gwas_1000g_bim_overlap import (
    SOURCES, GWAS_FTP, download, metadata_ld_source,
    canonical_chrom, pair, classify_variant, EXPECTED_SOURCE_MODEL,
)

AUTOSOMES=set(range(1,23))
EXPECTED_GWAS_ROWS={"EUR":7482032,"EAS":6725216}
PREVIOUS_CHR12={
    "EUR":{
        "selected_chr_GWAS_variant_rows":1203049,
        "UNORDERED_ALLELE_SET_EXACT":1012938,
        "PALINDROMIC_AMBIGUOUS_NO_EAF_ORIENTATION":176831,
        "COORDINATE_NOT_IN_LD_BIM":13279,
        "COORDINATE_PRESENT_ALLELES_DISCORDANT":1,
    },
    "EAS":{
        "selected_chr_GWAS_variant_rows":1087989,
        "UNORDERED_ALLELE_SET_EXACT":903835,
        "PALINDROMIC_AMBIGUOUS_NO_EAF_ORIENTATION":160761,
        "COORDINATE_NOT_IN_LD_BIM":23393,
        "COORDINATE_PRESENT_ALLELES_DISCORDANT":0,
    },
}
CATEGORIES=(
    "UNORDERED_ALLELE_SET_EXACT",
    "PALINDROMIC_AMBIGUOUS_NO_EAF_ORIENTATION",
    "STRAND_COMPLEMENT_UNRESOLVED",
    "COORDINATE_NOT_IN_LD_BIM",
    "MULTI_ALLELIC_REFERENCE_POSITION_REVIEW",
    "COORDINATE_PRESENT_ALLELES_DISCORDANT",
    "UNSUPPORTED_NON_SNP_OR_ALLELE",
)
COORDINATE_FOUND=(
    "UNORDERED_ALLELE_SET_EXACT",
    "PALINDROMIC_AMBIGUOUS_NO_EAF_ORIENTATION",
    "STRAND_COMPLEMENT_UNRESOLVED",
    "MULTI_ALLELIC_REFERENCE_POSITION_REVIEW",
    "COORDINATE_PRESENT_ALLELES_DISCORDANT",
)



def metadata_with_retry(ancestry,attempts=4):
    """Transient Zenodo API timeouts are not scientific source absence."""
    for attempt in range(1,attempts+1):
        try:
            print(f"SOURCE_METADATA attempt={attempt}/{attempts} ancestry={ancestry}",flush=True)
            return metadata_ld_source(ancestry)
        except (TimeoutError,urllib.error.URLError,OSError) as exc:
            print(f"SOURCE_METADATA_RETRY error={type(exc).__name__}: {str(exc)[:150]}",flush=True)
            if attempt>=attempts:
                raise
            time.sleep(attempt*4)
    raise RuntimeError("Unreachable metadata retry state")


def streamed_file_checksums(file_path,expected_size,expected_md5):
    """Publisher identity must match exact physical file bytes."""
    if not file_path.exists() or file_path.stat().st_size!=expected_size:
        return None
    md5=hashlib.md5()
    sha256=hashlib.sha256()
    with file_path.open("rb") as stream:
        for chunk in iter(lambda:stream.read(4*1024*1024),b""):
            md5.update(chunk)
            sha256.update(chunk)
    if md5.hexdigest().lower()!=expected_md5.lower():
        return None
    return {
        "publisher_MD5_verified":True,
        "downloaded_bytes":expected_size,
        "MD5":md5.hexdigest(),
        "SHA256_local":sha256.hexdigest(),
    }


def resumable_source_download(url,target,expected_md5,expected_size,attempts=4):
    """GitHub runner curl partial resume + full publisher MD5; fail closed.

    Handles interrupted original file reads but never treats curl success alone
    as an integrity PASS. Limited to public source domains and 650MB cap.
    """
    parsed=urlparse(url)
    if parsed.scheme!="https" or parsed.hostname not in {"zenodo.org","ftp.ebi.ac.uk"}:
        raise ValueError("Unexpected source URL host or scheme")
    if not isinstance(expected_size,int) or not 0<expected_size<650*1024*1024:
        raise ValueError("Unexpected original source size")
    for attempt in range(1,attempts+1):
        proven=streamed_file_checksums(target,expected_size,expected_md5)
        if proven:
            return proven
        if target.exists() and target.stat().st_size>=expected_size:
            # Corrupt/oversized "complete" download: retry from zero, never bypass MD5.
            target.unlink()
        print(f"ORIGINAL_DOWNLOAD attempt={attempt}/{attempts} "
              f"source={parsed.hostname} current_bytes={target.stat().st_size if target.exists() else 0}",
              flush=True)
        cmd=[
            "curl","--fail","--location","--silent","--show-error",
            "--retry","2","--retry-all-errors","--retry-delay","3",
            "--connect-timeout","30","--max-time","800",
            "--speed-time","45","--speed-limit","1024",
            "--continue-at","-",
            "--proto","=https","--proto-redir","=https",
            "--output",str(target),url,
        ]
        try:
            response=subprocess.run(cmd,capture_output=True,text=True,timeout=840,check=False)
        except (OSError,subprocess.TimeoutExpired) as exc:
            print(f"DOWNLOAD_TOOL_ERROR {type(exc).__name__}: {str(exc)[:160]}",flush=True)
            response=None
        if response and response.returncode!=0:
            print(f"DOWNLOAD_NETWORK_RETRY curl_exit={response.returncode}: "
                  f"{response.stderr[-240:]}",flush=True)
        # Even if curl exited with a nonzero network return, the final output can
        # be whole and publisher-identical when the connection closes late.
        proven=streamed_file_checksums(target,expected_size,expected_md5)
        if proven:
            print(f"ORIGINAL_SOURCE_MD5_PASS bytes={expected_size}",flush=True)
            return proven
        if attempt<attempts:
            time.sleep(attempt*3)
    raise ValueError(f"Publisher MD5/size not verified after {attempts} attempts: {parsed.hostname}")

def stage_bim_by_chr(archive,root):
    """Do not extract tar paths. Stream all .bim members to fixed-safe per-chr shards."""
    meta=Counter()
    observed_files=[]
    shards={c:root/f"reference_chr{c}.tsv" for c in sorted(AUTOSOMES)}
    with ExitStack() as stack:
        outputs={c:stack.enter_context(shards[c].open("w",encoding="utf-8"))
                 for c in shards}
        with tarfile.open(archive,"r:gz") as tar:
            for member in tar:
                if not member.isfile() or not member.name.lower().endswith(".bim"):
                    continue
                observed_files.append(member.name)
                inp=tar.extractfile(member)
                if inp is None:
                    raise ValueError("BIM cannot be read: "+member.name)
                for raw in inp:
                    parts=raw.decode("utf-8","strict").split()
                    if len(parts)!=6:
                        raise ValueError("Malformed BIM row")
                    chrom=canonical_chrom(parts[0])
                    if chrom not in AUTOSOMES:
                        meta["bim_non_autosome_rows"]+=1
                        continue
                    meta["all_autosomal_bim_rows"]+=1
                    try:
                        pos=int(parts[3])
                    except ValueError as exc:
                        raise ValueError("Noninteger BIM position") from exc
                    if pos<=0 or pair(parts[4],parts[5]) is None:
                        meta["bim_non_biallelic_or_invalid_rows"]+=1
                        continue
                    outputs[chrom].write(f"{pos}\t{parts[4]}\t{parts[5]}\n")
                    meta[f"bim_rows_chr{chrom}"]+=1
    if not observed_files:
        raise ValueError("No BIM in original LD archive")
    if any(meta[f"bim_rows_chr{c}"]<100 for c in AUTOSOMES):
        raise ValueError("At least one of 22 autosomes has unexpectedly few BIM SNPs")
    return shards,dict(meta),observed_files


def stage_gwas_by_chr(gwas_gzip,root,ancestry):
    """Read entire original gzip so EOF/CRC is verified and every SNP accounted for."""
    meta=Counter()
    shards={c:root/f"gwas_chr{c}.tsv" for c in sorted(AUTOSOMES)}
    with ExitStack() as stack:
        outputs={c:stack.enter_context(shards[c].open("w",encoding="utf-8"))
                 for c in shards}
        with gzip.open(gwas_gzip,"rt",encoding="utf-8",newline="") as source:
            names=source.readline().rstrip("\n").split("\t")
            necessary=("chromosome","base_pair_location","effect_allele","other_allele")
            if any(s not in names for s in necessary):
                raise ValueError("Original GWAS missing required harmonization columns")
            ix=[names.index(x) for x in necessary]
            for line in source:
                if not line.strip():
                    continue
                meta["all_original_gwas_rows"]+=1
                cols=line.rstrip("\n").split("\t")
                if len(cols)<len(names):
                    raise ValueError("Truncated GWAS row")
                chrom=canonical_chrom(cols[ix[0]])
                if chrom not in AUTOSOMES:
                    meta["gwas_non_autosome_rows"]+=1
                    continue
                try:
                    pos=int(cols[ix[1]])
                except ValueError as exc:
                    raise ValueError("GWAS position noninteger") from exc
                if pos<=0:
                    raise ValueError("GWAS coordinate nonpositive")
                outputs[chrom].write(f"{pos}\t{cols[ix[2]]}\t{cols[ix[3]]}\n")
                meta[f"gwas_rows_chr{chrom}"]+=1
    if meta["all_original_gwas_rows"]!=EXPECTED_GWAS_ROWS[ancestry]:
        raise ValueError("Full GWAS row count drift from previously publisher-verified evidence")
    if any(meta[f"gwas_rows_chr{c}"]<100 for c in AUTOSOMES):
        raise ValueError("At least one autosome has unexpectedly few GWAS rows")
    return shards,dict(meta)


def single_chr_overlap(chrom,bim_file,gwas_file):
    """Load one chromosome at a time, not all 22 coordinates into RAM."""
    ref=defaultdict(set)
    counts=Counter()
    with bim_file.open("rt",encoding="utf-8") as inp:
        for row in inp:
            pos,a,b=row.rstrip("\n").split("\t")
            site=int(pos)
            alleles=pair(a,b)
            if alleles is None:
                raise ValueError("BIM staging contained invalid pair")
            ref[(chrom,site)].add(alleles)
            counts["bim_rows"]+=1
    counts["bim_unique_positions"]=len(ref)
    counts["bim_multiallelic_positions"]=sum(len(x)>1 for x in ref.values())
    if not ref:
        raise ValueError("Empty reference chromosome")
    with gwas_file.open("rt",encoding="utf-8") as inp:
        for row in inp:
            pos,effect,other=row.rstrip("\n").split("\t")
            category=classify_variant(chrom,int(pos),effect,other,ref)
            if category not in CATEGORIES:
                raise ValueError("Unexpected allele classification "+category)
            counts["gwas_rows"]+=1
            counts[category]+=1
    if counts["gwas_rows"]<100:
        raise ValueError("Empty or incomplete GWAS chromosome")
    counts["chr_pos_overlap"]=sum(counts[key] for key in COORDINATE_FOUND)
    counts["nonpalindromic_unordered_allele_set_match"]=counts["UNORDERED_ALLELE_SET_EXACT"]
    if counts["chr_pos_overlap"]+counts["COORDINATE_NOT_IN_LD_BIM"]+counts["UNSUPPORTED_NON_SNP_OR_ALLELE"]!=counts["gwas_rows"]:
        raise ValueError("Coordinate-level classification denominator mismatch")
    counts["chr_pos_overlap_pct"]=round(100*counts["chr_pos_overlap"]/counts["gwas_rows"],4)
    counts["nonpalindromic_match_pct_of_all_gwas"]=round(
        100*counts["UNORDERED_ALLELE_SET_EXACT"]/counts["gwas_rows"],4)
    counts["palindromic_pct_of_all_gwas"]=round(
        100*counts["PALINDROMIC_AMBIGUOUS_NO_EAF_ORIENTATION"]/counts["gwas_rows"],4)
    return dict(counts)


def check_previous_chr12(ancestry,chromosome_data):
    summary=Counter()
    for c in (1,2):
        counts=chromosome_data[str(c)]
        summary["selected_chr_GWAS_variant_rows"]+=counts["gwas_rows"]
        for key in CATEGORIES:
            summary[key]+=counts.get(key,0)
    previous=PREVIOUS_CHR12[ancestry]
    wrong={k:{"previous":v,"new":summary[k]} for k,v in previous.items() if summary[k]!=v}
    if wrong:
        raise ValueError("GIGASTROKE 1/2 chr previous verified CI mismatch: "+str(wrong))
    return {"status":"PASS","previous_run_id":43021975612,
            "verified_counts":dict(previous)}


def finish_counts(chromosome_data,stage_gwas):
    aggregate=Counter()
    for c in sorted(AUTOSOMES):
        row=chromosome_data[str(c)]
        aggregate["all_22_autosome_gwas_rows"]+=row["gwas_rows"]
        aggregate["all_22_autosome_BIM_rows"]+=row["bim_rows"]
        aggregate["all_22_chr_pos_overlap"]+=row["chr_pos_overlap"]
        for key in CATEGORIES:
            aggregate[key]+=row.get(key,0)
    if aggregate["all_22_autosome_gwas_rows"]+stage_gwas.get("gwas_non_autosome_rows",0)!=stage_gwas["all_original_gwas_rows"]:
        raise ValueError("Whole GWAS GWAS-to-22-chromosome accounting failure")
    aggregate["22_autosome_coordinate_overlap_pct"]=round(
        100*aggregate["all_22_chr_pos_overlap"]/aggregate["all_22_autosome_gwas_rows"],4)
    aggregate["22_autosome_nonpalindromic_allele_match_pct"]=round(
        100*aggregate["UNORDERED_ALLELE_SET_EXACT"]/aggregate["all_22_autosome_gwas_rows"],4)
    return dict(aggregate)


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--ancestry",required=True,choices=sorted(SOURCES))
    p.add_argument("--out-dir",default="artifacts/is_gwas_ld_all22")
    args=p.parse_args()
    out=Path(args.out_dir)
    out.mkdir(parents=True,exist_ok=True)
    report={
        "ancestry":args.ancestry,
        "analysis_timestamp_utc":datetime.now(timezone.utc).isoformat(),
        "source":"ORIGINAL_MD5_VALIDATED_GIGASTROKE_AND_1000G_PHASE3",
        "status":"SOURCE_OR_QC_HOLD",
        "LD_reference_build":EXPECTED_SOURCE_MODEL["1000g_reference_build"],
        "gwas_build":EXPECTED_SOURCE_MODEL["gwas_genome_build"],
        "all_22_complete":False,
        "allele_effect_direction_harmonized":False,
        "reference_allele_frequencies_verified":False,
        "coloc_ready":False,
        "MAGMA_ready":False,
        "LDSC_SEG_ready":False,
        "LD_matrix_computed":False,
        "eQTL_source_adequate":False,
    }
    try:
        source=SOURCES[args.ancestry]
        ld=metadata_with_retry(args.ancestry)
        report["LD_publisher_source"]=ld
        accession=source["gwas_accession"]
        gwas_url=GWAS_FTP+accession+"/"+accession+"_buildGRCh37.tsv.gz"
        with tempfile.TemporaryDirectory(prefix="is_ld_all22_") as temp:
            root=Path(temp)
            ldf=root/"ld_original.tgz"
            gwasf=root/"gwas_original.tsv.gz"
            report["LD_checksum"]=resumable_source_download(ld["download_url"],ldf,ld["expected_publisher_md5"],ld["published_bytes"])
            print("REFERENCE_BIM_CHR_STAGING_BEGIN",flush=True)
            bim_shards,bim_tally,members=stage_bim_by_chr(ldf,root)
            print("REFERENCE_BIM_CHR_STAGING_PASS",flush=True)
            report["LD_BIM_source_counts"]=bim_tally
            report["LD_BIM_member_count"]=len(members)
            # Free compressed archive space; each chromosome reference BIM shard is enough for positional QC.
            ldf.unlink()
            report["GWAS_checksum"]=resumable_source_download(gwas_url,gwasf,source["gwas_publisher_md5"],source["gwas_publisher_size"])
            print("GWAS_FULL_GZIP_CHR_STAGING_BEGIN",flush=True)
            gwas_shards,gwas_tally=stage_gwas_by_chr(gwasf,root,args.ancestry)
            print("GWAS_FULL_GZIP_CHR_STAGING_PASS",flush=True)
            report["GWAS_source_counts"]=gwas_tally
            gwasf.unlink()
            by_chr={}
            for chrom in sorted(AUTOSOMES):
                by_chr[str(chrom)]=single_chr_overlap(chrom,bim_shards[chrom],gwas_shards[chrom])
                print(f"CHR_MATCH_PASS {chrom}/22 matched={by_chr[str(chrom)]['chr_pos_overlap_pct']}%",flush=True)
                # Fixed residual disk footprint
                bim_shards[chrom].unlink()
                gwas_shards[chrom].unlink()
            report["by_chr"]=by_chr
            report["previous_chr1_2_regression"]=check_previous_chr12(args.ancestry,by_chr)
            report["overall"]=finish_counts(by_chr,gwas_tally)
            print("ALL_22_AND_PRIOR_CHR12_REGRESSION_PASS",flush=True)
            report["status"]="ALL_22_AUTOSOME_ORIGINAL_MD5_AND_BIM_OVERLAP_PASS"
            report["all_22_complete"]=True
    except (urllib.error.URLError,OSError,TimeoutError,ValueError,tarfile.TarError,
            UnicodeDecodeError,KeyError,EOFError) as exc:
        report["error_type"]=type(exc).__name__
        report["error_short"]=str(exc)[:300]
    (out/"all22_audit.json").write_text(json.dumps(report,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    lines=[
        "# All-22 GIGASTROKE EUR/EAS × ancestry-matched 1000G Phase3 BIM SNP QC",
        "",
        f"**{args.ancestry}: {report['status']}**",
        "",
        "Complete original compressed sources are publisher MD5 verified. Per-chromosome BIM and GWAS staged as tiny allele-position records; no genotypes or patient IDs exported.",
        "",
        "**No** effect-allele harmonization, LD statistic or MAGMA/LDSC-SEG/coloc calculation.",
        "",
    ]
    if report.get("overall"):
        ov=report["overall"]
        lines.extend([
            f"- 22 autosome GWAS rows: {ov['all_22_autosome_gwas_rows']:,}",
            f"- LD reference BIM rows: {ov['all_22_autosome_BIM_rows']:,}",
            f"- coordinate overlap: {ov['22_autosome_coordinate_overlap_pct']:.4f}%",
            f"- nonpalindromic unordered allele set exact: {ov['22_autosome_nonpalindromic_allele_match_pct']:.4f}%",
            f"- palindromic ambiguous rows: {ov['PALINDROMIC_AMBIGUOUS_NO_EAF_ORIENTATION']:,}",
            f"- coordinate absent from reference: {ov['COORDINATE_NOT_IN_LD_BIM']:,}",
            f"- previous chr1+2 reproducibility: {report['previous_chr1_2_regression']['status']}",
            "",
            "| Chr | GWAS rows | Coordinate overlap % | Nonpalindromic allele match % |",
            "|---:|---:|---:|---:|",
        ])
        for c,row in sorted(report["by_chr"].items(),key=lambda x:int(x[0])):
            lines.append(f"| {c} | {row['gwas_rows']} | {row['chr_pos_overlap_pct']:.4f}% | {row['nonpalindromic_match_pct_of_all_gwas']:.4f}% |")
    if report.get("error_short"):
        lines.append(f"\nQC HOLD: {report['error_short']}")
    lines.extend([
        "",
        "BIM allele A1/A2 is not automatically the GWAS effect/non-effect allele. Frequency and strand orientation remain unharmonized.",
        "A/T and C/G palindromic SNPs are not allowed to pass effect-direction QC merely by allele-set equality.",
        "GWAS per-SNP N is not supplied by source; matched full-cis EAS eQTL still unavailable, so colocalization remains BLOCKED.",
    ])
    txt="\n".join(lines)+"\n"
    (out/"summary.md").write_text(txt,encoding="utf-8")
    print(txt)
    if os.environ.get("GITHUB_STEP_SUMMARY"):
        with open(os.environ["GITHUB_STEP_SUMMARY"],"a",encoding="utf-8") as fd:
            fd.write(txt)
    return 0 if report["all_22_complete"] else 1


if __name__=="__main__":
    raise SystemExit(main())
