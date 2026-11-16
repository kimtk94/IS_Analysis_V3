#!/usr/bin/env python3
"""Check all pinned BMI candidate chr:BP against canonical 22 EUR LD score records.

No raw reference or GWAS archives retained. Scientific LDSC GO remains HOLD.
The published Zenodo mirror MD5 and GIANT BMI original SHA-256 are reverified.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import tarfile
import urllib.request

from scripts.ckm_stroke.audit_bmi_duplicate_zero_p import iter_original_rows
from scripts.ckm_stroke.audit_eur_ldscore_zenodo import (
    FILE_URL, SIZE, EXPECTED_MD5, DigestReader, safe_member,
    LD_SCORE_PATTERN, M_5_50_PATTERN
)
from scripts.ckm_stroke.verify_bmi_ldsc_candidate import inspect

CHRS=set(range(1,23))
MAX_CANDIDATES=1_500_000
MAX_LINE=8192
MAX_INFLATED_PER_CHR=50_000_000
MAX_ARCHIVE_MEMBERS=130
MAX_ARCHIVE=40_000_000
EXPECTED_BMI_CANDIDATES=1_007_899
EXPECTED_BMI_DUPLICATES=9
EXPECTED_BMI_ROWS=2_336_269
MAX_EXAMPLES=18


def read_candidate_ids(path:Path):
    ids={}
    with gzip.open(path,"rt",encoding="ascii",newline="") as handle:
        if handle.readline().strip()!="SNP\tN\tZ\tA1\tA2":
            raise ValueError("BMI candidate columns inconsistent")
        for line in handle:
            snp=line.split("\t",1)[0]
            if snp in ids:
                raise ValueError("candidate rsID duplicate")
            if len(snp)>128 or not snp:
                raise ValueError("candidate SNP identifier invalid")
            ids[snp]=None
            if len(ids)>MAX_CANDIDATES:
                raise ValueError("candidate exceeds maximum rsIDs")
    if not ids:
        raise ValueError("no candidate rsIDs")
    return ids


def original_candidate_coordinates(positions, original_rows):
    counts=Counter()
    for index,row in original_rows:
        if index==0:
            counts["original_sha256_verified"]+=1
            counts["compressed_original_bytes"]=int(row[1])
            break
        counts["source_rows"]+=1
        if len(row)!=10:
            raise ValueError("BMI original field count")
        ch,bp,snp,*_=row
        if snp not in positions:
            continue
        chrom=int(ch)
        pos=int(bp)
        if not(1<=chrom<=22 and 0<pos<=300_000_000):
            raise ValueError("candidate original coordinate invalid")
        p=(chrom,pos)
        if positions[snp] is None:
            positions[snp]=p
            counts["distinct_candidate_rsids_found_original"]+=1
        elif positions[snp]!=p:
            raise ValueError("BMI candidate rsID maps to multiple coordinates")
        else:
            counts["candidate_extra_original_rows_identical_coordinate"]+=1
    else:
        raise ValueError("BMI original source verification sentinel absent")
    if len(positions)!=counts["distinct_candidate_rsids_found_original"]:
        raise ValueError("some candidate SNP IDs absent from source")
    if counts["original_sha256_verified"]!=1:
        raise ValueError("original source not verified")
    return counts


def match_ldscore_stream(f, chromosome, positions, seen, counts, mismatch_examples):
    with gzip.GzipFile(fileobj=f,mode="rb") as stream:
        head=stream.readline(MAX_LINE)
        if not head or len(head)>=MAX_LINE:
            raise ValueError("unbounded LDscore header")
        headers=head.decode("utf-8-sig").strip().split()
        for k in ("CHR","SNP","BP","L2"):
            if k not in headers:raise ValueError("LDscore missing column "+k)
        ci,si,pi=(headers.index(k) for k in ("CHR","SNP","BP"))
        inflated=len(head)
        for raw in stream:
            inflated+=len(raw)
            if inflated>MAX_INFLATED_PER_CHR or len(raw)>MAX_LINE:
                raise ValueError("LDscore uncompressed size or line cap exceeded")
            p=raw.decode("utf-8").split()
            if len(p)!=len(headers):
                raise ValueError("LDscore malformed row")
            if p[ci]!=str(chromosome):
                raise ValueError("canonical LDscore row has wrong chromosome")
            snp=p[si]
            counts["ldscore_total_rows"]+=1
            if snp in positions:
                counts["candidate_ldscore_matching_rsid_rows"]+=1
                if snp in seen:
                    counts["candidate_repeated_ldscore_rsid_rows"]+=1
                    if len(mismatch_examples)<MAX_EXAMPLES:
                        mismatch_examples.append({"SNP":snp,"type":"duplicate_LDscore_SNP"})
                    continue
                seen.add(snp)
                coord=positions[snp]
                try:bp=int(p[pi])
                except ValueError as e:raise ValueError("LDscore BP not integer") from e
                if coord==(chromosome,bp):
                    counts["candidate_coord_concordant"]+=1
                else:
                    counts["candidate_coord_discordant"]+=1
                    if len(mismatch_examples)<MAX_EXAMPLES:
                        mismatch_examples.append({
                            "SNP":snp,"BMI_CHR":coord[0],"BMI_BP":coord[1],
                            "LD_CHR":chromosome,"LD_BP":bp
                        })


def scan_zenodo_coordinates(positions, opener=urllib.request.urlopen,
                            *, published_size=SIZE, published_md5=EXPECTED_MD5):
    req=urllib.request.Request(FILE_URL,headers={"User-Agent":"CKM-Stroke-coordinate-reference-QC/1.0"})
    result_counts=Counter()
    examples=[]
    seen=set()
    seen_chr=set()
    seen_M=set()
    with opener(req,timeout=60) as response:
        from urllib.parse import urlsplit
        if response.status!=200 or urlsplit(response.geturl()).hostname!="zenodo.org":
            raise ValueError("redirect/status inconsistent with pinned Zenodo")
        length=response.headers.get("Content-Length")
        if not length or not length.isdigit() or int(length)!=published_size:
            raise ValueError("pinned Zenodo Content-Length differs")
        reader=DigestReader(response,MAX_ARCHIVE)
        with tarfile.open(fileobj=reader,mode="r|gz") as archive:
            members=0
            for member in archive:
                members+=1
                if members>MAX_ARCHIVE_MEMBERS:
                    raise ValueError("archive member count exceeds cap")
                safe_member(member.name)
                if not(member.isfile() or member.isdir()):
                    raise ValueError("special archive member not accepted")
                if member.isdir():continue
                hit=LD_SCORE_PATTERN.fullmatch(member.name)
                m=M_5_50_PATTERN.fullmatch(member.name)
                if hit:
                    chrom=int(hit.group(1))
                    if chrom in seen_chr:raise ValueError("duplicate chromosome LDscore member")
                    seen_chr.add(chrom)
                    if member.size>MAX_INFLATED_PER_CHR:
                        raise ValueError("LDscore zipped member oversize")
                    f=archive.extractfile(member)
                    if f is None:raise ValueError("missing LDscore content")
                    match_ldscore_stream(f,chrom,positions,seen,result_counts,examples)
                elif m:
                    chrom=int(m.group(1))
                    if chrom in seen_M:raise ValueError("duplicate M member")
                    seen_M.add(chrom)
                else:
                    result_counts["ignored_noncanonical_members"]+=1
        while reader.read(1_048_576):pass
        if reader.nbytes!=published_size or reader.md5.hexdigest()!=published_md5:
            raise ValueError("original Zenodo compressed MD5/length changed")
        digest=reader.sha256.hexdigest()
    if seen_chr!=CHRS or seen_M!=CHRS:
        raise ValueError("canonical 22 EUR LDscore/M_5_50 missing")
    result_counts["candidate_missing_from_canonical_LDscore"]=len(positions)-len(seen)
    result_counts["candidate_total"]=len(positions)
    result_counts["candidate_original_chr6_MHC_25_34Mb"]=sum(
        chr==6 and 25_000_000<=bp<=34_000_000 for chr,bp in positions.values()
    )
    result_counts["candidate_mhc_exclusion_candidate_only"]=result_counts["candidate_original_chr6_MHC_25_34Mb"]
    return {
        "counts":dict(result_counts),
        "mismatch_examples":examples,
        "reference_archive_sha256":digest,
        "reference_archive_md5":published_md5,
        "canonical_chromosomes":sorted(seen_chr),
        "ldscore_chr6_old_excluded":True,
    }


def audit(candidate,provenance, opener=urllib.request.urlopen):
    out={
        "schema":"ckm-stroke-bmi-candidate-vs-zenodo-eur-coordinate-qc-v1",
        "checked_at_utc":datetime.now(timezone.utc).isoformat(),
        "scientific_status":"HOLD",
        "ready_for_ldsc":False,
        "ldsc_h2_rg_performed":False,
        "full_gwas_retained":False,
        "technical_status":"HOLD_COORDINATE_CHECK",
    }
    try:
        prov=json.loads(Path(provenance).read_text(encoding="utf-8"))
        candidate_qc=inspect(Path(candidate),prov)
        ids=read_candidate_ids(Path(candidate))
        c=original_candidate_coordinates(ids,iter_original_rows(opener))
        if c["source_rows"]!=EXPECTED_BMI_ROWS or len(ids)!=EXPECTED_BMI_CANDIDATES:
            raise ValueError("original BMI candidate source profile differed")
        mapped=scan_zenodo_coordinates(ids,opener)
        s=mapped["counts"]
        if (s["candidate_total"]!=len(ids) or s.get("candidate_coord_discordant",0)>0
            or s.get("candidate_missing_from_canonical_LDscore",0)>0
            or s.get("candidate_repeated_ldscore_rsid_rows",0)>0):
            status="HOLD_COORDINATE_DISCORDANCE_OR_MISSING_SNP"
        else:
            status="PASS_SNP_ID_AND_GRCH37_COORDINATES_FOR_BMI_CANDIDATE"
        out.update({
            "technical_status":status,
            "candidate_integrity":candidate_qc,
            "original":dict(c),
            "coordinate_result":mapped,
            "candidate_expected_plain_sha256":prov["candidate_plaintext_sha256"],
            "limitations":"Reference co-locating rsIDs does not verify reference alleles in LDscore, outcome data allele orientation, per-cohort ancestry or sample overlap.",
        })
    except Exception as exc:
        out["technical_status"]="HOLD_SOURCE_OR_QC_ERROR"
        out["error_type"]=type(exc).__name__
        out["error_detail"]=str(exc)[:250]
    return out


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--candidate",type=Path,required=True)
    p.add_argument("--provenance",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    report=audit(args.candidate,args.provenance)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(report,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({"status":report["technical_status"],
                      "counts":report.get("coordinate_result",{}).get("counts"),
                      "scientific_status":report["scientific_status"],
                      "error":report.get("error_detail")},ensure_ascii=False))
    if report["technical_status"]!="PASS_SNP_ID_AND_GRCH37_COORDINATES_FOR_BMI_CANDIDATE":
        raise SystemExit(1)


if __name__=="__main__":
    main()
