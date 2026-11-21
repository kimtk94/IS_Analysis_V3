#!/usr/bin/env python3
"""Compare MRCIEU GRCh37 HapMap3 to bundled Zenodo EUR LDSC HM3 and 22 LD scores.

Reads pinned original public URLs, verifies full SHA-256/MD5, counts SNP/allele
set differences and canonical LDSC score SNPs, and writes aggregate-only JSON.
Neither reference is silently promoted to an LDSC-approved sumstats input.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
import re
import tarfile
import urllib.request

from scripts.ckm_stroke.audit_bmi_hm3_overlap import (
    hm3_catalog, HM3_REF_COMMIT, HM3_URL, orientation, allele_pair
)
from scripts.ckm_stroke.audit_eur_ldscore_zenodo import (
    FILE_URL, EXPECTED_MD5, SIZE, DigestReader, safe_member, read_record,
    LD_SCORE_PATTERN, M_5_50_PATTERN, MAX_INNER_BYTES
)

MRCIEU_SHA = "4ed5a031f88edff7f94663841160d0fff750c45d4d55a630515b4d6de082c49e"
BUNDLED_NAME = "eur_w_ld_chr/w_hm3.snplist"
BUNDLED_SIZE = 17_264_312
ARCHIVE_LIMIT = 40_000_000
SNP_LIMIT = 1_700_000
LD_SCORE_LIMIT = 3_000_000
LINES_CAP = 4096
SAMPLE_CAP = 12
CHR_IDS = set(range(1, 23))


def read_bundle_hm3(fileobj, expected_bytes=None):
    """Read one tar member without extraction, retaining only 1 integer/rsID."""
    if expected_bytes is None:
        expected_bytes = BUNDLED_SIZE
    h = hashlib.sha256()
    records = 0
    bytes_read = 0
    values = {}
    bad_ids = set()
    duplicates = 0
    conflicts = 0
    missing_label_rows = 0
    for raw in fileobj:
        bytes_read += len(raw)
        if bytes_read > expected_bytes or len(raw) > LINES_CAP:
            raise ValueError("bundled HapMap3 exceeded expected byte or row cap")
        h.update(raw)
        line = raw.rstrip(b"\r\n")
        if records == 0 and line == b"SNP\tA1\tA2":
            records += 1
            continue
        if records == 0:
            raise ValueError("bundled HapMap3 SNP A1 A2 header mismatch")
        records += 1
        if records > SNP_LIMIT:
            raise ValueError("bundled HM3 exceeds row cap")
        fields = line.decode("utf-8", "strict").split("\t")
        if len(fields) != 3:
            raise ValueError("bundled HapMap3 row width mismatch")
        snp, a, b = fields
        if not snp or len(snp) > 128:
            raise ValueError("invalid bundled rsID")
        pair = allele_pair(a, b)
        if snp == "NA":
            missing_label_rows += 1
            continue
        if snp in bad_ids:
            duplicates += 1
        elif snp in values:
            duplicates += 1
            if values[snp] != pair:
                bad_ids.add(snp)
                del values[snp]
                conflicts += 1
        else:
            values[snp] = pair
    if bytes_read != expected_bytes or records < 2:
        raise ValueError("bundled HapMap3 truncated member")
    if not values:
        raise ValueError("empty usable bundled HapMap3")
    return values, {
        "bytes": bytes_read,
        "sha256": h.hexdigest(),
        "lines_including_header": records,
        "usable_unique_rsids": len(values),
        "duplicated_rows_non_NA": duplicates,
        "conflicting_distinct_rsids_excluded": len(bad_ids),
        "first_conflicts": sorted(bad_ids)[:SAMPLE_CAP],
        "missing_NA_label_rows": missing_label_rows,
    }


def parse_ldscore_member(fileobj, chromosome: int, *, min_rows=1000):
    """Validate canonical LDscore row chromosome and collect its rsIDs only."""
    ld = set()
    nrows = 0
    inflated = 0
    dup_rows = 0
    with gzip.GzipFile(fileobj=fileobj, mode="rb") as stream:
        header = stream.readline(8192)
        if not header or len(header) > 8192:
            raise ValueError("missing LD score header")
        fields = header.decode("utf-8-sig", "strict").strip().split()
        if not {"CHR", "SNP", "BP", "L2"}.issubset(fields):
            raise ValueError("canonical LD score required header missing")
        ci, si, pi = (fields.index(s) for s in ("CHR", "SNP", "BP"))
        inflated += len(header)
        for line in stream:
            inflated += len(line)
            if inflated > MAX_INNER_BYTES or len(line) > LINES_CAP:
                raise ValueError("LD score member inflated beyond bound")
            parts = line.decode("utf-8", "strict").split()
            if len(parts) != len(fields):
                raise ValueError("LD score column count mismatch")
            if parts[ci] != str(chromosome):
                raise ValueError("LD score chromosome mismatch")
            if not parts[pi].isdigit() or int(parts[pi]) <= 0:
                raise ValueError("invalid LD score base-pair coordinate")
            snp=parts[si]
            if not snp or len(snp)>128:
                raise ValueError("invalid LD score marker")
            nrows += 1
            if nrows > LD_SCORE_LIMIT:
                raise ValueError("LD score rows exceed cap")
            if snp in ld:
                dup_rows += 1
            ld.add(snp)
    if nrows < min_rows:
        raise ValueError("too few canonical LD score SNP rows")
    return ld, {"rows": nrows, "unique": len(ld), "duplicates": dup_rows}


def parse_archive_response(response, *, expected_size=SIZE,
                           expected_md5=EXPECTED_MD5, min_ld_rows=1000):
    if getattr(response, "status", 200) != 200:
        raise ValueError("reference archive HTTP status unexpected")
    if not response.geturl().startswith("https://zenodo.org/"):
        raise ValueError("archive redirected away from pinned Zenodo host")
    length=response.headers.get("Content-Length")
    if not length or not length.isdecimal() or int(length)!=expected_size:
        raise ValueError("archive Content-Length differs from pinned Zenodo")
    hashed=DigestReader(response, ARCHIVE_LIMIT)
    bundle=None
    bundle_meta=None
    ldscore_ids=set()
    per_chr={}
    seen_m=set()
    archive_member_count=0
    seen_noncanon=[]
    with tarfile.open(fileobj=hashed,mode="r|gz") as tar:
        for member in tar:
            archive_member_count += 1
            if archive_member_count > 130:
                raise ValueError("too many archive members")
            safe_member(member.name)
            if not(member.isfile() or member.isdir()):
                raise ValueError("special archive member")
            if member.isdir():
                continue
            if member.name == BUNDLED_NAME:
                if bundle is not None or member.size != BUNDLED_SIZE:
                    raise ValueError("duplicate or mutated bundled HM3 member")
                f=tar.extractfile(member)
                if f is None:raise ValueError("cannot read bundled HM3")
                bundle,bundle_meta=read_bundle_hm3(f)
            elif LD_SCORE_PATTERN.fullmatch(member.name):
                ch=int(LD_SCORE_PATTERN.fullmatch(member.name).group(1))
                if ch in per_chr:
                    raise ValueError("duplicate canonical chr LD score")
                f=tar.extractfile(member)
                if f is None: raise ValueError("cannot open LD score file")
                ids,per_chr[ch]=parse_ldscore_member(f,ch,min_rows=min_ld_rows)
                ldscore_ids.update(ids)
            elif M_5_50_PATTERN.fullmatch(member.name):
                ch=int(M_5_50_PATTERN.fullmatch(member.name).group(1))
                if ch in seen_m:raise ValueError("duplicated canonical M_5_50")
                seen_m.add(ch)
                f=tar.extractfile(member)
                if f is None:raise ValueError("invalid M_5_50 file")
                content=f.read(200).strip().split()
                if len(content)!=1 or float(content[0])<=0:
                    raise ValueError("invalid M_5_50 count")
            else:
                seen_noncanon.append(member.name[:100])
    while hashed.read(1_048_576):
        pass
    if hashed.nbytes!=expected_size or hashed.md5.hexdigest()!=expected_md5:
        raise ValueError("full reference archive compressed MD5/size mismatch")
    if bundle is None or bundle_meta is None:
        raise ValueError("missing bundled HM3")
    if set(per_chr)!=CHR_IDS or seen_m!=CHR_IDS:
        raise ValueError("canonical 22 LD score or M_5_50 files incomplete")
    return bundle,bundle_meta,ldscore_ids,{
        "archive_bytes":hashed.nbytes,
        "archive_md5":hashed.md5.hexdigest(),
        "archive_sha256":hashed.sha256.hexdigest(),
        "members":archive_member_count,
        "canonical_ld_chromosomes":sorted(per_chr),
        "canonical_M_5_50_chromosomes":sorted(seen_m),
        "canonical_ldscore_per_chr":{str(k):v for k,v in sorted(per_chr.items())},
        "unique_ldscore_rsids":len(ldscore_ids),
        "noncanonical_members":seen_noncanon[:20],
    }


def summarize(mrcieu: dict, bundled: dict, ld_ids: set, mrc_meta: dict, bundle_meta: dict,
              ld_meta: dict):
    mrc=set(mrcieu)
    zen=set(bundled)
    both=mrc & zen
    orient=Counter()
    examples={}
    for snp in both:
        # allele_pair int coding is unordered; complement pair is computed
        # by explicit numeric XOR with base complement mapping (A<->T,C<->G).
        x=mrcieu[snp]
        y=bundled[snp]
        if x==y:
            category="pair_equal_unordered"
        else:
            # For unordered pair, complement of ACGT integers [0,1,2,3] = 3-i
            a,b=x>>2,x & 3  # unpack unordered numeric allele pair
            comp= (min(3-a,3-b)<<2) | max(3-a,3-b)
            category=("strand_complement_possible" if comp==y else "allele_incompatible")
        orient[category]+=1
        if category!="pair_equal_unordered":
            sample=examples.setdefault(category,[])
            if len(sample)<SAMPLE_CAP:
                sample.append(snp)
    mrc_only=mrc-zen
    zen_only=zen-mrc
    mrc_ld=mrc & ld_ids
    zen_ld=zen & ld_ids
    triple=both & ld_ids
    direct_triple=sum(1 for snp in triple if mrcieu[snp]==bundled[snp])
    result={
        "mrcieu":mrc_meta,
        "zenodo_hm3":bundle_meta,
        "zenodo_ldscore":ld_meta,
        "comparison":{
            "mrcieu_unique":len(mrc),
            "zenodo_unique":len(zen),
            "shared_rsids":len(both),
            "mrcieu_only":len(mrc_only),
            "zenodo_only":len(zen_only),
            "mrcieu_only_examples":sorted(mrc_only)[:SAMPLE_CAP],
            "zenodo_only_examples":sorted(zen_only)[:SAMPLE_CAP],
            "common_rsid_allele_categories":dict(orient),
            "discordance_examples":examples,
            "ldscore_union_unique":len(ld_ids),
            "mrcieu_intersection_with_ldscore":len(mrc_ld),
            "zenodo_intersection_with_ldscore":len(zen_ld),
            "mrcieu_zenodo_ldscore_triple_intersection":len(triple),
            "triple_direct_unordered_pair_agreement":direct_triple,
            "strict_catalog_interchangeable":(mrc==zen and orient.get("pair_equal_unordered",0)==len(both)),
            "three_way_ldsc_reference_technical_gate":("PASS_THREE_WAY_SNP_ID_AND_PAIR_COUNTED" if triple else "HOLD_NO_THREE_WAY_OVERLAP")
        },
        "source_approved_for_ldsc":False,
        "ldsc_h2_rg_performed":False,
        "ready_for_ldsc":False,
        "scientific_status":"HOLD",
        "no_raw_catalogs_retained":True
    }
    return result


def audit(opener=urllib.request.urlopen):
    out={
        "schema":"ckm-stroke-mrcieu-vs-zenodo-hm3-alelles-v1",
        "checked_at_utc":datetime.now(timezone.utc).isoformat(),
        "scientific_status":"HOLD",
        "ready_for_ldsc":False,
        "technical_status":"HOLD_UNVERIFIED_REFERENCE"
    }
    try:
        mrc,mrc_meta=hm3_catalog(opener)
        if mrc_meta["ref_sha256"]!=MRCIEU_SHA:
            raise ValueError("MRCIEU source original SHA256 changed")
        provenance=read_record(opener)
        req=urllib.request.Request(FILE_URL,headers={"User-Agent":"CKM-Stroke-HM3-Compare/1.0"})
        with opener(req,timeout=60) as response:
            z,zm,ld_ids,ld_meta=parse_archive_response(response)
        result=summarize(mrc,z,ld_ids,mrc_meta,zm,ld_meta)
        out.update(result)
        out["provenance"]=provenance
        out["technical_status"]="PASS_REFERENCE_CATALOG_DIFFERENCES_QUANTIFIED"
    except Exception as exc:
        out["technical_status"]="HOLD_REFERENCE_OR_SCHEMA"
        out["error_type"]=type(exc).__name__
        out["error_detail"]=str(exc)[:250]
    return out


def main():
    p=argparse.ArgumentParser()
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    report=audit()
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps({
        "status":report["technical_status"],
        "comparison":report.get("comparison"),
        "scientific_status":report["scientific_status"],
    },ensure_ascii=False))
    if report["technical_status"]!="PASS_REFERENCE_CATALOG_DIFFERENCES_QUANTIFIED":
        raise SystemExit(1)


if __name__=="__main__":
    main()
