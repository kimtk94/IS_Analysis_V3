#!/usr/bin/env python3
"""Prepare a fail-closed *candidate* BMI LDSC SNP N Z A1 A2 file.

Uses only exact pinned source bytes and the strict intersection of MRCIEU HM3,
Zenodo HM3, and canonical Zenodo 22-chromosome EUR LD score marker IDs.
This is not LDSC h2/rg, and the output is NOT an approved analysis input.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import gzip
import hashlib
import io
import json
import math
import os
from pathlib import Path
import urllib.request

from scripts.ckm_stroke.audit_bmi_hm3_overlap import hm3_catalog, allele_pair
from scripts.ckm_stroke.audit_bmi_duplicate_zero_p import iter_original_rows, signed_beta_se_z
from scripts.ckm_stroke.audit_hm3_zenodo_catalog_diff import (
    parse_archive_response, read_record, MRCIEU_SHA
)
from scripts.ckm_stroke.audit_eur_ldscore_zenodo import FILE_URL

ORIGINAL_BMI_ROWS = 2_336_269
EXPECTED_ZERO_P = 44
EXPECTED_DUPLICATE_EXTRAS = 9
EXPECTED_STRICT_REFERENCE = 1_172_984
DUPLICATE_RSIDS = frozenset((
    "rs11576885", "rs13598", "rs2160421", "rs3780697",
    "rs3819262", "rs573071", "rs627141", "rs7591686", "rs969206"
))
CANDIDATE_COLS = b"SNP\tN\tZ\tA1\tA2\n"
PALINDROMIC = frozenset((frozenset(("A", "T")), frozenset(("C", "G"))))


class OutputWriter:
    """Optional deterministic gzip writer; always hash canonical plaintext rows."""
    def __init__(self, dest: Path | None):
        self.dest=dest
        self.temp=None
        self.handle=None
        self.gzip_file=None
        self.sha256=hashlib.sha256()
        self.bytes=0

    def __enter__(self):
        if self.dest is not None:
            self.dest.parent.mkdir(parents=True, exist_ok=True)
            self.temp=self.dest.with_name(self.dest.name+".partial")
            self.handle=self.temp.open("wb")
            self.gzip_file=gzip.GzipFile(fileobj=self.handle,filename="",mode="wb",mtime=0,compresslevel=6)
        return self

    def write(self, payload: bytes):
        self.sha256.update(payload)
        self.bytes+=len(payload)
        if self.gzip_file is not None:
            self.gzip_file.write(payload)

    def __exit__(self,typ,value,tb):
        try:
            if self.gzip_file is not None:self.gzip_file.close()
        finally:
            if self.handle is not None:self.handle.close()
        if self.temp is not None:
            if typ is None:
                os.replace(self.temp,self.dest)
            else:
                self.temp.unlink(missing_ok=True)


def prepare_rows(original, mrcieu, bundled, ldscore, writer, *,
                 assert_source_profile=True):
    ref_ids=set(mrcieu)&set(bundled)&ldscore
    agree={s for s in ref_ids if mrcieu[s]==bundled[s]}
    if assert_source_profile and (len(ref_ids)!=EXPECTED_STRICT_REFERENCE or
                                  len(agree)!=EXPECTED_STRICT_REFERENCE):
        raise ValueError("strict EUR reference panel differs from audited pinned evidence")
    counts=Counter()
    prev_duplicate={}
    candidate_written=set()
    head_written=False
    writer.write(CANDIDATE_COLS)
    for index,fields in original:
        if index==0:
            counts["completed_source_passes"]+=1
            counts["original_compressed_bytes"]=int(fields[1])
            break
        counts["rows"]+=1
        if len(fields)!=10:
            raise ValueError("changed original BMI columns")
        ch,bp,snp,a1,a2,hrs,beta,se,p,n=fields
        if snp in DUPLICATE_RSIDS:
            raw=tuple(fields)
            if snp in prev_duplicate:
                counts["duplicate_extra_rows"]+=1
                if raw!=prev_duplicate[snp]:
                    raise ValueError("pinned duplicate differs in raw association fields")
                continue
            prev_duplicate[snp]=raw
        pv=float(p)
        if not math.isfinite(pv) or not 0 <= pv <= 1:
            raise ValueError("original P outside 0..1")
        if pv==0.0:
            counts["source_zero_p_nonduplicate"]+=1
        if snp not in agree:
            counts["excluded_not_triple_direct_pair"]+=1
            continue
        if frozenset((a1,a2)) in PALINDROMIC:
            counts["excluded_palindromic"]+=1
            continue
        # The reference panel is an unordered allele-pair concordance list;
        # preserve GIANT Tested_Allele as A1 and its beta sign unmodified.
        if allele_pair(a1,a2)!=mrcieu[snp]:
            counts["excluded_bmi_allele_mismatch"]+=1
            continue
        if snp in candidate_written:
            raise ValueError("unexpected duplicate rsID after pinned deduplication")
        candidate_written.add(snp)
        z=signed_beta_se_z(beta,se)
        sample=float(n)
        if not(math.isfinite(sample) and sample>0 and sample.is_integer()):
            raise ValueError("invalid per-variant integer sample size")
        chrom=int(ch)
        if not(1<=chrom<=22 and bp.isdigit() and 0<int(bp)<=300_000_000):
            raise ValueError("invalid chr/pos in emitted candidate")
        if pv==0.0:counts["candidate_zero_p_signed_z"]+=1
        # Decimal rendering stable across platforms, no synthetic P is written.
        line=f"{snp}\t{int(sample)}\t{format(z,'.15g')}\t{a1}\t{a2}\n".encode("ascii")
        writer.write(line)
        counts["candidate_rows"]+=1
    else:
        raise ValueError("source stream did not reach pinned SHA256 verification sentinel")
    if assert_source_profile:
        if (counts["rows"]!=ORIGINAL_BMI_ROWS or counts["duplicate_extra_rows"]!=EXPECTED_DUPLICATE_EXTRAS
            or len(prev_duplicate)!=len(DUPLICATE_RSIDS)
            or counts["completed_source_passes"]!=1
            or counts["source_zero_p_nonduplicate"]!=EXPECTED_ZERO_P
            or counts["excluded_bmi_allele_mismatch"]!=0):
            raise ValueError("pinned original BMI profile has changed")
    if not counts["candidate_rows"]:
        raise ValueError("no harmonized candidate variant rows")
    return {
        "counts":dict(counts),
        "reference_triple_rsids":len(ref_ids),
        "reference_equal_pair_rsids":len(agree),
        "candidate_plaintext_sha256":writer.sha256.hexdigest(),
        "candidate_plaintext_bytes":writer.bytes,
        "candidate_cols":CANDIDATE_COLS.decode("ascii").strip().split("\t"),
        "candidate_file_created":writer.dest is not None,
        "ready_for_ldsc":False,
        "scientific_status":"HOLD",
        "zero_p_policy":"SIGNED_Z_FROM_BETA_SE_NO_P_IMPUTATION",
        "allele_policy":"GIANT_TESTED_ALLELE_IS_A1; STRAND_COMPLEMENTS_NOT_AUTO_FLIPPED; PALINDROMIC_EXCLUDED",
        "duplicate_policy":"DROP_SECOND_ONLY_IF_IDENTICAL_FOR_9_VERIFIED_RSIDS",
    }


def load_references(opener=urllib.request.urlopen):
    mrc,meta=hm3_catalog(opener)
    if meta["ref_sha256"]!=MRCIEU_SHA:
        raise ValueError("MRCIEU original reference hash changed")
    read_record(opener)
    req=urllib.request.Request(FILE_URL,headers={"User-Agent":"CKM-Stroke-BMI-candidate-munge/1.0"})
    with opener(req,timeout=60) as response:
        bundle,_,ld_ids,archive=parse_archive_response(response)
    return mrc,bundle,ld_ids,{
        "mrcieu_sha256":meta["ref_sha256"],
        "zenodo_archive_md5":archive["archive_md5"],
        "zenodo_archive_sha256":archive["archive_sha256"],
    }


def audit(*,out_candidate: Path | None=None,opener=urllib.request.urlopen):
    result={
        "schema":"ckm-stroke-bmi-3panel-ldsc-candidate-v1",
        "at_utc":datetime.now(timezone.utc).isoformat(),
        "output_type":"DERIVED_NOT_SCIENTIFICALLY_APPROVED",
        "run_ldsc":False,"ready_for_ldsc":False,
        "scientific_status":"HOLD",
    }
    try:
        if out_candidate is not None:
            target=out_candidate.expanduser().resolve()
            repo_root=Path(__file__).resolve().parents[2]
            ignored_artifacts=(repo_root/"artifacts").resolve()
            if repo_root in target.parents and ignored_artifacts not in target.parents:
                raise ValueError("candidate output inside Git must be in ignored artifacts/")
        mrc,bundle,ld,meta=load_references(opener)
        with OutputWriter(out_candidate) as writer:
            # iter_original_rows verifies original published-source compressed SHA256.
            converted=prepare_rows(iter_original_rows(opener),mrc,bundle,ld,writer)
        result.update(converted)
        result["sources"]=meta
        result["technical_status"]="PASS_CONSERVATIVE_CANDIDATE_CREATED_NOT_LDSC_APPROVED"
    except Exception as exc:
        result["technical_status"]="HOLD_CANDIDATE_INPUT_OR_TRANSFORM_ERROR"
        result["error_type"]=type(exc).__name__
        result["error_detail"]=str(exc)[:220]
    return result


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out",type=Path,required=True)
    ap.add_argument("--emit-candidate",type=Path,default=None,
                    help="Path for local-derived .sumstats.gz (never committed)")
    args=ap.parse_args()
    if args.emit_candidate is not None and not args.emit_candidate.name.endswith(".sumstats.gz"):
        ap.error("candidate filename must end in .sumstats.gz")
    result=audit(out_candidate=args.emit_candidate)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps({
        "technical_status":result["technical_status"],
        "counts":result.get("counts"),
        "reference_triple_rsids":result.get("reference_triple_rsids"),
        "candidate_plaintext_sha256":result.get("candidate_plaintext_sha256"),
        "scientific_status":result["scientific_status"],
    }))
    if result["technical_status"]!="PASS_CONSERVATIVE_CANDIDATE_CREATED_NOT_LDSC_APPROVED":
        raise SystemExit(1)


if __name__=="__main__":
    main()
