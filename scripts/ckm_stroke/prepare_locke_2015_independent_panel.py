#!/usr/bin/env python3
"""Independently select Locke 2015 EUR BMI variants from Locke + pinned HM3/LD refs.

This NEVER opens or selects variants using Yengo 2018 BMI candidate/output.
The resulting file is a NOT-APPROVED intermediate, NOT LDSC munged sumstats.
Preserve decimal per-SNP N and distinguish beta/SE-derived Z from published GC-P.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
from statistics import median
import urllib.request

from scripts.ckm_stroke.audit_bmi_hm3_overlap import hm3_catalog, allele_pair
from scripts.ckm_stroke.audit_hm3_zenodo_catalog_diff import (
    parse_archive_response, read_record, MRCIEU_SHA
)
from scripts.ckm_stroke.audit_eur_ldscore_zenodo import FILE_URL
from scripts.ckm_stroke.audit_locke_2015_bmi_overlap import (
    SOURCE_URL, SOURCE_SIZE, EXPECTED_SOURCE_SHA256, DigestReader,
    EXPECTED_ORIGINAL_ROWS, EXPECTED_DUPLICATE_EXTRA_ROWS, MAX_ORIGINAL_ROWS,
    MAX_LINE_BYTES, HEAD, MAX_BYTES,
)

# Do not infer duplicates from filename 'uniq': 13 groups with 14 extra source rows.
DUPLICATE_SOURCE_RSIDS=frozenset((
    "rs11576885","rs11771665","rs13598","rs2160421","rs372268",
    "rs3752976","rs3780697","rs3819262","rs573071","rs5820007",
    "rs627141","rs7591686","rs969206",
))
EXPECTED_REF_TRIPLE=1_172_984
PALINDROMES={frozenset(("A","T")),frozenset(("C","G"))}
OUTPUT_HEADER=b"SNP\tN_REPORTED\tZ_FROM_BETA_SE\tA1\tA2\n"
MAX_CANDIDATE_ROWS=1_500_000


class CandidateWriter:
    """Reproducible gzip output under ignored artifacts only; SHA of plaintext."""
    def __init__(self,out:Path|None):
        self.out=out
        self.tmp=None
        self.file=None
        self.gz=None
        self.sha=hashlib.sha256()
        self.length=0
    def __enter__(self):
        if self.out is not None:
            self.out.parent.mkdir(parents=True,exist_ok=True)
            self.tmp=self.out.with_name(self.out.name+".partial")
            self.file=self.tmp.open("wb")
            self.gz=gzip.GzipFile(fileobj=self.file,mode="wb",filename="",mtime=0,compresslevel=6)
        return self
    def write(self,raw:bytes):
        self.sha.update(raw)
        self.length+=len(raw)
        if self.gz is not None:self.gz.write(raw)
    def __exit__(self,typ,value,tb):
        try:
            if self.gz is not None:self.gz.close()
        finally:
            if self.file is not None:self.file.close()
        if self.tmp is not None:
            if typ is None:os.replace(self.tmp,self.out)
            else:self.tmp.unlink(missing_ok=True)


def source_rows(fetcher=urllib.request.urlopen):
    req=urllib.request.Request(SOURCE_URL,headers={"User-Agent":"CKM-Stroke-Locke2015-Independently-Selected/1.0"})
    with fetcher(req,timeout=65) as response:
        if (getattr(response,"status",200)!=200 or
            response.geturl()!=SOURCE_URL or
            response.headers.get("Content-Length")!=str(SOURCE_SIZE)):
            raise ValueError("official source URL, size or status changed")
        reader=DigestReader(response)
        with gzip.GzipFile(fileobj=reader,mode="rb") as gz:
            header=gz.readline(MAX_LINE_BYTES+1).decode("ascii").rstrip("\r\n").split("\t")
            if tuple(header)!=HEAD:raise ValueError("Locke source schema changed")
            for idx,raw in enumerate(gz,start=1):
                if idx>MAX_ORIGINAL_ROWS or len(raw)>MAX_LINE_BYTES:
                    raise ValueError("original rows/line safety cap exceeded")
                fields=raw.decode("ascii","strict").rstrip("\r\n").split("\t")
                if len(fields)!=len(HEAD):raise ValueError("original row width changed")
                yield idx,fields
        while reader.read(1_048_576):pass
        if reader.n!=SOURCE_SIZE or reader.sha.hexdigest()!=EXPECTED_SOURCE_SHA256:
            raise ValueError("full compressed original source SHA256 or length changed")
        yield 0,["__LOCKE_FULL_SOURCE_SHA256_PASS__",str(reader.n)]


def load_strict_references(fetcher=urllib.request.urlopen):
    # Independent reference libraries; NO 2018 BMI variant calls or candidate files.
    mrc,meta=hm3_catalog(fetcher)
    if meta["ref_sha256"]!=MRCIEU_SHA:
        raise ValueError("MRCIEU HM3 raw reference SHA256 mutated")
    record=read_record(fetcher)
    req=urllib.request.Request(FILE_URL,headers={"User-Agent":"CKM-Stroke-Locke2015-Independent-HM3/1.0"})
    with fetcher(req,timeout=70) as response:
        bundled,bundle_metadata,ld_ids,archive=parse_archive_response(response)
    eligible=set(mrc).intersection(bundled).intersection(ld_ids)
    agreement={s for s in eligible if mrc[s]==bundled[s]}
    if len(eligible)!=EXPECTED_REF_TRIPLE or len(agreement)!=EXPECTED_REF_TRIPLE:
        raise ValueError("pinned strict triple HM3/EUR LDscore reference changed")
    return mrc,agreement,{
        "original_mrcieu_ref_sha256":meta["ref_sha256"],
        "zenodo_archive_md5":archive["archive_md5"],
        "zenodo_archive_sha256":archive["archive_sha256"],
        "strict_three_reference_rsid_count":len(eligible),
        "direct_unordered_allele_pair_agreement":len(agreement),
        "2018_BMI_candidate_used":False,
        "SNP_selection":"2015 original rsID ∩ MRCIEU HM3 ∩ Zenodo HM3 ∩ canonical EUR LDscore",
    }


def prepare_rows(rows,ref_pairs:dict,allowed:set,writer:CandidateWriter, *,
                 strict=True):
    """
    Reference contains unordered allele pair codes. Source order is retained.
    Known duplicated source IDs are quarantined entirely, not preferentially chosen.
    """
    stats=Counter()
    seen=set()
    dup_ids=set()
    selected_ids=set()
    ns=[]
    original_finished=False
    writer.write(OUTPUT_HEADER)
    for idx,fields in rows:
        if idx==0:
            if fields[0]!="__LOCKE_FULL_SOURCE_SHA256_PASS__":
                raise ValueError("full source sha sentinel differs")
            stats["source_sha_passes"]+=1
            stats["original_compressed_bytes"]=int(fields[1])
            original_finished=True
            break
        stats["original_rows"]+=1
        if len(fields)!=8:raise ValueError("original row shape changed")
        snp,a1,a2,freq,beta,se,p,n=fields
        if not snp:raise ValueError("missing original rsID")
        if snp in seen:
            stats["duplicate_extra_rows"]+=1
            dup_ids.add(snp)
            if snp not in DUPLICATE_SOURCE_RSIDS:
                raise ValueError("unrecognized new duplicated SNP ID")
        else:
            seen.add(snp)
        if snp in DUPLICATE_SOURCE_RSIDS:
            stats["source_duplicate_quarantine_rows"]+=1
            if snp in allowed:
                stats["reference_duplicate_source_rows_discarded"]+=1
            continue
        if snp not in allowed:
            stats["excluded_outside_strict_triple_reference"]+=1
            continue
        try:
            code=allele_pair(a1,a2)
        except ValueError:
            stats["excluded_invalid_source_alleles"]+=1
            continue
        if frozenset((a1,a2)) in PALINDROMES:
            stats["excluded_palindromic"]+=1
            continue
        if code!=ref_pairs[snp]:
            stats["excluded_allele_pair_discordance"]+=1
            continue
        try:
            n_value=float(n)
            b_value=float(beta)
            se_value=float(se)
            p_value=float(p)
            z=b_value/se_value
        except (ValueError,ZeroDivisionError):
            stats["excluded_bad_numeric"]+=1
            continue
        if not all(map(math.isfinite,(n_value,b_value,se_value,p_value,z))) or n_value<=0 or se_value<=0 or not 0<=p_value<=1:
            stats["excluded_bad_numeric"]+=1
            continue
        if snp in selected_ids:
            raise ValueError("unexpected repeated candidate SNP ID")
        selected_ids.add(snp)
        if not n_value.is_integer():stats["selected_fractional_N"]+=1
        if p_value==0.0:stats["selected_original_P_zero"]+=1
        if p_value<=5e-8:stats["selected_original_P_le_5e8"]+=1
        # Harmonic/imputed cohorts can report noninteger per-variant N.
        # Do not change N or use HapMap CEU as if study-wide EAF.
        text=f"{snp}\t{n}\t{format(z,'.15g')}\t{a1}\t{a2}\n".encode("ascii")
        writer.write(text)
        stats["selected_rows"]+=1
        ns.append(n_value)
        if len(selected_ids)>MAX_CANDIDATE_ROWS:
            raise ValueError("candidate size cap exceeded")
    if not original_finished or stats["source_sha_passes"]!=1:
        raise ValueError("original source did not reach full SHA256 validation")
    if strict and (stats["original_rows"]!=EXPECTED_ORIGINAL_ROWS or
       stats["duplicate_extra_rows"]!=EXPECTED_DUPLICATE_EXTRA_ROWS or
       dup_ids!=DUPLICATE_SOURCE_RSIDS):
        raise ValueError("Locke original row count, duplicate set or checksum changed")
    if not stats["selected_rows"]:raise ValueError("independent 2015 candidate empty")
    stats["original_unique_ids"]=len(seen)
    stats["distinct_source_duplicate_ids"]=len(dup_ids)
    ns.sort()
    def quantile(f):return ns[int((len(ns)-1)*f)]
    return {
        "counts":dict(stats),
        "2015_candidate_plaintext_sha256":writer.sha.hexdigest(),
        "2015_candidate_plaintext_bytes":writer.length,
        "2015_candidate_header":OUTPUT_HEADER.decode("ascii").strip().split("\t"),
        "candidate_file_written":writer.out is not None,
        "n_distribution":{"min":ns[0],"p10":quantile(0.1),"median":median(ns),
                          "p90":quantile(0.9),"max":ns[-1]},
        "effect_allele_policy":"LOCKE_2015_A1_EFFECT_HG19_FORWARD; no 2018 effects",
        "significance_policy":"Z_FROM_ROUNDED_BETA_SE_NOT_PUBLISHED_GENOMIC_CONTROL_P",
        "decimal_sample_N_policy":"PRESERVE_RAW_FRACTIONAL_N_NOT_AUTO_LDSC_INPUT",
        "scientific_status":"HOLD","ready_for_ldsc":False,"ldsc_h2_rg_MR_run":False,
    }


def audit(candidate:Path|None=None, *,fetcher=urllib.request.urlopen,
          source_iter=source_rows, reference_loader=load_strict_references,
          strict=True):
    report={
        "schema":"ckm-locke-2015-independent-three-reference-eur-bmi-candidate-v1",
        "checked_at_utc":datetime.now(timezone.utc).isoformat(),
        "technical_status":"HOLD_SOURCE_OR_REFERENCE_OR_MUNGE",
        "2018_exposure_candidates_used":False,
        "source":"Locke 2015 GIANT European BMI",
        "not_an_approved_ldsc_input":True,
        "scientific_status":"HOLD","ready_for_ldsc":False,
        "ldsc_h2_rg_MR_run":False,
        "gwas_participant_independence_established":False
    }
    try:
        if candidate is not None:
            root=Path(__file__).resolve().parents[2]
            target=candidate.expanduser().resolve()
            if root in target.parents and (root/"artifacts").resolve() not in target.parents:
                raise ValueError("candidate must live under ignored artifacts/ within repo")
        ref,selected,provenance=reference_loader(fetcher)
        with CandidateWriter(candidate) as writer:
            result=prepare_rows(source_iter(fetcher),ref,selected,writer,strict=strict)
        report.update(result)
        report["reference_sources"]=provenance
        report["technical_status"]="PASS_2015_ORIGINAL_INDEPENDENT_SNP_PANEL_ONLY"
    except Exception as exc:
        report["technical_status"]="HOLD_SOURCE_OR_REFERENCE_OR_MUNGE"
        report["error_type"]=type(exc).__name__
        report["error_detail"]=str(exc)[:260]
    return report


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out",required=True,type=Path)
    ap.add_argument("--emit-candidate",type=Path,default=None,
                    help="Ignored intermediate .candidate.tsv.gz; not LDScore-ready")
    args=ap.parse_args()
    if args.emit_candidate is not None and not args.emit_candidate.name.endswith(".candidate.tsv.gz"):
        ap.error("--emit-candidate must have .candidate.tsv.gz suffix")
    report=audit(args.emit_candidate)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(report,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({
        "status":report["technical_status"],
        "2015_independent_rows":report.get("counts",{}).get("selected_rows"),
        "plaintext_sha256":report.get("2015_candidate_plaintext_sha256"),
        "exclusions":report.get("counts"),
        "scientific_status":report["scientific_status"],
        "error":report.get("error_detail"),
    },ensure_ascii=False))
    if report["technical_status"]!="PASS_2015_ORIGINAL_INDEPENDENT_SNP_PANEL_ONLY":
        raise SystemExit(1)


if __name__=="__main__":main()
