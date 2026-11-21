#!/usr/bin/env python3
"""Independent 2015-derived BMI candidate vs 22 EUR LDscore CHR:BP and Ensembl GRCh37.

Does not query, import SNP IDs, or select any variants from Yengo 2018 BMI data.
Full source for this stage is SHA-pinned 2015-only intermediate and the canonical
Zenodo EUR LDscore reference. 44 forward reference sequence sentinels are sampled
independently from the 2015 panel. No LDSC h2/rg or MR performed.
"""
from __future__ import annotations
import argparse
from collections import Counter
from datetime import datetime,timezone
import gzip
import hashlib
import json
import math
from pathlib import Path
import tarfile
from urllib.parse import urlsplit
import urllib.request

from scripts.ckm_stroke.audit_eur_ldscore_zenodo import (
    FILE_URL,SIZE,EXPECTED_MD5,DigestReader,safe_member,LD_SCORE_PATTERN,
    M_5_50_PATTERN
)
from scripts.ckm_stroke.audit_bmi_grch37_reference_sentinels import inspect_batch

CANDIDATE_HEADER="SNP\tN_REPORTED\tZ_FROM_BETA_SE\tA1\tA2"
MAX_LINE=8192
MAX_ARCHIVE_MEMBERS=130
MAX_INFLATED_PER_CHR=50_000_000
MAX_CANDIDATE=1_500_000
MAX_ARCHIVE=40_000_000
EXPECTED_CANDIDATE_N=1_070_753
EXPECTED_CANDIDATE_SHA="a1cc5ec06cc81d2bb0a583eb76c0f7e81784cd5dc3f622270de61023f95de4bf"
SALT="locke-2015-EUR-LDscore-Ensembl44-original-only-v1|"
CHR_SET=set(range(1,23))


def read_candidate(path:Path,provenance:Path,*,pin=True):
    previous=json.loads(provenance.read_text(encoding="utf-8"))
    if (previous.get("technical_status")!="PASS_2015_ORIGINAL_INDEPENDENT_SNP_PANEL_ONLY" or
        previous.get("2018_exposure_candidates_used") is not False or
        previous.get("ready_for_ldsc") is not False):
        raise ValueError("2015 independent source candidate not authenticated")
    wanted=previous.get("counts",{}).get("selected_rows")
    if not isinstance(wanted,int) or wanted<100 or wanted>MAX_CANDIDATE:
        raise ValueError("2015 original candidate row count invalid")
    if path.stat().st_size>50_000_000:
        raise ValueError("intermediate candidate unusually large")
    sha=hashlib.sha256()
    candidates={}
    with gzip.open(path,"rb") as f:
        first=f.readline(MAX_LINE)
        if first.rstrip(b"\r\n").decode("ascii")!=CANDIDATE_HEADER:
            raise ValueError("2015 independent candidate header mutated")
        sha.update(first)
        for raw in f:
            if len(raw)>MAX_LINE:raise ValueError("candidate row exceeds cap")
            sha.update(raw)
            row=raw.rstrip(b"\r\n").decode("ascii").split("\t")
            if len(row)!=5:raise ValueError("candidate row field count")
            snp,n,z,a1,a2=row
            if snp in candidates:raise ValueError("duplicate 2015 candidate SNP")
            try:
                fv=float(n); zv=float(z)
            except ValueError as exc:raise ValueError("candidate numeric failure") from exc
            if not(math.isfinite(fv) and fv>0 and math.isfinite(zv)):
                raise ValueError("candidate missing or invalid N/Z")
            if len(a1)!=1 or len(a2)!=1 or a1==a2 or any(k not in "ACGT" for k in (a1,a2)):
                raise ValueError("candidate alleles invalid")
            if frozenset((a1,a2)) in (frozenset(("A","T")),frozenset(("C","G"))):
                raise ValueError("unexpected palindromic 2015 candidate")
            candidates[snp]=(a1,a2)
            if len(candidates)>MAX_CANDIDATE:
                raise ValueError("too many 2015 candidates")
    digest=sha.hexdigest()
    if len(candidates)!=wanted or digest!=previous.get("2015_candidate_plaintext_sha256"):
        raise ValueError("candidate content or count differs from pinned independent provenance")
    if pin and (len(candidates)!=EXPECTED_CANDIDATE_N or digest!=EXPECTED_CANDIDATE_SHA):
        raise ValueError("expected 2015 independent source-derived panel changed")
    return candidates,{"count":len(candidates),"plaintext_sha256":digest,"gzip_bytes":path.stat().st_size}


def update_sentinels(chrom:int,bp:int,snp:str,alleles:tuple,by_chr:dict):
    score=hashlib.sha256((SALT+snp).encode("ascii")).hexdigest()
    v={"chr":chrom,"pos":bp,"rsid":snp,"A1":alleles[0],
       "A2":alleles[1],"rank_sha256":score}
    slot=by_chr[chrom]
    slot.append(v)
    slot.sort(key=lambda x:(x["rank_sha256"],x["rsid"]))
    if len(slot)>2:del slot[2:]


def parse_ldscore_member(fileobj,chrom:int,panel:dict,seen:set,coords:set,by_chr:dict,
                         counts:Counter,*,allow_small=False):
    with gzip.GzipFile(fileobj=fileobj,mode="rb") as gz:
        header=gz.readline(MAX_LINE)
        if not header or len(header)>=MAX_LINE:
            raise ValueError("LDscore gzip header invalid")
        columns=header.decode("utf-8-sig").strip().split()
        if not {"CHR","SNP","BP","L2"}.issubset(columns):
            raise ValueError("LDscore required fields not present")
        chrom_i,snp_i,bp_i=(columns.index(x) for x in ("CHR","SNP","BP"))
        inflated=len(header)
        for raw in gz:
            inflated+=len(raw)
            if inflated>MAX_INFLATED_PER_CHR or len(raw)>MAX_LINE:
                raise ValueError("LDscore row or decompressed cap exceeded")
            row=raw.decode("ascii").split()
            if len(row)!=len(columns) or row[chrom_i]!=str(chrom):
                raise ValueError("LDscore canonical chromosome/schema changed")
            snp=row[snp_i]
            if not row[bp_i].isdecimal():raise ValueError("LDscore genomic BP invalid")
            bp=int(row[bp_i])
            if bp<=0 or bp>300_000_000:raise ValueError("LDscore BP out of GRCh37 autosomal bounds")
            counts["canonical_ldscore_rows"]+=1
            alleles=panel.get(snp)
            if alleles is None:continue
            counts["candidate_ref_id_rows"]+=1
            if snp in seen:
                counts["candidate_duplicate_ldscore_rsid"]+=1
                continue
            seen.add(snp)
            counts["candidate_ref_unique_rsids"]+=1
            coord=(chrom,bp)
            if coord in coords:counts["candidate_duplicate_ref_coordinates"]+=1
            coords.add(coord)
            counts["candidate_ldscore_chr_"+str(chrom)]+=1
            update_sentinels(chrom,bp,snp,alleles,by_chr)
    if not allow_small and counts["canonical_ldscore_rows"]<1000:
        raise ValueError("empty canonical LDscore")


def scan_reference(panel:dict,*,opener=urllib.request.urlopen,sequence_probe=inspect_batch):
    req=urllib.request.Request(FILE_URL,headers={"User-Agent":"CKM-Locke2015-Independent-Positions-EUR/1.0"})
    seen=set()
    coords=set()
    chromosomes=set()
    m_files=set()
    totals=Counter()
    by_chr={k:[] for k in range(1,23)}
    with opener(req,timeout=70) as response:
        u=urlsplit(response.geturl())
        if getattr(response,"status",200)!=200 or u.hostname!="zenodo.org" or u.scheme!="https":
            raise ValueError("reference redirected outside verified Zenodo")
        length=response.headers.get("Content-Length","")
        if not length.isdecimal() or int(length)!=SIZE:
            raise ValueError("reference size differs from Zenodo published")
        reader=DigestReader(response,MAX_ARCHIVE)
        with tarfile.open(fileobj=reader,mode="r|gz") as tf:
            for entry in tf:
                totals["archive_entries"]+=1
                if totals["archive_entries"]>MAX_ARCHIVE_MEMBERS:raise ValueError("too many tar entries")
                safe_member(entry.name)
                if not(entry.isdir() or entry.isfile()):raise ValueError("reference tar special file")
                if entry.isdir():continue
                ld=LD_SCORE_PATTERN.fullmatch(entry.name)
                m=M_5_50_PATTERN.fullmatch(entry.name)
                if ld:
                    ch=int(ld.group(1))
                    if ch in chromosomes:raise ValueError("duplicate canonical LDscore chr")
                    chromosomes.add(ch)
                    f=tf.extractfile(entry)
                    if f is None:raise ValueError("missing canonical LDscore member")
                    parse_ldscore_member(f,ch,panel,seen,coords,by_chr,totals)
                elif m:
                    ch=int(m.group(1))
                    if ch in m_files:raise ValueError("duplicate M_5_50 reference")
                    m_files.add(ch)
                else:totals["ignored_noncanonical_reference_members"]+=1
        while reader.read(1_048_576):pass
        if reader.nbytes!=SIZE or reader.md5.hexdigest()!=EXPECTED_MD5:
            raise ValueError("published Zenodo reference original compressed MD5 mismatch")
        reference_sha=reader.sha256.hexdigest()
    if chromosomes!=CHR_SET or m_files!=CHR_SET:
        raise ValueError("not 22 canonical chr LDscores + M files")
    totals["candidate_missing_from_eur_ldscore"]=len(panel)-len(seen)
    totals["candidate_total"]=len(panel)
    totals["candidate_unique_ref_coordinates"]=len(coords)
    if (totals["candidate_missing_from_eur_ldscore"]!=0 or
        totals.get("candidate_duplicate_ldscore_rsid",0)!=0 or
        totals.get("candidate_duplicate_ref_coordinates",0)!=0):
        raise ValueError("2015 candidate missing, repeated rsID or nonunique ref coordinate")
    if any(len(by_chr[k])!=2 for k in sorted(CHR_SET)):
        raise ValueError("22-chromosome reference sentinel sampling incomplete")
    chosen=[p for ch in range(1,23) for p in by_chr[ch]]
    samples,reference_counts=sequence_probe(chosen)
    good=reference_counts.get("reference_base_in_forward_allele_pair",0)
    return {
        "technical_status":"PASS_2015_INDEPENDENT_CANDIDATE_22_CHR_REFERENCE_AND_44_SENTINELS_ONLY" if good==44 and
            len(samples)==44 and sum(reference_counts.values())==44 else
            "HOLD_GRCH37_FORWARD_REFERENCE_SENTINEL",
        "reference_zenodo_md5":EXPECTED_MD5,
        "reference_zenodo_sha256":reference_sha,
        "canonical_chr_count":len(chromosomes),
        "canonical_M_5_50_count":len(m_files),
        "counts":dict(totals),
        "sentinels":samples,
        "sentinel_orientation_counts":reference_counts,
        "sentinel_count":len(samples),
        "scientific_status":"HOLD","ready_for_ldsc":False,
        "all_2015_candidate_positions_from_independent_LDscore":len(seen)==len(panel),
        "enrolled_cohort_independence_established":False,
    }


def audit(candidate:Path,provenance:Path,opener=urllib.request.urlopen,
          sequence_probe=inspect_batch):
    report={
        "schema":"ckm-locke-2015-independent-ldscore-grch37-reference-v1",
        "checked_utc":datetime.now(timezone.utc).isoformat(),
        "technical_status":"HOLD_REFERENCE_OR_POSITION",
        "2018_BMI_input_used":False,
        "scientific_status":"HOLD","ready_for_ldsc":False,
        "h2_rg_MR_run":False
    }
    try:
        ids,details=read_candidate(candidate,provenance)
        results=scan_reference(ids,opener=opener,sequence_probe=sequence_probe)
        report.update(results)
        report["2015_candidate"]=details
    except Exception as exc:
        report["technical_status"]="HOLD_REFERENCE_OR_POSITION"
        report["error_type"]=type(exc).__name__
        report["error_detail"]=str(exc)[:240]
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--candidate",required=True,type=Path)
    p.add_argument("--provenance",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    a=p.parse_args()
    report=audit(a.candidate,a.provenance)
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(report,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps({"status":report["technical_status"],
                      "ref_counts":report.get("counts"),
                      "sentinel":report.get("sentinel_orientation_counts"),
                      "error":report.get("error_detail"),
                      "scientific":report["scientific_status"]},ensure_ascii=False))
    if report["technical_status"]!="PASS_2015_INDEPENDENT_CANDIDATE_22_CHR_REFERENCE_AND_44_SENTINELS_ONLY":
        raise SystemExit(1)


if __name__=="__main__":main()
