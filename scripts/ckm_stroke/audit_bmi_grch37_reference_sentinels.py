#!/usr/bin/env python3
"""Stratified 44-SNP independent GRCh37 reference-sequence sentinel gate.

Samples TWO fully verified GIANT BMI candidate rsIDs per autosome by stable
SHA256 rank, resolves original GRCh37 chr:BP through pinned GIANT full gzip,
and requests one independent 44-region Ensembl GRCh37 FASTA sequence batch.
This only tests allele strand representation at 44 loci; it cannot certify
1,007,899 outcomes, cohort overlap, LDSC N, h2, rg, or causal effects.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import gzip
import hashlib
import json
from pathlib import Path
from urllib.parse import urlsplit
import urllib.request

from scripts.ckm_stroke.audit_bmi_duplicate_zero_p import iter_original_rows
from scripts.ckm_stroke.verify_bmi_ldsc_candidate import inspect

PINNED_N=1_007_899
PINNED_ROWS=2_336_269
PINNED_CANDIDATE_SHA="ed8f713827cd3281bea5f88a9cb42a3eabfecb1ee67c912ca7e86f151e4399bf"
ENSEMBL_URL="https://grch37.rest.ensembl.org/sequence/region/human"
SAMPLING_SALT="ckm-stroke-grch37-ref-forward-v1|"
N_PER_CHR=2
AUTOSOMES=tuple(range(1,23))
COMP={"A":"T","T":"A","C":"G","G":"C"}
MAX_API_BYTES=128_000


def candidate_alleles(candidate_path:Path,provenance_path:Path):
    evidence=json.loads(provenance_path.read_text(encoding="utf-8"))
    validated=inspect(candidate_path,evidence)
    if (validated["candidate_rows_verified"]!=PINNED_N or
        validated["candidate_plaintext_sha256"]!=PINNED_CANDIDATE_SHA):
        raise ValueError("candidate original pinned digest/size changed")
    panel={}
    with gzip.open(candidate_path,"rt",encoding="ascii") as gz:
        header=next(gz,"").strip().split("\t")
        if header!=["SNP","N","Z","A1","A2"]:
            raise ValueError("candidate column order mismatch")
        for raw in gz:
            snp,n,z,a1,a2=raw.strip().split("\t")
            if snp in panel:
                raise ValueError("nonunique candidate SNP")
            panel[snp]=(a1,a2)
    if len(panel)!=PINNED_N:
        raise ValueError("candidate complete SNP count changed")
    return panel,validated


def pick_sentinels(candidates:dict, source_rows, *, chromosomes=AUTOSOMES, per_chr=N_PER_CHR,
                   strict_original=True):
    if not chromosomes or per_chr<=0:
        raise ValueError("invalid sampling grid")
    by_chr={int(k):[] for k in chromosomes}
    seen=set()
    source_count=0
    full_verified=False
    duplicate_exact_repeats=0
    for ix,fields in source_rows:
        if ix==0:
            full_verified=fields[0]=="__SOURCE_COMPLETE_SHA256_PASS__"
            break
        source_count+=1
        if len(fields)!=10:
            raise ValueError("GIANT original row width changed")
        chrom,bp,snp,a1,a2=fields[:5]
        if snp not in candidates:
            continue
        if (a1,a2)!=candidates[snp]:
            raise ValueError("original BMI allele orientation differs from candidate")
        if not chrom.isdecimal() or not bp.isdecimal():
            raise ValueError("invalid BMI position")
        ch=int(chrom)
        if ch not in by_chr:
            if strict_original:
                raise ValueError("BMI candidate chromosome outside autosomes")
            continue
        pos=int(bp)
        if not 0<pos<=300_000_000:
            raise ValueError("invalid original GRCh37 coordinate")
        if snp in seen:
            duplicate_exact_repeats+=1
            continue
        seen.add(snp)
        key=hashlib.sha256((SAMPLING_SALT+snp).encode("ascii")).digest()
        score=key.hex()  # stable across operating systems
        record={"chr":ch,"pos":pos,"rsid":snp,"A1":a1,"A2":a2,"rank_sha256":score}
        slot=by_chr[ch]
        slot.append(record)
        slot.sort(key=lambda x:(x["rank_sha256"],x["rsid"]))
        if len(slot)>per_chr:
            del slot[per_chr:]
    else:
        raise ValueError("full original source did not reach SHA256 verification sentinel")
    if not full_verified:
        raise ValueError("source compressed SHA256 verification missing")
    if strict_original and (source_count!=PINNED_ROWS or len(seen)!=PINNED_N):
        raise ValueError("BMI full original source/candidate mismatch")
    if any(len(slot)!=per_chr for slot in by_chr.values()):
        raise ValueError("insufficient SNP coverage for one or more autosomes")
    return ([item for ch in sorted(by_chr) for item in by_chr[ch]], {
        "source_original_rows":source_count,
        "matched_candidate_rsids_in_original":len(seen),
        "verified_source_sha256":full_verified,
        "seen_exact_duplicate_snp_records":duplicate_exact_repeats,
        "sampling_rule":"Two smallest SHA256(salt+rsid) per chromosome; no significance bias",
        "sampling_salt":SAMPLING_SALT,
        "chromosomes":len(by_chr),"per_chromosome":per_chr
    })


def classify_reference_base(a1:str,a2:str,ref:str):
    if len(ref)!=1 or ref not in COMP:
        return "unknown_or_nonacgt_reference"
    if ref in (a1,a2):
        return "reference_base_in_forward_allele_pair"
    if COMP[ref] in (a1,a2):
        return "reverse_complement_only"
    return "reference_base_incompatible"


def inspect_batch(sentinels,provider=urllib.request.urlopen, *, url=ENSEMBL_URL):
    if not 1<=len(sentinels)<=50:
        raise ValueError("Ensembl per-batch maximum is 50")
    regions=[f'{s["chr"]}:{s["pos"]}..{s["pos"]}:1' for s in sentinels]
    if len(set(regions))!=len(regions):
        raise ValueError("duplicate sentinel genomic positions")
    request=urllib.request.Request(
        url,
        data=json.dumps({"regions":regions},separators=(",",":")).encode("ascii"),
        headers={
            "Accept":"application/json","Content-Type":"application/json",
            "User-Agent":"CKM-Stroke-GRCh37-forward-strand-sentinel-audit/1.0",
        },
        method="POST"
    )
    with provider(request,timeout=60) as response:
        u=urlsplit(response.geturl())
        if (getattr(response,"status",200)!=200 or u.scheme!="https" or
            u.hostname!="grch37.rest.ensembl.org"):
            raise ValueError("unexpected Ensembl GRCh37 source redirect")
        raw=response.read(MAX_API_BYTES+1)
        if len(raw)>MAX_API_BYTES:
            raise ValueError("unexpectedly large Ensembl batch response")
        results=json.loads(raw)
    if not isinstance(results,list) or len(results)!=len(regions):
        raise ValueError("Ensembl batch response count wrong")
    reference_by_query={}
    for rec in results:
        if not isinstance(rec,dict):
            raise ValueError("Ensembl response element not an object")
        query=rec.get("query")
        identifier=rec.get("id","")
        seq=rec.get("seq","").upper()
        if query not in regions or query in reference_by_query:
            raise ValueError("unexpected or duplicate Ensembl region")
        if not identifier.startswith("chromosome:GRCh37:") or len(seq)!=1:
            raise ValueError("Ensembl response is not GRCh37 1bp sequence")
        reference_by_query[query]=seq
    checked=[]
    counts=Counter()
    for sentinel,query in zip(sentinels,regions):
        ref=reference_by_query[query]
        group=classify_reference_base(sentinel["A1"],sentinel["A2"],ref)
        counts[group]+=1
        checked.append({
            "chr":sentinel["chr"],"pos":sentinel["pos"],"rsid":sentinel["rsid"],
            "A1":sentinel["A1"],"A2":sentinel["A2"],"forward_GRCh37_base":ref,
            "orientation_classification":group,"sampling_sha256":sentinel["rank_sha256"]
        })
    return checked,dict(counts)


def audit(candidate:Path,provenance:Path,*,provider=urllib.request.urlopen,
          source_rows=iter_original_rows):
    out={
        "schema":"ckm-stroke-bmi-ensembl-grch37-stratified-44-sentinels-v1",
        "checked_at_utc":datetime.now(timezone.utc).isoformat(),
        "technical_status":"HOLD_REFERENCE_NOT_VERIFIED",
        "reference":"Ensembl GRCh37 REST genome FASTA, forward strand",
        "reference_url":ENSEMBL_URL,
        "scientific_status":"HOLD","ready_for_ldsc":False,
        "full_bmi_and_outcome_allele_orientation_verified":False,
        "outcome_effect_beta_sign_recalculated":False,
        "h2_rg_MR_performed":False,
        "no_raw_genotypes_persisted":True,
    }
    try:
        panel,checked= candidate_alleles(candidate,provenance)
        sentinels,source_meta=pick_sentinels(panel,source_rows())
        examined,counts=inspect_batch(sentinels,provider)
        good=counts.get("reference_base_in_forward_allele_pair",0)
        out.update({
            "candidate_plaintext_sha256":checked["candidate_plaintext_sha256"],
            "sampling":source_meta,
            "sentinel_count":len(examined),
            "reference_base_categories":counts,
            "sentinels":examined,
            "technical_status":(
                "PASS_STRATIFIED_44_SENTINELS_ONLY"
                if good==44 and not any(v for k,v in counts.items() if k!="reference_base_in_forward_allele_pair")
                else "HOLD_SENTINEL_ALLELE_OR_REFERENCE_DISCORDANCE"),
        })
    except Exception as exc:
        out["technical_status"]="HOLD_REFERENCE_SOURCE_OR_SCHEMA"
        out["error_type"]=type(exc).__name__
        out["error_detail"]=str(exc)[:250]
    return out


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--candidate",required=True,type=Path)
    p.add_argument("--provenance",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    args=p.parse_args()
    result=audit(args.candidate,args.provenance)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({
        "technical_status":result["technical_status"],
        "sentinel_count":result.get("sentinel_count"),
        "orientation_counts":result.get("reference_base_categories"),
        "scientific_status":result["scientific_status"],
    }))
    if result["technical_status"]!="PASS_STRATIFIED_44_SENTINELS_ONLY":
        raise SystemExit(1)


if __name__=="__main__":
    main()
