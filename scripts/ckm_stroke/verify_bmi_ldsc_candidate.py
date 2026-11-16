#!/usr/bin/env python3
"""Independently validate derived BMI candidate .sumstats.gz without promoting it."""
import argparse
import gzip
import hashlib
import json
import math
from pathlib import Path

COLUMNS=b"SNP\tN\tZ\tA1\tA2"
PALINDROMES=(frozenset(("A","T")),frozenset(("C","G")))
MAX_CANDIDATE_ROWS=2_000_000
MAX_COMPRESSED_BYTES=60_000_000


def inspect(candidate: Path, provenance: dict):
    if provenance.get("technical_status")!="PASS_CONSERVATIVE_CANDIDATE_CREATED_NOT_LDSC_APPROVED":
        raise ValueError("candidate provenance did not pass technical transform")
    if provenance.get("scientific_status")!="HOLD" or provenance.get("ready_for_ldsc") is not False:
        raise ValueError("candidate cannot be marked as LDSC-approved")
    if candidate.stat().st_size>MAX_COMPRESSED_BYTES:
        raise ValueError("candidate gzip unexpectedly large")
    digest=hashlib.sha256()
    rsids=set()
    count=0
    z_min=None
    z_max=None
    with gzip.open(candidate,"rb") as input:
        for line in input:
            digest.update(line)
            if count==0:
                if line.rstrip(b"\r\n")!=COLUMNS:
                    raise ValueError("candidate columns/order incorrect")
                count+=1
                continue
            if len(line)>4096:
                raise ValueError("unexpected candidate row length")
            parts=line.rstrip(b"\r\n").decode("ascii").split("\t")
            if len(parts)!=5:
                raise ValueError("candidate row field count invalid")
            snp,n,z,a1,a2=parts
            if not snp or snp in rsids:
                raise ValueError("candidate has missing/duplicated rsID")
            rsids.add(snp)
            if not n.isdecimal() or int(n)<=0:
                raise ValueError("non-positive/non-integer candidate N")
            zs=float(z)
            if not math.isfinite(zs):
                raise ValueError("non-finite signed Z")
            if a1 not in "ACGT" or a2 not in "ACGT" or len(a1)!=1 or len(a2)!=1:
                raise ValueError("invalid candidate alleles")
            if a1==a2 or frozenset((a1,a2)) in PALINDROMES:
                raise ValueError("ambiguous or equal candidate alleles")
            z_min=zs if z_min is None else min(z_min,zs)
            z_max=zs if z_max is None else max(z_max,zs)
            count+=1
            if count>MAX_CANDIDATE_ROWS+1:
                raise ValueError("candidate SNP count exceeds cap")
    nrows=count-1
    if nrows<=0 or nrows!=provenance.get("counts",{}).get("candidate_rows"):
        raise ValueError("candidate row count differs from provenance")
    if digest.hexdigest()!=provenance.get("candidate_plaintext_sha256"):
        raise ValueError("candidate plaintext hash differs from provenance")
    return {
        "technical_status":"PASS_LDSC_CANDIDATE_FILE_INTEGRITY_ONLY",
        "candidate_rows_verified":nrows,
        "candidate_plaintext_sha256":digest.hexdigest(),
        "candidate_compressed_bytes":candidate.stat().st_size,
        "candidate_z_range":[z_min,z_max],
        "candidate_unique_rsids":len(rsids),
        "scientific_status":"HOLD",
        "ready_for_ldsc":False,
    }


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--candidate",type=Path,required=True)
    ap.add_argument("--provenance",type=Path,required=True)
    ap.add_argument("--out",type=Path,required=True)
    args=ap.parse_args()
    try:
        evidence=json.loads(args.provenance.read_text(encoding="utf-8"))
        result=inspect(args.candidate,evidence)
    except Exception as exc:
        result={"technical_status":"HOLD_CANDIDATE_FILE_INVALID",
                "error_type":type(exc).__name__,"error_detail":str(exc)[:200],
                "scientific_status":"HOLD","ready_for_ldsc":False}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(result))
    if result["technical_status"]!="PASS_LDSC_CANDIDATE_FILE_INTEGRITY_ONLY":
        raise SystemExit(1)


if __name__=="__main__":
    main()
