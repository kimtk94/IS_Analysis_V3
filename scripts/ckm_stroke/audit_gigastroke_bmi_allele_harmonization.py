#!/usr/bin/env python3
"""GIGASTROKE EUR original GWAS raw MD5 + BMI GRCh37 allele harmonization audit.

Streams each official original gzip ONCE, hashing publisher compressed bytes while
checking every row and matching the 1,007,899 validated GIANT BMI candidate loci.
No original GWAS or matched genetic effect data are persisted, only aggregate JSON.
A valid source MD5 and counted allele categories are NEVER scientific LDSC GO.
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import gzip
import hashlib
import json
import math
from pathlib import Path
from urllib.parse import urlsplit
import urllib.request

from scripts.ckm_stroke.audit_bmi_duplicate_zero_p import iter_original_rows
from scripts.ckm_stroke.verify_bmi_ldsc_candidate import inspect
from scripts.ckm_stroke.verify_gigastroke_raw import (
    MAX_BYTES, get_text, source_urls, load_manifest, check_sidecars,
)

COLUMNS=(
    "chromosome","base_pair_location","effect_allele_frequency","beta",
    "standard_error","p_value","odds_ratio","ci_lower","ci_upper",
    "effect_allele","other_allele"
)
EXPECTED_CANDIDATES=1_007_899
EXPECTED_GIANT_ROWS=2_336_269
MAX_ROWS=35_000_000
MAX_GZIP_LINE=8192
MAX_ERROR_EXAMPLES=6
EXPECTED_GZIP_SIZES={
    "AIS":283_482_604, "LAA":230_156_887,
    "CE":263_057_516, "SVD":237_118_943
}
NUCLEOTIDES={"A","C","G","T"}
COMPLEMENT={"A":"T","T":"A","C":"G","G":"C"}
GOOD={"same","swapped","complement_same","complement_swapped"}


class HashCompressedReader:
    def __init__(self,source,cap=MAX_BYTES):
        self.source=source
        self.cap=cap
        self.read_bytes=0
        self.md5=hashlib.md5(usedforsecurity=False)
        self.sha256=hashlib.sha256()

    def read(self,n=-1):
        n=1_048_576 if n is None or n < 0 else min(n,1_048_576)
        data=self.source.read(n)
        self.read_bytes+=len(data)
        if self.read_bytes>self.cap:
            raise ValueError("original compressed bytes exceeded cap")
        self.md5.update(data)
        self.sha256.update(data)
        return data


def build_bmi_coordinate_panel(candidate:Path,provenance:Path,
                               original=iter_original_rows):
    prior=json.loads(provenance.read_text(encoding="utf-8"))
    verified=inspect(candidate,prior)
    if verified["candidate_rows_verified"]!=EXPECTED_CANDIDATES:
        raise ValueError("unexpected pinned BMI candidate count")
    known={}
    with gzip.open(candidate,"rt",encoding="ascii",newline="") as handle:
        if handle.readline().strip()!="SNP\tN\tZ\tA1\tA2":
            raise ValueError("candidate header changed")
        for line in handle:
            snp,n,z,a1,a2=line.strip().split("\t")
            if snp in known:
                raise ValueError("duplicate BMI SNP")
            known[snp]=(a1,a2)
    if len(known)!=EXPECTED_CANDIDATES:
        raise ValueError("candidate count changed")
    positions={}
    stats=Counter()
    for lineno,row in original():
        if lineno==0:
            stats["verified_source_passes"]+=1
            stats["GIANT_compressed_bytes"]=int(row[1])
            break
        stats["GIANT_source_rows"]+=1
        snp=row[2]
        if snp not in known:
            continue
        chrom,pos=int(row[0]),int(row[1])
        if (row[3],row[4])!=known[snp]:
            raise ValueError("candidate allele orientation inconsistent with GIANT original")
        key=(chrom,pos)
        if key in positions:
            if positions[key]!=(snp,*known[snp]):
                raise ValueError("ambiguous GIANT BMI coordinates with distinct rsIDs")
            stats["repeated_BMI_source_coordinate"]+=1
        else:
            positions[key]=(snp,*known[snp])
            stats["BMI_candidates_mapped"]+=1
    else:
        raise ValueError("original GIANT source verification sentinel absent")
    if stats["verified_source_passes"]!=1 or stats["GIANT_source_rows"]!=EXPECTED_GIANT_ROWS:
        raise ValueError("GIANT original source size changed")
    if len(positions)!=EXPECTED_CANDIDATES:
        raise ValueError("missing BMI candidate coordinates")
    return positions,{
        "candidate_plaintext_sha256":verified["candidate_plaintext_sha256"],
        "candidate_total":len(positions),
        "GIANT_original_rows":stats["GIANT_source_rows"],
        "GIANT_compressed_bytes":stats["GIANT_compressed_bytes"],
    }


def orientation(bmi_a1:str,bmi_a2:str,out_ea:str,out_oa:str)->str:
    """Classify effect orientation relative to BMI Tested_Allele=A1, no guessing."""
    alleles=(bmi_a1,bmi_a2,out_ea,out_oa)
    if any(len(a)!=1 or a not in NUCLEOTIDES for a in alleles):
        return "non_ACGT_alleles"
    if bmi_a1==bmi_a2 or out_ea==out_oa:
        return "same_allele_invalid"
    if frozenset((bmi_a1,bmi_a2)) in (frozenset(("A","T")),frozenset(("C","G"))):
        return "palindromic_ambiguous"
    if (out_ea,out_oa)==(bmi_a1,bmi_a2):
        return "same"
    if (out_ea,out_oa)==(bmi_a2,bmi_a1):
        return "swapped"
    swapped=(COMPLEMENT[out_ea],COMPLEMENT[out_oa])
    if swapped==(bmi_a1,bmi_a2):
        return "complement_same"
    if swapped==(bmi_a2,bmi_a1):
        return "complement_swapped"
    return "allele_incompatible"


def aligned_effect(beta:float,se:float,label:str)->tuple[float,float]:
    """Translate an outcome log-odds effect to BMI Tested_Allele direction."""
    if label not in GOOD or not math.isfinite(beta) or not math.isfinite(se) or se<=0:
        raise ValueError("cannot align this outcome effect")
    oriented=(-beta if label in ("swapped","complement_swapped") else beta)
    return oriented,oriented/se


def examine_outcome_rows(rows,panel,*,max_rows=MAX_ROWS):
    """Streaming full rows, bounded error summaries. Does not retain SNP effects."""
    c=Counter()
    by_coord={}
    duplicate_coords=set()
    samples={}
    header=next(rows,None)
    if header is None or tuple(header.rstrip("\r\n").split("\t"))!=COLUMNS:
        raise ValueError("outcome original header changed")
    for number,line in enumerate(rows,start=1):
        if number>max_rows:
            raise ValueError("outcome row cap exceeded")
        c["original_outcome_rows"]+=1
        if len(line)>MAX_GZIP_LINE:
            raise ValueError("outcome row length cap exceeded")
        fields=line.rstrip("\r\n").split("\t")
        if len(fields)!=len(COLUMNS):
            c["invalid_row_width"]+=1
            continue
        chrom,bp=fields[0],fields[1]
        if not chrom.isdecimal() or not bp.isdecimal():
            c["invalid_or_nonautosomal_chrpos"]+=1
            continue
        chrnum,pos=int(chrom),int(bp)
        if not(1<=chrnum<=22 and 0<pos<=300_000_000):
            c["nonautosomal_or_invalid_chrpos"]+=1
            continue
        c["autosomal_rows_with_coordinates"]+=1
        coord=(chrnum,pos)
        bmi=panel.get(coord)
        if bmi is None:
            c["outcome_rows_absent_from_BMI_candidate"]+=1
            continue
        c["outcome_rows_with_BMI_coordinate"]+=1
        c["candidate_hits_chr"+str(chrnum)]+=1
        label=orientation(bmi[1],bmi[2],fields[9].upper(),fields[10].upper())
        if label in GOOD:
            try:
                beta=float(fields[3]);se=float(fields[4])
                eaf=float(fields[2]);pv=float(fields[5])
                if not all(map(math.isfinite,(beta,se,eaf,pv))) or se<=0 or not 0<=eaf<=1 or not 0<=pv<=1:
                    raise ValueError("matched numerical source invalid")
                signed_beta,signed_z=aligned_effect(beta,se,label)
            except ValueError:
                label="matched_numeric_invalid"
        if coord in duplicate_coords:
            c["additional_collision_outcome_rows"]+=1
            continue
        if coord in by_coord:
            c["outcome_multiple_rows_same_BMI_coordinate"]+=1
            duplicate_coords.add(coord)
            previous=by_coord.pop(coord)
            c["unique_"+previous]-=1
            continue
        by_coord[coord]=label
        c["unique_"+label]+=1
        if label in GOOD and len(samples.setdefault("directional_examples",[]))<MAX_ERROR_EXAMPLES:
            samples["directional_examples"].append({
                "chr":chrnum,"pos":pos,"rsid":bmi[0],
                "orientation":label,"raw_outcome_beta":beta,
                "beta_aligned_to_BMI_Tested_Allele":signed_beta,
                "z_aligned_to_BMI_Tested_Allele":signed_z,
                "outcome_eaf_raw":eaf,
                "transformed_eaf_not_claimed":True,
            })
        if label not in GOOD and len(samples.setdefault(label,[]))<MAX_ERROR_EXAMPLES:
            samples[label].append({"chr":chrnum,"pos":pos,"rsid":bmi[0]})
    c["candidate_coordinates_not_in_outcome"]=len(panel)-len(by_coord)-len(duplicate_coords)
    c["outcome_collision_coordinates_excluded"]=len(duplicate_coords)
    c["unique_outcome_coordinate_matches"]=len(by_coord)
    c["unique_harmonizable_BMI_positions"]=sum(c["unique_"+label] for label in GOOD)
    c["unique_orientation_sign_positive"]=c["unique_same"]+c["unique_complement_same"]
    c["unique_orientation_sign_negative"]=c["unique_swapped"]+c["unique_complement_swapped"]
    if c["candidate_coordinates_not_in_outcome"]<0:
        raise ValueError("invalid candidate missingness count")
    return {"counts":dict(c),"exclusion_examples":samples}


def scan_full(row,panel,*,fetcher=urllib.request.urlopen):
    raw,meta,sums=source_urls(row)
    check_sidecars(row,get_text(meta,fetcher),get_text(sums,fetcher))
    request=urllib.request.Request(raw,headers={"User-Agent":"CKM-Stroke-IS-EUR-harmonization/1.0"})
    with fetcher(request,timeout=90) as response:
        dest=urlsplit(response.geturl())
        if response.status!=200 or dest.scheme!="https" or dest.hostname!="ftp.ebi.ac.uk":
            raise ValueError("original publisher source redirected/status invalid")
        head_length=response.headers.get("Content-Length")
        expected_size=EXPECTED_GZIP_SIZES[row["phenotype"]]
        if not head_length or not head_length.isdecimal() or int(head_length)!=expected_size:
            raise ValueError("publisher raw compressed Content-Length differs pinned original")
        counter=HashCompressedReader(response,MAX_BYTES)
        with gzip.GzipFile(fileobj=counter,mode="rb") as gz:
            def lines():
                while True:
                    rawline=gz.readline(MAX_GZIP_LINE+1)
                    if not rawline:break
                    if len(rawline)>MAX_GZIP_LINE:
                        raise ValueError("outcome row longer than safe bound")
                    yield rawline.decode("utf-8-sig","strict")
            matched=examine_outcome_rows(lines(),panel)
        # Count/MD5 full original stream, not only the first requested prefix.
        while counter.read(1_048_576):
            pass
        if counter.read_bytes!=expected_size:
            raise ValueError("publisher raw stream incomplete")
        expected_md5=row["published_raw_md5"].lower()
        if counter.md5.hexdigest()!=expected_md5:
            raise ValueError("full source gzip MD5 differs from official publisher sidecar")
        if matched["counts"].get("invalid_row_width",0):
            raise ValueError("source row column count invalid")
        if matched["counts"].get("unique_harmonizable_BMI_positions",0)==0:
            raise ValueError("no uniquely harmonizable BMI positions")
        return {
            "phenotype":row["phenotype"],"accession":row["accession"],
            "trait":row["reported_trait"],"build":"GRCh37","ancestry":"EUR",
            "original_compressed_bytes":counter.read_bytes,
            "raw_md5_verified":True,
            "published_md5":expected_md5,
            "source_sha256":counter.sha256.hexdigest(),
            "publisher_yaml_md5_sidecar_verified":True,
            "source_header_columns":list(COLUMNS),
            "outcome_has_rsid_column":False,
            "outcome_per_variant_N_available":False,
            "publisher_total_N":int(row["sample_total"]),
            "directional_alignment":"Allele direction classified; signed beta/Z checked for examples only; no per-SNP harmonized effects exported",
            "technical_status":"PASS_RAW_MD5_AND_ALLELE_CLASSIFICATION_ONLY",
            **matched,
            "scientific_status":"HOLD","ready_for_ldsc":False,
        }


def audit(candidate,provenance,manifest,*,phenotype="ALL",fetcher=urllib.request.urlopen):
    report={
        "schema":"ckm-stroke-gigastroke-eur-full-md5-bmi-harmonization-v1",
        "checked_utc":datetime.now(timezone.utc).isoformat(),
        "source_archives_retained":False,
        "candidate_upload_to_GitHub":False,
        "LDSC_h2_or_rg_executed":False,
        "scientific_status":"HOLD",
        "ready_for_ldsc":False,
        "outcomes":[],
    }
    records=load_manifest(manifest)
    try:
        panel,meta=build_bmi_coordinate_panel(candidate,provenance)
        report["BMI_reference_panel"]=meta
    except Exception as exc:
        report["technical_status"]="HOLD_BMI_PANEL_PREPARATION"
        report["error_type"]=type(exc).__name__
        report["error_detail"]=str(exc)[:240]
        return report
    for row in records:
        if phenotype!="ALL" and row["phenotype"]!=phenotype:
            continue
        try:
            report["outcomes"].append(scan_full(row,panel,fetcher=fetcher))
        except Exception as exc:
            report["outcomes"].append({
                "phenotype":row["phenotype"],"accession":row["accession"],
                "technical_status":"HOLD_ORIGINAL_SOURCE_OR_HARMONIZATION",
                "error_type":type(exc).__name__,
                "error_detail":str(exc)[:240],
                "raw_md5_verified":False,
                "ready_for_ldsc":False,
                "scientific_status":"HOLD",
            })
    report["technical_status"]=(
        "PASS_ALL_REQUESTED_ORIGINAL_MD5_AND_ALLELE_COUNTS"
        if report["outcomes"] and all(
            x["technical_status"]=="PASS_RAW_MD5_AND_ALLELE_CLASSIFICATION_ONLY"
            for x in report["outcomes"])
        else "HOLD_ONE_OR_MORE_ORIGINAL_OUTCOMES"
    )
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--candidate",type=Path,required=True)
    p.add_argument("--provenance",type=Path,required=True)
    p.add_argument("--manifest",type=Path,default=Path("config/ckm_stroke/gigastroke_eur_official.tsv"))
    p.add_argument("--phenotype",choices=("ALL","AIS","LAA","CE","SVD"),default="ALL")
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    result=audit(args.candidate,args.provenance,args.manifest,phenotype=args.phenotype)
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({
        "technical_status":result.get("technical_status"),
        "BMI_candidate":result.get("BMI_reference_panel",{}).get("candidate_total"),
        "outcomes":[{"phenotype":r["phenotype"],"status":r["technical_status"],
                     "md5_verified":r.get("raw_md5_verified",False),
                     "harmonizable":r.get("counts",{}).get("unique_harmonizable_BMI_positions"),
                     "error":r.get("error_detail")} for r in result["outcomes"]],
        "scientific_status":result["scientific_status"]},ensure_ascii=False))
    if result.get("technical_status")!="PASS_ALL_REQUESTED_ORIGINAL_MD5_AND_ALLELE_COUNTS":
        raise SystemExit(1)


if __name__=="__main__":
    main()
