#!/usr/bin/env python3
"""Official UCSC hg19/GRCh37 -> hg38/GRCh38 chr10 coordinate and REF screen.

Scientific HOLD: 5 uniquely-mapped positions and actual reference alleles do
NOT prove that every full-cis QTL/GWAS SNP can be safely converted or joined.
No use of these mapped bounds for Coloc eligibility until SNP-by-SNP checking.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import math
import tempfile
import urllib.request
from collections import Counter
from pathlib import Path

BASE = "https://hgdownload.soe.ucsc.edu/goldenPath"
FILES = {
    "chain": {
        "url":f"{BASE}/hg19/liftOver/hg19ToHg38.over.chain.gz",
        "md5":"35887f73de5e2231656504d1f6430900",
        "limit":2_000_000,
    },
    "hg19": {
        "url":f"{BASE}/hg19/chromosomes/chr10.fa.gz",
        "md5":"508a1e86cc0f64257dd0c9775462eebf",
        "limit":56_000_000,
    },
    "hg38": {
        "url":f"{BASE}/hg38/chromosomes/chr10.fa.gz",
        "md5":"b7d478cb06a84400f74a8932cd5fad0d",
        "limit":56_000_000,
    },
}
REGION = {
    "id":"IS_XDATA_G0015",
    "chromosome":"chr10",
    "GRCh37_start_1based":104663676,
    "GRCh37_end_1based":105720351,
    "context_GRCh37_start_1based":104163676,
    "context_GRCh37_end_1based":106220351,
    "historical_lead_GRCh37_1based":105596641,
    "historical_lead_allele_pair":["G","C"],
}
ANCHORS = {
    "context_start":104163676,
    "cis_region_start":104663676,
    "cis_region_midpoint":(104663676+105720351)//2,
    "palindromic_lead":105596641,
    "cis_region_end":105720351,
    "context_end":106220351,
}
COMP = str.maketrans("ACGTN","TGCAN")


def require(ok:bool,msg:str)->None:
    if not ok:raise ValueError(msg)


def fetch_checked(name:str,output:Path,open_url=urllib.request.urlopen)->dict:
    info=FILES[name]
    req=urllib.request.Request(info["url"],headers={"User-Agent":"IS-P3A-Masters-UCSC-hg19-hg38-source-validation/1.0"})
    md5=hashlib.md5(usedforsecurity=False)
    sha=hashlib.sha256()
    n=0
    with open_url(req,timeout=95) as response:
        require(getattr(response,"status",200)==200,"UCSC source must return HTTP 200")
        with output.open("wb") as f:
            while True:
                chunk=response.read(1<<20)
                if not chunk:break
                n+=len(chunk)
                require(n<=info["limit"],"oversized UCSC source: "+name)
                md5.update(chunk)
                sha.update(chunk)
                f.write(chunk)
    require(n>100000 and md5.hexdigest()==info["md5"],
            "original UCSC MD5 mismatch/missing file: "+name+
            " bytes="+str(n)+" observed_md5="+md5.hexdigest()+
            " expected_md5="+info["md5"])
    return {"url":info["url"],"bytes":n,"published_UCSC_MD5":info["md5"],
            "observed_SHA256":sha.hexdigest(),"verified":True}


def chromosome_fasta_sequence(file:Path)->bytes:
    """Uncompress only a single human chr10 FASTA; no patient data."""
    seq=bytearray()
    with gzip.open(file,"rb") as inp:
        header=inp.readline().strip().split()
        require(header and header[0]==b">chr10","wrong UCSC chromosome FASTA")
        for line in inp:
            if line.startswith(b">"):
                raise ValueError("unexpected second chromosome sequence")
            s=line.strip().upper()
            require(set(s).issubset(set(b"ACGTN")),"non-DNA FASTA bases")
            seq.extend(s)
            require(len(seq)<160_000_000,"unbounded chromosome")
    require(len(seq)>125_000_000,"incomplete chr10 FASTA")
    return bytes(seq)


def single_coordinate(lifter,chrom:str,pos1:int)->dict:
    out=lifter.convert_coordinate(chrom,pos1-1) or []
    if len(out)!=1:
        return {"source_position_1based":pos1,"status":"UNMAPPED_OR_AMBIGUOUS",
                "mapping_count":len(out)}
    target,target0,strand,score=out[0]
    if not isinstance(target0,(float,int)) or not math.isfinite(target0):
        return {"source_position_1based":pos1,"status":"INVALID_COORDINATE"}
    return {"source_position_1based":pos1,"target_chromosome":str(target),
            "target_position_1based":int(target0)+1,"strand":str(strand),
            "alignment_score":int(score),
            "mapping_count":1,
            "status":"UNIQUELY_MAPPED"}


def anchors_audit(lifter,hg19:bytes,hg38:bytes)->dict:
    records={}
    for label,pos in ANCHORS.items():
        a=single_coordinate(lifter,"chr10",pos)
        if a["status"]=="UNIQUELY_MAPPED":
            mapped=a["target_position_1based"]
            if a["target_chromosome"]!="chr10" or a["strand"] not in ("+","-") \
               or not (1<=pos<=len(hg19)) or not (1<=mapped<=len(hg38)):
                a["status"]="WRONG_CHROM_OR_OUT_OF_BOUNDS"
            else:
                old=chr(hg19[pos-1])
                new=chr(hg38[mapped-1])
                expected=old if a["strand"]=="+" else old.translate(COMP)
                a.update({"hg19_REF":old,"hg38_REF":new,
                          "expected_hg38_REF_after_strand":expected,
                          "reference_alleles_concordant":new==expected and new in "ACGT"})
                if not a["reference_alleles_concordant"]:
                    a["status"]="REFERENCE_MISMATCH"
                if label=="palindromic_lead":
                    a["archived_allele_pair_ref_member"]=old in REGION["historical_lead_allele_pair"]
                    if not a["archived_allele_pair_ref_member"]:
                        a["status"]="LEAD_ARCHIVE_ALLELE_PAIR_VS_FASTA_MISMATCH"
        records[label]=a
    positions=list(ANCHORS)
    points=[records[k] for k in positions]
    unique=all(x["status"]=="UNIQUELY_MAPPED" for x in points)
    same_strand=len({x.get("strand") for x in points})==1
    plus_or_minus=points[0].get("strand") if unique else None
    coords=[x.get("target_position_1based") for x in points]
    ordering=unique and (
        (plus_or_minus=="+" and coords==sorted(coords)) or
        (plus_or_minus=="-" and coords==sorted(coords,reverse=True))
    )
    boundary=(
       records["cis_region_start"].get("target_position_1based"),
       records["cis_region_end"].get("target_position_1based"))
    context=(
       records["context_start"].get("target_position_1based"),
       records["context_end"].get("target_position_1based"))
    approved=bool(unique and same_strand and ordering and plus_or_minus=="+" and
                  all(x.get("reference_alleles_concordant") for x in points) and
                  records["palindromic_lead"].get("archived_allele_pair_ref_member") is True)
    return {
        "anchor_status":"PASS_SIX_ANCHORS_AND_TWO_FASTAS" if approved else
                        "HOLD_COORDINATE_OR_REFERENCE_DISCORDANCE",
        "six_positions_only_not_every_variant":True,
        "mapped_anchors":records,
        "all_unique_hg38_chr10_plus_monotonic_REF_concordant":approved,
        "candidate_query_region_GRCh38_1based":
            {"chrom":"10","start":min(context),"end":max(context),
             "is_lower_bound_of_full_cis_with_guaranteed_coverage":False}
                if approved else None,
        "mapped_cis_endpoints_GRCh38_1based":
            list(boundary) if approved else None,
        "all_variant_liftover_verified":False,
        "full_cis_region_equivalence_verified":False,
        "liftover_is_not_snp_harmonization":True,
    }


def audit_results(manifest:dict,sources:dict,reference_sha:dict,
                  anchor_report:dict)->dict:
    groups=manifest.get("groups",[])
    require(manifest.get("schema")=="is-p3a-36-real-group-keys-from-legacy-393-v1"
            and manifest.get("origin",{}).get("commit_sha")==
                "0bf646fd6f29322d17da878253f37def24c01c74"
            and manifest.get("scientific_status")=="HOLD"
            and len(groups)==36
            and len({(x["region_id"],x["QTL_dataset"],x["gene_id"]) for x in groups})==36
            and len({x["gene_id"] for x in groups})==20
            and Counter(x["QTL_dataset"] for x in groups)==
                {"QTD000131":10,"QTD000171":11,"QTD000434":15}
            and all(x["region_id"]=="IS_XDATA_G0015" and
                    x["approved_for_new_coloc"] is False for x in groups),
            "36 real QTL groups changed")
    return {
       "schema":"is-p3a-chr10-official-ucsc-grch37-to-grch38-anchor-ref-audit-v1",
       "technical_status":("PASS_SOURCE_MD5_AND_ANCHOR_REF_VALIDATION" if
                  anchor_report["all_unique_hg38_chr10_plus_monotonic_REF_concordant"]
                   else "HOLD_ALIGNMENT_OR_REFERENCE_INTEGRITY"),
       "scientific_status":"HOLD",
       "source_provenance":sources,
       "chr10_uncompressed_fasta_base_SHA256":reference_sha,
       "declared_GRCh37_region":REGION,
       "liftover":anchor_report,
       "n_real_group_keys":36,
       "n_distinct_genes":20,
       "formal_coloc_abf_eligible":0,"formal_coloc_abf_executed":0,
       "formal_coloc_susie_eligible":0,"formal_coloc_susie_executed":0,
       "groups":[
           {"region_id":x["region_id"],"QTL_dataset":x["QTL_dataset"],
            "gene_id":x["gene_id"],"gene_symbol":x["gene_symbol"],
            "candidate_region_38_locator":anchor_report.get("candidate_query_region_GRCh38_1based"),
            "per_SNP_source_QTL_GWAS_build_and_REF_checked":False,
            "effective_gene_QTL_N_approved":False,
            "GWAS_case_control_N_approved":False,
            "coloc_abf_eligible":False,
            "status":"HOLD_FULL_CIS_VARIANT_LIFTOVER_SIGNED_ALLELES_SAMPLE_N"}
           for x in groups],
    }


def main()->int:
    parser=argparse.ArgumentParser()
    parser.add_argument("--manifest",type=Path,required=True)
    parser.add_argument("--out",type=Path,required=True)
    args=parser.parse_args()
    try:
        manifest=json.loads(args.manifest.read_text(encoding="utf8"))
        with tempfile.TemporaryDirectory(prefix="IS_UCSC_liftover_") as temp:
            folder=Path(temp)
            provenance={}
            for name in FILES:
                provenance[name]=fetch_checked(name,folder/(name+".gz"))
            chain=folder/"hg19ToHg38.over.chain"
            with gzip.open(folder/"chain.gz","rb") as inp,chain.open("wb") as out:
                while (chunk:=inp.read(1<<20)):out.write(chunk)
            from pyliftover import LiftOver
            converter=LiftOver(str(chain))
            hg19=chromosome_fasta_sequence(folder/"hg19.gz")
            hg38=chromosome_fasta_sequence(folder/"hg38.gz")
            seqhash={"GRCh37_chr10":hashlib.sha256(hg19).hexdigest(),
                     "GRCh38_chr10":hashlib.sha256(hg38).hexdigest()}
            anchor_report=anchors_audit(converter,hg19,hg38)
            result=audit_results(manifest,provenance,seqhash,anchor_report)
    except (OSError,ValueError,TypeError,KeyError,OverflowError,ImportError) as e:
        result={"schema":"is-p3a-chr10-official-ucsc-grch37-to-grch38-anchor-ref-audit-v1",
                "technical_status":"HOLD_REFERENCE_NOT_CHECKED",
                "scientific_status":"HOLD","formal_coloc_abf_eligible":0,
                "formal_coloc_abf_executed":0,"error":str(e)[:300]}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n",encoding="utf8")
    print(json.dumps({k:result.get(k) for k in (
          "technical_status","scientific_status","n_real_group_keys",
          "formal_coloc_abf_executed","error")},ensure_ascii=False))
    if result.get("liftover"):
        print(json.dumps({"mapped_GRCh38":result["liftover"].get(
             "candidate_query_region_GRCh38_1based"),
             "anchor_status":result["liftover"]["anchor_status"]}))
    return 0 if result["technical_status"].startswith("PASS_") else 1


if __name__=="__main__":
    raise SystemExit(main())
