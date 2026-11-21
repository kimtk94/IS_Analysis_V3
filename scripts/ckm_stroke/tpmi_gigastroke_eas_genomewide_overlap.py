#!/usr/bin/env python3
"""Measured whole-genome TPMI 433.21 -> GIGASTROKE EAS AIS SNP intersection.

Public aggregate QC only, not a GWAS release or inferred independent replication.
Validate full original TPMI SHA256, official UCSC hg38ToHg19 chain MD5,
known Ensembl control locus, full EAS AIS source SHA256, allele orientation,
and exclude palindromic or ambiguous variants. Never approve N_eff, LDSC,
MR or colocalization from this descriptive overlap alone.
"""
from __future__ import annotations

import argparse
from collections import Counter
import gzip
import hashlib
import json
import math
from pathlib import Path
import tempfile
import urllib.parse
import urllib.request

from scripts.ckm_stroke.audit_tpmi_43321_full_gwas_qc import (
    EXPECTED_HEADER, EXPECTED_SHA256, EXPECTED_BYTES, URL as TPMI_URL,
)
from scripts.ckm_stroke.qc_gigastroke_eas_ais_original import (
    URL as AIS_URL, SCHEMA as AIS_HEADER,
)

CHAIN_URL=("https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/"
           "hg38ToHg19.over.chain.gz")
CHAIN_MD5="ff3031d93792f4cbb86af44055efd903"
AIS_SHA256="c521d5695414bdedcc84d2ea21d588e330c75a672f36a87c3fd092b4c55687a3"
TPMI_ROWS=1739048
AIS_ROWS=6725216
TARGET_CHR="19"
TARGET_POS38=15462568
TARGET_POS37=15573379
DNA=set("ACGT")
COMPLEMENT=str.maketrans("ACGT","TGCA")

def download_pinned(url,hostname,path,sha,byte_limit,exact_size=None):
    parsed=urllib.parse.urlsplit(url)
    if parsed.scheme!="https" or parsed.hostname!=hostname:
        raise ValueError("disallowed original source host")
    req=urllib.request.Request(
        url,headers={"User-Agent":"CKM-TPMI43321-EAS-genomewide-overlap/1.0"})
    digest=hashlib.sha256()
    size=0
    with urllib.request.urlopen(req,timeout=180) as response:
        resolved=urllib.parse.urlsplit(response.geturl())
        if (response.status!=200 or resolved.scheme!="https"
                or resolved.hostname!=hostname):
            raise ValueError("original download changed host or status")
        declared=response.headers.get("Content-Length")
        if exact_size is not None and declared and int(declared)!=exact_size:
            raise ValueError("original source declared size drift")
        with path.open("wb") as dest:
            while chunk:=response.read(1048576):
                size+=len(chunk)
                if size>byte_limit:
                    raise ValueError("source exceeded bounded download size")
                digest.update(chunk)
                dest.write(chunk)
    if (exact_size is not None and size!=exact_size) or digest.hexdigest()!=sha:
        raise ValueError("pinned original source size/SHA256 failed")
    return {"sha256":digest.hexdigest(),"bytes":size}

def download_chain(path):
    url=CHAIN_URL
    req=urllib.request.Request(
        url,headers={"User-Agent":"CKM-TPMI43321-UCSC-chain-verified/1.0"})
    h=hashlib.md5(usedforsecurity=False)
    size=0
    with urllib.request.urlopen(req,timeout=90) as resp:
        parsed=urllib.parse.urlsplit(resp.geturl())
        if resp.status!=200 or parsed.hostname!="hgdownload.soe.ucsc.edu":
            raise ValueError("UCSC chain URL/status changed")
        with path.open("wb") as out:
            while chunk:=resp.read(1048576):
                size+=len(chunk)
                if size>10000000: raise ValueError("UCSC chain file too large")
                out.write(chunk);h.update(chunk)
    if h.hexdigest()!=CHAIN_MD5:
        raise ValueError("UCSC published hg38ToHg19 chain MD5 mismatch")
    return {"md5":h.hexdigest(),"bytes":size}

def is_snv(ref,alt,a1):
    return (len(ref)==len(alt)==len(a1)==1
            and ref in DNA and alt in DNA and a1 in (ref,alt)
            and ref!=alt)

def translate_snv(ref,alt,a1,strand):
    if strand not in ("+","-"):
        raise ValueError("invalid chain strand")
    if strand=="-":
        return (ref.translate(COMPLEMENT),alt.translate(COMPLEMENT),
                a1.translate(COMPLEMENT))
    return ref,alt,a1

def allele_align(source_ref,source_alt,source_a1,source_freq,ea,oa,eaf):
    """Align two non-palindromic alleles. Sign is outcome beta for TPMI A1."""
    if not is_snv(source_ref,source_alt,source_a1):
        return "BAD_SOURCE_ALLELES",None,None
    if any(len(a)!=1 or a not in DNA for a in (ea,oa)) or ea==oa:
        return "BAD_OUTCOME_ALLELES",None,None
    if {source_ref,source_alt} in ({"A","T"},{"C","G"}):
        return "PALINDROMIC_UNRESOLVED",None,None
    if {ea,oa}!={source_ref,source_alt}:
        return "ALLELES_MISMATCH",None,None
    sign=1 if ea==source_a1 else -1
    expected_freq=eaf if sign==1 else 1-eaf
    if not (math.isfinite(source_freq) and 0<=source_freq<=1
            and math.isfinite(eaf) and 0<=eaf<=1):
        return "BAD_ALLELE_FREQUENCY",None,None
    return "MATCHED",sign,abs(source_freq-expected_freq)

def association_numeric_qc(or_value,tpmi_se,tpmi_p,tpmi_n,ais_beta,ais_se,ais_p,ais_or):
    """Numerical eligibility for allele-matched pairs, no N_eff or causality claims."""
    try:
        values=(float(or_value),float(tpmi_se),float(tpmi_p),float(tpmi_n),
                float(ais_beta),float(ais_se),float(ais_p),float(ais_or))
    except (TypeError, ValueError, OverflowError):
        return "NON_NUMERIC",None
    odds,s1,p1,n,beta,s2,p2,odds2=values
    if not all(math.isfinite(x) for x in values):
        return "NONFINITE",None
    if not (odds>0 and s1>0 and 0<=p1<=1 and n>0 and
            s2>0 and 0<=p2<=1 and odds2>0):
        return "OUT_OF_RANGE",None
    logor1=math.log(odds)
    # Some providers use SPA/Firth-adjusted p; flag large drift but do not
    # silently revise or discard source P based on a Wald approximation.
    odds_drift=abs(math.log(odds2)-beta)
    return "PASS_FINITE_RANGE",odds_drift

def unique_lift(lifter,chrom,pos):
    matches=lifter.convert_coordinate("chr"+chrom,pos-1)
    if not matches:return "UNMAPPED",None
    if len(matches)!=1:return "MULTIMAPPED",None
    target,zero,strand,_score=matches[0]
    if (target!="chr"+chrom or strand not in ("+","-")
            or not isinstance(zero,(int,float)) or int(zero)!=zero
            or zero<0):
        return "TRANSCHROM_OR_INVALID",None
    return "UNIQUE",(chrom,int(zero)+1,strand)

def tpmi_index(path,lifter):
    mapstat=Counter()
    idx={}
    significant=0
    with gzip.open(path,"rt",encoding="utf-8-sig") as f:
        header=f.readline().rstrip("\r\n").split("\t")
        if header!=EXPECTED_HEADER:raise ValueError("TPMI original header changed")
        rows=0
        for line in f:
            col=line.rstrip("\r\n").split("\t")
            if len(col)!=17:raise ValueError("TPMI source row malformed")
            rows+=1
            chrom,pos,variant,ref,alt,prov,a1,omitted,af,model,obs,odd,se,z,p,err,rsid=col
            if chrom not in {str(i) for i in range(1,23)}:
                mapstat["NON_AUTOSOMAL"]+=1;continue
            try:
                pval=float(p)
                if pval<=5e-8:significant+=1
            except ValueError:raise ValueError("TPMI P-value not numeric")
            if not is_snv(ref,alt,a1):
                mapstat["NON_BIALLELIC_SNV"]+=1;continue
            if model!="ADD" or err not in (".","NA"):
                mapstat["INVALID_TEST_OR_ERROR"]+=1;continue
            try:
                p38=int(pos)
                freq=float(af);or_value=float(odd);se_val=float(se)
                nobs=int(obs)
            except ValueError:
                mapstat["INVALID_NUMERIC"]+=1;continue
            if not(all(math.isfinite(v) for v in (freq,or_value,se_val,pval))
                   and 0<=freq<=1 and or_value>0 and se_val>0
                   and 0<=pval<=1 and nobs>0):
                mapstat["INVALID_NUMERIC"]+=1;continue
            condition,coord=unique_lift(lifter,chrom,p38)
            if condition!="UNIQUE":
                mapstat[condition]+=1;continue
            c37,p37,strand=coord
            r37,a37,eff37=translate_snv(ref,alt,a1,strand)
            if (chrom=="19" and p38==TARGET_POS38
                    and (p37!=TARGET_POS37 or (r37,a37,eff37)!=("T","C","C"))):
                raise ValueError("UCSC liftover disagrees with independently proven TPMI lead")
            if {r37,a37} in ({"A","T"},{"C","G"}):
                mapstat["PALINDROMIC"]+=1;continue
            key=(c37,p37)
            idx.setdefault(key,[]).append((r37,a37,eff37,freq,or_value,se_val,pval,nobs))
            mapstat["UNIQUE_MAPPED_VALID_SNV"]+=1
        if rows!=TPMI_ROWS or significant!=1:
            raise ValueError("TPMI full row count or genome-wide lead drift")
    return idx,dict(mapstat),rows,significant

def scan_eas(path,idx):
    counts=Counter()
    chr_counts=Counter()
    frequency_diff=Counter()
    matched_per_chrom=Counter()
    used_keys=set()
    lead_saw=False
    with gzip.open(path,"rt",encoding="utf-8-sig") as fp:
        header=fp.readline().rstrip("\r\n").split("\t")
        if tuple(header)!=AIS_HEADER:
            raise ValueError("AIS original file column schema drift")
        rows=0
        for line in fp:
            cols=line.rstrip("\r\n").split("\t")
            if len(cols)!=len(AIS_HEADER):
                raise ValueError("AIS original association width drift")
            rows+=1
            chrom,p=cols[0],cols[1]
            chr_counts[chrom]+=1
            if chrom==TARGET_CHR and p==str(TARGET_POS37):lead_saw=True
            if not p.isdecimal():raise ValueError("AIS original position not numeric")
            key=(chrom,int(p))
            options=idx.get(key)
            if not options:continue
            counts["COORDINATE_ROW_INTERSECTIONS"]+=1
            if len(options)!=1:
                counts["TPMI_MULTIALLELIC_COORDINATE_AMBIGUOUS"]+=1;continue
            ref,alt,a1,tpmi_freq,tpmi_or,tpmi_se,tpmi_p,tpmi_n=options[0]
            try:eaf=float(cols[2])
            except ValueError:
                counts["OUTCOME_INVALID_EAF"]+=1;continue
            status,sign,delta=allele_align(ref,alt,a1,tpmi_freq,
                                          cols[9].upper(),cols[10].upper(),eaf)
            counts[status]+=1
            if status!="MATCHED":continue
            if key in used_keys:
                raise ValueError("duplicate matched EAS coordinates")
            used_keys.add(key)
            matched_per_chrom[chrom]+=1
            num_status,or_drift=association_numeric_qc(
                tpmi_or,tpmi_se,tpmi_p,tpmi_n,
                cols[3],cols[4],cols[5],cols[6])
            counts["NUMERIC_"+num_status]+=1
            if or_drift is not None:
                if or_drift>0.05:counts["EAS_LOGOR_BETA_DIFF_GT_0.05"]+=1
                if or_drift>0.005:counts["EAS_LOGOR_BETA_DIFF_GT_0.005"]+=1
            if num_status=="PASS_FINITE_RANGE":
                counts["NUMERIC_ELIGIBLE_MATCHES"]+=1
                if delta>0.20:counts["NUMERIC_AF_DIFF_GT_0.20"]+=1
                if delta>0.10:counts["NUMERIC_AF_DIFF_GT_0.10"]+=1
                if delta>0.05:counts["NUMERIC_AF_DIFF_GT_0.05"]+=1
            if delta>0.20:frequency_diff["ABS_DIFF_GT_0.20"]+=1
            elif delta>0.10:frequency_diff["ABS_DIFF_GT_0.10_LE_0.20"]+=1
            elif delta>0.05:frequency_diff["ABS_DIFF_GT_0.05_LE_0.10"]+=1
            else:frequency_diff["ABS_DIFF_LE_0.05"]+=1
        if rows!=AIS_ROWS:raise ValueError("AIS original full source row count drift")
    return dict(counts),dict(matched_per_chrom),dict(chr_counts),dict(frequency_diff),lead_saw,rows

def execute():
    try:
        from pyliftover import LiftOver
    except ImportError as e:
        raise RuntimeError("pyliftover 0.4.1 not installed") from e
    with tempfile.TemporaryDirectory(prefix="tpmi_eas_aggregate_") as td:
        root=Path(td); tpmi=root/"tpmi.gz";outcome=root/"eas.gz"
        chainzip=root/"hg38ToHg19.over.chain.gz"; chainplain=root/"chain.over.chain"
        src=download_pinned(TPMI_URL,urllib.parse.urlsplit(TPMI_URL).hostname,
                            tpmi,EXPECTED_SHA256,100000000,EXPECTED_BYTES)
        chain=download_chain(chainzip)
        with gzip.open(chainzip,"rb") as zin,chainplain.open("wb") as zout:
            while chunk:=zin.read(1048576):
                zout.write(chunk)
        lifter=LiftOver(str(chainplain))
        control,coordinate=unique_lift(lifter,"19",TARGET_POS38)
        if control!="UNIQUE" or coordinate!=("19",TARGET_POS37,"+"):
            raise ValueError("UCSC chain lead mapping conflicts with independent Ensembl validation")
        idx,mapstat,tpmi_rows,significant=tpmi_index(tpmi,lifter)
        del lifter
        eas=download_pinned(AIS_URL,"ftp.ebi.ac.uk",outcome,
                            AIS_SHA256,1500000000)
        counted,matched_by_chr,outcome_by_chr,freqdiff,lead_saw,eas_rows=scan_eas(outcome,idx)
        result={
            "schema":"TPMI43321-EAS-AIS-genomewide-overlap-aggregate-v1",
            "technical_status":"PASS_ORIGINAL_BUILD_AND_ALLELE_AGGREGATE_QC",
            "TPMI_source_sha256":src["sha256"],
            "TPMI_original_rows":tpmi_rows,
            "TPMI_original_GW_significant_count":significant,
            "TPMI_build":"GRCh38",
            "EAS_AIS_source_sha256":eas["sha256"],
            "EAS_AIS_original_rows":eas_rows,
            "EAS_build":"GRCh37",
            "UCSC_chain_md5_official":chain["md5"],
            "UCSC_chain_source":CHAIN_URL,
            "Ensembl_lead_crosscheck":"PASS_GRCh38_19_15462568_TO_GRCh37_19_15573379",
            "TPMI_liftover_QC_counts":mapstat,
            "mapped_TPMI_unique_positions":len(idx),
            "source_allele_match_QC_counts":counted,
            "EAS_by_chromosome_original_rows":outcome_by_chr,
            "matched_by_chromosome":matched_by_chr,
            "frequency_delta_buckets":freqdiff,
            "numeric_qc_is_necessary_not_sufficient_for_MR_LDSC":True,
            "official_binary_Neff_still_unapproved":True,
            "source_rare_lead_GRCh37_in_EAS_AIS":lead_saw,
            "source_rare_lead_eligible_for_independent_replication":False,
            "cohort_independence_verified":False,
            "TPMI_phenotype_effective_N_approved":False,
            "EAS_AIS_effective_N_approved":False,
            "LDSC_h2_rg_MR_coloc_performed":False,
            "per_variant_original_effects_exported":False,
            "scientific_status":"HOLD",
        }
        if counted.get("MATCHED",0)==0:
            result["technical_status"]="HOLD_NO_MATCHED_NONPALINDROMIC_SNVS"
        return result

def main():
    p=argparse.ArgumentParser();p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    try:result=execute()
    except Exception as e:
        result={"schema":"TPMI43321-EAS-AIS-genomewide-overlap-aggregate-v1",
                "technical_status":"HOLD_SOURCE_LIFTOVER_OR_HARMONIZATION",
                "error_type":type(e).__name__,"error":str(e)[:400],
                "scientific_status":"HOLD"}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n")
    print(json.dumps({"status":result["technical_status"],
                      "TPMI_rows":result.get("TPMI_original_rows"),
                      "EAS_rows":result.get("EAS_AIS_original_rows"),
                      "mapped_positions":result.get("mapped_TPMI_unique_positions"),
                      "overlap":result.get("source_allele_match_QC_counts"),
                      "allele_freq_bins":result.get("frequency_delta_buckets"),
                      "error":result.get("error"),
                      "scientific_status":"HOLD"}))
if __name__=="__main__":main()
