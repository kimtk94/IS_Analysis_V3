#!/usr/bin/env python3
"""Stream GIGASTROKE EAS AIS original sumstats and verify source schema and coverage.

No genome-wide overlap with TPMI is claimed: TPMI source was GRCh38 and must be
lifted, and cohort independence/effective N remain unverified.
"""
from __future__ import annotations
import argparse,gzip,hashlib,json,math,urllib.parse,urllib.request
from pathlib import Path

ACCESSION="GCST90104545"
URL=("https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
     "GCST90104001-GCST90105000/GCST90104545/"
     "GCST90104545_buildGRCh37.tsv.gz")
TARGET_CHROM="19";TARGET_POS="15573379";TARGET_REF="T";TARGET_ALT="C"
TPMI_A1="C";TPMI_LOGOR=0.5215696989267228
SCHEMA=("chromosome","base_pair_location","effect_allele_frequency","beta",
        "standard_error","p_value","odds_ratio","ci_lower","ci_upper",
        "effect_allele","other_allele")
CAP_BYTES=1_500_000_000
class StreamHash:
    def __init__(self,source):
        self.source=source;self.n=0;self.sha=hashlib.sha256()
        self.md5=hashlib.md5(usedforsecurity=False)
    def read(self,size=-1):
        n=1048576 if size<0 else min(size,1048576)
        b=self.source.read(n)
        self.n+=len(b)
        if self.n>CAP_BYTES:raise ValueError("source exceeds compressed size cap")
        self.sha.update(b);self.md5.update(b)
        return b

def verify():
    req=urllib.request.Request(URL,headers={
        "User-Agent":"TPMI-43321-GIGASTROKE-EAS-AIS-full-original-QC/1.0",
        "Accept":"application/octet-stream"})
    with urllib.request.urlopen(req,timeout=180) as response:
        host=urllib.parse.urlsplit(response.geturl())
        if (response.status!=200 or host.scheme!="https"
            or host.hostname!="ftp.ebi.ac.uk"):
            raise ValueError("source provenance has changed")
        expected=response.headers.get("Content-Length")
        if expected is None or not expected.isdecimal():
            raise ValueError("source Content-Length unavailable")
        expected=int(expected)
        if not 1000000<expected<=CAP_BYTES:
            raise ValueError("source compressed size outside permitted range")
        reader=StreamHash(response)
        count=0; chr19=0; exact=[]
        with gzip.GzipFile(fileobj=reader,mode="rb") as fp:
            header=fp.readline().decode("utf-8-sig").rstrip("\r\n").split("\t")
            if tuple(header)!=SCHEMA:
                raise ValueError("source header mismatch "+repr(header[:15]))
            for line in fp:
                if len(line)>8192:raise ValueError("oversized original data record")
                count+=1
                x=line.rstrip(b"\r\n").decode("utf-8").split("\t")
                if len(x)!=len(SCHEMA):
                    raise ValueError("row source schema drift")
                if x[0]!="19":continue
                chr19+=1
                if x[1]!=TARGET_POS:continue
                exact.append(x)
                if len(exact)>1:raise ValueError("duplicate exact genomic coordinate")
            # Validate full compressed stream, including any trailing bytes.
            while reader.read(1048576):pass
        if reader.n!=expected:
            raise ValueError("published compressed byte count and read bytes mismatch")
        result={
           "schema":"GIGASTROKE-EAS-AIS-original-coverage-QC-v1",
           "study_accession":ACCESSION,
           "official_source_url":URL,
           "ancestry":"East Asian",
           "trait":"Any ischemic stroke",
           "source_build":"GRCh37",
           "catalog_cases":19032,
           "catalog_controls":237242,
           "source_compressed_bytes":reader.n,
           "computed_full_source_sha256":reader.sha.hexdigest(),
           "computed_full_source_md5":reader.md5.hexdigest(),
           "publisher_independent_checksum_matched":False,
           "full_gzip_CRC_integrity_pass":True,
           "source_header_verified":True,
           "source_rows_scanned":count,
           "chr19_rows":chr19,
           "tested_target":"19:15573379:T:C",
           "exact_locus_rows":len(exact),
           "technical_status":"PASS_SOURCE_VERIFIED_NO_EXACT_LOCUS" if not exact
                              else "PASS_SOURCE_VERIFIED_EXACT_LOCUS_REVIEW",
           "scientific_status":"HOLD",
           "source_has_approved_effective_sample_size":False,
           "TPMI_cohort_overlap_verified":False,
           "genomewide_TPMI_GIGASTROKE_overlap_analysed":False,
           "pQTL_MR_LDSC_coloc_approved":False}
        if exact:
            x=dict(zip(SCHEMA,exact[0]))
            try:
                ea=x["effect_allele"].upper();oa=x["other_allele"].upper()
                beta=float(x["beta"]);se=float(x["standard_error"])
                p=float(x["p_value"]);eaf=float(x["effect_allele_frequency"])
                if (not all(math.isfinite(v) for v in (beta,se,p,eaf))
                    or se<=0 or not 0<=p<=1 or not 0<=eaf<=1):
                    raise ValueError("invalid association statistics")
                orientation=1 if (ea,oa)==(TARGET_ALT,TARGET_REF) else (
                            -1 if (ea,oa)==(TARGET_REF,TARGET_ALT) else None)
                result["AIS_raw_effect"]={
                    "effect_allele":ea,"other_allele":oa,
                    "beta":beta,"SE":se,"p_value":p,"EAF":eaf}
                result["TPMI_A1_aligned_beta"]=beta*orientation if orientation else None
                result["allele_alignment_status"]=(
                    "PASS_ALLELES_ORIENTED" if orientation
                    else "HOLD_MISMATCHED_ALLELES")
                result["direction_concordant"]=(
                    beta*orientation*TPMI_LOGOR>0 if orientation else None)
                result["direction_is_independent_replication"]=False
            except ValueError as exc:
                result["technical_status"]="HOLD_INVALID_SINGLE_LOCUS_EFFECT"
                result["single_locus_error"]=str(exc)
        return result

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--out",type=Path,required=True)
    a=parser.parse_args()
    try:
        result=verify()
    except Exception as exc:
        result={"schema":"GIGASTROKE-EAS-AIS-original-coverage-QC-v1",
                "study_accession":ACCESSION,
                "technical_status":"HOLD_SOURCE_OR_SCHEMA",
                "error_type":type(exc).__name__,"error":str(exc)[:340],
                "scientific_status":"HOLD"}
    a.out.parent.mkdir(parents=True,exist_ok=True)
    a.out.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n")
    print(json.dumps({
        "status":result["technical_status"],
        "source_rows":result.get("source_rows_scanned"),
        "chr19_rows":result.get("chr19_rows"),
        "exact_target_rows":result.get("exact_locus_rows"),
        "observed_effect":result.get("AIS_raw_effect"),
        "source_sha256":result.get("computed_full_source_sha256"),
        "error":result.get("error"),
        "scientific_status":"HOLD"},ensure_ascii=False))

if __name__=="__main__":main()
