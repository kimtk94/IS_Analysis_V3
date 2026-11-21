#!/usr/bin/env python3
"""One-locus GIGASTROKE AIS EUR source verification, independent of TPMI summary file.
Original source is streamed once, full publisher MD5 is mandatory. Output only one
public association, never an inferred sample size or independent replication claim.
"""
from __future__ import annotations
import argparse, gzip, hashlib, json, math, urllib.request
from pathlib import Path
from urllib.parse import urlsplit
from scripts.ckm_stroke.prepare_tpmi_43321_gigastroke_lookup import EXPECTED_SHA

URL=("https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/"
     "GCST90104001-GCST90105000/GCST90104540/"
     "GCST90104540_buildGRCh37.tsv.gz")
PUBLISHER_MD5="1bfe4ae8a5042fb24cb0a562b05b0d2f"
SIZE=283_482_604
COLUMNS=("chromosome","base_pair_location","effect_allele_frequency","beta",
         "standard_error","p_value","odds_ratio","ci_lower","ci_upper",
         "effect_allele","other_allele")
COMPLEMENT=str.maketrans("ACGT","TGCA")
def orient(a1,a2,ea,oa,eaf):
    if any(x not in ("A","C","G","T") for x in (a1,a2,ea,oa)) or a1==a2 or ea==oa:
        return "INVALID_ALLELES",None
    if {a1,a2} in ({"A","T"},{"C","G"}):
        return "PALINDROMIC_FREQUENCY_UNRESOLVED",None
    for transformed,alleles in (("same",(ea,oa)),("complement",(ea.translate(COMPLEMENT),oa.translate(COMPLEMENT)))):
        if alleles==(a1,a2):return transformed+"_same",1
        if alleles==(a2,a1):return transformed+"_swapped",-1
    return "ALLELE_MISMATCH",None

def target_gate(t):
    if (t.get("schema")!="ckm-TPMI43321-GIGASTROKE-GRCh37-target-request-v1"
        or t.get("technical_status")!="PASS_VERIFIED_TARGET_ONLY"
        or t.get("target_build")!="GRCh37"
        or t.get("scientific_status")!="HOLD"
        or t.get("GIGASTROKE_original_row_matched") is not False
        or t.get("binary_case_control_Neff_approved") is not False):
        raise ValueError("TPMI GRCh37 verified target is missing")
    ch=t.get("target_chr"); pos=t.get("target_pos")
    ref=t.get("target_REF");alt=t.get("target_ALT");a1=t.get("target_A1")
    if not(str(ch).isdecimal() and 1<=int(ch)<=22 and type(pos) is int and 0<pos<=300000000
        and all(isinstance(x,str) and x in ("A","C","G","T") for x in (ref,alt,a1))
        and ref!=alt and a1 in (ref,alt)
        and t.get("target_variant_key")==f"{ch}:{pos}:{ref}:{alt}"):
        raise ValueError("invalid source target locus or allele")
    tpmi_beta=t.get("TPMI_A1_logOR")
    tpmi_se=t.get("TPMI_A1_logOR_SE")
    tpmi_p=t.get("TPMI_original_P")
    if (not all(type(v) in (float,int) and math.isfinite(v) for v in (tpmi_beta,tpmi_se,tpmi_p))
        or tpmi_se<=0 or not 0<=tpmi_p<=5e-8):
        raise ValueError("missing valid TPMI original A1 logOR, SE or P")
    return str(ch),pos,ref,alt,a1

class HashReader:
    def __init__(self,source):
        self.source=source;self.n=0;self.md5=hashlib.md5(usedforsecurity=False)
    def read(self,size=-1):
        chunk=self.source.read(min(size,1048576) if size>=0 else 1048576)
        self.n+=len(chunk)
        if self.n>SIZE:raise ValueError("original gzip exceeds pinned byte cap")
        self.md5.update(chunk)
        return chunk

def scan(source,target):
    ch,pos,ref,alt,a1=target_gate(target)
    reader=HashReader(source)
    hits=[];nrows=0
    with gzip.GzipFile(fileobj=reader,mode="rb") as gz:
        head=gz.readline().decode("utf-8-sig").rstrip("\r\n").split("\t")
        if tuple(head)!=COLUMNS:raise ValueError("GIGASTROKE original columns changed")
        for row in gz:
            if len(row)>8192:raise ValueError("overlong GWAS source record")
            parts=row.decode("utf-8").rstrip("\r\n").split("\t")
            if len(parts)!=len(COLUMNS):raise ValueError("GIGASTROKE source row schema drift")
            nrows+=1
            if parts[0]!=ch or parts[1]!=str(pos):continue
            hits.append(parts)
            if len(hits)>1:raise ValueError("duplicate AIS outcome locus")
    # gzip may stop at member end; consume trailing bytes for strict full compressed checksum
    while reader.read(1048576):pass
    if reader.n!=SIZE or reader.md5.hexdigest()!=PUBLISHER_MD5:
        raise ValueError("GIGASTROKE original compressed size or published MD5 mismatch")
    out={"schema":"ckm-TPMI43321-GIGASTROKE-AIS-one-locus-v1",
         "original_source_accession":"GCST90104540","original_source_publisher_MD5":PUBLISHER_MD5,
         "original_source_bytes":reader.n,"original_row_count_scanned":nrows,
         "target_variant_key":target["target_variant_key"],
         "source_target_rows":len(hits),"original_source_complete_integrity_pass":True,
         "TPMI_Neff_approved":False,"GIGASTROKE_Neff_approved":False,
         "independent_replication_claimed":False,"scientific_status":"HOLD"}
    if len(hits)==0:
        out["technical_status"]="HOLD_NO_EXACT_GIGASTROKE_COORDINATE";return out
    x=hits[0]
    ea,oa=x[9].upper(),x[10].upper()
    try:
        eaf,beta,se,p,odds=map(float,(x[2],x[3],x[4],x[5],x[6]))
    except ValueError:raise ValueError("non-numeric original GIGASTROKE association")
    if (not all(map(math.isfinite,(eaf,beta,se,p,odds))) or
        not 0<=eaf<=1 or se<=0 or not 0<=p<=1 or odds<=0):
        raise ValueError("invalid outcome effect or frequency")
    orientation,sign=orient(a1,alt if a1==ref else ref,ea,oa,eaf)
    signed=beta*sign if sign is not None else None
    tpmi_beta=target["TPMI_A1_logOR"]
    same_direction=(signed*tpmi_beta>0) if signed is not None and tpmi_beta!=0 and signed!=0 else None
    out.update({"TPMI_original_A1_logOR":tpmi_beta,
                "TPMI_original_A1_SE":target["TPMI_A1_logOR_SE"],
                "observed_effect_direction_concordant":same_direction,
                "effect_direction_concordance_is_not_independent_replication":True,
                "GIGASTROKE_effect_allele":ea,"GIGASTROKE_other_allele":oa,
                "GIGASTROKE_raw_beta":beta,"GIGASTROKE_raw_SE":se,
                "GIGASTROKE_raw_P":p,"GIGASTROKE_raw_EAF":eaf,
                "effect_allele_orientation":orientation,
                "signed_beta_to_TPMI_A1":signed,
                "technical_status":"PASS_EXACT_LOCUS_AND_ALLELE_MATCH" if sign is not None else "HOLD_ALLELE_ORIENTATION"})
    return out

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--target",required=True,type=Path)
    parser.add_argument("--out",required=True,type=Path)
    args=parser.parse_args()
    try:
        target=json.loads(args.target.read_text("utf-8"))
        target_gate(target)
        req=urllib.request.Request(URL,headers={"User-Agent":"CKM-TPMI43321-one-locus-validated-research/1.0"})
        with urllib.request.urlopen(req,timeout=120) as response:
            u=urlsplit(response.geturl())
            if (response.status!=200 or u.scheme!="https" or u.hostname!="ftp.ebi.ac.uk"
                or response.headers.get("Content-Length")!=str(SIZE)):
                raise ValueError("GIGASTROKE official source URL or length mismatch")
            result=scan(response,target)
    except Exception as exc:
        result={"technical_status":"HOLD_GIGASTROKE_SOURCE_OR_TARGET",
                "error_type":type(exc).__name__,"error_detail":str(exc)[:260],
                "scientific_status":"HOLD"}
    args.out.parent.mkdir(parents=True,exist_ok=True)
    args.out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
    print(json.dumps({"status":result["technical_status"],"locus":result.get("target_variant_key"),
                      "orientation":result.get("effect_allele_orientation"),"science":"HOLD"}))
    if result["technical_status"]!="PASS_EXACT_LOCUS_AND_ALLELE_MATCH":raise SystemExit(1)

if __name__=="__main__":main()
