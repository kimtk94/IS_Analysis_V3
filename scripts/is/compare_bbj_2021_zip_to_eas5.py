#!/usr/bin/env python3
"""Read-only BBJ 2021 source-vs-archive comparison, streaming ZIP/GZIP.
No extraction, no sample-size adoption, no colocalization, no source edits.
"""
import argparse,csv,gzip,io,json,math,re,zipfile
from collections import Counter
from pathlib import Path

REQUIRED_SOURCE={"CHR","POS","Allele1","Allele2","BETA","SE","p.value"}
ARCHIVE_KEYS={"mapped_genomic_variant_GRCh37","BBJ_beta_ALT","BBJ_p",
              "BBJ_allele_exact_Fasta37"}
SNP=re.compile(r"^(?:chr)?([0-9]{1,2}|X|Y|MT):([1-9][0-9]*):([ACGT]+):([ACGT]+)$",re.I)

def parsed_variant(s):
    m=SNP.fullmatch(str(s).strip())
    if not m:return None
    chrom,pos,ref,alt=m.groups()
    chrom=chrom.upper()
    if chrom.isdigit():chrom=str(int(chrom))
    ref,alt=ref.upper(),alt.upper()
    if ref==alt:return None
    return (chrom,int(pos),ref,alt)

def numeric(s,is_p=False):
    try:
        f=float(s)
        if math.isfinite(f) and (not is_p or 0<=f<=1):return f
    except (TypeError,ValueError,OverflowError):pass
    return None

def close_p(x,y):
    if x is None or y is None:return False
    if x==0 or y==0:return x==y
    return abs(math.log10(x)-math.log10(y))<=0.05

def close_beta(x,y):
    return x is not None and y is not None and math.isclose(x,y,rel_tol=.015,abs_tol=.00005)

def audit(archive:Path, zip_path:Path, *, max_rows=30_000_000):
    candidates={}
    counts=Counter()
    with gzip.open(archive,"rt",encoding="utf-8-sig",newline="") as f:
        r=csv.DictReader(f,delimiter="\t")
        if not ARCHIVE_KEYS.issubset(r.fieldnames or []):
            raise ValueError("Archive required fields absent")
        for row in r:
            if row["BBJ_allele_exact_Fasta37"].lower() not in ("true","1","yes"):
                continue
            key=parsed_variant(row["mapped_genomic_variant_GRCh37"])
            b=numeric(row["BBJ_beta_ALT"]);p=numeric(row["BBJ_p"],True)
            if not key or b is None or p is None:continue
            if key in candidates and candidates[key]!=(b,p):
                candidates[key]=None
                counts["archive_ambiguous_duplicate_rows"]+=1
            elif key not in candidates:candidates[key]=(b,p)
    counts["archive_candidate_variant_ids"]=len(candidates)
    seen=set()
    with zipfile.ZipFile(zip_path) as z:
        auto=[n for n in z.namelist() if n.endswith(".auto.txt.gz")]
        if len(auto)!=1:raise ValueError("Expected single autosomal BBJ GWAS member")
        with z.open(auto[0]) as packed:
            with gzip.GzipFile(fileobj=packed) as raw:
                with io.TextIOWrapper(raw,encoding="utf-8-sig",newline="") as handle:
                    reader=csv.DictReader(handle,delimiter="\t")
                    if not REQUIRED_SOURCE.issubset(reader.fieldnames or []):
                        raise ValueError("Source GWAS columns absent")
                    for row in reader:
                        counts["official_source_rows"]+=1
                        if counts["official_source_rows"]>max_rows:
                            raise ValueError("Max source scan bound exceeded")
                        chrom=row["CHR"].strip().upper().removeprefix("CHR")
                        if chrom.isdigit():chrom=str(int(chrom))
                        try:pos=int(row["POS"])
                        except (ValueError,TypeError):continue
                        a1=row["Allele1"].strip().upper()
                        a2=row["Allele2"].strip().upper()
                        if len(a1)!=1 or len(a2)!=1 or a1 not in "ACGT" or a2 not in "ACGT" or a1==a2:continue
                        # Two possible orientations; *not* proof of the true effect allele.
                        direct=(chrom,pos,a1,a2)
                        flipped=(chrom,pos,a2,a1)
                        if direct in candidates:
                            key=direct;orientation="A1_REF_A2_ALT"
                        elif flipped in candidates:
                            key=flipped;orientation="A2_REF_A1_ALT"
                        else:continue
                        if key in seen:
                            counts["duplicate_source_overlap_rows"]+=1
                            continue
                        seen.add(key)
                        counts["same_coordinate_allele_pair"]+=1
                        expected=candidates[key]
                        if expected is None:
                            counts["ambiguous_archived_beta_p"]+=1
                            continue
                        old_beta,old_p=expected
                        b=numeric(row["BETA"]);p=numeric(row["p.value"],True)
                        if b is None or p is None:
                            counts["invalid_source_beta_p"]+=1
                            continue
                        counts["comparable_effect_p"]+=1
                        if close_p(old_p,p):counts["same_p_with_log_tolerance"]+=1
                        else:counts["different_p"]+=1
                        if close_beta(old_beta,b):counts["same_beta_sign_and_magnitude"]+=1
                        if close_beta(old_beta,-b):counts["opposite_beta_sign_same_magnitude"]+=1
                        if close_p(old_p,p) and (close_beta(old_beta,b) or close_beta(old_beta,-b)):
                            counts["p_and_beta_match_either_sign"]+=1
                        if orientation=="A1_REF_A2_ALT":
                            if close_beta(old_beta,b) and close_p(old_p,p):
                                counts["hypothesis_effect_A2_agree"]+=1
                            if close_beta(old_beta,-b) and close_p(old_p,p):
                                counts["hypothesis_effect_A1_agree"]+=1
                        else:
                            if close_beta(old_beta,-b) and close_p(old_p,p):
                                counts["hypothesis_effect_A2_agree"]+=1
                            if close_beta(old_beta,b) and close_p(old_p,p):
                                counts["hypothesis_effect_A1_agree"]+=1
    comparable=counts["comparable_effect_p"]
    either=counts["p_and_beta_match_either_sign"]
    return {
        "schema":"is-bbj-2021-source-effect-comparison-v1",
        "summary":dict(counts),
        "p_beta_concordance_fraction_either_sign":round(either/comparable,6) if comparable else None,
        "study_release": "Sakaue2021_BBJ_IS_candidate",
        "genomic_build_verified_independently":False,
        "allele_effect_orientation_verified_independently":False,
        "gwas_source_accepted":False,
        "sample_size_adopted":False,
        "coloc_executed":False,
        "interpretation":"Observed source-stat agreement only; incompatible release or allele semantics must be resolved before any GWAS metadata adoption."
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--archive",type=Path,required=True)
    p.add_argument("--zip",type=Path,required=True)
    a=p.parse_args()
    print("BBJ_SOURCE_COMPARISON_JSON="+json.dumps(audit(a.archive,a.zip),separators=(",",":"),sort_keys=True))
if __name__=="__main__":main()
