#!/usr/bin/env python3
"""Read-only IS GWAS source matching; writes aggregate diagnostics only.

Accepts a locally curated, provenance-recorded GRCh37 ALT-normalized official
TSV with columns variant_GRCh37, beta_ALT, p_value. Never trusts labels,
study names or N values without independent audit and source-file hashes.
"""
import argparse
import collections
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path

ARCHIVE_BETA={"BBJ":"BBJ_beta_ALT","EAS_AIS":"EAS_AIS_beta_ALT"}
ARCHIVE_P={"BBJ":"BBJ_p","EAS_AIS":"EAS_AIS_p"}
MAX_OFFICIAL_ROWS=30000000

def open_table(path):
    return gzip.open(path,"rt",encoding="utf-8-sig",newline="") if path.suffix==".gz" else path.open(encoding="utf-8-sig",newline="")

def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda:f.read(1024*1024),b""):h.update(chunk)
    return h.hexdigest()

def variant(s):
    import re
    m=re.fullmatch(r"(?:chr)?([1-9]|1[0-9]|2[0-2]|X|Y|MT):([1-9][0-9]*):([ACGT]+):([ACGT]+)",str(s or "").strip(),re.I)
    if not m:return None
    chrom,pos,ref,alt=m.groups()
    chrom=chrom.upper()
    ref,alt=ref.upper(),alt.upper()
    if ref==alt:return None
    return chrom+":"+pos+":"+ref+":"+alt

def number(x,p=False):
    try:
        v=float(x)
        if not math.isfinite(v) or (p and not (0<=v<=1)):return None
        return v
    except (ValueError,TypeError):return None

def same_p(a,b):
    if a==0 or b==0:return a==b
    return abs(math.log10(a)-math.log10(b))<=0.05

def match(archive,official,source,*,max_rows=MAX_OFFICIAL_ROWS):
    if source not in ARCHIVE_BETA:raise ValueError("Unsupported source")
    old={}
    counts=collections.Counter()
    with open_table(archive) as f:
        reader=csv.DictReader(f,delimiter="\t")
        needed={"mapped_genomic_variant_GRCh37",ARCHIVE_BETA[source],ARCHIVE_P[source]}
        if not needed.issubset(reader.fieldnames or []):raise ValueError("Archive schema missing")
        for row in reader:
            v=variant(row["mapped_genomic_variant_GRCh37"])
            beta=number(row[ARCHIVE_BETA[source]])
            p=number(row[ARCHIVE_P[source]],p=True)
            if v is None or beta is None or p is None:continue
            if v in old and old[v]!=(beta,p):
                old[v]=None
                counts["archive_conflicting_duplicate_rows"]+=1
            elif v not in old:old[v]=(beta,p)
    seen=set()
    with open_table(official) as f:
        reader=csv.DictReader(f,delimiter="\t")
        if not {"variant_GRCh37","beta_ALT","p_value"}.issubset(reader.fieldnames or []):
            raise ValueError("Official file must be independently normalized to GRCh37 ALT")
        for i,row in enumerate(reader):
            if i>=max_rows:raise ValueError("Bounded official scan exceeded")
            counts["official_rows"]+=1
            v=variant(row["variant_GRCh37"])
            if v is None or v not in old:continue
            if v in seen:
                counts["duplicate_official_variant"]+=1
                continue
            seen.add(v)
            if old[v] is None:
                counts["archive_ambiguous_variant"]+=1
                continue
            beta=number(row["beta_ALT"])
            p=number(row["p_value"],p=True)
            if beta is None or p is None:
                counts["official_invalid_effect_or_p"]+=1
                continue
            a,b=old[v]
            match_p=same_p(p,b)
            if math.isclose(beta,a,rel_tol=0.015,abs_tol=0.00005) and match_p:
                counts["same_ALT_beta_and_p"]+=1
            elif math.isclose(beta,-a,rel_tol=0.015,abs_tol=0.00005) and match_p:
                counts["opposite_beta_orientation"]+=1
            else:
                counts["discordant_beta_or_p"]+=1
    n=counts["same_ALT_beta_and_p"]+counts["opposite_beta_orientation"]+counts["discordant_beta_or_p"]
    f=counts["same_ALT_beta_and_p"]/n if n else 0
    # This audit is intentionally not authorized to adopt cohort/sample metadata.
    return {
      "schema":"is-p3a-independent-source-match-v1",
      "source":source,"archive_sha256":sha(archive),"official_sha256":sha(official),
      "counts":dict(counts),"concordance_fraction":round(f,6),
      "source_statistics_screen":"PASS_REVIEW_REQUIRED" if n>=20 and counts["same_ALT_beta_and_p"]>=20 and f>=0.95 and counts["duplicate_official_variant"]==0 and counts["archive_ambiguous_variant"]==0 else "BLOCKED",
      "approved_source":False,"sample_size_adopted":False,"coloc_executed":False,
      "limitations":"Requires external review of exact source release, variant coordinates, allele normalization and phenotype provenance. This result alone never authorizes sample size adoption."
    }

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--archive",type=Path,required=True)
    p.add_argument("--official-normalized",type=Path,required=True)
    p.add_argument("--source",choices=sorted(ARCHIVE_BETA),required=True)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    result=match(args.archive,args.official_normalized,args.source)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n")
    print(json.dumps({k:v for k,v in result.items() if k not in ("archive_sha256","official_sha256")}))
if __name__=="__main__":main()
