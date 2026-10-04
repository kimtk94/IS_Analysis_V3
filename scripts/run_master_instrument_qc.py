#!/usr/bin/env python3
"""Generic MASTER instrument QC for normalized exposure summaries.

Expected normalized columns:
protein_id,gene_symbol,rsid,effect_allele,other_allele,beta,se
Optional:
eaf,p,n,chrom,pos,genome_build,cis_trans

QC:
- valid beta/se and F statistic
- optional P threshold
- optional cis-only filter
- optional MHC exclusion
- optional deduplication within protein by rsID/variant
- palindromic flag (not automatically removed unless requested)

This module does not perform LD clumping. Source-specific independent-signal
sets may pass directly; raw GWAS/QTL inputs must be LD-pruned upstream or in a
separate ancestry-matched LD step.
"""
from __future__ import annotations
import argparse,csv,gzip,math
from pathlib import Path

PAL={frozenset(("A","T")),frozenset(("C","G"))}

def read_tsv(path):
    p=Path(path);op=gzip.open if p.suffix==".gz" else open
    with op(p,"rt",encoding="utf-8",newline="") as f:yield from csv.DictReader(f,delimiter="	")

def fnum(x):
    try:v=float(str(x).strip())
    except Exception:return None
    return v if math.isfinite(v) else None

def chrom_num(x):
    s=str(x or "").strip().lower().replace("chr","")
    try:return int(s)
    except Exception:return None

def is_mhc(chrom,pos):
    c=chrom_num(chrom);p=fnum(pos)
    return c==6 and p is not None and 25_000_000<=p<=34_000_000

def is_pal(a,b):
    a=str(a or "").upper();b=str(b or "").upper()
    return frozenset((a,b)) in PAL

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--input",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    ap.add_argument("--min-fstat",type=float,default=10.0)
    ap.add_argument("--p-threshold",type=float)
    ap.add_argument("--cis-only",action="store_true")
    ap.add_argument("--exclude-mhc",action="store_true")
    ap.add_argument("--drop-palindromic",action="store_true")
    a=ap.parse_args()

    out=[];seen=set();counts={}
    def bump(k):counts[k]=counts.get(k,0)+1

    for r in read_tsv(a.input):
        bx,se=fnum(r.get("beta")),fnum(r.get("se"))
        if bx is None or se in (None,0):
            bump("invalid_beta_se");continue
        f=(bx/se)**2
        if f<a.min_fstat:
            bump("weak_f");continue
        p=fnum(r.get("p"))
        if a.p_threshold is not None and (p is None or p>a.p_threshold):
            bump("p_fail");continue
        if a.cis_only and str(r.get("cis_trans","")).strip().lower() not in {"cis",""}:
            bump("non_cis");continue
        if a.exclude_mhc and is_mhc(r.get("chrom"),r.get("pos")):
            bump("mhc");continue
        pal=is_pal(r.get("effect_allele"),r.get("other_allele"))
        if a.drop_palindromic and pal:
            bump("palindromic");continue
        key=(str(r.get("protein_id","")),str(r.get("rsid","") or r.get("variant_id","")))
        if key in seen:
            bump("duplicate");continue
        seen.add(key)
        x=dict(r)
        x["f_stat"]=f
        x["palindromic"]=int(pal)
        x["instrument_qc"]="ok"
        out.append(x)
        bump("kept")

    fields=list(out[0]) if out else [
      "protein_id","gene_symbol","rsid","variant_id","chrom","pos","effect_allele","other_allele",
      "beta","se","eaf","p","n","ancestry","genome_build","platform","cis_trans","f_stat","palindromic","instrument_qc"
    ]
    a.output.parent.mkdir(parents=True,exist_ok=True)
    op=gzip.open if a.output.suffix==".gz" else open
    with op(a.output,"wt",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter="	",lineterminator="
",extrasaction="ignore")
        w.writeheader();w.writerows(out)
    print("MASTER_INSTRUMENT_QC_PASS "+" ".join(f"{k}={v}" for k,v in sorted(counts.items()))+f" output={a.output}")

if __name__=="__main__":main()
