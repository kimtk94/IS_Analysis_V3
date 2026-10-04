#!/usr/bin/env python3
"""Normalize exposure summary statistics into MASTER exposure schema.

This is source-agnostic and uses an explicit column-map TSV:
canonical_column    source_column

Required canonical columns:
protein_id,gene_symbol,rsid,effect_allele,other_allele,beta,se
Optional:
eaf,p,n,chrom,pos,variant_id,ancestry,genome_build,platform,cis_trans

If a canonical field is supplied as a constant CLI option, it can be omitted
from the source table/map.
"""
from __future__ import annotations
import argparse,csv,gzip,bz2,lzma,math
from pathlib import Path

REQ={"protein_id","gene_symbol","rsid","effect_allele","other_allele","beta","se"}
OUT=["protein_id","gene_symbol","rsid","variant_id","chrom","pos",
     "effect_allele","other_allele","beta","se","eaf","p","n",
     "ancestry","genome_build","platform","cis_trans"]

def open_text(path):
    p=Path(path);n=p.name.lower()
    if n.endswith(".gz"):return gzip.open(p,"rt",encoding="utf-8",errors="replace",newline="")
    if n.endswith(".bz2"):return bz2.open(p,"rt",encoding="utf-8",errors="replace",newline="")
    if n.endswith(".xz"):return lzma.open(p,"rt",encoding="utf-8",errors="replace",newline="")
    return p.open("r",encoding="utf-8",errors="replace",newline="")

def read_tsv(path):
    with open_text(path) as f:yield from csv.DictReader(f,delimiter="	")

def fmap(path):
    out={}
    with Path(path).open("r",encoding="utf-8",newline="") as f:
        for r in csv.DictReader(f,delimiter="	"):
            c=str(r.get("canonical_column","")).strip()
            s=str(r.get("source_column","")).strip()
            if c and s:out[c]=s
    return out

def fnum(x):
    try:v=float(str(x).strip())
    except Exception:return None
    return v if math.isfinite(v) else None

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--input",type=Path,required=True)
    ap.add_argument("--column-map",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    ap.add_argument("--ancestry",default="")
    ap.add_argument("--genome-build",default="")
    ap.add_argument("--platform",default="")
    ap.add_argument("--cis-only",action="store_true")
    a=ap.parse_args()

    mp=fmap(a.column_map)
    constants={"ancestry":a.ancestry,"genome_build":a.genome_build,"platform":a.platform}
    missing=[x for x in REQ if x not in mp]
    if missing:raise SystemExit("column map missing required canonical fields: "+",".join(sorted(missing)))

    rows=[];dropped=0
    for r in read_tsv(a.input):
        x={k:(r.get(mp[k],"") if k in mp else constants.get(k,"")) for k in OUT}
        x["protein_id"]=str(x["protein_id"]).strip()
        x["gene_symbol"]=str(x["gene_symbol"]).strip().upper()
        x["rsid"]=str(x["rsid"]).strip()
        x["effect_allele"]=str(x["effect_allele"]).strip().upper()
        x["other_allele"]=str(x["other_allele"]).strip().upper()
        for k in ("beta","se","eaf","p","n"):
            v=fnum(x[k]);x[k]="" if v is None else v
        if x["beta"]=="" or x["se"]=="" or x["se"]<=0:
            dropped+=1;continue
        if a.cis_only and str(x.get("cis_trans","")).strip().lower() not in {"cis",""}:
            dropped+=1;continue
        rows.append(x)

    a.output.parent.mkdir(parents=True,exist_ok=True)
    op=gzip.open if a.output.suffix==".gz" else open
    with op(a.output,"wt",encoding="utf-8",newline="") as f:
        w=csv.DictWriter(f,fieldnames=OUT,delimiter="	",lineterminator="
")
        w.writeheader();w.writerows(rows)
    print(f"MASTER_EXPOSURE_NORMALIZE_PASS rows={len(rows)} dropped={dropped} output={a.output}")

if __name__=="__main__":main()
