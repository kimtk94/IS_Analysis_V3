#!/usr/bin/env python3
"""CKD manuscript Figure 3 from a SHA-pinned *existing* nine-gene evidence matrix.

No inference, imputation or live data access. Outputs an SVG figure plus a
machine-readable TSV with source SHA256. EUR ABF vs EUR SuSiE only; does not
misrepresent EAS dual-LD as a directly comparable same-population estimate.
"""
import argparse
import csv
import hashlib
import math
from pathlib import Path

GENES=("ACP1","CPVL","F12","GSTA1","GSTA3","HLA-E","INHBC","SDCCAG8","UMOD")
FIELDS=("gene_symbol","eur_abf_h4","eur_susie_max_h4","eur_stage2_class","eur_comparison_status")
def parse_h4(raw,name):
    if str(raw or "").strip().upper() in ("","NA","N/A","NULL","NONE","NAN","MISSING"):
        return None
    try: x=float(raw)
    except (ValueError,TypeError): raise ValueError(name+" invalid probability")
    if not math.isfinite(x) or not 0<=x<=1: raise ValueError(name+" outside [0,1]")
    return x
def load(path):
    with path.open(newline="",encoding="utf-8-sig") as f:
        reader=csv.DictReader(f,delimiter="\t")
        if not reader.fieldnames or not set(FIELDS).issubset(reader.fieldnames):
            raise ValueError("Missing required evidence columns")
        rows=list(reader)
    names=[x.get("gene_symbol","").strip().upper() for x in rows]
    if len(names)!=9 or set(names)!=set(GENES) or len(set(names))!=9:
        raise ValueError("Expected exactly 9 unique frozen CKD genes")
    return {x["gene_symbol"].strip().upper():x for x in rows}
def create(source,out,expected_sha):
    if not source.is_file(): raise ValueError("Missing original CKD evidence matrix")
    source_sha=hashlib.sha256(source.read_bytes()).hexdigest()
    if expected_sha and source_sha.lower()!=expected_sha.lower():
        raise ValueError("Source SHA256 mismatch - refusing to plot")
    if out.exists() and any(out.iterdir()):
        raise ValueError("Output must be new/empty; refusing overwrite")
    data=load(source)
    prepared=[]
    for gene in GENES:
        r=data[gene]
        abf=parse_h4(r.get("eur_abf_h4"),gene+" ABF")
        susie=parse_h4(r.get("eur_susie_max_h4"),gene+" SuSiE")
        status=" ".join(str(r.get(k,"")) for k in ("eur_stage2_class","eur_comparison_status")).upper()
        flag=("UNRESOLVED" if "UNRESOLVED" in status or "NONCONVERGED" in status or "FAILED" in status else
              "DISCORDANT" if "DISCORDANT" in status or "NOT_CONFIRMED" in status else
              "SOURCE_STATUS_OTHER")
        prepared.append({"gene_symbol":gene,"eur_abf_pph4":abf,"eur_susie_pph4":susie,
                         "eur_signal_status":flag,"source_sha256":source_sha})
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(10.8,6.5))
    ys=list(range(9))
    a=[(row["eur_abf_pph4"],i-.13) for i,row in enumerate(prepared) if row["eur_abf_pph4"] is not None]
    b=[(row["eur_susie_pph4"],i+.13) for i,row in enumerate(prepared) if row["eur_susie_pph4"] is not None]
    if a:ax.scatter([x for x,_ in a],[y for _,y in a],marker="o",s=65,label="EUR coloc.abf")
    if b:ax.scatter([x for x,_ in b],[y for _,y in b],marker="D",s=55,label="EUR coloc.susie")
    for i,row in enumerate(prepared):
        if row["eur_susie_pph4"] is None:
            ax.annotate("NA / unresolved",(0.025,i+.13),fontsize=9,va="center")
    ax.axvline(.8,linestyle="--",linewidth=1,label="Interpretation reference PP.H4 = 0.8")
    ax.set_xlim(-.03,1.05)
    ax.set_yticks(ys,list(GENES));ax.invert_yaxis()
    ax.set_xlabel("PP.H4 (posterior probability of shared variant hypothesis)")
    ax.set_title("CKD: European GWAS–plasma pQTL regional colocalization",loc="left")
    ax.grid(axis="x",alpha=.25)
    ax.legend(loc="lower right")
    fig.text(.10,.025,"Source SHA256: "+source_sha[:18]+"...  |  Summary outputs only; not experimental causality.",fontsize=8)
    fig.subplots_adjust(left=.16,bottom=.14,right=.98,top=.90)
    out.mkdir(parents=True,exist_ok=True)
    fig.savefig(out/"CKD_FIG3_EUR_ABF_VS_SUSIE.svg")
    plt.close(fig)
    with (out/"CKD_FIG3_EVIDENCE.tsv").open("w",newline="",encoding="utf-8") as f:
        writer=csv.DictWriter(f,fieldnames=list(prepared[0]),delimiter="\t");writer.writeheader()
        for row in prepared:
            writer.writerow({**row,"eur_abf_pph4":row["eur_abf_pph4"] if row["eur_abf_pph4"] is not None else "NA",
                             "eur_susie_pph4":row["eur_susie_pph4"] if row["eur_susie_pph4"] is not None else "NA"})
    return {"genes":len(prepared),"source_sha256":source_sha,"abf_nonmissing":len(a),"susie_nonmissing":len(b)}
if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--input",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    p.add_argument("--sha256",required=True)
    a=p.parse_args()
    print(create(a.input,a.out,a.sha256))
