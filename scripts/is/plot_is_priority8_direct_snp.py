#!/usr/bin/env python3
"""Generate exact-data SVG for eight SNP-level, not summary-reweighted, coloc curves."""
import argparse
import csv
import html
import math
import xml.etree.ElementTree as ET
from pathlib import Path

GENES=("FGF5","CALHM2","NEURL1","C4orf22","INA","SH3PXD2A","COL4A2","COL4A1")
PALETTE=("#2b789c","#349276","#ac752b","#805fa2","#ba5079","#7c8733","#be5845","#63768e")
P12=(1e-6,3e-6,1e-5,3e-5,1e-4)

def run(source,out):
    if out.exists(): raise FileExistsError("NO_OVERWRITE")
    with source.open(newline="",encoding="utf8") as handle:
        lines=list(csv.DictReader(handle,delimiter="\t"))
    if len(lines)!=120: raise ValueError("SOURCE_MUST_BE_120_REPLAYS")
    curves={}
    for gene in GENES:
        curve=[x for x in lines if x["gene"]==gene and x["qtl_n_mode"]=="FIRST"]
        if len(curve)!=5 or {float(r["p12"]) for r in curve}!=set(P12):
            raise ValueError("MISSING_GENE_P12_"+gene)
        curves[gene]={float(x["p12"]):float(x["PP_H4"]) for x in curve}
    W,H=1120,690;x0,y0,x1,y1=118,148,780,570
    X=lambda p:x0+(math.log10(p)+6)/2*(x1-x0)
    Y=lambda h:y1-h*(y1-y0)
    s=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W} {H}" width="{W}" height="{H}">',
       '<rect width="100%" height="100%" fill="#ffffff"/>',
       '<text x="34" y="43" font-family="Arial,sans-serif" font-size="24" font-weight="bold" fill="#1f3545">Ischemic stroke: direct SNP-level coloc prior sensitivity</text>',
       '<text x="34" y="74" font-family="Arial,sans-serif" font-size="14" fill="#516a79">Eight selected GTEx gene–tissue pairs • coloc 5.2.3 • direct beta/SE/MAF SNP rerun</text>',
       '<text x="34" y="98" font-family="Arial,sans-serif" font-size="14" fill="#a04734">Original 646-test screen, not expanded 80-locus discovery or causal gene proof</text>']
    for v in [0,.2,.4,.6,.8,1]:
        y=Y(v)
        s.extend([f'<line x1="{x0}" x2="{x1}" y1="{y:.1f}" y2="{y:.1f}" stroke="#e2e8ed"/>',
                  f'<text x="{x0-10}" y="{y+4:.1f}" text-anchor="end" font-size="12" fill="#475968">{v:.1f}</text>'])
    for p,label in [(1e-6,"10^-6"),(3e-6,"3e-6"),(1e-5,"10^-5"),(3e-5,"3e-5"),(1e-4,"10^-4")]:
        x=X(p)
        s.extend([f'<line x1="{x:.1f}" x2="{x:.1f}" y1="{y0}" y2="{y1}" stroke="{"#9bbcc2" if p==1e-5 else "#eef0f3"}"/>',
                  f'<text x="{x:.1f}" y="{y1+26}" text-anchor="middle" font-size="13" fill="#34495e">{label}</text>'])
    for i,gene in enumerate(GENES):
        color=PALETTE[i];curve=curves[gene]
        positions=[(X(p),Y(curve[p])) for p in P12]
        s.append('<polyline points="'+' '.join(f'{x:.2f},{y:.2f}' for x,y in positions)+f'" fill="none" stroke="{color}" stroke-width="3"/>')
        for x,y in positions:s.append(f'<circle cx="{x:.2f}" cy="{y:.2f}" r="4.2" fill="#fff" stroke="{color}" stroke-width="2.5"/>')
        y=167+48*i
        s.extend([f'<line x1="821" x2="847" y1="{y}" y2="{y}" stroke="{color}" stroke-width="3.5"/>',
                  f'<text x="859" y="{y+5}" font-family="Arial,sans-serif" font-size="16" fill="#2e4657">{html.escape(gene)} ({curve[1e-5]:.3f})</text>'])
    s+=['<text x="445" y="652" text-anchor="middle" font-size="14" fill="#465e6b">p12 prior (logarithmic axis, original baseline p12=1e-5)</text>',
        '<text x="28" y="357" text-anchor="middle" transform="rotate(-90 28 357)" font-size="15" fill="#465e6b">PP.H4</text>',
        '</svg>']
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text("\n".join(s)+"\n",encoding="utf8")
    ET.parse(out)
    print("SOURCE_SNP_FIGURE_VERIFIED",out)

if __name__=="__main__":
    p=argparse.ArgumentParser()
    p.add_argument("--source",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    a=p.parse_args();run(a.source,a.out)
