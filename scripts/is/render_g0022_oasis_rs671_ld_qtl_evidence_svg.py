#!/usr/bin/env python3
"""Descriptive SNP tagging r2 vs cis-eQTL p chart, never interpreted as coloc."""
import argparse,csv,html,math,xml.etree.ElementTree as ET
from collections import defaultdict
from pathlib import Path

def figure(rows,path):
    groups=defaultdict(list)
    for r in rows:
        k=(r["gene"],r["cell"])
        if r.get("coloc_claim")!="NOT_VALID":
            raise ValueError("Causal claim gate must remain closed")
        groups[k].append(r)
    if len(groups)!=5:
        raise ValueError("Expected five genuine OASIS gene-cell evidence groups")
    cell_labels={"Mono-L1":"Monocytes","B_Activated-L2":"Activated B cells"}
    width,height=1300,1100
    s=[f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width} {height}" width="{width}" height="{height}">']
    s.append('<rect width="100%" height="100%" fill="#fff"/>')
    s.append('<style>text{font-family:Arial,Helvetica,sans-serif;fill:#263343} .ttl{font-size:22px;font-weight:bold} .sub{font-size:13px;fill:#516175} .hdr{font-size:16px;font-weight:bold} .ax{font-size:11px;fill:#65748a}</style>')
    s.append('<text class="ttl" x="55" y="44">IS G0022 — OASIS eQTL signals stratified by rs671 LD</text>')
    s.append('<text class="sub" x="55" y="72">1000 Genomes EAS n=504 LD reference; Y: marginal QTL -log10(p), X: r² to ALDH2 rs671.</text>')
    s.append('<text class="sub" x="55" y="92">Plot includes only genome-build and allele-exact matched GWAS–OASIS SNPs with callable reference LD.</text>')
    sortedkeys=sorted(groups, key=lambda k:({"ALDH2":0,"BRAP":1,"RPH3A":2}[k[0]],k[1]))
    for i,key in enumerate(sortedkeys):
        gene,cell=key
        col=i%2;row=i//2
        x0=74+col*610;y0=157+row*280;pw=480;ph=181
        s.append(f'<text class="hdr" x="{x0}" y="{y0-22}">{html.escape(gene)} / {html.escape(cell_labels.get(cell,cell))} (n={len(groups[key]):,})</text>')
        s.append(f'<rect x="{x0}" y="{y0}" width="{pw}" height="{ph}" fill="#fbfcfe" stroke="#b5c3d1"/>')
        raw=groups[key]
        ycap=max(10., math.ceil(max(-math.log10(max(float(r["qtl_p"]),1e-290)) for r in raw)/10.)*10.)
        for xr in (0,.2,.4,.6,.8,1):
            xx=x0+pw*xr
            s.append(f'<line x1="{xx:.2f}" y1="{y0}" x2="{xx:.2f}" y2="{y0+ph}" stroke="#e6ecf2"/>')
            s.append(f'<text class="ax" text-anchor="middle" x="{xx:.2f}" y="{y0+ph+17}">{xr:g}</text>')
        xx=x0+pw*.8
        s.append(f'<line x1="{xx:.2f}" y1="{y0}" x2="{xx:.2f}" y2="{y0+ph}" stroke="#db9a23" stroke-width="1.3" stroke-dasharray="5,4"/>')
        for yr in (0,.25,.50,.75,1):
            yy=y0+ph*(1-yr)
            s.append(f'<line x1="{x0}" y1="{yy:.2f}" x2="{x0+pw}" y2="{yy:.2f}" stroke="#e6ecf2"/>')
            s.append(f'<text class="ax" text-anchor="end" x="{x0-7}" y="{yy+4:.2f}">{ycap*yr:g}</text>')
        for r in raw:
            v=float(r["EAS504_r2_to_rs671"]);p=float(r["qtl_p"])
            if not (0<=v<=1 and 0<p<=1):raise ValueError("Bad SNP LD/p")
            xx=x0+pw*v
            yy=y0+ph*(1-min(-math.log10(max(p,1e-290))/ycap,1))
            s.append(f'<circle cx="{xx:.2f}" cy="{yy:.2f}" r="1.75" fill="#2772ad" fill-opacity=".55"/>')
        lead=min(raw,key=lambda r:float(r["qtl_p"]))
        xx=x0+pw*float(lead["EAS504_r2_to_rs671"])
        yy=y0+ph*(1-min(-math.log10(max(float(lead["qtl_p"]),1e-290))/ycap,1))
        s.append(f'<circle cx="{xx:.2f}" cy="{yy:.2f}" r="4" fill="#bd443b" stroke="white" stroke-width=".8"/>')
        s.append(f'<text class="sub" x="{x0+pw/2}" y="{y0+ph+39}" text-anchor="middle">EAS reference LD r² with rs671</text>')
    s.append('<text class="sub" x="74" y="1004">Gold line: high-LD reference threshold r²=0.8. Red dot: strongest QTL among SNPs with callable EAS reference LD.</text>')
    s.append('<text class="sub" x="74" y="1032">JPT (n=104) sensitivity in accompanying TSV; it is a subset of EAS reference, not independent validation.</text>')
    s.append('<text class="sub" x="74" y="1060">NOT colocalization: OASIS full tested cis variant coverage, in-study LD and variant effect direction remain incomplete.</text>')
    s.append('</svg>')
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text("\n".join(s)+"\n")
    ET.parse(path)
    return dict(groups=len(groups),plotted_rows=sum(map(len,groups.values())),svg_bytes=path.stat().st_size)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--input",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    args=p.parse_args()
    with args.input.open(newline="") as f:
        rows=list(csv.DictReader(f,delimiter="\t"))
    print(figure(rows,args.out))
if __name__=="__main__":
    main()
