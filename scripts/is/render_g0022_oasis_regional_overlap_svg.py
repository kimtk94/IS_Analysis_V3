#!/usr/bin/env python3
"""Render defensible descriptive G0022 GIGASTROKE vs OASIS regional signal panels.

Data are exact positional/REF/ALT SNV matches only. The panels are *not* coloc
plots or independent replication. SVG output requires Python standard library.
"""
import argparse
import csv
import html
import math
from pathlib import Path

GENE_COLORS={"ALDH2":"#1873a6","BRAP":"#1873a6","RPH3A":"#1873a6"}
GWAS_COL="#a7afb9"
OASIS_COL="#2174a5"


def esc(value):
    return html.escape(str(value),quote=True)


def y_for_logp(lp, ytop, plot_h, ymax):
    return ytop+plot_h-min(max(lp,0),ymax)/ymax*plot_h


def format_svg(path, rows):
    observed=[r for r in rows if r["cell"]=="Mono-L1"]
    if not observed:
        raise ValueError("No monocyte rows found")
    W=1320;H=1130
    xmin=110.5e6;xmax=113.9e6
    margin_left=88;plotW=525;gap=95
    leftW=margin_left; rightW=margin_left+plotW+gap
    plotH=218; rowtops=[157,435,713]
    mx=lambda pos,start:start+(pos-xmin)/(xmax-xmin)*plotW
    s=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}" role="img" aria-label="IS G0022 East Asian AIS and Japanese OASIS monocyte cis eQTL shared-variant association evidence">']
    s.append('<rect width="100%" height="100%" fill="white"/>')
    s.append('<style>text{font-family:Arial,Helvetica,sans-serif;fill:#243040} .title{font-size:23px;font-weight:bold} .subtitle{font-size:13px;fill:#49586a} .axis{font-size:12px;fill:#596575} .panel{font-size:17px;font-weight:bold}</style>')
    s.append('<text class="title" x="70" y="46">IS G0022 — shared-SNP regional association comparison</text>')
    s.append('<text class="subtitle" x="70" y="72">GIGASTROKE East Asian AIS (GRCh37 → GRCh38) vs OASIS Japanese monocytes; strict rsID-independent coordinate + REF/ALT matching</text>')
    s.append('<text class="subtitle" x="70" y="92">Signals are marginal and correlated; vertical markers indicate AIS credible-set sites, not tested OASIS associations.</text>')
    s.append(f'<rect x="{leftW}" y="113" width="12" height="12" fill="{GWAS_COL}"/>')
    s.append(f'<text x="{leftW+20}" y="124" class="axis">AIS GWAS −log₁₀(p)</text>')
    s.append(f'<rect x="{rightW}" y="113" width="12" height="12" fill="{OASIS_COL}"/>')
    s.append(f'<text x="{rightW+20}" y="124" class="axis">OASIS monocyte eQTL −log₁₀(p)</text>')
    markers=[(111730205,"rs11066015","#ab7e1a"),(111803962,"rs671","#c24a3f")]
    for idx,gene in enumerate(("ALDH2","BRAP","RPH3A")):
        sub=[r for r in observed if r["gene"]==gene]
        if not sub:continue
        ymax_gwas=22
        ymax_qtl=max(9.0,math.ceil(max(-math.log10(max(float(r["oasis_p"]),1e-290)) for r in sub)/3)*3)
        top=rowtops[idx]
        s.append(f'<text x="24" y="{top+plotH/2:.0f}" class="panel" transform="rotate(-90 24 {top+plotH/2:.0f})">{gene}</text>')
        for x0,ymax in ((leftW,ymax_gwas),(rightW,ymax_qtl)):
            s.append(f'<rect x="{x0}" y="{top}" width="{plotW}" height="{plotH}" fill="#fbfcfe" stroke="#c4ccd7"/>')
            for val in range(0,int(ymax)+1,3 if ymax<=24 else 6):
                yy=y_for_logp(val,top,plotH,ymax)
                s.append(f'<line x1="{x0}" y1="{yy:.1f}" x2="{x0+plotW}" y2="{yy:.1f}" stroke="#ebeff4"/>')
                s.append(f'<text class="axis" text-anchor="end" x="{x0-7}" y="{yy+4:.1f}">{val}</text>')
            for v in (111,112,113):
                xx=mx(v*1e6,x0)
                s.append(f'<line x1="{xx:.1f}" y1="{top}" x2="{xx:.1f}" y2="{top+plotH}" stroke="#ebeff4"/>')
                s.append(f'<text class="axis" text-anchor="middle" x="{xx:.1f}" y="{top+plotH+17}">{v} Mb</text>')
            for pos,label,color in markers:
                xx=mx(pos,x0)
                s.append(f'<line x1="{xx:.1f}" y1="{top}" x2="{xx:.1f}" y2="{top+plotH}" stroke="{color}" stroke-dasharray="4,5" opacity=".75"/>')
        for r in sub:
            pos=float(r["mapped_grch38_pos"])
            xx=mx(pos,leftW)
            gwp=-math.log10(max(float(r["gwas_p"]),1e-290))
            qp=-math.log10(max(float(r["oasis_p"]),1e-290))
            s.append(f'<circle cx="{xx:.1f}" cy="{y_for_logp(gwp,top,plotH,ymax_gwas):.1f}" r="1.45" fill="{GWAS_COL}" fill-opacity=".48"/>')
            xq=mx(pos,rightW)
            s.append(f'<circle cx="{xq:.1f}" cy="{y_for_logp(qp,top,plotH,ymax_qtl):.1f}" r="1.55" fill="{OASIS_COL}" fill-opacity=".55"/>')
        s.append(f'<text x="{leftW+5}" y="{top+plotH+40}" class="subtitle">Shared source-matched SNVs: {len(sub):,}; association data only</text>')
        s.append(f'<text x="{rightW+5}" y="{top+plotH+40}" class="subtitle">QTL values from OASIS Mono-L1 source</text>')
    s.append('<text x="88" y="1009" class="subtitle">Dashed gold: rs11066015 (observed in both sources); red: rs671 (NOT observed in these OASIS browser panels).</text>')
    s.append('<text x="88" y="1034" class="subtitle">QTL effect-allele semantics, complete cis SNP universe, matched QTL LD and effective per-cell N not verified.</text>')
    s.append('<text x="88" y="1059" class="subtitle">P-value profiles and LD-tagging are not evidence of shared causal variant, directional mediation or therapeutic target.</text>')
    s.append('<text x="88" y="1084" class="subtitle">Source: GIGASTROKE GCST90104545; OASIS Japan Omics Browser; UCSC hg19ToHg38 chain. Analysis 2026-10-10.</text>')
    s.append('</svg>')
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text("\n".join(s)+"\n",encoding="utf-8")
    return {"path":str(path),"gene_mono_matched_counts":{g:sum(r["gene"]==g for r in observed) for g in ("ALDH2","BRAP","RPH3A")}}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--matched-tsv",type=Path,required=True)
    p.add_argument("--out",type=Path,required=True)
    a=p.parse_args()
    with a.matched_tsv.open(newline="") as f:
        rows=list(csv.DictReader(f,delimiter="\t"))
    print(format_svg(a.out,rows))


if __name__=="__main__":
    main()
