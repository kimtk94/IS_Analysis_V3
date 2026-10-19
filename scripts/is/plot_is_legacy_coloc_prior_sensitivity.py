#!/usr/bin/env python3
"""Draw model-conditioned ABF p12 prior sensitivity, from frozen focal grid."""
import argparse
import csv
import html
import math
from pathlib import Path

FOCUS=("FGF5","CALHM2","SH3PXD2A","COL4A2","ALDH2")
COLORS={"FGF5":"#216a90","CALHM2":"#2f8b69","SH3PXD2A":"#ae6020",
        "COL4A2":"#7354a1","ALDH2":"#ad5162"}


def run(src,out):
    if out.exists():raise FileExistsError(f"Refusing overwrite: {out}")
    with src.open(newline="",encoding="utf8") as h:
        rows=list(csv.DictReader(h,delimiter="\t"))
    by={}
    for r in rows:
        if r["gene_symbol"] in FOCUS:
            by.setdefault(r["gene_symbol"],[]).append(
                (float(r["conditional_p12"]),float(r["PP_H4"])))
    if set(by)!=set(FOCUS) or any(len(v)!=5 for v in by.values()):
        raise ValueError("Expected five prespecified genes x five priors")
    W,H=1120,695
    x0,y0,x1,y1=136,126,813,578
    X=lambda p:x0+(math.log10(p)+6)/2*(x1-x0)
    Y=lambda pp:y1-pp*(y1-y0)
    s=[f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">',
       '<rect width="100%" height="100%" fill="white"/>',
       '<style>text {font-family:Arial,Helvetica,sans-serif;}</style>',
       '<text x="54" y="46" font-size="26" fill="#172f42" font-weight="bold">Legacy IS coloc ABF: prior sensitivity</text>',
       '<text x="54" y="76" font-size="15" fill="#4d6779">Selected four-locus GTEx gene–tissue tests • summary-only conditional recalculation</text>',
       '<text x="54" y="102" font-size="13" fill="#637b87">Assumed baseline p12 = 10⁻⁵; p1 = p2 = 10⁻⁴ | NOT an independent SNP-level rerun</text>']
    for v in (0,.2,.4,.6,.8,1):
        y=Y(v)
        s.extend([f'<line x1="{x0}" x2="{x1}" y1="{y:.2f}" y2="{y:.2f}" stroke="{"#bb6363" if v==.8 else "#e4eaf0"}" stroke-dasharray="{"6 6" if v==.8 else "none"}"/>',
                  f'<text x="{x0-13}" y="{y+5:.2f}" text-anchor="end" font-size="13" fill="#455f70">{v:.1f}</text>'])
    for p,lab in ((1e-6,"10⁻⁶"),(3e-6,"3×10⁻⁶"),(1e-5,"10⁻⁵"),(3e-5,"3×10⁻⁵"),(1e-4,"10⁻⁴")):
        x=X(p)
        s.extend([f'<line x1="{x:.2f}" y1="{y0}" x2="{x:.2f}" y2="{y1}" stroke="{"#b0ccd4" if p==1e-5 else "#ecf0f2"}"/>',
                  f'<text x="{x:.2f}" y="{y1+28}" text-anchor="middle" font-size="13" fill="#455f70">{lab}</text>'])
    s.append(f'<text x="{(x0+x1)/2}" y="661" text-anchor="middle" font-size="15" fill="#364f62">Assumed shared-variant prior probability per SNP (p12, logarithmic axis)</text>')
    s.append(f'<text x="33" y="{(y0+y1)/2}" transform="rotate(-90 33 {(y0+y1)/2})" text-anchor="middle" font-size="15" fill="#364f62">Conditional PP.H4</text>')
    for i,g in enumerate(FOCUS):
        data=sorted(by[g])
        if sorted(p for p,v in data)!=[1e-6,3e-6,1e-5,3e-5,1e-4]:
            raise ValueError("Unexpected prior grid")
        color=COLORS[g]
        coords=" ".join(f'{X(p):.2f},{Y(v):.2f}' for p,v in data)
        s.append(f'<polyline points="{coords}" fill="none" stroke="{color}" stroke-width="3.2"/>')
        for p,v in data:
            s.append(f'<circle cx="{X(p):.2f}" cy="{Y(v):.2f}" r="4.5" fill="white" stroke="{color}" stroke-width="2.6"/>')
        legend_y=208+i*45
        s.extend([f'<line x1="849" x2="879" y1="{legend_y}" y2="{legend_y}" stroke="{color}" stroke-width="4"/>',
                  f'<text x="891" y="{legend_y+5}" font-size="17" fill="{color}" font-weight="bold">{html.escape(g)}</text>'])
    s.extend(['<text x="840" y="481" font-size="14" fill="#465d6e">Red dashed line: PP.H4 = 0.8</text>',
              '<text x="840" y="506" font-size="13" fill="#546c7d">Not a validated causal threshold</text>',
              '<text x="55" y="687" font-size="12" fill="#6c7e8a">Selected baseline-best tissues; input prior settings and SNP-level QTL/LD provenance still unconfirmed</text>',
              '</svg>'])
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text("\n".join(s)+"\n",encoding="utf8")
    print("IS_P12_FIGURE_PASS",out.stat().st_size,out)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--top-grid",required=True,type=Path)
    p.add_argument("--out",required=True,type=Path)
    a=p.parse_args()
    run(a.top_grid,a.out)

if __name__=="__main__":
    main()
