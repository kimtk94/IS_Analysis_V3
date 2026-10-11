#!/usr/bin/env python3
"""Render GRCh37 association tracks from fixed exact-allele BBJ/GIGA/QTD joint set."""
import csv,math,argparse,collections,html
from pathlib import Path
LEADS={"BBJ_IS_L001":[81182554,81181072],"BBJ_IS_L002":[105459834,105599770]}
def run(src,out):
 with src.open() as f: allrows=list(csv.DictReader(f,delimiter="\t"))
 if len(allrows)!=20689:raise ValueError("Frozen evidence input drift")
 out.mkdir(parents=True,exist_ok=True)
 summary=[]
 for loc,genes in (("BBJ_IS_L001",("FGF5","C4orf22")),("BBJ_IS_L002",("CALHM2","NEURL"))):
  rows=[r for r in allrows if r["locus"]==loc]
  # By-gene repeated genomic SNP effects from stroke are not independent.
  snps={}
  for r in rows:
   if r["BBJ_exact_Fasta37"]=="True" and r["GIGA_EAS_AIS_exact_Fasta37"]=="True" and r["palindromic_review"]=="False":
    snps[r["SNP_GRCh37"]]=r
  if len(snps)<2500:raise ValueError("GWAS common SNP set unexpected")
  xmin=min(int(k.split(":")[1]) for k in snps);xmax=max(int(k.split(":")[1]) for k in snps)
  W,H=1160,880
  def xcoord(pos):return 108+(pos-xmin)/(xmax-xmin)*1000
  lines=['<?xml version="1.0" encoding="utf-8"?>',
   f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" viewBox="0 0 {W} {H}">',
   '<rect width="100%" height="100%" fill="white"/>',
   f'<text x="60" y="44" font-size="24" font-family="sans-serif" font-weight="bold" fill="#213953">IS {html.escape(loc)} · GRCh37 multi-track evidence</text>',
   '<text x="60" y="70" font-size="13" font-family="sans-serif" fill="#607386">Fixed EBI QTD000161 cis SNP intersection; genomic positions shared by BBJ and GIGASTROKE</text>']
  tracks=[("BBJ IS GWAS","BBJ_p","#245c94",132),("GIGASTROKE EAS AIS","GIGA_EAS_AIS_p","#c4773e",364)]
  for title,key,color,top in tracks:
   bottom=top+194
   lines.append(f'<text x="60" y="{top-14}" font-size="17" font-family="sans-serif" font-weight="bold" fill="{color}">{title}</text>')
   vals=[-math.log10(max(float(r[key]),1e-300)) for r in snps.values() if float(r[key])>0]
   maxy=max(10,math.ceil(max(vals)/2)*2)
   def yy(z):return bottom-min(z,maxy)/maxy*182
   for tick in (0,3,6,9):
    if tick>maxy:continue
    yt=yy(tick);lines.append(f'<line x1="108" x2="1108" y1="{yt:.1f}" y2="{yt:.1f}" stroke="#e6ebef"/>')
    lines.append(f'<text x="96" y="{yt+4:.1f}" text-anchor="end" font-size="12" font-family="sans-serif" fill="#67798a">{tick}</text>')
   for k,r in snps.items():
    z=-math.log10(max(float(r[key]),1e-300))
    lines.append(f'<circle cx="{xcoord(int(k.split(":")[1])):.2f}" cy="{yy(z):.2f}" r="1.8" fill="{color}" opacity=".5"/>')
   thr=yy(-math.log10(5e-8));lines.append(f'<line x1="108" x2="1108" y1="{thr:.1f}" y2="{thr:.1f}" stroke="#a34854" stroke-dasharray="5 5"/>')
  for g_i,gene in enumerate(genes):
   top=595+g_i*100;bottom=top+76
   q=[r for r in rows if r["symbol"]==gene and r["BBJ_exact_Fasta37"]=="True" and r["GIGA_EAS_AIS_exact_Fasta37"]=="True" and r["palindromic_review"]=="False"]
   lines.append(f'<text x="60" y="{top-9}" font-size="15" font-family="sans-serif" fill="#586d51">EBI QTD000161 · {gene} nominal eQTL</text>')
   zmax=max(4,max(-math.log10(max(float(r["QTD000161_nominal_p"]),1e-300)) for r in q))
   for r in q:
    z=-math.log10(max(float(r["QTD000161_nominal_p"]),1e-300))
    lines.append(f'<circle cx="{xcoord(int(r["SNP_GRCh37"].split(":")[1])):.2f}" cy="{bottom-min(z,zmax)/zmax*67:.2f}" r="1.7" fill="#6c9564" opacity=".55"/>')
  for pos in LEADS[loc]:
   if not xmin<=pos<=xmax:continue
   xx=xcoord(pos)
   lines.append(f'<line x1="{xx:.2f}" x2="{xx:.2f}" y1="115" y2="786" stroke="#9b5063" stroke-dasharray="3 7" opacity=".8"/>')
   lines.append(f'<text x="{xx:.2f}" y="112" font-size="11" text-anchor="middle" font-family="sans-serif" fill="#9b5063">{pos}</text>')
  for j in range(5):
   px=xmin+round((xmax-xmin)*j/4);xx=xcoord(px)
   lines.append(f'<text x="{xx:.1f}" y="814" font-size="11" font-family="sans-serif" text-anchor="middle" fill="#536d81">{px:,}</text>')
  lines.append('<text x="570" y="841" text-anchor="middle" font-family="sans-serif" font-size="13" fill="#293e51">Chromosomal position (GRCh37, bp)</text>')
  lines.append('<text x="570" y="862" text-anchor="middle" font-family="sans-serif" font-size="12" fill="#a44250">Nominal p values; no SuSiE/conditional proof · Lead markers from unrestricted GWAS may be excluded from plotted SNPs</text>')
  lines.append('</svg>')
  dest=out/f"IS_{loc}_BBJ_GIGA_EBI_QTL_MULTITRACK_20261011.svg";dest.write_text("\n".join(lines)+"\n")
  summary.append(dict(locus=loc,n_unique_exact_nonpal_snp=len(snps),qtl_gene_labels=",".join(genes),svg=str(dest),p_values_nominal=True))
 with (out/"IS_PRIORITY4_LOCUS_MULTITRACK_SUMMARY.tsv").open("w") as f:
  w=csv.DictWriter(f,fieldnames=list(summary[0]),delimiter="\t");w.writeheader();w.writerows(summary)
 return summary
if __name__=="__main__":
 p=argparse.ArgumentParser();p.add_argument("--source",type=Path,required=True);p.add_argument("--out",type=Path,required=True)
 a=p.parse_args();print(run(a.source,a.out))
