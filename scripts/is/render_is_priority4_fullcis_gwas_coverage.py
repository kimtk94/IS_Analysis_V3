#!/usr/bin/env python3
"""Standalone SVG of independently computed full-cis QTL and GWAS exact coverage."""
import argparse,csv,xml.etree.ElementTree as ET
from pathlib import Path
NS="http://www.w3.org/2000/svg";ET.register_namespace("",NS)
METRICS=(
 ("biallelic_snp_FASTA37_verified","QTL SNP valid liftover","#7284a5"),
 ("exact_BBJ_GWAS_overlap","BBJ IS exact","#247baa"),
 ("exact_GIGASTROKE_EAS_AIS_overlap","GIGASTROKE EAS AIS exact","#db8b37"),
 ("exact_BOTH_nonpal_overlap","Both GWAS; non-palindromic","#268a75"),
 ("EBI_only_with_BOTH_GWAS","NEW vs archive; both GWAS","#9b569a"))
def element(root,tag,attrs={},txt=None):
 e=ET.SubElement(root,"{"+NS+"}"+tag,{k:str(v) for k,v in attrs.items()})
 if txt is not None:e.text=str(txt)
 return e
def render(table,output):
 with table.open(newline="") as f:rows=list(csv.DictReader(f,delimiter="\t"))
 if {r["symbol"] for r in rows}!={"FGF5","C4orf22","CALHM2","NEURL"}:raise ValueError("Unexpected four-gene source")
 W,H=1220,755
 root=ET.Element("{"+NS+"}svg",{"width":str(W),"height":str(H),"viewBox":f"0 0 {W} {H}","role":"img",
     "aria-label":"Four ischemic stroke priority genes EBI reprocessed full cis source BBJ and GIGASTROKE EAS AIS allele-exact GWAS overlap"})
 element(root,"rect",{"width":W,"height":H,"fill":"#ffffff"})
 element(root,"text",{"x":45,"y":48,"font-family":"sans-serif","font-size":23,"font-weight":"bold","fill":"#192d42"},
    "IS priority genes · full-cis eQTL versus GWAS variant coverage")
 element(root,"text",{"x":45,"y":76,"font-family":"sans-serif","font-size":13,"fill":"#4d637e"},
    "Numerator: FASTA-anchored SNP overlap; denominator: all tested QTD000161 cis variants per gene")
 basex=230;barwidth=570
 for i,row in enumerate(rows):
  y=108+i*150
  sym=row["symbol"];full=int(row["full_ebi_reprocessed_cis_snps"])
  element(root,"text",{"x":43,"y":y+14,"font-family":"sans-serif","font-size":17,"font-weight":"bold","fill":"#203b56"},sym)
  element(root,"text",{"x":43,"y":y+34,"font-family":"sans-serif","font-size":11,"fill":"#6d7d8e"},f"Full cis = {full:,}")
  for j,(col,label,color) in enumerate(METRICS):
   val=int(row[col]);pct=100*val/full
   yy=y+j*24
   element(root,"rect",{"x":basex,"y":yy,"width":barwidth,"height":16,"fill":"#edf1f5","rx":2})
   element(root,"rect",{"x":basex,"y":yy,"width":round(barwidth*val/full,2),"height":16,"fill":color,"rx":2})
   element(root,"text",{"x":basex+barwidth+9,"y":yy+13,"font-family":"sans-serif","font-size":12,"fill":"#32465c"},f"{val:,} ({pct:.1f}%)")
   if i==0:
    element(root,"rect",{"x":970,"y":113+j*24,"width":10,"height":10,"fill":color})
    element(root,"text",{"x":989,"y":123+j*24,"font-family":"sans-serif","font-size":11,"fill":"#52647a"},label)
 element(root,"text",{"x":610,"y":722,"text-anchor":"middle","font-family":"sans-serif","font-size":12,"fill":"#994b4b"},
   "Reference overlap and marginal association only · NOT GTEx donor LD, formal colocalization, or causal proof")
 output.parent.mkdir(exist_ok=True,parents=True)
 ET.ElementTree(root).write(output,encoding="utf-8",xml_declaration=True)
 return len(rows)
if __name__=="__main__":
 p=argparse.ArgumentParser()
 for k in ("table","output"):p.add_argument("--"+k,type=Path,required=True)
 a=p.parse_args();print("RENDERED",render(a.table,a.output))
