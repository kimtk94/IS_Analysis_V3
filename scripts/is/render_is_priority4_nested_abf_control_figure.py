#!/usr/bin/env python3
"""Descriptive SVG of the controlled same-source narrow vs broad cis-SNP ABF contrast."""
import argparse,csv,xml.etree.ElementTree as ET
from pathlib import Path
NS="http://www.w3.org/2000/svg";ET.register_namespace("",NS)
def add(svg,name,attrs=None,label=None):
 x=ET.SubElement(svg,"{"+NS+"}"+name,{k:str(v) for k,v in (attrs or {}).items()})
 if label is not None:x.text=label
 return x
def render(table,out):
 with table.open(newline="") as f:r=list(csv.DictReader(f,delimiter="\t"))
 if {x["gene_symbol"] for x in r}!={"FGF5","C4orf22","CALHM2","NEURL"}:raise ValueError("Four genes required")
 root=ET.Element("{"+NS+"}svg",{"width":"1120","height":"580","viewBox":"0 0 1120 580","role":"img",
  "aria-label":"Four ischemic stroke genes held-constant EBI BBJ ABF nested vs expanded SNP sets"})
 add(root,"rect",{"width":1120,"height":580,"fill":"#fff"})
 add(root,"text",{"x":38,"y":43,"font-family":"sans-serif","font-size":22,"font-weight":"bold","fill":"#163449"},
  "Controlled SNP-universe test: same EBI QTL and BBJ GWAS")
 add(root,"text",{"x":38,"y":68,"font-family":"sans-serif","font-size":13,"fill":"#4f6980"},
  "Fixed QTL N=175; p1=p2=1e-4, p12=1e-5 · only the matched cis-SNP inclusion set differs")
 left=250;right=868;width=right-left
 for tick in (0,.2,.4,.6,.8,1):
  X=left+width*tick
  add(root,"line",{"x1":X,"x2":X,"y1":128,"y2":475,"stroke":"#e8eef3"})
  add(root,"text",{"x":X,"y":493,"font-size":12,"text-anchor":"middle","fill":"#52677e","font-family":"sans-serif"},f"{tick:.1f}")
 add(root,"text",{"x":(left+right)/2,"y":525,"font-size":13,"text-anchor":"middle","fill":"#52677e","font-family":"sans-serif"},"Exploratory PP.H4")
 for i,x in enumerate(r):
  y=134+i*86
  symbol=x["gene_symbol"];n=float(x["baseline_nested_PP_H4"]);v=float(x["baseline_full_PP_H4"])
  add(root,"text",{"x":38,"y":y+17,"font-size":16,"font-weight":"bold","fill":"#294258","font-family":"sans-serif"},symbol)
  add(root,"text",{"x":38,"y":y+39,"font-size":11,"fill":"#8190a2","font-family":"sans-serif"},f"n {x['n_nested']} → {x['n_full']} variants")
  for offset,value,fill in ((0,n,"#8e9eae"),(23,v,"#197da4")):
   add(root,"rect",{"x":left,"y":y+offset,"width":width*value,"height":18,"rx":2,"fill":fill})
   add(root,"text",{"x":left+width*value+7,"y":y+offset+14,"font-size":12,"font-family":"sans-serif","fill":"#40546a"},f"{value:.3f}")
  add(root,"text",{"x":955,"y":y+30,"text-anchor":"middle","font-size":14,"font-weight":"bold","font-family":"sans-serif","fill":"#9b464b"},
   f"Δ {float(x['baseline_PP_H4_full_minus_nested']):+.3f}")
 add(root,"rect",{"x":60,"y":542,"width":16,"height":11,"fill":"#8e9eae"})
 add(root,"text",{"x":84,"y":552,"font-size":11,"font-family":"sans-serif","fill":"#53677c"},"Nested archived-overlap SNP set, same QTD stats")
 add(root,"rect",{"x":470,"y":542,"width":16,"height":11,"fill":"#197da4"})
 add(root,"text",{"x":494,"y":552,"font-size":11,"font-family":"sans-serif","fill":"#53677c"},"Expanded EBI full cis ∩ BBJ SNP set")
 add(root,"text",{"x":1100,"y":569,"text-anchor":"end","font-size":10,"font-family":"sans-serif","fill":"#aa5050"},
   "Posterior is prior/N sensitive; no validated SuSiE or causal sharing")
 out.parent.mkdir(parents=True,exist_ok=True)
 ET.ElementTree(root).write(out,encoding="utf-8",xml_declaration=True)
 return len(r)
if __name__=="__main__":
 p=argparse.ArgumentParser()
 for k in ("table","out"):p.add_argument("--"+k,type=Path,required=True)
 a=p.parse_args();print("FIG_N_GENES",render(a.table,a.out))
