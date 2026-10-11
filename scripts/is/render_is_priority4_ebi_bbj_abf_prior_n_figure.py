#!/usr/bin/env python3
"""Publication-style descriptive SVG of 4×3×5 coloc.abf sensitivity scenarios."""
import argparse,csv,math,xml.etree.ElementTree as ET
from pathlib import Path
NS="http://www.w3.org/2000/svg";ET.register_namespace("",NS)
COL={125:"#96a8b7",150:"#e2a165",175:"#197cad"}
def node(p,name,attrs=None,value=None):
 x=ET.SubElement(p,"{"+NS+"}"+name,{k:str(v) for k,v in (attrs or {}).items()})
 if value is not None:x.text=str(value)
 return x
def render(inp,out):
 with inp.open(newline="") as f:r=list(csv.DictReader(f,delimiter="\t"))
 if len(r)!=60:raise ValueError("Requires full 60 scenarios")
 v={}
 for x in r:
  if x["status"]!="EXPLORATORY_ABF_SUCCESS":raise ValueError("Missing ABF scenario")
  key=(x["gene_symbol"],int(x["QTL_N_assumption"]),float(x["assumed_p12"]))
  if key in v:raise ValueError("Repeated scenario")
  v[key]=float(x["PP_H4"])
 W,H=1200,810
 svg=ET.Element("{"+NS+"}svg",{"width":str(W),"height":str(H),"viewBox":f"0 0 {W} {H}",
 "role":"img","aria-label":"Coloc ABF PP H4 across assumed p12 and QTL effective N for four IS priority genes, exploratory only"})
 node(svg,"rect",{"width":W,"height":H,"fill":"#fff"})
 node(svg,"text",{"x":30,"y":43,"font-family":"sans-serif","font-size":22,"font-weight":"bold","fill":"#203950"},
  "Ischemic stroke · reprocessed EBI × BBJ ABF sensitivity")
 node(svg,"text",{"x":30,"y":67,"font-family":"sans-serif","font-size":13,"fill":"#607288"},
  "PP.H4 against assumed shared-causal prior (p12); separate curves assume QTL N = 125, 150 or 175")
 P=[1e-6,3e-6,1e-5,3e-5,1e-4]
 for idx,g in enumerate(("FGF5","C4orf22","CALHM2","NEURL")):
  col=idx%2;row=idx//2;x0=95+col*580;y0=145+row*307;w=440;h=205
  node(svg,"text",{"x":x0,"y":y0-20,"font-size":18,"font-weight":"bold","fill":"#243d57","font-family":"sans-serif"},g)
  def xx(p):return x0+w*(math.log10(p)+6)/2
  def yy(p):return y0+h*(1-p)
  for y in (0,.2,.4,.6,.8,1):
   Y=yy(y);node(svg,"line",{"x1":x0,"x2":x0+w,"y1":Y,"y2":Y,"stroke":"#e6eaf0"})
   node(svg,"text",{"x":x0-13,"y":Y+4,"text-anchor":"end","fill":"#586b80","font-size":11,"font-family":"sans-serif"},f"{y:.1f}")
  for p in P:
   X=xx(p)
   node(svg,"line",{"x1":X,"x2":X,"y1":y0,"y2":y0+h,"stroke":"#e6edf2"})
   node(svg,"text",{"x":X,"y":y0+h+22,"text-anchor":"middle","fill":"#586b80","font-size":11,"font-family":"sans-serif"},f"{p:g}")
  node(svg,"line",{"x1":x0,"x2":x0+w,"y1":yy(.8),"y2":yy(.8),"stroke":"#ba4d52","stroke-dasharray":"5 4","stroke-width":1.3})
  node(svg,"line",{"x1":xx(1e-5),"x2":xx(1e-5),"y1":y0,"y2":y0+h,"stroke":"#63869b","stroke-dasharray":"4 5"})
  for N,color in COL.items():
   points=" ".join(f"{xx(p):.1f},{yy(v[g,N,p]):.1f}" for p in P)
   node(svg,"polyline",{"points":points,"stroke":color,"fill":"none","stroke-width":2.8})
   for p in P:
    node(svg,"circle",{"cx":xx(p),"cy":yy(v[g,N,p]),"r":3.5,"fill":color})
  base=v[g,175,1e-5]
  node(svg,"text",{"x":x0+w-8,"y":y0+19,"text-anchor":"end","fill":"#1d6e9e","font-size":12,"font-family":"sans-serif"},f"N175, p12=1e-5: H4={base:.3f}")
  node(svg,"text",{"x":x0+w/2,"y":y0+h+40,"text-anchor":"middle","fill":"#566a80","font-size":11,"font-family":"sans-serif"},"Assumed p12 (log scale)")
  node(svg,"text",{"x":x0-52,"y":y0+h/2,"transform":f"rotate(-90 {x0-52} {y0+h/2})","text-anchor":"middle","font-family":"sans-serif","font-size":11,"fill":"#566a80"},"PP.H4")
 for i,(N,color) in enumerate(COL.items()):
  xx1=185+i*170
  node(svg,"line",{"x1":xx1,"x2":xx1+32,"y1":775,"y2":775,"stroke":color,"stroke-width":3})
  node(svg,"text",{"x":xx1+41,"y":780,"font-family":"sans-serif","font-size":12,"fill":"#51667b"},f"Assumed QTL N={N}")
 node(svg,"text",{"x":1190,"y":780,"text-anchor":"end","font-family":"sans-serif","font-size":12,"fill":"#a63e41"},
  "Exploratory single-signal ABF · not SuSiE or validated causal colocalization")
 out.parent.mkdir(parents=True,exist_ok=True);ET.ElementTree(svg).write(out,encoding="utf-8",xml_declaration=True)
 return 60
if __name__=="__main__":
 p=argparse.ArgumentParser()
 for k in ("input","out"):p.add_argument("--"+k,type=Path,required=True)
 a=p.parse_args()
 print("SCENARIOS_RENDERED",render(a.input,a.out))
