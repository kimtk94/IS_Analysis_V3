#!/usr/bin/env python3
"""SVG of reference-only population LD with 2000-sample bootstrap intervals."""
import csv,sys,xml.etree.ElementTree as E
from pathlib import Path
NS="http://www.w3.org/2000/svg";E.register_namespace("",NS)
data={x["population"]:x for x in csv.DictReader(open(sys.argv[1]),delimiter="\t")}
names=["JPT104","otherEAS400","EAS504","EUR503"]
if not all(x in data for x in names):raise ValueError("Bootstrap data missing")
root=E.Element("{"+NS+"}svg",{"width":"1050","height":"470","viewBox":"0 0 1050 470"})
def add(t,a,txt=None):
 e=E.SubElement(root,"{"+NS+"}"+t,{k:str(v) for k,v in a.items()})
 if txt is not None:e.text=str(txt)
add("rect",{"width":1050,"height":470,"fill":"#fff"})
add("text",{"x":32,"y":44,"font-size":22,"font-family":"sans-serif","fill":"#20384e"}, "IS chr10 · JPT/EAS versus EUR reference LD")
add("text",{"x":32,"y":74,"font-size":13,"font-family":"sans-serif","fill":"#52677c"}, "Candidate 10:105694487:G:T vs GWAS lead 10:105459834:T:C; 2,000 bootstrap resamples")
def xcoord(r):return 230+(r+.5)*685
for a in (-.5,-.4,-.3,-.2,-.1,0,.1,.2,.3):
 x=xcoord(a);add("line",{"x1":x,"x2":x,"y1":114,"y2":365,"stroke":"#dce5eb"})
 add("text",{"x":x,"y":391,"text-anchor":"middle","font-size":12,"font-family":"sans-serif"},str(a))
colors=["#187e9e","#699cb1","#4b68a1","#db8343"]
for i,p in enumerate(names):
 d=data[p];y=157+i*55;lo=float(d["bootstrap_r_CI_low"]);hi=float(d["bootstrap_r_CI_high"]);v=float(d["signed_r_observed"])
 add("text",{"x":32,"y":y+5,"font-family":"sans-serif","font-size":15,"fill":colors[i]},p)
 add("line",{"x1":xcoord(lo),"x2":xcoord(hi),"y1":y,"y2":y,"stroke":colors[i],"stroke-width":6,"stroke-linecap":"round"})
 add("circle",{"cx":xcoord(v),"cy":y,"r":6,"fill":colors[i]})
 add("text",{"x":950,"y":y+5,"font-family":"sans-serif","font-size":13,"fill":colors[i]},f"r={v:+.3f}")
add("text",{"x":520,"y":432,"text-anchor":"middle","font-family":"sans-serif","font-size":13,"fill":"#943e45"},"Reference genotype sampling only; NOT conditional GWAS or causal colocalization")
out=Path(sys.argv[2]);out.parent.mkdir(parents=True,exist_ok=True)
E.ElementTree(root).write(out,encoding="utf-8",xml_declaration=True)
print("FIGURE_COMPLETE",out)
