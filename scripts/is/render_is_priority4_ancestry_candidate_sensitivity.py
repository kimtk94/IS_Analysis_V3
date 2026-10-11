#!/usr/bin/env python3
"""Render directly observed EAS/EUR LD comparison; no inferred causality."""
import argparse,csv,xml.etree.ElementTree as ET
from pathlib import Path
NS="http://www.w3.org/2000/svg";ET.register_namespace("",NS)
def elem(parent,name,attrs={},txt=None):
 e=ET.SubElement(parent,"{"+NS+"}"+name,{k:str(v) for k,v in attrs.items()})
 if txt is not None:e.text=str(txt)
 return e
def render(input,output):
 with input.open(newline="") as f:rr=list(csv.DictReader(f,delimiter="\t"))
 if len(rr)!=10:raise ValueError("5x2 ancestry-shrinkage rows required")
 dat={(r["ancestry"],float(r["alpha"])):r for r in rr}
 if len(dat)!=10:raise ValueError("Duplicate scenario")
 W,H=1090,520
 root=ET.Element("{"+NS+"}svg",{"width":str(W),"height":str(H),"viewBox":f"0 0 {W} {H}",
 "role":"img","aria-label":"EAS versus EUR reference signed LD and heuristic lead-conditional z for stroke candidate chr10 105694487 G T"})
 elem(root,"rect",{"width":W,"height":H,"fill":"white"})
 elem(root,"text",{"x":35,"y":46,"font-family":"sans-serif","font-size":23,"font-weight":"bold","fill":"#192f49"},
      "Stroke chr10 · ancestry changes the LD residual interpretation")
 elem(root,"text",{"x":35,"y":74,"font-family":"sans-serif","font-size":13,"fill":"#566c80"},
      "Candidate 10:105694487:G:T vs lead 10:105459834:T:C · matched ALT dosage")
 colors={"EAS504":"#16839a","EUR503":"#d57c3c"}
 # Left: signed r, baseline alpha 0.
 elem(root,"text",{"x":55,"y":121,"font-size":16,"font-weight":"bold","font-family":"sans-serif","fill":"#263f58"},
      "Signed LD with BBJ lead (r)")
 xcenter=287;y1=205;y2=281;scale=680
 elem(root,"line",{"x1":80,"x2":515,"y1":y1+30,"y2":y1+30,"stroke":"#d9e3ec"})
 elem(root,"line",{"x1":xcenter,"x2":xcenter,"y1":164,"y2":331,"stroke":"#9eb3c3"})
 for y,anc in ((y1,"EAS504"),(y2,"EUR503")):
  r=float(dat[anc,0]["lead_candidate_signed_r"])
  end=xcenter+scale*r
  elem(root,"line",{"x1":xcenter,"x2":end,"y1":y,"y2":y,"stroke":colors[anc],"stroke-width":11,"stroke-linecap":"round"})
  elem(root,"circle",{"cx":end,"cy":y,"r":6,"fill":colors[anc]})
  elem(root,"text",{"x":58,"y":y+29,"font-family":"sans-serif","font-size":13,"fill":"#344f67"},anc+f"  r={r:+.4f}")
 elem(root,"text",{"x":290,"y":354,"text-anchor":"middle","font-size":12,"font-family":"sans-serif","fill":"#62778a"},
      "Same SNPs, opposite-signed correlations")
 # Right: heuristic residual z, separate ancestor-specific trajectories.
 lx=622;rx=1030;top=154;bottom=352
 def xx(a):return lx+(rx-lx)*a/.20
 def yy(z):return bottom-(z-3.6)/2.3*(bottom-top)
 for v in (4,4.5,5,5.45):
  y=yy(v)
  elem(root,"line",{"x1":lx,"x2":rx,"y1":y,"y2":y,"stroke":"#e0e7ed"})
  elem(root,"text",{"x":lx-9,"y":y+4,"text-anchor":"end","font-size":11,"font-family":"sans-serif","fill":"#617486"},str(v))
 elem(root,"text",{"x":lx,"y":123,"font-family":"sans-serif","font-size":16,"font-weight":"bold","fill":"#263f58"},"Heuristic lead-conditioned z (NOT GWAS p)")
 y=yy(5.45)
 elem(root,"line",{"x1":lx,"x2":rx,"y1":y,"y2":y,"stroke":"#af424a","stroke-dasharray":"5 5","stroke-width":2})
 for anc in ("EAS504","EUR503"):
  row=[dat[anc,a] for a in (0,.01,.05,.1,.2)]
  coords=[(xx(float(r["alpha"])),yy(float(r["approx_lead_residual_z"]))) for r in row]
  elem(root,"polyline",{"points":" ".join(f"{x:.2f},{y:.2f}" for x,y in coords),"stroke":colors[anc],
       "fill":"none","stroke-width":3.5})
  for x,y in coords:elem(root,"circle",{"cx":x,"cy":y,"r":4,"fill":colors[anc]})
 for a in (0,.05,.1,.2):
  elem(root,"text",{"x":xx(a),"y":bottom+22,"text-anchor":"middle","font-family":"sans-serif","font-size":12,"fill":"#63768b"},f"{a:g}")
 elem(root,"text",{"x":(lx+rx)/2,"y":bottom+43,"text-anchor":"middle","font-family":"sans-serif","font-size":12,"fill":"#5b7084"},"Reference LD shrinkage alpha")
 for i,anc in enumerate(("EAS504","EUR503")):
  elem(root,"line",{"x1":670+i*165,"x2":696+i*165,"y1":426,"y2":426,"stroke":colors[anc],"stroke-width":4})
  elem(root,"text",{"x":702+i*165,"y":431,"font-family":"sans-serif","font-size":12,"fill":"#465e72"},anc)
 elem(root,"text",{"x":540,"y":479,"text-anchor":"middle","font-size":13,"font-family":"sans-serif","fill":"#9b3e45"},
      "Reference-LD approximation only: NOT a formally calibrated conditional test, independence or causal colocalization")
 output.parent.mkdir(parents=True,exist_ok=True)
 ET.ElementTree(root).write(output,encoding="utf-8",xml_declaration=True)
 return 10
if __name__=="__main__":
 p=argparse.ArgumentParser()
 p.add_argument("--input",required=True,type=Path);p.add_argument("--out",required=True,type=Path)
 a=p.parse_args();print("DRAWN",render(a.input,a.out))
