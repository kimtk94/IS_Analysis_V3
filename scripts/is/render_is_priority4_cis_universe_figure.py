#!/usr/bin/env python3
"""Standard-library SVG: GTEx v8 archived source vs reprocessed eQTL Catalogue."""
import argparse,csv,xml.etree.ElementTree as ET
from pathlib import Path
NS="http://www.w3.org/2000/svg";ET.register_namespace("",NS)
def el(root,name,attrs=None,text=None):
 e=ET.SubElement(root,"{"+NS+"}"+name,{k:str(v) for k,v in (attrs or {}).items()})
 if text is not None:e.text=str(text)
 return e
def render(table,output):
 with table.open(newline="") as f:rows=list(csv.DictReader(f,delimiter="\t"))
 if {x["symbol"] for x in rows}!={"FGF5","C4orf22","CALHM2","NEURL"}:raise ValueError("Four genes required")
 root=ET.Element("{"+NS+"}svg",{"width":"1160","height":"680","viewBox":"0 0 1160 680","role":"img","aria-label":"Four IS priority genes original archived SNP overlap versus eQTL Catalogue reprocessed full cis universe"})
 el(root,"rect",{"width":1160,"height":680,"fill":"#fff"})
 el(root,"text",{"x":32,"y":43,"font-family":"sans-serif","font-size":22,"font-weight":"bold","fill":"#15263c"},
    "GTEx brain cerebellar hemisphere — cis variant coverage")
 el(root,"text",{"x":32,"y":69,"font-family":"sans-serif","font-size":13,"fill":"#52657b"},
    "Four archived IS gene–tissue assays compared with eQTL Catalogue QTD000161 reprocessing")
 label=[("catalogue_full_cis_tested_snp_count","Full cis variants (reprocessed)","#8796ab"),
        ("catalogue_unique_snp_in_archived_coordinate_region","Region variants (reprocessed)","#66acc0"),
        ("archived_common_snp_count","Archived BBJ × GTEx input","#4c7aa8"),
        ("archived_overlap_with_catalogue_region","Exact shared variants","#258976")]
 scale=620/6000
 for i,row in enumerate(sorted(rows,key=lambda x:(x["locus"],x["symbol"]))):
  y=115+i*129
  el(root,"text",{"x":38,"y":y+18,"font-family":"sans-serif","font-size":17,"font-weight":"bold","fill":"#203b56"},row["symbol"])
  el(root,"text",{"x":38,"y":y+40,"font-family":"sans-serif","font-size":11,"fill":"#7b8999"},row["locus"])
  for j,(field,title,color) in enumerate(label):
   value=int(row[field]);yy=y+j*25
   el(root,"rect",{"x":205,"y":yy,"width":value*scale,"height":18,"rx":2,"fill":color})
   el(root,"text",{"x":205+value*scale+8,"y":yy+14,"fill":"#34455c","font-family":"sans-serif","font-size":12},f"{value:,}")
   if i==0:
    el(root,"rect",{"x":870,"y":118+j*24,"width":13,"height":13,"fill":color})
    el(root,"text",{"x":891,"y":130+j*24,"font-family":"sans-serif","font-size":12,"fill":"#415572"},title)
  fraction=float(row["shared_SNP_count_div_full_catalogue_cis_count"])
  el(root,"text",{"x":870,"y":y+114,"font-family":"sans-serif","font-size":12,"fill":"#5a6272"},f"Shared / full cis: {100*fraction:.1f}%")
 el(root,"text",{"x":580,"y":661,"text-anchor":"middle","font-family":"sans-serif","font-size":12,"fill":"#8c3636"},
    "Reprocessed tested variant denominators ≠ original GTEx v8 cohort model; no formal coloc inference.")
 output.parent.mkdir(parents=True,exist_ok=True)
 ET.ElementTree(root).write(output,encoding="utf-8",xml_declaration=True)
 return rows
if __name__=="__main__":
 p=argparse.ArgumentParser()
 for x in ("table","output"):p.add_argument("--"+x,type=Path,required=True)
 a=p.parse_args();print("FIGURE_GENES",len(render(a.table,a.output)))
