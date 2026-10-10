#!/usr/bin/env python3
"""Standalone stdlib SVG of ancestry-stratified 1000G reference LD (noncausal)."""
import argparse,csv,math,xml.etree.ElementTree as ET
from pathlib import Path
NS="http://www.w3.org/2000/svg"
ET.register_namespace("",NS)
def tag(parent,name,attrs=None,text=None):
    e=ET.SubElement(parent,"{"+NS+"}"+name,{k:str(v) for k,v in (attrs or {}).items()})
    if text is not None:e.text=text
    return e
def read(path):
    with path.open(newline="") as f:return list(csv.DictReader(f,delimiter="\t"))
def render(snp,main,output):
    variants=read(snp);genes=read(main)
    root=ET.Element("{"+NS+"}svg",{"width":"1200","height":"585","viewBox":"0 0 1200 585","role":"img","aria-label":"EAS vs EUR reference LD r-squared scatter by ischemic stroke locus"})
    tag(root,"rect",{"x":0,"y":0,"width":1200,"height":585,"fill":"#ffffff"})
    tag(root,"text",{"x":600,"y":38,"fill":"#15263c","font-family":"sans-serif","font-size":20,"font-weight":"bold","text-anchor":"middle"},"Ischemic stroke: ancestry-stratified reference LD")
    groups=[("BBJ_IS_L001","chr4 · FGF5 / C4orf22",75),("BBJ_IS_L002","chr10 · CALHM2 / NEURL",650)]
    for locus,title,left in groups:
        top=120;w=470;h=370
        tag(root,"text",{"x":left+w/2,"y":85,"font-family":"sans-serif","fill":"#253b55","font-size":17,"text-anchor":"middle"},title)
        def xx(v):return left+w*v
        def yy(v):return top+h*(1-v)
        for tick in (0,.25,.5,.75,1):
            y=yy(tick);x=xx(tick)
            tag(root,"line",{"x1":left,"x2":left+w,"y1":y,"y2":y,"stroke":"#e8edf3","stroke-width":"1"})
            tag(root,"line",{"x1":x,"x2":x,"y1":top,"y2":top+h,"stroke":"#edf1f5","stroke-width":"1"})
            tag(root,"text",{"x":left-11,"y":y+4,"text-anchor":"end","fill":"#526579","font-size":11,"font-family":"sans-serif"},f"{tick:g}")
            tag(root,"text",{"x":x,"y":top+h+20,"text-anchor":"middle","fill":"#526579","font-size":11,"font-family":"sans-serif"},f"{tick:g}")
        tag(root,"rect",{"x":left,"y":top,"width":w,"height":h,"fill":"none","stroke":"#a1adbb"})
        tag(root,"line",{"x1":left,"y1":top+h,"x2":left+w,"y2":top,"stroke":"#a8b1bf","stroke-dasharray":"5,4"})
        ss=[r for r in variants if r["locus"]==locus];shown=0
        for r in ss:
            try:
                a=float(r["EAS_r2_with_GWAS_lead"]);b=float(r["EUR_r2_with_GWAS_lead"])
            except (ValueError,TypeError):continue
            if not (math.isfinite(a) and math.isfinite(b)):continue
            tag(root,"circle",{"cx":round(xx(a),2),"cy":round(yy(b),2),"r":2,"fill":"#4987b3","fill-opacity":.38})
            shown+=1
        tag(root,"text",{"x":left+9,"y":top+22,"font-size":12,"fill":"#31495f","font-family":"sans-serif"},f"Informative SNPs: {shown:,}/{len(ss):,}")
        for g in genes:
            if g["locus"]!=locus:continue
            # Duplicate marker at FGF5/C4orf22 pair is labeled jointly.
            if g["symbol"]=="C4orf22":continue
            x=float(g["EAS_reference_r2"]);y=float(g["EUR_reference_r2"])
            px=xx(x);py=yy(y)
            tag(root,"line",{"x1":px-6,"x2":px+6,"y1":py-6,"y2":py+6,"stroke":"#c44a42","stroke-width":3})
            tag(root,"line",{"x1":px-6,"x2":px+6,"y1":py+6,"y2":py-6,"stroke":"#c44a42","stroke-width":3})
            label="FGF5 / C4orf22" if g["symbol"]=="FGF5" else g["symbol"]
            tag(root,"text",{"x":px+9,"y":py+17,"fill":"#9a2d29","font-family":"sans-serif","font-size":12,"font-weight":"bold"},label)
        tag(root,"text",{"x":left+w/2,"y":top+h+48,"fill":"#42546a","font-family":"sans-serif","font-size":13,"text-anchor":"middle"},"1000G EAS504  r² to GWAS lead")
    tag(root,"text",{"x":17,"y":304,"transform":"rotate(-90 17 304)","fill":"#42546a","font-family":"sans-serif","font-size":13,"text-anchor":"middle"},"1000G EUR503  r² to GWAS lead")
    tag(root,"text",{"x":600,"y":568,"fill":"#627184","font-size":11,"font-family":"sans-serif","text-anchor":"middle"},"Descriptive reference LD only. Not original BBJ/GTEx in-study LD, formal colocalization or a causal-gene result.")
    output.parent.mkdir(parents=True,exist_ok=True)
    ET.ElementTree(root).write(output,encoding="utf-8",xml_declaration=True)
    return sum(1 for r in variants)
if __name__=="__main__":
    p=argparse.ArgumentParser()
    for k in ("snp","main","output"):p.add_argument("--"+k,type=Path,required=True)
    a=p.parse_args()
    print("VARIANT_ROWS=",render(a.snp,a.main,a.output))
