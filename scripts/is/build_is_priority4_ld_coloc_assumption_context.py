#!/usr/bin/env python3
"""Couple reference LD proxy neighborhoods to archived ABF/prior sensitivity.

Scientific boundary: PPH4 IS NOT recalculated based on LD differences.
The QTL source is the original matched SNP subset, not a complete tested universe.
"""
import argparse,csv,hashlib,json
from pathlib import Path
from xml.etree import ElementTree as ET

NS="http://www.w3.org/2000/svg"
ET.register_namespace("",NS)

def get_rows(path):
 with path.open(newline="") as f:return list(csv.DictReader(f,delimiter="\t"))
def svg_node(root,name,attrib=None,content=None):
 e=ET.SubElement(root,"{"+NS+"}"+name,{k:str(v) for k,v in (attrib or {}).items()})
 if content is not None:e.text=str(content)
 return e
def render(rows,path):
 root=ET.Element("{"+NS+"}svg",{"width":"1160","height":"460","viewBox":"0 0 1160 460","role":"img",
      "aria-label":"Reference LD r squared and historical ABF prior sensitivity comparison"})
 svg_node(root,"rect",{"width":1160,"height":460,"fill":"#fff"})
 svg_node(root,"text",{"x":36,"y":39,"font-family":"sans-serif","font-size":21,"font-weight":"bold","fill":"#1e3450"},"IS locus evidence · reference LD versus archived coloc assumptions")
 svg_node(root,"text",{"x":36,"y":64,"font-family":"sans-serif","font-size":12,"fill":"#607086"},"Separate scales with same 0–1 range; numeric comparison does not infer causal sharing.")
 cols=((230,480,"LD tagging of strongest eQTL","r², GWAS lead → eQTL lead"),
       (650,900,"Original ABF / p12 sweep","PP.H4; original prior not independently attested"))
 for a,b,title,subtitle in cols:
  svg_node(root,"text",{"x":a,"y":109,"font-size":14,"font-weight":"bold","font-family":"sans-serif","fill":"#223d58"},title)
  svg_node(root,"text",{"x":a,"y":126,"font-size":11,"font-family":"sans-serif","fill":"#698098"},subtitle)
  for x in (0,.25,.5,.75,1):
   xp=a+(b-a)*x
   svg_node(root,"line",{"x1":xp,"x2":xp,"y1":142,"y2":340,"stroke":"#e4eaf0"})
   svg_node(root,"text",{"x":xp,"y":357,"text-anchor":"middle","fill":"#688099","font-size":11,"font-family":"sans-serif"},f"{x:g}")
 for i,g in enumerate(rows):
  y=161+i*56
  svg_node(root,"text",{"x":34,"y":y+5,"font-size":15,"font-weight":"bold","font-family":"sans-serif","fill":"#223d58"},g["symbol"])
  svg_node(root,"text",{"x":34,"y":y+20,"font-size":10,"font-family":"sans-serif","fill":"#8493a2"},g["locus"])
  for x,fill,label in ((float(g["EAS_reference_r2"]), "#1972a2","EAS"),(float(g["EUR_reference_r2"]), "#dd8435","EUR")):
   xp=230+250*x
   yp=y-5 if label=="EAS" else y+9
   svg_node(root,"circle",{"cx":xp,"cy":yp,"r":5.5,"fill":fill})
  # Posterior min/max under the assumed prior set; archive baseline black dot
  lo=float(g["prior_min_H4"]);hi=float(g["prior_max_H4"]);h=float(g["archived_PP_H4"])
  svg_node(root,"line",{"x1":650+250*lo,"x2":650+250*hi,"y1":y+2,"y2":y+2,
                         "stroke":"#8594aa","stroke-width":4,"stroke-linecap":"round"})
  svg_node(root,"circle",{"cx":650+250*h,"cy":y+2,"r":6,"fill":"#1e3450"})
  svg_node(root,"text",{"x":925,"y":y+5,"font-family":"sans-serif","font-size":12,"fill":"#506176"},
           f"H4 {h:.3f}; eQTL p<1e-4: {g['n_nominal_QTL_p_lt_1e4']}")
 svg_node(root,"circle",{"cx":232,"cy":383,"r":5,"fill":"#1972a2"})
 svg_node(root,"text",{"x":245,"y":387,"font-family":"sans-serif","font-size":12,"fill":"#4f6278"},"1000G EAS504")
 svg_node(root,"circle",{"cx":386,"cy":383,"r":5,"fill":"#dd8435"})
 svg_node(root,"text",{"x":399,"y":387,"font-family":"sans-serif","font-size":12,"fill":"#4f6278"},"1000G EUR503")
 svg_node(root,"line",{"x1":661,"x2":701,"y1":382,"y2":382,"stroke":"#8594aa","stroke-width":4})
 svg_node(root,"circle",{"cx":681,"cy":382,"r":5,"fill":"#1e3450"})
 svg_node(root,"text",{"x":712,"y":387,"font-family":"sans-serif","font-size":12,"fill":"#4f6278"},"Original H4 & prior-sweep range")
 svg_node(root,"text",{"x":580,"y":428,"font-family":"sans-serif","font-size":12,"text-anchor":"middle","fill":"#8d393a"},
   "NOT formal multi-signal colocalization. Reference LD does NOT update coloc.abf PP.H4.")
 path.parent.mkdir(exist_ok=True,parents=True)
 ET.ElementTree(root).write(path,encoding="utf-8",xml_declaration=True)

def run(gene,prior,replay,out,figure):
 g=get_rows(gene)
 j=get_rows(prior)
 a=get_rows(replay)
 if len(g)!=4 or len(j)!=30 or len(a)!=646:raise ValueError("Archive cohort count changed")
 one={(r["gene_id"],r["locus"],r["dataset_key"]):r for r in j}
 old={(r["gene_base"],r["locus"],r["dataset_key"]):r for r in a}
 result=[]
 for x in g:
  key=(x["gene_id"],x["locus"],x["tissue"])
  if key not in one or key not in old:raise ValueError("Gene/tissue/prior join not unique")
  pp=one[key];arch=old[key]
  if arch["status"]!="PASS":raise ValueError("Historical assay not replayed")
  h=float(arch["PP_H4"])
  if abs(h-float(pp["archived_H4"]))>1e-8:raise ValueError("Prior baseline cannot replicate old H4")
  r=dict(locus=x["locus"],symbol=x["symbol"],gene_id=x["gene_id"],tissue=x["tissue"],
         strongest_QTL_SNP=x["qtl_top_variant_b37"],nominal_QTL_p=x["qtl_top_nominal_p"],
         EAS_reference_r2=x["qtl_top_r2_eas"],EUR_reference_r2=x["qtl_top_r2_eur"],
         qtl_top_proxy_category=x["qtl_top_proxy_category"],
         n_nominal_QTL_p_lt_1e4=x["n_nominal_QTL_p_lt_1e4"],
         GWAS_proxy_QTL_p_lt_1e4_high_LD_in_both=x["nominal_QTL_high_ld_both"],
         archived_PP_H3=arch["PP_H3"],archived_PP_H4=h,
         prior_min_H4=pp["prior_H4_min"],prior_max_H4=pp["prior_H4_max"],
         assumption_sensitivity_only=True,
         posterior_reestimated_after_ancestry_LD=False,
         hypothesis="NO_CAUSAL_ASSERTION")
  result.append(r)
 result.sort(key=lambda x:("BBJ_IS_L002"!=x["locus"],-x["archived_PP_H4"],x["symbol"]))
 out.mkdir(parents=True,exist_ok=True)
 path=out/"IS_PRIORITY4_LD_COLOC_ASSUMPTION_CONTEXT.tsv"
 with path.open("w",newline="") as f:
  w=csv.DictWriter(f,fieldnames=list(result[0]),delimiter="\t");w.writeheader();w.writerows(result)
 render(result,figure)
 summary=dict(source_genes=4,distinct_loci=2,original_candidate_gene_universe=2225,
   genes_top_QTL_lead_r2_ge_0p8_in_both=sum(r["qtl_top_proxy_category"]=="HIGH_LD_BOTH" for r in result),
   current_analysis_recomputes_PP_H4_using_LD=False,
   original_coloc_ancestry_LD_input_provided=False,
   original_p12_execution_provenance_verified=False,
   figure_path=str(figure),rows=result,
   source_SHA256={"gene":hashlib.sha256(gene.read_bytes()).hexdigest(),"prior":hashlib.sha256(prior.read_bytes()).hexdigest(),"replay":hashlib.sha256(replay.read_bytes()).hexdigest()})
 (out/"IS_PRIORITY4_LD_COLOC_ASSUMPTION_CONTEXT.json").write_text(json.dumps(summary,indent=2)+"\n")
 return summary

if __name__=="__main__":
 p=argparse.ArgumentParser()
 for k in ("gene","prior","replay","out","figure"):p.add_argument("--"+k,type=Path,required=True)
 a=p.parse_args()
 x=run(a.gene,a.prior,a.replay,a.out,a.figure)
 print(json.dumps({k:v for k,v in x.items() if k not in ("rows","source_SHA256")},indent=2))
