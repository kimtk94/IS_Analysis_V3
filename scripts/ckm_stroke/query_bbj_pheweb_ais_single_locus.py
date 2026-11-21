#!/usr/bin/env python3
"""Query exact GRCh37 TPMI locus against independent BBJ PheWeb AIS phenocode.
Fail closed: appearance in BBJ PheWAS for other diseases is NOT an AIS effect.
"""
import argparse,json,math,urllib.request
from pathlib import Path
URL="https://pheweb.jp/api/variant/19:15573379:T:C"
TARGET={"chrom":"19","pos":15573379,"ref":"T","alt":"C"}
def run():
 p=argparse.ArgumentParser();p.add_argument("--out",required=True,type=Path);a=p.parse_args()
 out={"schema":"TPMI43321-BBJ-PheWeb-IS-one-locus-v1","target":"19:15573379:T:C",
     "provider":"PheWeb.jp","endpoint":URL,
     "dataset":"Sakaue/Kanai et al 2021 BBJ 22664 IS cases 152022 controls",
     "not_equivalent_to":"OpenGWAS bbj-a-129 Ishigaki 2019 17671 cases 192383 controls",
     "independence_from_TPMI":"NOT_PROVEN",
     "scientific_status":"HOLD"}
 try:
  with urllib.request.urlopen(urllib.request.Request(URL,headers={"Accept":"application/json","User-Agent":"TPMI-43321-BBJ-IS-audit/1.0"}),timeout=45) as r:
   if r.status!=200 or r.geturl()!=URL:raise ValueError("API source was redirected")
   body=r.read(1000001)
  if len(body)>1000000:raise ValueError("response too large")
  d=json.loads(body)
  if not all(str(d.get(k))==str(v) for k,v in TARGET.items()):raise ValueError("API variant coordinate or alleles inconsistent")
  phenos=d.get("phenos")
  if not isinstance(phenos,list):raise ValueError("PheWeb API has no phenos list")
  isrows=[v for v in phenos if isinstance(v,dict) and v.get("phenocode")=="IS"]
  out.update({"technical_status":"PASS_EXACT_VARIANT_API_RESPONSE","returned_phenotypes":len(phenos),
    "source_rsid":d.get("rsids",d.get("rsid")),
    "bbj_ischemic_stroke_row_count":len(isrows),
    "other_phenotypes_not_treated_as_IS_replication":True})
  if len(isrows)>1:raise ValueError("ambiguous duplicate IS phenotype rows")
  if len(isrows)==0:
   out["technical_status"]="HOLD_NO_BBJ_PHEWEB_IS_ASSOCIATION_FOR_TARGET"
  else:
   row=isrows[0]
   attrs={k:row.get(k) for k in ("phenocode","phenostring","num_cases","num_controls","beta","sebeta","pval","af","af_alt","maf","nearest_genes")}
   beta,se,pval=attrs["beta"],attrs["sebeta"],attrs["pval"]
   if not all(isinstance(x,(int,float)) and math.isfinite(x) for x in (beta,se,pval)) or se<=0 or not 0<=pval<=1:
    raise ValueError("Invalid numeric BBJ IS effect")
   out["bbj_IS_original"]=attrs
   out["effect_allele_orientation"]="HOLD_UNLESS_PROVIDER_METADATA_CONFIRMS_ALT_BETA"
   out["technical_status"]="PASS_IS_ASSOCIATION_FOUND_ORIENTATION_REVIEW"
 except Exception as exc:
  out.update({"technical_status":"HOLD_SOURCE_OR_VARIANT", "error_type":type(exc).__name__,"error":str(exc)[:240]})
 a.out.parent.mkdir(parents=True,exist_ok=True)
 a.out.write_text(json.dumps(out,indent=2,ensure_ascii=False)+"\n",encoding="utf-8")
 print(json.dumps(out,ensure_ascii=False))
if __name__=="__main__":run()
