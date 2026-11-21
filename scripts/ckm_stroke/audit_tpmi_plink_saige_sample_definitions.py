#!/usr/bin/env python3
"""Trace published TPMI 433.21 counts to exact PLINK/SAIGE methodological provenance.

No guess for actual observed cases and controls. Source PheWeb 433.21 reports
9,249 cases and 304,660 controls; independent PLINK2 GWAS OBS_CT=204,099.
Extract 2025 primary JATS methods and link labels from the actual PheWeb
phenotype page. Report evidence with paragraph snippets and remain HOLD.
"""
from __future__ import annotations
import argparse,html,json,re,urllib.request,xml.etree.ElementTree as ET
from html.parser import HTMLParser
from pathlib import Path

PHENO_URL="https://pheweb.ibms.sinica.edu.tw/pheno/433.21"
FULL_TEXT={
 "Chen_TPMI_PRS_2025":"PMC12675292",
 "Yang_TPMI_cohort_2025":"PMC12675286",
}
def grab(url,maxbytes=11000000):
 req=urllib.request.Request(url,headers={"User-Agent":"IS_Analysis_V3-TPMI-N-eff-methods/1.0"})
 with urllib.request.urlopen(req,timeout=55) as r:
  if r.status!=200 or r.geturl()!=url:raise ValueError("source host or redirect")
  buf=r.read(maxbytes+1)
 if len(buf)>maxbytes:raise ValueError("source exceeds metadata cap")
 return buf
class Page(HTMLParser):
 def __init__(self):super().__init__();self.links=[];self.text=[];self.anchor=None
 def handle_starttag(self,tag,attrs):
  if tag=="a":self.anchor=dict(attrs).get("href");self.text=[]
 def handle_endtag(self,tag):
  if tag=="a" and self.anchor:
   self.links.append({"label":" ".join(" ".join(self.text).split())[:140],"href":self.anchor[:600]})
   self.anchor=None
 def handle_data(self,s):
  if self.anchor is not None:self.text.append(s)
def paper(pmcid):
 url=f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML"
 root=ET.fromstring(grab(url))
 matches={}
 terms={
  "PLINK_independent":["plink","independent"],
  "Firth":["firth","logistic"],
  "SAIGE":["saige"],
  "case_control_methods":["case","control","gwas"],
  "sample_exclusions":["remove","exclude","sample"],
 }
 paragraphs=[" ".join("".join(node.itertext()).split()) for node in root.iter() if node.tag.split("}")[-1]=="p"]
 for key,words in terms.items():
  seen=[s for s in paragraphs if all(w.lower() in s.lower() for w in words)]
  matches[key]={"paragraph_count":len(seen),"excerpts":[s[:1900] for s in seen[:3]]}
 return {"source":url,"status":"PASS_PRIMARY_XML","method_evidence":matches}
def main():
 ap=argparse.ArgumentParser();ap.add_argument("--out",required=True,type=Path);a=ap.parse_args()
 out={"schema":"TPMI43321-PLINK2-SAIGE-sample-basis-audit-v1",
      "TPMI_pheweb_cases":9249,"TPMI_pheweb_controls":304660,
      "TPMI_pheweb_total":313909,"source_PLINK2_OBS_CT":204099,
      "N_discrepancy":109810,
      "regression_case_count":None,"regression_control_count":None,
      "per_SNP_case_control_effective_N":None,
      "page_total_equal_PLINK2_effective_N":False,
      "scientific_status":"HOLD_NO_SOURCE_SUPPORTED_CASE_CONTROL_BREAKDOWN"}
 try:
  page=Page();page.feed(grab(PHENO_URL,400000).decode("utf-8","replace"))
  links=[x for x in page.links if any(k in (x["label"]+" "+x["href"]).lower() for k in ("saige","plink","summary","scor","statistics","download"))]
  out["pheweb"]={"status":"PASS_PHEWEB_LINKS","source":PHENO_URL,"data_links":links[:40]}
 except Exception as e:
  out["pheweb"]={"status":"HOLD_PHEWEB","error":str(e)[:230]}
 out["primary_papers"]={}
 for name,pmcid in FULL_TEXT.items():
  try:out["primary_papers"][name]=paper(pmcid)
  except Exception as e:out["primary_papers"][name]={"status":"HOLD_SOURCE","error":str(e)[:240]}
 a.out.parent.mkdir(parents=True,exist_ok=True)
 a.out.write_text(json.dumps(out,indent=2,ensure_ascii=False)+"\n")
 print(json.dumps({"pheweb":out["pheweb"],
                   "primary_papers":out["primary_papers"],
                   "N_discrepancy":109810,
                   "regression_case_control_N_verified":False,
                   "scientific_status":out["scientific_status"]},ensure_ascii=False))
if __name__=="__main__":main()
