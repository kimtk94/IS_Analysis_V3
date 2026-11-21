#!/usr/bin/env python3
"""Audit exact TPMI PheWeb Phecode 433.21 phenotype and case/control description.

Phecode != same-number ICD-9-CM diagnosis. Broad GIGASTROKE AIS is a
related but non-identical endpoint unless the case definitions are reconciled.
"""
from __future__ import annotations
import argparse,html,json,re,urllib.request
from html.parser import HTMLParser
from pathlib import Path

URLS={
 "TPMI_phenotype":"https://pheweb.ibms.sinica.edu.tw/pheno/433.21",
 "TPMI_phenotype_filter":"https://pheweb.ibms.sinica.edu.tw/pheno-filter/433.21",
 "TPMI_phenotypes_index":"https://pheweb.ibms.sinica.edu.tw/phenotypes",
}
class Visible(HTMLParser):
 def __init__(self):super().__init__();self.parts=[];self.skip=0
 def handle_starttag(self,tag,attrs):
  if tag in ("script","style"):self.skip+=1
 def handle_endtag(self,tag):
  if tag in ("script","style") and self.skip:self.skip-=1
 def handle_data(self,data):
  if not self.skip:self.parts.append(data)
def fetch(url):
 req=urllib.request.Request(url,headers={"User-Agent":"TPMI-43321-Phecode-source-definition-audit/1.0"})
 with urllib.request.urlopen(req,timeout=35) as r:
  if r.status!=200 or r.geturl()!=url:raise ValueError("PheWeb page changed or redirected")
  body=r.read(2500000)
 if len(body)>=2500000:raise ValueError("PheWeb response truncated")
 parser=Visible();parser.feed(body.decode("utf-8","replace"))
 text=" ".join(" ".join(parser.parts).split())
 if "phenotypes" in url:
  match=re.search(r"(?:^|\s)433\.21\s+(.{0,400})",text)
  piece=match.group(0) if match else None
 else:
  piece=text[:900]
 return {"status":"PASS_HTTP","evidence_text":piece,
         "verified_source_contains_code":("433.21" in text),
         "source_mentions_cerebral_infarction":("infarction" in text.lower())}
def main():
 p=argparse.ArgumentParser();p.add_argument("--out",type=Path,required=True)
 a=p.parse_args()
 docs={}
 for name,url in URLS.items():
  try:docs[name]={"url":url,**fetch(url)}
  except Exception as e:docs[name]={"url":url,"status":"HOLD_HTTP",
                 "error_type":type(e).__name__,"error":str(e)[:220]}
 d={"schema":"TPMI-Phecode-43321-phenotype-ontology-audit-v1",
    "phecode":"433.21",
    "external_standard_phecode_label":"Cerebral artery occlusion, with cerebral infarction",
    "external_standard_ICD9CM_43321_label":"Occlusion and stenosis of vertebral artery with cerebral infarction",
    "numeric_code_across_systems_not_semantically_identical":True,
    "TPMI_evidence":docs,
    "GIGASTROKE_AIS_identity_assumed":False,
    "TOAST_subtype_equivalence_assumed":False,
    "diagnosis_definition_final_status":"HOLD_REVIEW_SOURCE_MAPPING",
    "scientific_status":"HOLD"}
 a.out.parent.mkdir(parents=True,exist_ok=True)
 a.out.write_text(json.dumps(d,indent=2,ensure_ascii=False)+"\n")
 print(json.dumps({"phecode":d["phecode"],"TPMI_evidence":docs,
                   "scientific_status":"HOLD"},ensure_ascii=False))
if __name__=="__main__":main()
