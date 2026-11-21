#!/usr/bin/env python3
"""Inventory original GIGASTROKE 2022 supplementary material with provenance.

Do not infer absence of TPMI from a table without confirming its EAS AIS
discovery composition; supplementary correction was published Dec 2022.
"""
import argparse,json,re,urllib.request,xml.etree.ElementTree as ET
from pathlib import Path
URL="https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9524349/fullTextXML"
XLINK="{http://www.w3.org/1999/xlink}href"
def main():
 p=argparse.ArgumentParser();p.add_argument("--out",type=Path,required=True);a=p.parse_args()
 result={"schema":"GIGASTROKE-2022-supplement-inventory-v1","source":URL,
         "scientific_status":"HOLD","EAS_AIS_discovery_cohort_list_verified":False}
 try:
  req=urllib.request.Request(URL,headers={"User-Agent":"IS_Analysis_V3_GIGASTROKE_supplement-audit/1.0"})
  with urllib.request.urlopen(req,timeout=60) as r:
   if r.status!=200 or r.geturl()!=URL:raise RuntimeError("unapproved XML source")
   body=r.read(13000000)
  if len(body)>=13000000:raise ValueError("primary source XML size cap")
  root=ET.fromstring(body)
  media=[]
  for elem in root.iter():
   if elem.tag.split("}")[-1] in ("supplementary-material","media","ext-link"):
    href=elem.attrib.get(XLINK) or elem.attrib.get("href")
    label=" ".join("".join(elem.itertext()).split())[:240]
    if elem.tag.split("}")[-1]=="ext-link" and not any(x in (href or "").lower() for x in ("supp","41586","mmedia")):continue
    media.append({"tag":elem.tag.split("}")[-1],"id":elem.attrib.get("id"),
                  "href":href,"mime":elem.attrib.get("mimetype"),
                  "subtype":elem.attrib.get("mime-subtype"),"label":label})
  result.update({"status":"PASS_PRIMARY_XML","supplement_refs":media,"XML_bytes":len(body),
                 "is_publisher_correction_2022_checked":False})
 except Exception as exc:
  result.update({"status":"HOLD_SOURCE","error":str(exc)[:300]})
 a.out.parent.mkdir(parents=True,exist_ok=True)
 a.out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n")
 print(json.dumps({"status":result.get("status"),"references":result.get("supplement_refs",[]),
                   "error":result.get("error"),"science":"HOLD"},ensure_ascii=False))
if __name__=="__main__":main()
