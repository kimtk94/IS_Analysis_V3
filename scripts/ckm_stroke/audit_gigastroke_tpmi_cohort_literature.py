#!/usr/bin/env python3
"""Read primary full-text public evidence relevant to TPMI/GIGASTROKE overlap.

Do not equate TWB PGS validation with GWAS discovery. Approximate TPMI-TWB
overlap is NOT GIGASTROKE AIS discovery overlap. Participant-level matching and
case/control among TPMI PLINK2 OBS_CT remain absent; science HOLD.
"""
import argparse,html,json,re,urllib.request,xml.etree.ElementTree as ET
from pathlib import Path

PAPERS={
 "gigastroke_2022":{
  "pmcid":"PMC9524349","doi":"10.1038/s41586-022-05165-3",
  "checks":{"Taiwan Biobank":["Taiwan Biobank","TWB"],
            "BBJ":["Biobank Japan","BBJ"],
            "EAS":["East Asian","East-Asian"]}},
 "taiwan_early_stroke_2026":{
  "pmcid":"PMC13239414",
  "checks":{"TPMI_TWB":["TPMI","TWB"],
            "overlap":["overlap"],
            "rough_one_percent":["approximately 1%","approximately 1 %","about 1%"]}}
}
def fetch_text(pmcid):
 url=f"https://www.ebi.ac.uk/europepmc/webservices/rest/{pmcid}/fullTextXML"
 req=urllib.request.Request(url,headers={"User-Agent":"IS-Analysis-V3-stroke-literature-evidence/1.0"})
 with urllib.request.urlopen(req,timeout=45) as r:
  if r.status!=200 or r.geturl()!=url:raise RuntimeError("source redirect/HTTP unexpected")
  blob=r.read(12500001)
 if len(blob)>12500000:raise RuntimeError("oversized primary full-text")
 tree=ET.fromstring(blob)
 # Use paragraphs; distinguish primary article and any supplementary narrative.
 paras=[]
 for elem in tree.iter():
  if elem.tag=="p":
   s=" ".join("".join(elem.itertext()).split())
   if s:paras.append(s)
 return url,paras
def main():
 p=argparse.ArgumentParser();p.add_argument("--out",required=True,type=Path)
 args=p.parse_args()
 papers={}
 for name,paper in PAPERS.items():
  try:
   url,paragraphs=fetch_text(paper["pmcid"])
   checks={}
   for key,tokens in paper["checks"].items():
    matches=[s for s in paragraphs if any(t.lower() in s.lower() for t in tokens)]
    checks[key]={"paragraph_count":len(matches),"evidence_excerpt":matches[0][:1500] if matches else None}
   papers[name]={"status":"PASS_PRIMARY_FULLTEXT","source":url,
                 "doi":paper.get("doi"),"checks":checks}
  except Exception as exc:
   papers[name]={"status":"HOLD_PRIMARY_SOURCE_UNAVAILABLE","pmcid":paper["pmcid"],
                 "error_type":type(exc).__name__,"error":str(exc)[:240]}
 result={"schema":"TPMI-GIGASTROKE-source-cohort-evidence-audit-v1",
         "papers":papers,
         "TPMI_TWB_approximate_overlap_is_GIGASTROKE_AIS_discovery_overlap":False,
         "TWB_iPGS_validation_must_not_be_equated_to_AIS_GWAS_discovery":True,
         "GIGASTROKE_EAS_discovery_full_study_specific_cohort_manifest_verified":False,
         "TPMI_GIGASTROKE_participant_level_nonoverlap_verified":False,
         "TPMI_PLINK2_OBS_CT":204099,
         "TPMI_PheWeb_cases":9249,"TPMI_PheWeb_controls":304660,
         "regression_observed_cases_and_controls_verified":False,
         "LDSC_MR_coloc_approved":False,
         "scientific_status":"HOLD_INDEPENDENCE_AND_N_EFF"}
 args.out.parent.mkdir(parents=True,exist_ok=True)
 args.out.write_text(json.dumps(result,indent=2,ensure_ascii=False)+"\n")
 print(json.dumps(result,ensure_ascii=False))
if __name__=="__main__":main()
