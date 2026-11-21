#!/usr/bin/env python3
"""Verify Supplementary Table 1 EAS AIS source cohorts and totals against Catalog.

Source: GIGASTROKE 2022 paper supplementary table 1 XLSX.
J = any stroke cases, K = AIS cases, P = shared discovery controls.
CKBB cases are NA for AIS, but CKBB controls are in total. Do not treat
iPGS-validation Taiwan Biobank as AIS discovery. A public cohort list does
NOT certify participant-level nonoverlap between studies.
"""
import argparse,json,re
from pathlib import Path
from scripts.ckm_stroke.read_gigastroke_eas_cohort_supplement import URLS,fetch,xml_book

EAS_CASES={"KCPS2_ASN**":46,"BIOBANK_JAPAN*":17651,
           "HISAYAMA_FSR":1113,"INTERSTROKE_ASN":222}
EAS_CASE_CONTROLS={"KCPS2_ASN**":12366,"BIOBANK_JAPAN*":159610,
                   "HISAYAMA_FSR":901,"INTERSTROKE_ASN":698}
EXPECTED_CASES=19032;EXPECTED_CONTROLS=237242
def value_at(row,col):
 vals={re.match(r"[A-Z]+",cell).group():str(value).strip() for cell,value in row}
 return vals.get(col)
def original_table_eas():
 from io import BytesIO
 from zipfile import is_zipfile
 for url in URLS[:2]:
  try:
   body,source=fetch(url)
   if not is_zipfile(BytesIO(body)):continue
   sheet=next(x for x in xml_book(body) if x["sheet"]=="Suppl.Table1")
   return source,sheet["supp_table1_rows_45_to_76"]
  except Exception:continue
 raise ValueError("original publisher supplement unavailable")
def audit():
 source,rows=original_table_eas()
 r50=next(r for r in rows if any(k=="A50" for k,v in r))
 if value_at(r50,"A")!="EAS":raise ValueError("missing EAS ancestry table anchor")
 region=[row for row in rows if any(re.search(r"\d+$",k) and 50<=int(re.search(r"\d+$",k).group())<=64 for k,v in row)]
 records=[]
 for row in region:
  rownum=int(re.search(r"\d+$",row[0][0]).group())
  name=value_at(row,"B")
  if not name:raise ValueError("missing EAS cohort label at "+str(rownum))
  case=value_at(row,"K");controls=value_at(row,"P")
  ncase=int(case) if case and case.isdecimal() else None
  ncontrols=int(controls) if controls and controls.isdecimal() else None
  records.append({"supp_table_row":rownum,"study":name,
    "EAS_AIS_cases":ncase,"EAS_discovery_controls":ncontrols,
    "AIS_case_contributor":ncase is not None and ncase>0,
    "classification":"incident_only_not_in_global_meta" if rownum==60
                     else "EAS_discovery"})
 usable=[r for r in records if r["classification"]=="EAS_discovery"]
 contributors={x["study"]:x["EAS_AIS_cases"] for x in usable if x["AIS_case_contributor"]}
 if contributors!=EAS_CASES:
  raise ValueError("AIS contributing cohort identity or n changed "+str(contributors))
 for row in usable:
  if row["study"] in EAS_CASE_CONTROLS and row["EAS_discovery_controls"]!=EAS_CASE_CONTROLS[row["study"]]:
   raise ValueError("published AIS contributor control count changed")
 ckbb=[x for x in usable if x["study"].startswith("CKBB_ASN")]
 if len(ckbb)!=10 or any(x["EAS_AIS_cases"] is not None for x in ckbb):
  raise ValueError("CKBB regional AI S cases not reported as NA")
 cases=sum(x["EAS_AIS_cases"] or 0 for x in usable)
 controls=sum(x["EAS_discovery_controls"] or 0 for x in usable)
 if (cases,controls)!=(EXPECTED_CASES,EXPECTED_CONTROLS):
  raise ValueError(f"EAS AIS source total differs from catalog: cases={cases} controls={controls}")
 total=next(r for r in rows if any(k=="A65" for k,v in r))
 if value_at(total,"A")!="TOTAL EAS" or int(value_at(total,"K"))!=cases or int(value_at(total,"P"))!=controls:
  raise ValueError("source published total row disagrees")
 candidates=[x["study"] for x in usable]
 if any("TPMI" in x.upper() or "TAIWAN" in x.upper() or re.search(r"\bTWB\b",x.upper()) for x in candidates):
  raise ValueError("Taiwan/TPMI listed: overlap assessment required")
 return {"schema":"GIGASTROKE-EAS-AIS-primary-discovery-cohort-v1",
   "official_supplement_url":source,
   "table":"Suppl.Table1, rows 50–65, combined incident and prevalent AIS K, discovery controls P",
   "AIS_discovery_cohorts":records,
   "GIGASTROKE_EAS_AIS_discovery_cases":cases,
   "GIGASTROKE_EAS_AIS_discovery_controls":controls,
   "CKBB_regional_controls_total":sum(x["EAS_discovery_controls"] for x in ckbb),
   "AIS_case_contributor_count":len(contributors),
   "AIS_case_contributors":contributors,
   "EAS_discovery_study_rows":len(usable),
   "TPMI_listed_in_discovery":False,
   "TWB_listed_in_discovery":False,
   "taiwan_biobank_only_reported_as_PGS_validation_in_primary_text":True,
   "source_based_cohort_list_status":"PASS_PRIMARY_SUPPLEMENT_COHORT_INVENTORY_AND_COUNTS",
   "participant_level_independence_verified":False,
   "case_control_N_eff_verified":False,
   "scientific_status":"HOLD_PARTICIPANT_OVERLAP_AND_N_EFF"}
def main():
 p=argparse.ArgumentParser();p.add_argument("--out",required=True,type=Path);a=p.parse_args()
 try:r=audit()
 except Exception as e:
  r={"schema":"GIGASTROKE-EAS-AIS-primary-discovery-cohort-v1",
     "source_based_cohort_list_status":"HOLD_SOURCE_OR_COHORT_TOTAL",
     "error":str(e)[:400],"scientific_status":"HOLD"}
 a.out.parent.mkdir(parents=True,exist_ok=True)
 a.out.write_text(json.dumps(r,indent=2,ensure_ascii=False)+"\n")
 print(json.dumps({"status":r["source_based_cohort_list_status"],
  "cases":r.get("GIGASTROKE_EAS_AIS_discovery_cases"),
  "controls":r.get("GIGASTROKE_EAS_AIS_discovery_controls"),
  "ckbb_controls":r.get("CKBB_regional_controls_total"),
  "case_contributors":r.get("AIS_case_contributors"),
  "error":r.get("error"),"scientific_status":r["scientific_status"]},ensure_ascii=False))
if __name__=="__main__":main()
