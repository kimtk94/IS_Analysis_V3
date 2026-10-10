#!/usr/bin/env python3
"""30 source lead exact-variant GWAS alleles: BBJ Japanese vs GWS EAS AIS/AS.

Sources may share cohorts. Identical variant alignment is NOT independent
replication, allelic association mechanism or colocalization. Source-specific
REF/ALT swapped allele pairs are review-only, no automatic complement rescue.
"""
import argparse,csv,json,hashlib,math
from collections import Counter,defaultdict
from pathlib import Path
def read(p):
 with p.open(newline="") as f:return list(csv.DictReader(f,delimiter="\t"))
def sha(p):
 import hashlib
 h=hashlib.sha256()
 with p.open("rb") as f:
  for b in iter(lambda:f.read(4194304),b""):h.update(b)
 return h.hexdigest()
def normalize(lead,r):
 if r["group_source_lead_present"]!="True":
  return {"status":"SOURCE_LEAD_NOT_REPORTED","p":"","beta":"","ALT_AF":""}
 src=r["group_source_lead_variant"]
 l=lead.split(":");s=src.split(":")
 if len(l)!=4 or len(s)!=4 or l[:2]!=s[:2] or set(l[2:])!=set(s[2:]):
  raise ValueError("Original-region source allele pair mismatch")
 if (s[2],s[3])==(l[2],l[3]):
  direction=1;status="EXACT_REF_ALT_PAIR"
 elif (s[2],s[3])==(l[3],l[2]):
  direction=-1;status="SWAPPED_REF_ALT_PAIR_REVIEW"
 else:raise ValueError("Ambiguous source alt")
 beta=float(r["group_source_lead_beta_ALT"])
 af=float(r["group_source_lead_ALT_AF"])
 p=float(r["group_source_lead_p"])
 if not all(math.isfinite(x) for x in (beta,af,p)) or not 0<=af<=1 or not 0<=p<=1:
  raise ValueError("Bad genomic GWAS effect parameters")
 return {"status":status,"p":p,"beta":direction*beta,
   "ALT_AF":af if direction==1 else 1-af}
def build(rows):
 groups=defaultdict(dict)
 for r in rows:
  g=groups[r["region_id"]];trait=r["source_dataset"]
  if trait in g:raise ValueError("Duplicate source/region")
  g[trait]=r
 if len(groups)!=30:raise ValueError("30 original regions required")
 out=[]
 for region in sorted(groups):
  g=groups[region]
  if not {"BBJ_IS","GWS_EAS_AIS","GWS_EAS_AS"}.issubset(g):
   raise ValueError("Source data missing")
  lead=g["BBJ_IS"]["source_region_lead"]
  if len({x["source_region_lead"] for x in g.values()})!=1:
   raise ValueError("Different region source leads")
  b=normalize(lead,g["BBJ_IS"])
  a=normalize(lead,g["GWS_EAS_AIS"])
  allstroke=normalize(lead,g["GWS_EAS_AS"])
  allowed={"EXACT_REF_ALT_PAIR","SWAPPED_REF_ALT_PAIR_REVIEW"}
  ok=b["status"] in allowed and a["status"] in allowed
  if ok:
   dir_agree=(b["beta"]*a["beta"]>0)
   afgap=abs(b["ALT_AF"]-a["ALT_AF"])
  else:dir_agree=None;afgap=None
  out.append({
    "region_id":region,"lead_GRCh37":lead,
    "original_region_phenotypes":g["BBJ_IS"]["original_source_phenotypes"],
    "lead_palindromic":set(lead.split(":")[2:]) in ({"A","T"},{"C","G"}),
    "BBJ_exact_lead_pair_status":b["status"],
    "EAS_AIS_exact_lead_pair_status":a["status"],
    "EAS_AS_exact_lead_pair_status":allstroke["status"],
    "BBJ_source_lead_p":b["p"],"EAS_AIS_source_lead_p":a["p"],
    "EAS_AS_source_lead_p":allstroke["p"],
    "BBJ_source_lead_beta_to_region_ALT":b["beta"],
    "EAS_AIS_source_lead_beta_to_region_ALT":a["beta"],
    "BBJ_ALT_AF":b["ALT_AF"],"EAS_AIS_ALT_AF":a["ALT_AF"],
    "BBJ_vs_AIS_same_ALT_beta_sign_when_exact_pair":dir_agree if dir_agree is not None else "",
    "BBJ_vs_AIS_ALT_AF_absolute_difference":afgap if afgap is not None else "",
    "independent_cohort_replication_proven":False,
    "association_only_not_causal":"YES"
  })
 return out
def main():
 p=argparse.ArgumentParser()
 p.add_argument("--gwas",required=True,type=Path)
 p.add_argument("--out-dir",required=True,type=Path)
 a=p.parse_args()
 if a.out_dir.exists() and any(a.out_dir.iterdir()):
  raise FileExistsError("No prior audit overwrite")
 rows=build(read(a.gwas))
 npaired=sum(x["BBJ_exact_lead_pair_status"]!="SOURCE_LEAD_NOT_REPORTED" and
              x["EAS_AIS_exact_lead_pair_status"]!="SOURCE_LEAD_NOT_REPORTED" for x in rows)
 if npaired!=28:raise ValueError("Original 28 dual-source lead matches changed")
 summary={
  "source":"SOURCE_LEAD_ALLELE_PAIR_DESCRIPTIVE_NO_REPLICATION",
  "n_regions":len(rows),
  "n_exact_source_lead_paired_BBJ_GWS_AIS":npaired,
  "n_paired_same_beta_sign":sum(x["BBJ_vs_AIS_same_ALT_beta_sign_when_exact_pair"] is True for x in rows),
  "n_paired_any_source_refalt_orientation_reversal":sum(x["BBJ_exact_lead_pair_status"]=="SWAPPED_REF_ALT_PAIR_REVIEW" or x["EAS_AIS_exact_lead_pair_status"]=="SWAPPED_REF_ALT_PAIR_REVIEW" for x in rows),
  "n_paired_opposite_beta_sign":sum(x["BBJ_vs_AIS_same_ALT_beta_sign_when_exact_pair"] is False for x in rows),
  "n_paired_AF_delta_ge_0p1":sum(isinstance(x["BBJ_vs_AIS_ALT_AF_absolute_difference"],float) and
      x["BBJ_vs_AIS_ALT_AF_absolute_difference"]>=.1 for x in rows),
  "n_sources_missing_EAS_AIS_lead":sum(x["EAS_AIS_exact_lead_pair_status"]=="SOURCE_LEAD_NOT_REPORTED" for x in rows),
  "GWAS_studies_may_have_shared_subjects":True,
  "same_beta_sign_not_replication":True,
  "source_sha256":sha(a.gwas)
 }
 a.out_dir.mkdir(parents=True,exist_ok=True)
 with (a.out_dir/"IS_EAS30_ORIGINAL_LEAD_BBJ_AIS_AS_ALLELE_ALIGNMENT.tsv").open("w",newline="") as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t");w.writeheader();w.writerows(rows)
 (a.out_dir/"IS_EAS30_SOURCE_LEAD_ALIGNMENT_MANIFEST.json").write_text(json.dumps(summary,indent=2)+"\n")
 print("IS_EAS30_EXACT_LEAD_ALLELES_AUDITED",json.dumps(summary),flush=True)
if __name__=="__main__":main()
