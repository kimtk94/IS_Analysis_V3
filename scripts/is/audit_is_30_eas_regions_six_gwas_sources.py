#!/usr/bin/env python3
"""Full-file read-only GWAS coverage for 30 provisional IS_XDATA GRCh37 groups.

Every source GWAS is streamed ONCE. Reports exact native group lead allele-set
match plus a contextual group +/-500kb window. Window overlap is NOT LD-defined,
regions are NOT independent causal loci, cross-study direction isn't replication.
No source/archived coloc data is modified.
"""
import argparse,csv,gzip,hashlib,json,math,sys
from collections import Counter
from pathlib import Path

TRAITS={
 "BBJ_IS":"japan/bbj/BBJ_IS_GRCh37.canonical.tsv.gz",
 "GWS_EAS_AS":"gigastroke/eas/GCST90104544_AS_GRCh37.canonical.tsv.gz",
 "GWS_EAS_AIS":"gigastroke/eas/GCST90104545_AIS_GRCh37.canonical.tsv.gz",
 "GWS_EAS_CES":"gigastroke/eas/GCST90104546_CES_GRCh37.canonical.tsv.gz",
 "GWS_EAS_LAS":"gigastroke/eas/GCST90104547_LAS_GRCh37.canonical.tsv.gz",
 "GWS_EAS_SVS":"gigastroke/eas/GCST90104548_SVS_GRCh37.canonical.tsv.gz",
}
def sha256(p):
 h=hashlib.sha256()
 with p.open("rb") as f:
  for buf in iter(lambda:f.read(4*1024*1024),b""):h.update(buf)
 return h.hexdigest()
def table(p):
 with p.open(newline="") as f:return list(csv.DictReader(f,delimiter="\t"))
def write(p,rows):
 if not rows:raise ValueError("Empty")
 with p.open("w",newline="") as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t")
  w.writeheader();w.writerows(rows)
def regions_read(path,flank):
 rows=table(path)
 if len(rows)!=30 or len({r["group_id"] for r in rows})!=30:
  raise ValueError("Must contain exactly 30 unique XDATA regions")
 bychr={}
 for r in rows:
  if not r["group_id"].startswith("IS_XDATA_"):
   raise ValueError("Non-XDATA region")
  chrom=r["chr"].replace("chr","")
  s,e=int(r["start"]),int(r["end"])
  if not(1<=int(chrom)<=22 and 1<=s<=e):raise ValueError("bad chrom/window")
  r["window_start_GRCh37"]=max(1,s-flank)
  r["window_end_GRCh37"]=e+flank
  bychr.setdefault(chrom,[]).append(r)
 return rows,bychr
def variant_same(a,b):
 p=a.split(":");q=b.split(":")
 return len(p)==4 and len(q)==4 and p[:2]==q[:2] and set(p[2:])==set(q[2:])
def empty_record(group,trait):
 return {"region_id":group["group_id"],"source_dataset":trait,
 "original_source_phenotypes":group["phenotypes"],
 "component_has_GWS":group["has_gws"],
 "chromosome_GRCh37":group["chr"],
 "original_region_start":group["start"],"original_region_end":group["end"],
 "context_window_start":group["window_start_GRCh37"],
 "context_window_end":group["window_end_GRCh37"],
 "source_region_lead":group["lead_variant"],
 "n_rows_in_context":0,"n_valid_variant_allele_records":0,
 "n_genomewide_p_lt_5e8":0,"n_suggestive_p_lt_1e6":0,
 "n_nominal_p_lt_0p05":0,"n_invalid_or_ambiguous_record":0,
 "best_p":"","best_variant":"","best_beta_ALT":"","best_effect_ALLELE":"",
 "group_source_lead_present":False,"group_source_lead_variant":"",
 "group_source_lead_p":"","group_source_lead_beta_ALT":"",
 "group_source_lead_ALT_AF":"","group_source_lead_allele_status":"NOT_PRESENT",
 "analysis_scope":"DESCRIPTIVE_GWAS_COVERAGE_NOT_CAUSAL_COLOC"}
def scan(source,trait,groups,bychr):
 ledger={r["group_id"]:empty_record(r,trait) for r in groups}
 processed=0;matching=0;dataset_allele_bad=0;invalidp=0
 with gzip.open(source,"rt",encoding="utf-8",errors="strict") as f:
  hdr=f.readline().rstrip("\r\n").split("\t")
  req=("chr","pos","variant_id","effect_allele","other_allele","beta","p","eaf")
  if not all(k in hdr for k in req):
   raise ValueError("Source GWAS schema missing "+str(source))
  j={k:hdr.index(k) for k in req}
  for line in f:
   processed+=1
   # chr/pos are before the effect and P columns in all canonical inputs
   pref=line.split("\t",max(j["chr"],j["pos"])+1)
   if len(pref)<=max(j["chr"],j["pos"]):continue
   chrom=pref[j["chr"]].replace("chr","")
   reg=bychr.get(chrom)
   if not reg:continue
   try:pos=int(pref[j["pos"]])
   except ValueError:continue
   matched=[r for r in reg if r["window_start_GRCh37"]<=pos<=r["window_end_GRCh37"]]
   if not matched:continue
   parts=line.rstrip("\r\n").split("\t")
   matching+=1
   try:
    var=parts[j["variant_id"]]
    tokens=var.split(":")
    ea=parts[j["effect_allele"]].upper()
    oa=parts[j["other_allele"]].upper()
    beta=float(parts[j["beta"]]);p=float(parts[j["p"]]);af=float(parts[j["eaf"]])
    ok=len(tokens)==4 and tokens[0]==chrom and int(tokens[1])==pos and (
      (ea==tokens[2] and oa==tokens[3]) or (ea==tokens[3] and oa==tokens[2])
    ) and math.isfinite(beta) and math.isfinite(p) and 0<=p<=1 and 0<=af<=1
   except (ValueError,IndexError):
    ok=False
   if not ok:
    dataset_allele_bad+=1
    for r in matched:ledger[r["group_id"]]["n_invalid_or_ambiguous_record"]+=1
    continue
   aligned_beta=beta if ea==tokens[3] else -beta
   aligned_af=af if ea==tokens[3] else 1-af
   for r in matched:
    d=ledger[r["group_id"]]
    d["n_rows_in_context"]+=1
    d["n_valid_variant_allele_records"]+=1
    d["n_genomewide_p_lt_5e8"]+=p<5e-8
    d["n_suggestive_p_lt_1e6"]+=p<1e-6
    d["n_nominal_p_lt_0p05"]+=p<.05
    if d["best_p"]=="" or p<float(d["best_p"]):
     d.update({"best_p":p,"best_variant":var,
       "best_beta_ALT":aligned_beta,"best_effect_ALLELE":ea})
    if variant_same(var,r["lead_variant"]):
     if d["group_source_lead_present"]:
      d["group_source_lead_allele_status"]="DUPLICATE_SOURCE_MATCH"
     else:
      d.update({"group_source_lead_present":True,
      "group_source_lead_variant":var,"group_source_lead_p":p,
      "group_source_lead_beta_ALT":aligned_beta,
      "group_source_lead_ALT_AF":aligned_af,
      "group_source_lead_allele_status":"SOURCE_VARIANT_UNORDERED_ALLELE_PAIR_MATCH"})
  if processed%5000000==0:print("PROGRESS",trait,processed,flush=True)
 result=list(ledger.values())
 if len(result)!=30:raise ValueError("Expected 30 coverage records")
 print("SOURCE_GWAS_SCAN_DONE",trait,"whole_rows",processed,
   "context_rows",sum(d["n_rows_in_context"] for d in result),
   "bad_records",dataset_allele_bad,flush=True)
 return result,{"trait":trait,"total_native_rows_scanned":processed,
  "matched_native_rows_in_any_context":matching,
  "invalid_records_in_any_context":dataset_allele_bad,
  "path":str(source),"sha256":sha256(source)}
def main():
 ap=argparse.ArgumentParser()
 ap.add_argument("--regions",required=True,type=Path)
 ap.add_argument("--out-dir",required=True,type=Path)
 ap.add_argument("--processed-root",required=True,type=Path)
 ap.add_argument("--flank",type=int,default=500000)
 a=ap.parse_args()
 if a.flank!=500000:raise ValueError("Window intentionally frozen +/-500kb")
 if a.out_dir.exists() and any(a.out_dir.iterdir()):raise FileExistsError("No overwrite")
 groups,bychr=regions_read(a.regions,a.flank)
 allrows=[];sources=[]
 for trait,relative in TRAITS.items():
  path=a.processed_root/relative
  if not path.is_file():raise FileNotFoundError(path)
  r,manifest=scan(path,trait,groups,bychr)
  allrows.extend(r);sources.append(manifest)
 a.out_dir.mkdir(parents=True,exist_ok=True)
 write(a.out_dir/"IS_30_REGIONS_SIX_OUTCOMES_GWAS_CONTEXT_500KB.tsv",allrows)
 cross=[]
 for r in groups:
  x=[v for v in allrows if v["region_id"]==r["group_id"]]
  cross.append({
   "region_id":r["group_id"],"chromosome_GRCh37":r["chr"],
   "region_original_phenotypes":r["phenotypes"],
   "region_has_GWS":r["has_gws"],"original_region_lead":r["lead_variant"],
   "n_source_datasets_with_any_context_variant":sum(v["n_valid_variant_allele_records"]>0 for v in x),
   "n_source_datasets_reporting_original_region_lead":sum(v["group_source_lead_present"] is True for v in x),
   "n_BBJ_GWS_EAS_original_region_lead_both":all(v["group_source_lead_present"] is True for v in x if v["source_dataset"] in ("BBJ_IS","GWS_EAS_AIS")),
   "n_outcomes_with_5e8_context_signal":sum(v["n_genomewide_p_lt_5e8"]>0 for v in x),
   "BBJ_500kb_min_p":next(v["best_p"] for v in x if v["source_dataset"]=="BBJ_IS"),
   "EAS_AIS_500kb_min_p":next(v["best_p"] for v in x if v["source_dataset"]=="GWS_EAS_AIS"),
   "interpretation":"REGIONAL_COVERAGE_NOT_INDEPENDENT_LOCI_OR_REPLICATION"
  })
 write(a.out_dir/"IS_30_REGIONS_OUTCOME_AVAILABILITY_MATRIX.tsv",cross)
 info={"status":"READ_ONLY_SIX_EXISTING_OUTCOME_COVERAGE_30_PROVISIONAL_GRCH37_INTERVALS",
  "n_provisional_regions":30,"n_outcome_sources":len(TRAITS),
  "n_region_dataset_rows":len(allrows),"flank_bp":a.flank,
  "regions_are_not_independent_LD_defined_loci":True,
  "no_genetic_causal_claim":True,"no_new_data_download":True,
  "original_region_source_sha256":sha256(a.regions),
  "source_scans":sources}
 (a.out_dir/"IS_30_REGION_GWAS_SOURCE_SCAN_MANIFEST.json").write_text(json.dumps(info,indent=2)+"\n")
 print("IS_EAS30_OUTCOME_GWAS_SOURCE_SCAN_PASS",json.dumps({
  "n_dataset_records":len(allrows),
  "region_leads_with_BBJ_and_EAS_AIS":sum(x["n_BBJ_GWS_EAS_original_region_lead_both"] for x in cross),
  "regions_with_EAS_AIS_min_p_5e8":sum(float(x["EAS_AIS_500kb_min_p"])<5e-8 for x in cross if x["EAS_AIS_500kb_min_p"]!=""),
 }),flush=True)
if __name__=="__main__":main()
