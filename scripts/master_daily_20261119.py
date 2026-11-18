#!/usr/bin/env python3
"""CI-only, fail-closed research input-contract inventory; no biological PASS."""
import csv, hashlib, json, os, sys
from pathlib import Path
from datetime import datetime, timezone

CONTRACTS = {
 "is": {"tokens": ["gwas", "harmon"], "required": [["chr","chromosome"],["pos","position"],["ref","other_allele"],["alt","effect_allele"],["beta","effect"],["se","standard_error"],["p","pvalue"]], "gate":"EUR/EAS full-file build/allele harmonization and matched LD"},
 "glaucoma": {"tokens":["dex","tm","count"],"required":[["donor","donor_id"],["sample","sample_id"],["condition","treatment"],["pool","batch"]],"gate":"donor-paired integer count design and leave-pool-out"},
 "hcc": {"tokens":["aliquot","tmt","channel"],"required":[["aliquot","aliquot_id"],["case","case_id","patient_id"],["channel","tmt_channel"]],"gate":"330 TMT channels mapped uniquely to official aliquots and cases"},
 "radiomics": {"tokens":["radiomics","patient","manifest"],"required":[["patient_id","case_id"],["image","image_path"],["rna","rna_sample"]],"gate":"patient-level NIfTI/RNA/clinical mapping and voxel spacing"},
 "ocular": {"tokens":["ena","ocular","run"],"required":[["run_accession"],["sample_accession"],["experiment_accession"],["fastq_md5"],["timepoint","time_point"],["strain"]],"gate":"46 runs actual FASTQ MD5 and timepoint/strain provenance"},
 "muscle": {"tokens":["muscle","loci","ld"],"required":[["rsid","snp","variant_id"],["ancestry","population"],["ld","r2"]],"gate":"20 loci and ancestry-matched LD reference provenance"},
 "oral": {"tokens":["oral","donor","clinical"],"required":[["donor","donor_id"],["sample","sample_id"],["endpoint","outcome"]],"gate":"donor-level expression to clinical endpoint uniqueness"},
 "skin": {"tokens":["skin","mr","instrument"],"required":[["snp","rsid"],["effect_allele","ea"],["other_allele","oa"],["beta"],["se"],["pval","p_value"],["build","genome_build"]],"gate":"phenotype-specific GWAS instrument/build/allele QC"}
}
def main():
 key=sys.argv[1]
 contract=CONTRACTS[key]
 # Only tracked candidate tabular metadata, never raw data or patient values in artifacts.
 paths=[]
 for root, dirs, files in os.walk("."):
  dirs[:]=[d for d in dirs if d not in (".git","node_modules",".venv","venv","results","data","raw","fastq")]
  for name in files:
   p=Path(root)/name
   if p.suffix.lower() not in (".tsv",".csv"):continue
   if not any(t in str(p).lower() for t in contract["tokens"]):continue
   try:
    size=p.stat().st_size
    if 0<size<=2_000_000:paths.append(p)
   except OSError:pass
 result={"research":key,"date":"2026-11-19","technical_status":"PASS","scientific_status":"HOLD","required_gate":contract["gate"],"candidate_files_checked":0,"schema_matched_files":[],"reason":"NO_VERIFIED_PRIMARY_DATA_MAPPING"}
 for p in sorted(paths)[:100]:
  try:
   with p.open(encoding="utf-8-sig",errors="replace",newline="") as f:
    header=next(csv.reader(f,delimiter="\t" if p.suffix.lower()==".tsv" else ","))
   names={h.strip().lower() for h in header}
   matched=all(any(alt.lower() in names for alt in group) for group in contract["required"])
   result["candidate_files_checked"]+=1
   if matched:
    result["schema_matched_files"].append({"path":str(p),"sha256":hashlib.sha256(p.read_bytes()).hexdigest(),"columns_verified":len(contract["required"])})
  except (OSError,UnicodeError,StopIteration,csv.Error):continue
 if result["schema_matched_files"]:
  result["reason"]="HEADER_ONLY_MATCH_PRIMARY_SOURCE_AND_ROW_LEVEL_INTEGRITY_NOT_PROVEN"
 Path("results").mkdir(exist_ok=True)
 Path("results",key+".json").write_text(json.dumps(result,ensure_ascii=False,indent=2)+"\n")
 print(json.dumps(result,ensure_ascii=False))
if __name__=="__main__":main()
