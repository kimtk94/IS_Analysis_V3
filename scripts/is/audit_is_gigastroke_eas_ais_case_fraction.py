#!/usr/bin/env python3
"""Cross-check published GIGASTROKE East Asian AS/AIS counts with local EBI YAML.

Values from independent published analyses of Mishra et al. 2022 GIGASTROKE.
The EBI YAML itself states only TOTAL sample_size, no cases/controls.
This is a secondary-source sample fraction, NOT original cohort-by-cohort audit.
Never fill BBJ and do not claim independent GWAS replication or new coloc.
"""
import argparse,csv,hashlib,json,re
from pathlib import Path
CASES={
 "GCST90104544":{
  "trait":"AS","cases":27413,"controls":237242,"total":264655,
  "evidence":"https://kp4cd.org/node/1138"},
 "GCST90104545":{
  "trait":"AIS","cases":19032,"controls":237242,"total":256274,
  "evidence":"https://pmc.ncbi.nlm.nih.gov/articles/PMC11354504/"},
}
SECONDARY_CORROBORATION="https://pmc.ncbi.nlm.nih.gov/articles/PMC13621555/"
def sha(p):
 return hashlib.sha256(p.read_bytes()).hexdigest()
def check(path,study,source=CASES):
 if not path.is_file():raise FileNotFoundError(path)
 text=path.read_text()
 match=re.search(r"(?m)^\s+sample_size:\s*(\d+)\s*$",text)
 if not match:raise ValueError("Source EBI YAML total sample size missing")
 local=int(match.group(1))
 row=source[study]
 if row["cases"]+row["controls"]!=row["total"] or row["total"]!=local:
  raise ValueError("Disagrees with official GWAS total")
 if not re.search(r'(?m)^genome_assembly:\s*GRCh37\s*$',text):
  raise ValueError("GWAS assembly source changed")
 return {
  "GIGASTROKE_GWAS_ID":study,"trait":row["trait"],
  "EBI_total_genotyped_sample_size":local,
  "secondary_published_cases":row["cases"],
  "secondary_published_controls":row["controls"],
  "secondary_fraction_cases":row["cases"]/local,
  "secondary_case_fraction_count_consistency":"MATCHES_LOCAL_EBI_TOTAL",
  "case_count_source":row["evidence"],
  "additional_2026_corrob":SECONDARY_CORROBORATION,
  "case_fraction_status":"SECONDARY_PUBLISHED_RECONCILED_NOT_ORIGINAL_COHORT_LEVEL",
  "per_variant_N_and_case_counts_available":False,
  "BBJ_is_independent_stroke_replication_proven":False,
  "cohort_membership_and_overlap_checked_individually":False,
  "coloc_casecontrol_input_from_original_cohort_table_attested":False
 }
def main():
 ap=argparse.ArgumentParser()
 ap.add_argument("--metadata-dir",type=Path,required=True)
 ap.add_argument("--out-dir",type=Path,required=True)
 a=ap.parse_args()
 if a.out_dir.exists() and any(a.out_dir.iterdir()):raise FileExistsError("No overwrite")
 rows=[]
 for study in ("GCST90104544","GCST90104545"):
  p=a.metadata_dir/f"{study}-meta.yaml"
  rows.append(check(p,study))
 a.out_dir.mkdir(parents=True,exist_ok=True)
 with (a.out_dir/"IS_GIGASTROKE_EAS_AS_AIS_CASE_FRACTION_SECONDARY_CROSSCHECK.tsv").open("w",newline="") as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter="\t")
  w.writeheader();w.writerows(rows)
 manifest={
  "status":"SECONDARY_CASE_CONTROL_FRACTION_VERIFIED_VS_LOCAL_EBI_TOTAL_NO_ORIGINAL_COHORT_BREAKDOWN",
  "AIS_n_cases":19032,"AIS_n_controls":237242,"AIS_total":256274,
  "AIS_secondary_case_fraction":19032/256274,
  "source_primary_GIGASTROKE_paper":"https://www.nature.com/articles/s41586-022-05165-3",
  "AIS_case_control_secondary_source":CASES["GCST90104545"]["evidence"],
  "AIS_case_control_recent_corrob":SECONDARY_CORROBORATION,
  "n_assessed_traits":len(rows),"case_fraction_source_type":"SECONDARY_PAPERS_RECONCILED_EBI_META",
  "sample_overlap_BBJ_not_excluded":True,
  "molecular_QTL_LD_not_verified":True,
  "new_coloc_H4_not_calculated":True,
  "local_YAML_SHA256":{str(a.metadata_dir/f"{sid}-meta.yaml"):sha(a.metadata_dir/f"{sid}-meta.yaml") for sid in CASES}
 }
 (a.out_dir/"IS_GIGASTROKE_EAS_AIS_CASE_FRACTION_PROVENANCE.json").write_text(json.dumps(manifest,indent=2)+"\n")
 print("IS_GIGA_EAS_AIS_19032_CASE_237242_CONTROL_META_CROSSCHECK_PASS",json.dumps({
  "n_case":manifest["AIS_n_cases"],"n_control":manifest["AIS_n_controls"],
  "fraction":manifest["AIS_secondary_case_fraction"]}),flush=True)
if __name__=="__main__":main()
