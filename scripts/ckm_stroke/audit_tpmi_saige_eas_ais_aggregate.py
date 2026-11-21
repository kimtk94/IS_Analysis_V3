#!/usr/bin/env python3
"""Independent fail-closed aggregate audit of a real SAIGE x EAS AIS run artifact.

Checks source provenance, original genomewide counts, cohort-model
independence semantics, SNP filtering conservation, aligned counts and allele
frequency buckets. This validates already-exported aggregate observations,
not allele/effect truth or patient-level independence. No GWAS row access.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
from collections import Counter

EXPECTED={
 "SAIGE_original_sha256_pinned":"9f235bfec6265442b7699accecbbd2a69c5acb369f7ee1bcc9e79ec94285a33e",
 "GIGASTROKE_EAS_AIS_source_sha256_pinned":"c521d5695414bdedcc84d2ea21d588e330c75a672f36a87c3fd092b4c55687a3",
 "UCSC_hg38ToHg19_chain_md5":"ff3031d93792f4cbb86af44055efd903",
 "SAIGE_original_rows":15588353,
 "SAIGE_full_original_QC_PASS":8042797,
 "effective_N_SAIGE_case":9249,
 "effective_N_SAIGE_control":304660,
 "scientific_status":"HOLD",
 "NOT_TPMI_PLINK2":True,
 "TPMI_source_model":"SAIGE",
 "LDSC_MR_coloc_performed":False,
 "EAS_sample_overlap_certified":False,
 "individual_SNP_effects_exported":False,
}
def check(condition,why):
 if not condition:raise ValueError(why)
def validate(d):
 check(isinstance(d,dict),"output artifact JSON not object")
 for k,v in EXPECTED.items():
  check(d.get(k)==v,f"source metadata/hold gate drift: {k}")
 check(d.get("technical_status")=="PASS_TECHNICAL_SAIGE_EAS_ALLELE_INTERSECTION",
       "original allele harmonization did not pass")
 s=d.get("SAIGE_source_filter_liftover_QC",{})
 e=d.get("GIGASTROKE_EAS_source_index_QC",{})
 m=d.get("allele_overlap_QC",{})
 f=d.get("AF_delta_bins",{})
 check(s.get("all_original_rows")==15588353,"SAIGE original rows drift")
 check(s.get("SOURCE_QC_FAIL")==7545556,"SAIGE original FAIL count drift")
 check(s.get("all_original_rows")-s.get("SOURCE_QC_FAIL")==8042797,
       "SAIGE original PASS conservation failed")
 check(s.get("CANDIDATE",0)+s.get("PALINDROMIC_UNRESOLVED",0)+
       s.get("NOT_BIALLELIC_SNV",0)+s.get("MALFORMED_MARKER_ID",0)==8042797,
       "source-QC-PASS eligible candidate partition changed")
 check(s.get("LIFT_UNIQUE",0)+s.get("LIFT_UNMAPPED",0)+
       s.get("LIFT_TRANSCHROM_OR_INVALID",0)==s.get("CANDIDATE",0),
       "liftover outcomes do not partition source candidate rows")
 check(s.get("unique_saige_lifted_coords",0)>0 and
       s.get("LIFT_UNIQUE",0)>=s.get("unique_saige_lifted_coords",0),
       "unique SAIGE coords exceeds successful lifted rows")
 check(s.get("LIFT_UNIQUE",0)-s.get("unique_saige_lifted_coords",0)==
       s.get("saige_ambiguous_lifted_coords",0),
       "expected exactly one duplicate record for each ambiguous coordinate")
 check(e.get("all_original_rows")==6725216 and
       e.get("candidate_outcome_rows")==6725216 and
       e.get("unique_outcome_coordinates")==6725216 and
       e.get("outcome_duplicate_coordinate_positions")==0,
       "GIGASTROKE exact source or duplicate count drift")
 check(m.get("COORDINATE_INTERSECTIONS",0)==
       m.get("MATCHED",0)+m.get("ALLELES_MISMATCH",0),
       "matched and mismatched source alleles do not exhaust coordinates")
 check(0<m.get("MATCHED",0)<=s.get("unique_saige_lifted_coords",0),
       "matched SNP number outside valid SAIGE coordinates")
 check(m.get("COORDINATE_INTERSECTIONS",0)<=e["unique_outcome_coordinates"],
       "coordinate intersections cannot exceed EAS outcome source")
 check(sum(f.values())==m["MATCHED"],"AF frequency bins do not sum to matched")
 check(m.get("AF_DELTA_GT_0.1",0)==
       f.get("DELTA_GT_0.1_LE_0.2",0)+f.get("DELTA_GT_0.2",0),
       "AF gap >0.1 disagrees with frequency buckets")
 check(m.get("INFO_GE_0.8_AND_SAIGE_MAF_GE_0.01",0)<=
       min(m["MATCHED"],m.get("INFO_GE_0.8",0),m.get("SAIGE_MAF_GE_0.01",0)),
       "joint INFO/MAF eligible count exceeds marginals")
 check(d.get("SAIGE_effect_allele")=="Allele2","wrong association allele")
 check(d.get("SAIGE_source_GRCh38_official")=="https://pheweb.ibms.sinica.edu.tw/about",
       "missing official build source")
 check(d.get("GIGASTROKE_EAS_AIS_source_build")=="GRCh37",
       "wrong GIGASTROKE build")
 out={
  "schema":"TPMI-SAIGE-EAS-AIS-original-overlap-independent-conservation-audit-v1",
  "audit_status":"PASS_AGGREGATE_CONSERVATION_AND_PROVENANCE",
  "source_harmonization_action":"https://github.com/kimtk94/IS_Analysis_V3/actions/runs/44920288867",
  "SAIGE_source_model_only":True,
  "SAIGE_source_QC_PASS":s["all_original_rows"]-s["SOURCE_QC_FAIL"],
  "SAIGE_unambiguous_original_candidate_count":s["CANDIDATE"],
  "SAIGE_lifted_unique_coordinate_count":s["unique_saige_lifted_coords"],
  "allele_matched":m["MATCHED"],
  "INFO_ge_0_8_MAF_ge_0_01":m.get("INFO_GE_0.8_AND_SAIGE_MAF_GE_0.01"),
  "AF_delta_gt_0_05":sum(n for k,n in f.items() if k!="DELTA_LE_0.05"),
  "AF_delta_gt_0_10":m["AF_DELTA_GT_0.1"],
  "source_allele_mismatches":m.get("ALLELES_MISMATCH"),
  "allele_matched_pct_of_source_QC_PASS":round(
    100*m["MATCHED"]/(s["all_original_rows"]-s["SOURCE_QC_FAIL"]),3),
  "AF_delta_gt_0_05_pct_of_allele_matched":round(
    100*sum(n for k,n in f.items() if k!="DELTA_LE_0.05")/m["MATCHED"],3),
  "AF_delta_gt_0_10_pct_of_allele_matched":round(
    100*m["AF_DELTA_GT_0.1"]/m["MATCHED"],3),
  "technical_only_QC":True,
  "phenotype_comparability_verified":False,
  "participant_nonoverlap_verified":False,
  "LDSC_MR_coloc_approved":False,
  "scientific_status":"HOLD"
 }
 return out
def main():
 p=argparse.ArgumentParser()
 p.add_argument("--input",required=True,type=Path)
 p.add_argument("--out",required=True,type=Path)
 args=p.parse_args()
 d=json.loads(args.input.read_text())
 out=validate(d)
 args.out.parent.mkdir(parents=True,exist_ok=True)
 args.out.write_text(json.dumps(out,indent=2,ensure_ascii=False)+"\n")
 print(json.dumps(out,ensure_ascii=False))
if __name__=="__main__":main()
