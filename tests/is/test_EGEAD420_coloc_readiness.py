#!/usr/bin/env python3
"""Only synthetic source-evidence fixtures; cannot produce false coloc PASS."""
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"scripts"/"is"))
from audit_EGEAD420_coloc_readiness import audit_inputs,run,SOURCE_STATUS,COHORT_STATUS


def fixture():
    g={
        "source_status":SOURCE_STATUS,
        "published_whole_ZIP_MD5_verified_by_full_download":False,
        "nominal_source_build":"GRCh38",
        "source_GWAS_build":"GIGASTROKE GRCh37",
        "first_actual_nominal_association_records":{
            "source_nominal_file":"CD16p_Mono_nominal.txt",
            "whitespace_delimited_columns_per_line":13,
            "actual_ImmuNexUT_13_column_source_data_QC":{
                "source_13col_header_verified":True,
                "n_actual_complete_source_association_rows":79,
                "source_nominal_non_significant_p_gt_0_05":75,
                "source_nominal_significant_p_le_0_05":4,
                "source_first_gene_id":"ENSG00000227232.5",
                "source_reported_cis_variants_for_first_gene":835,
                "slope_beta_present":True,
                "beta_effect_allele":"ALT",
                "ref_alt_alleles_present":True,
                "slope_se_column_present":False,
                "allele_frequency_MAF_EAF_present":False,
                "per_celltype_sample_n_column_present":False,
                "full_cis_association_completeness":"NOT_VERIFIED",
            }
        }
    }
    c={
        "status":COHORT_STATUS,
        "original_SDRF_publisher_MD5_verified":True,
        "aggregate_evidence":{
            "original_public_SDRF_distinct_source_names":416,
            "healthy_control_WGS_records":79,
            "immune_disease_WGS_records":337,
            "per_celltype_effective_eQTL_N_observed":False,
        },
    }
    return g,c


class ColocReadinessTest(unittest.TestCase):
    def test_genuine_source_conventions_yield_hold(self):
        source,cohort=fixture()
        result=audit_inputs(source,cohort)
        self.assertEqual(result["scientific_decision"],"SOURCE_ROW_QC_PASS_COLOC_HOLD")
        self.assertEqual(result["coloc_ABF_pvalue_MAF_N_route"],
                         "BLOCKED_MISSING_MAF_AND_CELLTYPE_N")
        self.assertEqual(result["original_nominal_association_rows_sampled"],79)

    def test_raw_partial_zip_not_full_source_md5(self):
        source,cohort=fixture()
        source["published_whole_ZIP_MD5_verified_by_full_download"]=True
        with self.assertRaisesRegex(ValueError,"43GB"):
            audit_inputs(source,cohort)

    def test_unverified_sample_n_is_not_allowed(self):
        source,cohort=fixture()
        cohort["aggregate_evidence"]["per_celltype_effective_eQTL_N_observed"]=True
        with self.assertRaisesRegex(ValueError,"per-cell n"):
            audit_inputs(source,cohort)

    def test_inconsistent_significance_groups_rejected(self):
        source,cohort=fixture()
        source["first_actual_nominal_association_records"]["actual_ImmuNexUT_13_column_source_data_QC"]["source_nominal_non_significant_p_gt_0_05"]=80
        with self.assertRaisesRegex(ValueError,"P-value"):
            audit_inputs(source,cohort)

    def test_invalid_qtl_source_build_rejected(self):
        source,cohort=fixture()
        source["nominal_source_build"]="GRCh37"
        with self.assertRaisesRegex(ValueError,"build"):
            audit_inputs(source,cohort)

    def test_missing_actual_source_artifact_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            with self.assertRaisesRegex(ValueError,"not downloaded"):
                run(root/"source",root/"out")

    def test_artifact_output_retains_rigor(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            source,cohort=fixture()
            a=root/"evidence"/"source";b=root/"evidence"/"cohort"
            a.mkdir(parents=True);b.mkdir()
            (a/"E_GEAD_420_source_zip64_evidence.json").write_text(json.dumps(source))
            (b/"E_GEAD_420_cohort_metadata.json").write_text(json.dumps(cohort))
            result=run(root/"evidence",root/"results")
            self.assertTrue((root/"results"/"EGEAD420_reproducible_go_hold.json").is_file())
            self.assertEqual(result["coloc_ABF_beta_varbeta_route"],"BLOCKED_MISSING_VARBETA")

if __name__=="__main__":
    unittest.main()
