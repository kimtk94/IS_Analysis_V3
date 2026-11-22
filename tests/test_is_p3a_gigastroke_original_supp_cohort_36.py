"""Fail-closed fake original Suppl.Table1 rows: no network or real participants."""
import importlib
import unittest

M=importlib.import_module("scripts.is.audit_p3a_gigastroke_original_supp_cohort_36")


def source_rows():
    """Construct exactly 15 EAS body rows incl 1 incident-only row 60."""
    cohort=[
       ("KCPS2_ASN**",46,12366),
       ("BIOBANK_JAPAN*",17651,159610),
       ("HISAYAMA_FSR",1113,901),
       ("INTERSTROKE_ASN",222,698),
    ]
    # 10 CKBB discovery-control regional strata, all AIS cases NA
    # Exact regional controls unknown for synthetic; sum to 63,667.
    controls=[5000,5000,5000,5000,5000,5000,5000,5000,5000,18667]
    cohort+= [(f"CKBB_ASN_{i+1}",None,n) for i,n in enumerate(controls)]
    assert sum(x[1] or 0 for x in cohort)==19032
    assert sum(x[2] or 0 for x in cohort)==237242
    d={3:{"A":"EUR"},4:{"B":"EUR_COHORT_A","K":"62000","P":"1200000"},
       5:{"B":"EUR_COHORT_B","K":"50","P":"10000"},
       6:{"B":"EUR_COHORT_C","K":"30","P":"10000"},
       7:{"B":"EUR_COHORT_D","K":"20","P":"14808"},
       8:{"A":"TOTAL EUR","K":"62100","P":"1234808"},
       50:{"A":"EAS"}}
    # EAS has row 50..64 body: 4 AIS contributor, 10 CKBB, and row60
    rows=[x for x in range(50,65) if x!=60]
    for line,(name,cases,control) in zip(rows,cohort):
        d[line]={"B":name,"K":str(cases) if cases is not None else "NA",
                 "P":str(control)}
    d[50]["A"]="EAS"
    d[60]={"B":"INCIDENT_ONLY_NOT_DISCOVERY","K":"4","P":"102"}
    d[65]={"A":"TOTAL EAS","K":"19032","P":"237242"}
    return d


class PrimaryCohortTests(unittest.TestCase):
    def test_EAS_anchor_original_K_P_and_reconcile(self):
        a=M.extract_ancestry(source_rows(),"EAS")
        self.assertEqual((a["original_total_from_source_K_AIS_cases"],
                          a["original_total_from_source_P_discovery_controls"]),
                         (19032,237242))
        self.assertTrue(a["source_individual_cohort_rows_reconcile_total"])
        self.assertTrue(a["primary_total_matches_target_GWAS_case_control_counts"])

    def test_EAS_original_BBJ_is_not_independent(self):
        a=M.validate_eas(M.extract_ancestry(source_rows(),"EAS"))
        self.assertEqual(a["original_BBJ_AIS_cases"],17651)
        self.assertEqual(a["original_BBJ_shared_discovery_controls"],159610)
        self.assertTrue(a["BBJ_is_EAS_discovery_contributing_cohort"])
        self.assertTrue(a["BBJ_reuse_as_independent_EAS_replication_prohibited"])
        self.assertFalse(a["sample_level_independence_of_other_datasets_proven"])

    def test_AIS_distinct_from_any_stroke_source_column_J(self):
        rows=source_rows()
        rows[51]["J"]="222222"
        a=M.extract_ancestry(rows,"EAS")
        self.assertEqual(a["original_total_from_source_K_AIS_cases"],19032)

    def test_EAS_incident_only_row60_excluded_from_discovery(self):
        a=M.extract_ancestry(source_rows(),"EAS")
        self.assertEqual(a["recomputed_usable_AIS_cases"],19032)
        self.assertEqual(len([x for x in a["cohort_rows"]
            if x["source_row_is_incident_only_separate_from_AIS_meta"]]),1)

    def test_BBJ_case_count_source_drift_not_approved(self):
        rows=source_rows()
        rows[52]["K"]="17652"
        with self.assertRaises(ValueError):
            M.validate_eas(M.extract_ancestry(rows,"EAS"))

    def test_BBJ_controls_source_drift_not_approved(self):
        rows=source_rows()
        rows[52]["P"]="159611"
        with self.assertRaises(ValueError):
            M.validate_eas(M.extract_ancestry(rows,"EAS"))

    def test_wrong_total_case_discovery_hold(self):
        rows=source_rows();rows[65]["K"]="20000"
        with self.assertRaises(ValueError):
            M.validate_eas(M.extract_ancestry(rows,"EAS"))

    def test_CKBB_na_not_imputed_as_cases(self):
        a=M.validate_eas(M.extract_ancestry(source_rows(),"EAS"))
        self.assertEqual(a["AIS_discovery_source_crosswalk_status"],
                         "PASS_EAS_AIS_ORIGINAL_SUPP_COORT_SUMS_AND_BBJ")
        self.assertEqual(len([x for x in a["cohort_rows"]
              if x["published_study_name"].startswith("CKBB_ASN")]),10)

    def test_EUR_independent_row_sum_reconciles(self):
        a=M.extract_ancestry(source_rows(),"EUR")
        self.assertEqual(a["original_total_from_source_K_AIS_cases"],62100)
        self.assertEqual(a["original_total_from_source_P_discovery_controls"],1234808)
        self.assertTrue(a["source_individual_cohort_rows_reconcile_total"])

    def test_EUR_source_row_sum_discrepancy_flagged(self):
        rows=source_rows();rows[5]["P"]="34800"
        a=M.extract_ancestry(rows,"EUR")
        self.assertFalse(a["source_individual_cohort_rows_reconcile_total"])

    def test_no_EUR_anchor_fail_closed(self):
        rows=source_rows();del rows[3]
        with self.assertRaises(ValueError):
            M.extract_ancestry(rows,"EUR")

    def test_missing_original_total_is_not_accepted(self):
        rows=source_rows();rows[65].pop("P")
        with self.assertRaises(ValueError):
            M.extract_ancestry(rows,"EAS")

    def test_header_case_control_not_merged_with_PGS_validation(self):
        self.assertEqual(M.EXPECTED["EAS"],{"cases":19032,"controls":237242})
        self.assertEqual(M.EXPECTED["EUR"],{"cases":62100,"controls":1234808})

    def test_original_publisher_url_stable(self):
        self.assertIn("41586_2022_5165_MOESM4_ESM.xlsx",M.URLS[0])
        self.assertTrue(M.URLS[0].startswith("https://static-content.springer.com/"))


if __name__=="__main__":
    unittest.main()
