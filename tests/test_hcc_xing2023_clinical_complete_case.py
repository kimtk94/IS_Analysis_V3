"""Synthetic tests for Xing2023 complete-case eligibility; no patient data."""
import pathlib
import sys
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts" / "hcc"))
from probe_xing2023_clinical_complete_case import (
    COLUMN, PROFILE_DEFS, valid_value, profile_audit
)

def make_mock():
    rows = {}
    for n in range(152):
        r = {
            COLUMN["age"]: str(40+n % 35),
            COLUMN["sex"]: "F" if n % 2 else "M",
            COLUMN["AFP"]: str(n+10),
            COLUMN["HBV_marker"]: "0" if n % 2 else "1",
            COLUMN["HCV_marker"]: "0" if n % 3 else "1",
            COLUMN["MVI"]: "0" if n % 3 else "1",
            COLUMN["tumor_size"]: "3.1",
            COLUMN["TNM_descriptor"]: "T1",
            COLUMN["TNM_stage"]: "I",
            COLUMN["BCLC"]: "B",
            COLUMN["RFS_candidate_code"]: "0" if n % 4 else "1",
            COLUMN["RFS_candidate_duration"]: "4.5",
            COLUMN["survival_candidate_code"]: "1" if n % 3 else "0",
            COLUMN["survival_candidate_duration"]: "24"
        }
        if n < 5:
            r[COLUMN["HCV_marker"]] = "."
        if 5 <= n < 12:
            r[COLUMN["AFP"]] = "positive"
        rows[f"synthetic-{n}"] = r
    return rows

class Tests(unittest.TestCase):
    def test_profiles_are_distinct_and_not_claimable(self):
        r = profile_audit(make_mock())
        p = r["profiles"]
        self.assertEqual(r["n_matched_structurally"],152)
        self.assertEqual(p["endpoint_only_RFS"]["candidate_complete_case_n"],152)
        self.assertEqual(p["baseline_RFS"]["candidate_complete_case_n"],152)
        self.assertEqual(p["clinical_extended_RFS"]["candidate_complete_case_n"],145)
        self.assertEqual(p["viral_marker_coverage"]["candidate_complete_case_n"],147)
        self.assertEqual(p["clinical_extended_viral_RFS"]["candidate_complete_case_n"],140)
        self.assertEqual(p["clinical_extended_viral_survival"]["candidate_complete_case_n"],140)
        self.assertEqual(r["fields"]["AFP"]["valid_n"],145)
        self.assertEqual(r["fields"]["HCV_marker"]["valid_n"],147)
        self.assertFalse(r["prognostic_coding_reviewed"])
        self.assertFalse(r["survival_model_allowed"])
        self.assertFalse(r["patient_ids_or_values_logged"])
        self.assertEqual(
            p["clinical_extended_viral_RFS"]["uninterpreted_event_code_0_n"]+
            p["clinical_extended_viral_RFS"]["uninterpreted_event_code_1_n"], 140)

    def test_binary_and_duration_fail_closed(self):
        r=make_mock()
        one=next(iter(r.values()))
        one[COLUMN["RFS_candidate_code"]]="2"
        self.assertFalse(valid_value(one,"RFS_candidate_code"))
        one[COLUMN["RFS_candidate_code"]]="0"
        one[COLUMN["RFS_candidate_duration"]]="-1"
        self.assertFalse(valid_value(one,"RFS_candidate_duration"))
        one[COLUMN["RFS_candidate_duration"]]="1"
        self.assertTrue(valid_value(one,"RFS_candidate_duration"))

    def test_missing_source_column_is_not_complete(self):
        r=make_mock()
        one=next(iter(r.values()))
        one.pop(COLUMN["TNM_stage"])
        self.assertFalse(valid_value(one,"TNM_stage"))
        rep=profile_audit(r)
        self.assertEqual(rep["profiles"]["clinical_extended_RFS"]["candidate_complete_case_n"],144)

    def test_profile_is_fixed_predeclared_schema(self):
        self.assertEqual(len(PROFILE_DEFS),9)
        self.assertIn("AFP",PROFILE_DEFS["clinical_extended_RFS"])
        self.assertNotIn("HBV_marker",PROFILE_DEFS["baseline_RFS"])

if __name__=="__main__":
    unittest.main()
