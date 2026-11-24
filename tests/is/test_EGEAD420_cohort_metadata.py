#!/usr/bin/env python3
"""Synthetic privacy and cohort-structure QC for public E-GEAD-420 SDRF."""
import csv
import io
from pathlib import Path
import sys
import unittest

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"scripts"/"is"))
from audit_EGEAD420_cohort_metadata import cohort_summary,EXPECTED,EXPECTED_HEALTHY,FIELD_KEYS


def fixture(healthy=79,missing=0):
    text=io.StringIO()
    writer=csv.DictWriter(text,fieldnames=list(FIELD_KEYS),delimiter="\t",lineterminator="\n")
    writer.writeheader()
    for i in range(416):
        d={
            "Source Name":f"SOURCE_TEST_{i}",
            "Characteristics[disease]":"healthy control" if i<healthy else "Systemic Lupus Erythematosus",
            "Characteristics[sex]":"female" if i%2==0 else "male",
            "Characteristics[age]":"50",
            "Characteristics[Phase]":"1" if i%2==0 else "2",
            "Derived Array Data Matrix File":"nominal_eQTL_2.tar.gz",
            "Comment[Derived Array Data Matrix File md5]":"41544fb671177aeee7fd2d0634b90ae2"
        }
        if missing and i==1:
            d["Characteristics[sex]"]="unknown"
        writer.writerow(d)
    return text.getvalue().encode()


class CohortAuditTests(unittest.TestCase):
    def test_all_416_distinct_and_79_healthy(self):
        result=cohort_summary(fixture())
        self.assertEqual(result["original_public_SDRF_wgs_sample_records"],416)
        self.assertEqual(result["healthy_control_WGS_records"],79)
        self.assertEqual(result["immune_disease_WGS_records"],337)
        self.assertTrue(result["no_participant_identifiers_exported"])
        self.assertFalse(result["per_celltype_effective_eQTL_N_observed"])
        self.assertNotIn("SOURCE_TEST_",repr(result))

    def test_missing_healthy_count_rejected(self):
        with self.assertRaisesRegex(ValueError,"healthy"):
            cohort_summary(fixture(healthy=78))

    def test_wrong_number_of_samples_rejected(self):
        rows=fixture().decode().splitlines()
        with self.assertRaisesRegex(ValueError,"SDRF WGS subject rows"):
            cohort_summary(("\n".join(rows[:-1])+"\n").encode())

    def test_invalid_mismatch_original_qtl_md5_fails(self):
        raw=fixture()
        raw=raw.replace(b"41544fb671177aeee7fd2d0634b90ae2",b"00000000000000000000000000000000",1)
        with self.assertRaisesRegex(ValueError,"MD5"):
            cohort_summary(raw)

    def test_missing_sex_is_not_false_complete(self):
        x=cohort_summary(fixture(missing=1))
        self.assertEqual(x["fields_with_explicit_missing_records"]["sex"],1)

    def test_source_labels_not_stroke_patient_samples(self):
        x=cohort_summary(fixture())
        self.assertFalse(x["source_study_is_stroke_cohort"])
        self.assertFalse(x["total_sample_size_416_valid_per_each_cell_type"]) if "total_sample_size_416_valid_per_each_cell_type" in x else None

if __name__=="__main__":
    unittest.main()
