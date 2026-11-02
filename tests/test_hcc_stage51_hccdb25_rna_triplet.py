#!/usr/bin/env python3
"""Synthetic-only tests. No real identifiers, expression or clinical records."""
import pathlib
import sys
import tempfile
import unittest
sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[1]/"scripts"/"hcc"))
from audit_stage51_hccdb25_rna_triplet import (
    CODE_URL,EXPECTED,evaluate_triplet,parse_expression,parse_sample_transposed,
    parse_patient,source_code_reference,review
)

def create_triplet(root,expr=None,sample=None,clinical=None):
    a=expr if expr is not None else "Index\tSymbol\tS01\tS02\n1\tA1BG\t1.5\t2.5\n2\tABO\t3.2\t4.5\n"
    b=sample if sample is not None else "SAMPLE_ID\tS01\tS02\nPATIENT_ID\tP01\tP02\nTYPE\tHCC\tHCC\n"
    c=clinical if clinical is not None else "PATIENT_ID\tOS\tOS_STATUS\nP01\t12\tAlive\nP02\t14\tDead\n"
    for key,value in (("expression",a),("sample",b),("patient",c)):
        (root/EXPECTED[key]).write_text(value,encoding="utf-8")

class MockResponse:
    status=200
    def __init__(self,data):self.data=data
    def geturl(self):return CODE_URL
    def read(self,n):return self.data[:n]
    def __enter__(self):return self
    def __exit__(self,*args):return False

class Stage51Tests(unittest.TestCase):
    def test_fixture_success_is_still_source_hold(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=pathlib.Path(tmp)
            create_triplet(root)
            x=evaluate_triplet(root)
            self.assertEqual(x["status"],"TRIPLET_INTERNAL_CONSISTENCY_PASS_SOURCE_HOLD")
            self.assertEqual(x["aggregate"]["genes_n"],2)
            self.assertEqual(x["aggregate"]["tumor_HCC_sample_metadata_n"],2)
            self.assertTrue(x["internal_RNA_patient_mapping_consistent"])
            self.assertFalse(review(x)["original_RNA_source_authorized"])
            self.assertTrue(review(x)["G30_D"].startswith("STOP"))
            self.assertEqual(len(x["source_file_metadata"]["expression"]["file_sha256"]),64)
    def test_missing_triplet(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.assertEqual(evaluate_triplet(tmp)["status"],"TRIPLET_SOURCE_FILES_MISSING")
    def test_duplicate_expression_sample_stops(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=pathlib.Path(tmp)
            create_triplet(root,expr="Symbol\tS01\tS01\nA1BG\t1\t2\nABO\t3\t4\n")
            x=evaluate_triplet(root)
            self.assertEqual(x["status"],"TRIPLET_VALIDATION_HOLD")
            self.assertEqual(x["validation_error_code"],"EXPRESSION_DUPLICATE_OR_MISSING_SAMPLE_IDS")
    def test_duplicate_patient_id_stops(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=pathlib.Path(tmp)
            create_triplet(root,clinical="PATIENT_ID\tOS\nP01\t12\nP01\t20\n")
            x=evaluate_triplet(root)
            self.assertEqual(x["status"],"TRIPLET_VALIDATION_HOLD")
            self.assertEqual(x["validation_error_code"],"PATIENT_DUPLICATE_PATIENT_ID")
    def test_clinical_join_missing_stops(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=pathlib.Path(tmp)
            create_triplet(root,clinical="PATIENT_ID\tOS\nP01\t12\nP03\t20\n")
            x=evaluate_triplet(root)
            self.assertEqual(x["status"],"TRIPLET_INTERNAL_CROSSWALK_HOLD")
            self.assertEqual(x["aggregate"]["sample_patient_ids_missing_clinical_n"],1)
    def test_tissue_HCC_not_assumed(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=pathlib.Path(tmp)
            create_triplet(root,sample="SAMPLE_ID\tS01\tS02\nPATIENT_ID\tP01\tP02\nTYPE\tNormal\tNormal\n")
            x=evaluate_triplet(root)
            self.assertEqual(x["status"],"TRIPLET_INTERNAL_CROSSWALK_HOLD")
    def test_missing_sample_metadata_stops(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=pathlib.Path(tmp)
            create_triplet(root,sample="SAMPLE_ID\tS01\tS03\nPATIENT_ID\tP01\tP02\nTYPE\tHCC\tHCC\n")
            x=evaluate_triplet(root)
            self.assertEqual(x["status"],"TRIPLET_INTERNAL_CROSSWALK_HOLD")
            self.assertEqual(x["aggregate"]["expression_samples_missing_metadata_n"],1)
    def test_non_numeric_and_nan_expression_stops(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=pathlib.Path(tmp)
            create_triplet(root,expr="Symbol\tS01\tS02\nA1BG\t1.0\t2.0\nABO\tNaN\tInfinity\n")
            x=evaluate_triplet(root)
            self.assertEqual(x["validation_error_code"],"EXPRESSION_VALUE_NONFINITE")
    def test_duplicate_metadata_key_stops(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=pathlib.Path(tmp)
            create_triplet(root,sample="SAMPLE_ID\tS01\tS02\nPATIENT_ID\tP01\tP02\nTYPE\tHCC\tHCC\nTYPE\tHCC\tHCC\n")
            x=evaluate_triplet(root)
            self.assertEqual(x["validation_error_code"],"SAMPLE_DUPLICATE_OR_EMPTY_METADATA_KEYS")
    def test_research_code_reference_not_source_file(self):
        source=b'# OEP000321\n' + b'\n'.join(k.encode() for k in EXPECTED.values())
        x=source_code_reference(opener=lambda req,timeout=None:MockResponse(source))
        self.assertEqual(x["status"],"RESEARCH_CODE_LOCAL_TRIPLET_REFERENCES_VERIFIED")
        self.assertFalse(x["actual_RNA_file_download_verified"])
        self.assertFalse(x["public_download_endpoint_verified"])
    def test_reference_without_required_file_is_hold(self):
        x=source_code_reference(opener=lambda req,timeout=None:MockResponse(b'OEP000321 HCCDB25.sample.txt'))
        self.assertEqual(x["status"],"REFERENCE_INCOMPLETE_HOLD")
    def test_biology_claim_upgrade_prohibited(self):
        with self.assertRaises(ValueError):
            review({"source_authorized_hccdb25_download":True,
                    "verified_PDC_official_aliquot_bridge":False})
    def test_missing_patient_field_stops(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=pathlib.Path(tmp)
            create_triplet(root,clinical="ID\tOS\nP01\t12\nP02\t19\n")
            x=evaluate_triplet(root)
            self.assertEqual(x["validation_error_code"],"PATIENT_TABLE_PATIENT_ID_REQUIRED")

if __name__=="__main__":unittest.main()
