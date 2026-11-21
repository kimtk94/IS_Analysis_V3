"""Offline reference-overlap and exact duplicate fixture tests (no internet)."""
import gzip
import hashlib
import io
import unittest
from unittest.mock import patch

from scripts.ckm_stroke.audit_bmi_hm3_overlap import (
    HM3_URL, BMI_URL, audit, hm3_catalog, inspect_bmi, orientation, allele_pair,
)


class Response(io.BytesIO):
    def __init__(self, data: bytes, url: str, length=None, status=200):
        super().__init__(data)
        self.url = url
        self.status = status
        self.headers = {"Content-Length": str(len(data) if length is None else length)}

    def geturl(self):
        return self.url


REF = (
    b"SNP\tA1\tA2\n"
    b"rs1\tC\tA\n"
    b"rs2\tA\tT\n"
    b"rs3\tT\tG\n"
    b"rs4\tA\tC\n"
)
BMI_ROWS = [
    "1\t100\trs1\tA\tC\t0.1\t0.1\t0.05\t0.1\t600000",
    "1\t100\trs1\tA\tC\t0.1\t0.1\t0.05\t0.1\t600000",
    "1\t101\trs2\tA\tT\t0.1\t0.1\t0.05\t0.1\t600000",
    "1\t102\trs3\tA\tC\t0.1\t0.1\t0.05\t0.1\t600000",
    "1\t103\trs4\tA\tG\t0.1\t0.1\t0.05\t0.1\t600000",
]
BMI_HEADER = "CHR\tPOS\tSNP\tTested_Allele\tOther_Allele\tFreq_Tested_Allele_in_HRS\tBETA\tSE\tP\tN\n"
BMI = gzip.compress((BMI_HEADER + "\n".join(BMI_ROWS) + "\n").encode())


def opener(req, timeout):
    if req.full_url == HM3_URL:
        return Response(REF, HM3_URL)
    if req.full_url == BMI_URL:
        return Response(BMI, BMI_URL)
    raise AssertionError("unexpected test fetch")


class OverlapTests(unittest.TestCase):
    def test_orientation_classes(self):
        self.assertEqual(orientation("A", "C", allele_pair("C", "A")), "direct_alleles")
        self.assertEqual(orientation("A", "T", allele_pair("A", "T")), "palindromic_ambiguous")
        self.assertEqual(orientation("C", "A", allele_pair("G", "T")), "strand_complement_possible")
        self.assertEqual(orientation("A", "G", allele_pair("A", "C")), "allele_mismatch")

    def test_reference_parser_and_digest(self):
        ref, metadata = hm3_catalog(opener)
        self.assertEqual(len(ref), 4)
        self.assertEqual(metadata["ref_sha256"], hashlib.sha256(REF).hexdigest())
        self.assertEqual(metadata["source_status"], "REF_BYTES_AND_SCHEMA_PASS")

    def test_duplicate_counts_and_hm3_compatibility(self):
        ref, _ = hm3_catalog(opener)
        with patch("scripts.ckm_stroke.audit_bmi_hm3_overlap.EXPECTED_BMI_SHA256",
                   hashlib.sha256(BMI).hexdigest()):
            output = inspect_bmi(ref, opener)
        counts = output["counts"]
        self.assertEqual(counts["rows"], 5)
        self.assertEqual(counts["duplicate_rsids"], 1)
        self.assertEqual(counts["duplicate_coordinates"], 1)
        self.assertEqual(counts["duplicate_allele_keys"], 1)
        self.assertEqual(counts["hm3_rsids_matched_rows"], 5)
        self.assertEqual(counts["hm3_unique_rsids_intersection"], 4)
        self.assertEqual(counts["direct_alleles"], 2)
        self.assertEqual(counts["palindromic_ambiguous"], 1)
        self.assertEqual(counts["strand_complement_possible"], 1)
        self.assertEqual(counts["allele_mismatch"], 1)

    def test_conflicted_reference_rsid_is_excluded_from_allele_matching(self):
        repeated = REF + b"rs1\tG\tT\n"
        def corrupted(req, timeout):
            if req.full_url == HM3_URL:
                return Response(repeated, HM3_URL)
            return opener(req, timeout)
        ref, meta = hm3_catalog(corrupted)
        self.assertNotIn("rs1", ref)
        self.assertEqual(meta["conflicted_unique_rsids_excluded"], 1)
        self.assertEqual(meta["source_status"], "PASS_REFERENCE_WITH_CONFLICT_EXCLUSIONS")
        with patch("scripts.ckm_stroke.audit_bmi_hm3_overlap.EXPECTED_BMI_SHA256",
                   hashlib.sha256(BMI).hexdigest()):
            result = inspect_bmi(ref, opener)
        self.assertEqual(result["counts"]["hm3_rsids_matched_rows"], 3)
        self.assertEqual(result["counts"]["hm3_unique_rsids_intersection"], 3)

    def test_bad_reference_header_rejected(self):
        def bad(req, timeout):
            if req.full_url == HM3_URL:
                return Response(b"RSID\tA1\tA2\nrs1\tA\tC\n", HM3_URL)
            return opener(req, timeout)
        with self.assertRaisesRegex(ValueError, "header"):
            hm3_catalog(bad)

    def test_fails_closed_to_scientific_hold(self):
        with patch("scripts.ckm_stroke.audit_bmi_hm3_overlap.EXPECTED_BMI_SHA256",
                   hashlib.sha256(BMI).hexdigest()):
            result = audit(opener)
        self.assertFalse(result["source_approved_for_ldsc"])
        self.assertEqual(result["scientific_status"], "HOLD")
        self.assertFalse(result["ldsc_h2_rg_performed"])

    def test_bmi_sha_mismatch_stops(self):
        ref, _ = hm3_catalog(opener)
        with self.assertRaisesRegex(ValueError, "SHA256 changed"):
            inspect_bmi(ref, opener)


if __name__ == "__main__":
    unittest.main()
