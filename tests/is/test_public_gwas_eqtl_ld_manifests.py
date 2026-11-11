#!/usr/bin/env python3
"""Offline tests for controlled external source metadata audit."""
import sys
import unittest
from pathlib import Path
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts" / "is"))
from audit_public_gwas_eqtl_ld_manifests import scan_yaml_text, GWAS, ZENODO


class ManifestAuditTests(unittest.TestCase):
    def test_exact_gwas_ischemic_stroke_ids(self):
        self.assertEqual(GWAS, {"EUR_IS": "GCST90104540", "EAS_IS": "GCST90104545"})

    def test_eas_eqtl_repository_record_is_distinct_from_ld_panel(self):
        self.assertIn("/21296030", ZENODO["MAEEA_EAS_CIS_EQTL"])
        self.assertIn("/7768714", ZENODO["ONE_KG_GRCH37_EUR_EAS_REFERENCE"])

    def test_yaml_extracts_source_checksums(self):
        sample = ("gwas_id: GCST90104540\n"
                  "genome_assembly: GRCh37\n"
                  "data_file_name: GCST90104540_buildGRCh37.tsv.gz\n"
                  "data_file_md5sum: 00000000000000000000000000000000\n"
                  "sample_size: 1296908\n"
                  "case_count: 62100\n"
                  "control_count: 1234808\n")
        d = scan_yaml_text(sample)
        self.assertEqual(d["gwas_id"], "GCST90104540")
        self.assertEqual(d["data_file_md5sum"], "00000000000000000000000000000000")
        self.assertEqual(d["genome_assembly"], "GRCh37")

    def test_missing_yaml_is_not_silent_pass(self):
        self.assertEqual(scan_yaml_text("malformed metadata"), {})

    def test_yml_scalar_quotes(self):
        d = scan_yaml_text('gwas_id: "GCST90104545"\n')
        self.assertEqual(d["gwas_id"], "GCST90104545")


if __name__ == "__main__":
    unittest.main()
