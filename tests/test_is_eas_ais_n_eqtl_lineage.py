import importlib.util
import unittest
from pathlib import Path

PATH=Path(__file__).resolve().parents[1]/"scripts/is/audit_is_eas_ais_n_eqtl_lineage.py"
spec=importlib.util.spec_from_file_location("is_n_lineage",PATH)
m=importlib.util.module_from_spec(spec)
spec.loader.exec_module(m)

YAML_BASE=b"""gwas_id: GCST90104545
data_file_name: GCST90104545_buildGRCh37.tsv.gz
data_file_md5sum: 904e12e33ec3c2eeefe4dc93bec4995c
genome_assembly: GRCh37
file_type: pre-GWAS-SSF
is_harmonised: false
"""
STUDIES=b"""<html><body>List of studies
QTS000015 GTEx (v8) 49 tissues 838 donors
Imported studies GTEx_v8 49 tissues 948 donors independent workflows
</body></html>"""
RELEASES=b"""<html><body>eQTL Catalogue release 6 April 2023
complete set of cis-eQTL summary statistics .all.tsv.gz
eQTL Catalogue r8 pre-release September 2026
QTS000015 GTEx_v10 dataset ids have remained the same
</body></html>"""

class LineageTests(unittest.TestCase):
    def fake(self,url):
        return {m.YAML_URL:YAML_BASE,m.EQTL_STUDIES:STUDIES,
                m.EQTL_RELEASE:RELEASES,m.QTL_URL:b"\x1f\x8b\x08\x00\x00\x00"}[url]

    def test_missing_case_control_fields_hold(self):
        x=m.audit(self.fake)
        self.assertEqual(x["gwas_yaml"]["source_file_identity_and_yaml_digest"],"PASS")
        self.assertIsNone(x["gwas_yaml"]["illustrative_s_from_yaml_not_adopted"])
        self.assertFalse(x["gwas_yaml"]["association_level_N_or_case_fraction_approved"])
        self.assertIn("HOLD_MISSING",x["gwas_yaml"]["N_gate"])

    def test_consistent_yaml_n_still_review_only(self):
        raw=YAML_BASE+b"sample_size: 256274\ncase_count: 19032\ncontrol_count: 237242\n"
        result=m.inspect_gwas_yaml(raw)
        self.assertEqual(result["N_gate"],"YAML_N_PRESENT_REQUIRES_VARIANT_AND_META_ANALYSIS_REVIEW")
        self.assertTrue(result["published_reference_matches_yaml"])
        self.assertAlmostEqual(result["illustrative_s_from_yaml_not_adopted"],19032/256274)
        self.assertFalse(result["coloc_approved"])

    def test_conflicting_N_does_not_promote(self):
        x=m.inspect_gwas_yaml(YAML_BASE+b"sample_size: 999\ncase_count: 19032\ncontrol_count: 237242\n")
        self.assertFalse(x["yaml_internal_N_consistency"])
        self.assertIsNone(x["illustrative_s_from_yaml_not_adopted"])
        self.assertFalse(x["association_level_N_or_case_fraction_approved"])

    def test_untrusted_gwas_md5_keeps_hold(self):
        x=m.inspect_gwas_yaml(YAML_BASE.replace(b"904e12e33ec3c2eeefe4dc93bec4995c",b"0"*32)+
            b"sample_size: 256274\ncase_count: 19032\ncontrol_count: 237242\n")
        self.assertEqual(x["source_file_identity_and_yaml_digest"],"HOLD")
        self.assertNotEqual(x["N_gate"],"YAML_N_PRESENT_REQUIRES_VARIANT_AND_META_ANALYSIS_REVIEW")

    def test_invalid_case_count_keeps_hold(self):
        x=m.inspect_gwas_yaml(YAML_BASE+b"case_count: N/A\ncontrol_count: 237242\n")
        self.assertFalse(x["file_specific_case_control_fields_available"])
        self.assertIsNone(x["yaml_case_count"])

    def test_conflicting_duplicate_yaml_key_raises(self):
        with self.assertRaises(ValueError):
            m.inspect_gwas_yaml(YAML_BASE+b"gwas_id: GCST90104544\n")

    def test_official_studies_distinguish_original_imported(self):
        x=m.inspect_qtl_metadata(STUDIES,RELEASES)
        self.assertTrue(x["official_reprocessed_GTEx_listing"])
        self.assertTrue(x["official_separate_imported_GTEx_v8_listing"])
        self.assertTrue(x["official_release6_complete_cis_documented"])
        self.assertTrue(x["GTEx_v10_ID_reuse_warning_observed"])
        self.assertFalse(x["current_QTD000161_is_original_GTExv8_FastQTL"])
        self.assertFalse(x["original_GTExv8_cis_results_verified"])

    def test_missing_imported_original_label_does_not_upgrade(self):
        r=m.inspect_qtl_metadata(STUDIES.replace(b"Imported studies GTEx_v8",b"Other imported files"),
                                 RELEASES)
        self.assertEqual(r["QTL_gate"],"HOLD_OFFICIAL_QTL_LINEAGE_INCOMPLETE")

    def test_bad_network_keeps_each_gate_hold(self):
        def bad(url):raise OSError("network")
        x=m.audit(bad)
        self.assertEqual(x["gwas_yaml"]["N_gate"],"HOLD_GWAS_YAML_UNAVAILABLE")
        self.assertEqual(x["qtl_provenance"]["QTL_gate"],"HOLD_OFFICIAL_QTL_METADATA_UNAVAILABLE")
        self.assertFalse(x["coloc_performed"])

    def test_qtl_prefix_unavailable_does_not_erase_documentation(self):
        def q(url):
            if url==m.QTL_URL:raise OSError("not found")
            return self.fake(url)
        x=m.audit(q)
        self.assertEqual(x["qtl_provenance"]["qtl_prefix_gate"],"HOLD_QTL_PREFIX_UNAVAILABLE")
        self.assertEqual(x["qtl_provenance"]["QTL_gate"],
                          "REPROCESSED_VS_IMPORTED_SOURCE_DISTINGUISHED_ORIGINAL_HOLD")

    def test_url_allowlist(self):
        with self.assertRaises(ValueError):
            m.metadata_fetch("https://untrusted.example/eqtl")
        with self.assertRaises(ValueError):
            m.metadata_fetch("https://ftp.ebi.ac.uk/elsewhere")

if __name__=="__main__":
    unittest.main()
