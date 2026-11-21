"""Fail-closed offline tests for alternative stroke GWAS sample N provenance."""
from __future__ import annotations
import unittest
from scripts.ckm_stroke.audit_alternative_stroke_gwas_neff_sources import (
    MEGA_COLS,FINN_COLS,FINN_R13_MANIFEST,
    check_mega_description,check_finn_schema,check_finn_access,
    classify_finn_manifest,review_sources,
)

def mega_doc():
    return (
      "<html><body>MEGASTROKE Europeans 40,585 cases 406,111 controls "
      "We present results for SNPs passing the MEGASTROKE filters "
      "(n_cases > 50% and oevar_imp > 0.5). "
      + " ".join(x+": value" for x in MEGA_COLS)+"</body></html>"
    )

def finn_schema():
    return (
      "# Summary association statistics\n"
      "GWAS summary statistics genome build GRCh38 estimated using regenie\n"
      "| Column name | Description |\n| --- | --- |\n"
      + "".join("| **\x60"+x+"\x60** | value |\n" for x in FINN_COLS)
      + "## Fine-mapping results\n"
      "| n_cases | a fine-mapping unrelated value |\n"
    )

def finn_access():
    return ("Release 13\nTo download summary statistics fill the online form.\n"
            + "https://storage.googleapis.com/finngen-public-data-r13/summary_stats/finngen_R13_manifest.tsv")

def manifest():
    return ("phenocode\tphenostring\tsample_info\n"
            "I9_STR\tStroke\ttext\n"
            "I9_CEREBINFARCT\tCerebral infarction\tmore\n"
            "I9_CHD\tCoronary heart disease\tmore\n")

class OriginalSourceGateTests(unittest.TestCase):
    def test_megastroke_seven_columns_does_not_promote_n_cases_filter(self):
        r=check_mega_description(mega_doc())
        self.assertEqual(len(r["documented_summary_columns"]),7)
        self.assertFalse(r["public_per_variant_Neff_column_documented"])
        self.assertTrue(r["variant_level_n_cases_used_for_filter_only"])
        self.assertEqual(r["EUR_cases"],40585)
    def test_missing_meg_hash_schema_fails(self):
        doc=mega_doc().replace("StdErr:","StdError:")
        with self.assertRaisesRegex(ValueError,"StdErr"):
            check_mega_description(doc)
    def test_finn_schema_is_grch38_and_no_neff(self):
        r=check_finn_schema(finn_schema())
        self.assertEqual(len(r["verified_summary_columns"]),13)
        self.assertFalse(r["variant_N_or_Neff_column_present"])
        self.assertEqual(r["reference_genome"],"GRCh38")
    def test_neff_is_not_allowed_if_it_appears_in_summary_header(self):
        doc=finn_schema().replace("## Fine-mapping results",
            "| **\x60n_eff\x60** | newly provided |\n## Fine-mapping results")
        with self.assertRaisesRegex(ValueError,"unexpected"):
            check_finn_schema(doc)
    def test_finn_access_must_contain_current_release(self):
        self.assertTrue(check_finn_access(finn_access())["official_download_requires_request_or_access_step"])
        with self.assertRaisesRegex(ValueError,"release"):
            check_finn_access(finn_access().replace("Release 13","Release 12"))
    def test_manifest_stroke_endpoint_names_are_not_validated_diagnoses(self):
        r=classify_finn_manifest(manifest())
        self.assertEqual(r["likely_stroke_endpoint_codes"],["I9_CEREBINFARCT","I9_STR"])
        self.assertFalse(r["per_variant_N_eff_verified_from_manifest"])
        self.assertTrue(r["endpoint_code_match_is_not_validated_stroke_case_definition"])
    def test_hemorrhagic_not_promoted_to_ischemic_replication(self):
        synthetic=("phenocode\tphenostring\tsample_info\n"
                   "I9_STR\tStroke unspecified\tna\n"
                   "I9_STR_EMBOLIC\tEmbolic stroke\tna\n"
                   "I9_ICH\tIntracerebral hemorrhage\tna\n"
                   "I9_SAH\tSubarachnoid hemorrhage\tna\n")
        # Interpret fixture escape sequences as TSV, preserving source parser.
        r=classify_finn_manifest(synthetic.replace("\\\\t","\\t").replace("\\\\n","\\n"))
        self.assertEqual(r["ischemic_stroke_candidate_codes_require_definition_review"],
                         ["I9_STR","I9_STR_EMBOLIC"])
        self.assertEqual(r["hemorrhagic_endpoint_codes_EXCLUDED_from_ischemic_stroke"],
                         ["I9_ICH","I9_SAH"])
        self.assertTrue(r["no_verified_LAA_or_SVD_specific_endpoint_claim"])
        self.assertTrue(all(x["clinical_case_definition_verified"] is False
                            for x in r["likely_stroke_endpoint_labels"]))

    def test_mutated_manifest_header_fails_closed(self):
        with self.assertRaisesRegex(ValueError,"phenotype code"):
            classify_finn_manifest("wrong\tlabel\nI9_STR\tStroke\n")
        with self.assertRaisesRegex(ValueError,"row width"):
            classify_finn_manifest("phenocode\tname\nI9_STR\n")
    def test_source_review_with_manifest_success_keeps_neff_hold(self):
        r=review_sources(mega_doc(),finn_schema(),finn_access(),manifest())
        self.assertEqual(r["technical_status"],"PASS_OFFICIAL_CONSORTIUM_SCHEMA_NEFF_GAP_REVIEWED")
        self.assertEqual(r["findings"]["verified_alternative_outcomes_satisfying_primary_LDSC_Neff_criterion"],0)
        self.assertFalse(r["findings"]["alternative_per_variant_N_eff_problem_solved"])
        self.assertEqual(r["FinnGen_R13"]["manifest"]["likely_stroke_endpoint_count"],2)
        self.assertFalse(r["ready_for_ldsc"])
        self.assertFalse(r["LDSC_h2_intercept_rg_MR_coloc_executed"])
    def test_manifest_unavailable_keeps_source_schema_pass_but_no_endpoint_claim(self):
        r=review_sources(mega_doc(),finn_schema(),finn_access(),
                         manifest_error="HTTPError: 403")
        self.assertEqual(r["FinnGen_R13"]["manifest"]["manifest_status"],
                         "NOT_RETRIEVED_OR_ACCESS_RESTRICTED")
        self.assertEqual(r["FinnGen_R13"]["manifest"]["likely_stroke_endpoint_codes"],[])
        self.assertEqual(r["scientific_status"],"HOLD")

if __name__=="__main__":
    unittest.main()
