"""Fail-closed readback of real 24 EAS GWAS original source summary artifact."""
from __future__ import annotations
import copy
import importlib
from pathlib import Path
import unittest

M=importlib.import_module("scripts.is.audit_gws_eas_24_leads_original")
R=importlib.import_module("scripts.is.review_gws_eas_24_original_leads")
OLD=Path("legacy/docs/data/is/IS_30_REGIONS_SIX_OUTCOMES_GWAS_CONTEXT_500KB.tsv")

def fake():
    a=M.read_legacy_180(OLD)
    rows=[]
    for r in a["selected"]:
        raw={"chrom":r["chromosome"],"pos":r["position"],
             "beta":r["historical_beta_ALT"],
             "se":.01,"p":r["historical_p"],
             "eaf":r["historical_ALT_EAF"],
             "effect_allele":r["historical_alt_allele"],
             "other_allele":r["historical_ref_allele"]}
        rows.append(M.compare(r,[raw]))
    return {
        "schema":"is-gws-eas-30-provisional-region-24-nonpal-original-variant-audit-v1",
        "technical_status":"PASS_FULL_PUBLISHER_FILE_AND_24_LEAD_SOURCE_SCREEN",
        "scientific_status":"HOLD",
        "historical_GWAS_table_blob_sha1":M.EXPECTED_LEGACY_BLOB,
        "verified_whole_file_original_source":{
            "accession":"GCST90104545",
            "compressed_original_file_bytes":M.GATE.SOURCES["EAS"]["bytes"],
            "official_publisher_md5_and_pinned_sha256":True},
        "original_GWAS_full_variant_rows_scanned":6725216,
        "original_nonpalindromic_archived_leads_tested":24,
        "independently_excluded_palindromic_source_leads":4,
        "source_missing_lead_regions":2,
        "original_same_position_allele_pair_match_count":24,
        "original_signed_beta_p_eaf_all_three_match_count":24,
        "unique_chromosomes_with_concordant_leads":len({x["chromosome"] for x in rows}),
        "technical_source_archive_signed_reproduction_screen":"PASS_REVIEW_REQUIRED",
        "published_original_study_Ncase_Ncontrol_independently_approved":False,
        "BBJ_stroke_source_provenance_verified":False,
        "matched_QTL_cis_effects_verified":False,
        "coloc_abf_executed":0,
        "coloc_susie_executed":0,
        "tested_leads":rows,
        "palindromic_leads_excluded":a["withheld"],
        "historical_missing_leads":a["missing"]
    }

class Original24ReadbackTests(unittest.TestCase):
    def test_authentic_24_reproduced_but_never_coloc(self):
        out=R.inspect(fake(),OLD)
        self.assertEqual(out["nonpalindromic_leads_tested"],24)
        self.assertEqual(out["signed_beta_p_eaf_reproduced"],24)
        self.assertEqual(out["palindromic_leads_scientifically_held"],4)
        self.assertEqual(out["original_leads_missing"],2)
        self.assertEqual(out["P3A_QTL_gene_assay_count"],36)
        self.assertEqual(out["P3A_unique_candidate_gene_count"],20)
        self.assertEqual(out["P3A_coloc_abf"],0)
        self.assertEqual(out["scientific_status"],"HOLD")

    def test_source_archive_not_approved_with_malformed_original_study_digest(self):
        d=fake()
        d["verified_whole_file_original_source"]["official_publisher_md5_and_pinned_sha256"]=False
        with self.assertRaisesRegex(ValueError,"publisher source"):
            R.inspect(d,OLD)

    def test_one_assay_deleted_fail_closed(self):
        d=fake()
        d["tested_leads"].pop()
        with self.assertRaisesRegex(ValueError,"24 true original"):
            R.inspect(d,OLD)

    def test_gene_cohort_sample_N_false_approval_denied(self):
        d=fake()
        d["published_original_study_Ncase_Ncontrol_independently_approved"]=True
        with self.assertRaisesRegex(ValueError,"scientific eligibility"):
            R.inspect(d,OLD)

    def test_mutated_lead_allele_key_denied(self):
        d=fake()
        d["tested_leads"][0]["GRCh37_variant_key"]="1:1:A:T"
        with self.assertRaisesRegex(ValueError,"variant key"):
            R.inspect(d,OLD)

    def test_mutated_effect_match_denied(self):
        d=fake()
        d["tested_leads"][0]["ALT_oriented_beta_match"]=False
        with self.assertRaisesRegex(ValueError,"joint agreement"):
            R.inspect(d,OLD)

    def test_status_gates_not_silently_promoted(self):
        d=fake()
        d["technical_source_archive_signed_reproduction_screen"]="PASS_CAUSAL"
        with self.assertRaisesRegex(ValueError,"acceptance screen"):
            R.inspect(d,OLD)

    def test_untrusted_partner_cohort_claim_denied(self):
        d=fake()
        d["BBJ_stroke_source_provenance_verified"]=True
        with self.assertRaisesRegex(ValueError,"scientific eligibility"):
            R.inspect(d,OLD)

if __name__=="__main__":
    unittest.main()
