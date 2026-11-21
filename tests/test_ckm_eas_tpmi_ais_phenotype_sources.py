"""Network-free original EAS stroke phenotype and source overlap safety gates."""
import unittest
from scripts.ckm_stroke import audit_eas_tpmi_ais_phenotype_sources as a

def page(text):
    return {"status":"HTTP_200","source_sha256":"a"*64,
            "text":"<html><body>"+text+"</body></html>"}
def cases(name,n,c):
    return page("<h1>"+name+"</h1><p>"+str(n)+" cases, "+str(c)+" controls</p>")
def source():
    return {
      "tpmi_about":page("About TPMI PheWeb provides GWAS and PRS. All positions GRCh38."),
      "tpmi_433_21":cases("433.21 Cerebral artery occlusion, with cerebral infarction",9249,304660),
      "tpmi_433_2":cases("433.2 Occlusion of cerebral arteries",12463,304660),
      "tpmi_433_3":cases("433.3 Cerebral ischemia",12468,304660),
      "ckb_i63":cases("i63",14302,67954),
      "pgs_019250":page("PGS Catalog Polygenic Score PGS019250 GRCh38"),
      "pgs_019234":page("PGS Catalog Polygenic Score PGS019234 GRCh38"),
    }

class TestEASPhenotypeSource(unittest.TestCase):
    def test_official_source_pages_keep_scientific_hold(self):
        r=a.evidence(source())
        self.assertEqual(r["technical_status"],
              "PASS_OFFICIAL_PHEWEB_SOURCE_PROVENANCE_WITH_PER_PHENOTYPE_GATES")
        self.assertEqual(r["original_public_endpoint_checks"]["tpmi_433_21"]["n_cases"],9249)
        self.assertEqual(r["original_public_endpoint_checks"]["tpmi_433_3"]["n_cases"],12468)
        self.assertFalse(r["research_decision"]["primary_EAS_AIS_replication_GO"])
        self.assertFalse(r["research_decision"]["independent_replication_claim"])
        self.assertFalse(r["ready_for_ldsc"])
        self.assertTrue(all(not x["is_original_GWAS_summary_statistics"]
            for x in r["PGS_Catalog_not_outcome_sumstats"].values()))
    def test_missing_endpoint_is_not_verified(self):
        x=source()["tpmi_433_21"]
        x["status"]="HTTP_404";x["text"]=None
        r=a.parse_gwas_page("tpmi_433_21",x)
        self.assertFalse(r["publicly_indexed_endpoint_with_source_verified"])
        self.assertIsNone(r["n_cases"])
    def test_mismatched_phenotype_name_cannot_pass(self):
        r=a.parse_gwas_page("tpmi_433_21",cases("433.21 Intracerebral hemorrhage",9249,304660))
        self.assertEqual(r["status"],"HOLD_SOURCE_PAGE_PARSING")
    def test_gigastroke_and_tpmi_not_declared_disjoint(self):
        r=a.evidence(source())
        self.assertFalse(r["study_overlap_risk"]["TPMI_no_GIGASTROKE_participant_ID_overlap_proven"])
        self.assertTrue(r["study_overlap_risk"]["BBJ_prior_GIGASTROKE_or_MEGASTROKE_contribution_known"])
        self.assertFalse(r["external_TPMI_2026_study"]["equivalence_to_2025_TPMI_PheWeb_phenotypes_verified"])
    def test_cerebral_ischemia_is_not_a_toast_case_definition(self):
        r=a.evidence(source())["original_public_endpoint_checks"]["tpmi_433_3"]
        self.assertFalse(r["not_proven_TOAST_or_imaging_adjudicated"] is False)
        self.assertEqual(r["clinical_category"],"BROAD_ISCHEMIA_NOT_PROVEN_ACUTE_INFARCTION")
    def test_absent_tpmi_methods_build_rejected(self):
        with self.assertRaisesRegex(ValueError,"GRCh38"):
            a.parse_about(page("TPMI GWAS summary tables in GRCh37"))
    def test_bad_pgs_metadata_cannot_be_labeled_original_sumstats(self):
        with self.assertRaisesRegex(ValueError,"identifier"):
            a.parse_pgs("pgs_019250",page("PGS019250 unrelated"))
    def test_url_host_fail_closed(self):
        with self.assertRaisesRegex(ValueError,"hostname"):
            a.safe_fetch("https://example.com/hack")
    def test_invalid_source_count_fails(self):
        with self.assertRaisesRegex(ValueError,"maximum"):
            a.parse_gwas_page("tpmi_433_3",cases("Cerebral ischemia",1200001,304660))
    def test_live_portal_absent_pgs_not_promoted(self):
        r=a.parse_pgs("pgs_019250",{"status":"HTTP_429","text":None})
        self.assertFalse(r["is_original_GWAS_summary_statistics"])
        self.assertEqual(r["status"],"HTTP_429")

if __name__=="__main__":
    unittest.main()
