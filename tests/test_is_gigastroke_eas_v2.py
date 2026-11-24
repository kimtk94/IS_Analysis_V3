import importlib.util
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PATH = ROOT / "scripts/is/audit_gigastroke_eas_v2.py"
spec = importlib.util.spec_from_file_location("giga_v2", PATH)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

TRAITS = {
    "GCST90104544": "Stroke",
    "GCST90104545": "Any ischemic stroke",
    "GCST90104546": "Cardioembolic stroke",
    "GCST90104547": "Large-artery atherosclerotic stroke",
    "GCST90104548": "Small vessel stroke",
}

def mocked(url):
    parts = url.split("/")
    study = parts[-2] if parts[-1] == "ancestries" else parts[-1]
    if parts[-1] == "ancestries":
        return [{"ancestral_group": "East Asian", "number_of_individuals": "Some metadata"}]
    return {"accessionId":study, "diseaseTrait":TRAITS[study]}


class GIGASTROKEAncestryTests(unittest.TestCase):
    def test_all_five_synthetic_studies_pass_metadata_not_science(self):
        result = mod.run(mocked)
        self.assertTrue(result["accession_level_screen_complete"])
        self.assertEqual(len(result["reports"]), 5)
        self.assertEqual(result["scientific_status"], "HOLD")
        self.assertFalse(result["MR_or_coloc_performed"])
        self.assertFalse(result["GWAS_original_files_downloaded"])
        self.assertFalse(result["independent_EAS_replication_verified"])
        self.assertFalse(result["reports"][1]["coloc_approved"])

    def test_44_and_45_must_not_be_swapped(self):
        def swapped(url):
            val=mocked(url)
            if url.endswith("GCST90104544"):
                val["diseaseTrait"] = TRAITS["GCST90104545"]
            if url.endswith("GCST90104545"):
                val["diseaseTrait"] = TRAITS["GCST90104544"]
            return val
        report=mod.run(swapped)
        self.assertFalse(report["accession_level_screen_complete"])
        self.assertFalse(report["reports"][0]["official_study_trait_match"])
        self.assertFalse(report["reports"][1]["official_study_trait_match"])

    def test_study_id_mismatch_is_hold(self):
        x=mod.inspect("GCST90104545", {"accessionId":"GCST90104544",
                      "diseaseTrait":"Any ischemic stroke"},
                      [{"ancestral_group":"East Asian"}])
        self.assertFalse(x["exact_study_id"])
        self.assertEqual(x["metadata_gate"],"HOLD_METADATA_INCOMPLETE_OR_MISMATCH")

    def test_ancestry_mismatch_is_hold(self):
        x=mod.inspect("GCST90104545",mocked(mod.official_url("GCST90104545")),
                      [{"ancestral_group":"European"}])
        self.assertFalse(x["official_EAS_ancestry_match"])
        self.assertEqual(x["metadata_gate"],"HOLD_METADATA_INCOMPLETE_OR_MISMATCH")

    def test_unavailable_source_is_hold_not_false_pass(self):
        def fail(url):
            raise OSError("Unavailable")
        result=mod.run(fail)
        self.assertEqual([r["failure_class"] for r in result["reports"]],["OSError"]*5)
        self.assertFalse(result["accession_level_screen_complete"])

    def test_reject_unapproved_url(self):
        with self.assertRaises(ValueError):
            mod.fetch_json("https://untrusted.example/studies/GCST90104545")
        with self.assertRaises(ValueError):
            mod.official_url("GCST90104540")

    def test_case_control_sample_numbers_only_reference(self):
        result=mod.run(mocked)
        assert [r["published_reference_cases_not_adopted"] for r in result["reports"]] == [
            27413,19032,926,1735,5532]
        self.assertEqual(result["reports"][1]["published_reference_total_not_adopted"],256274)
        self.assertEqual(result["reports"][0]["published_reference_total_not_adopted"],264655)
        self.assertEqual(result["reports"][2]["published_reference_total_not_adopted"],238168)
        self.assertEqual(result["reports"][3]["published_reference_total_not_adopted"],238977)
        self.assertEqual(result["reports"][4]["published_reference_total_not_adopted"],242774)

    def test_nested_results_supported(self):
        m={"data":{"studies":[{"accession_id":"GCST90104545",
                               "reported_trait":"Ischaemic stroke"}]}}
        x=mod.inspect("GCST90104545",m,
                      {"data":[{"ancestralGroup":{"name":"East Asian"}}]})
        self.assertTrue(x["exact_study_id"])
        self.assertTrue(x["official_EAS_ancestry_match"])
        self.assertTrue(x["official_study_trait_match"])

    def test_study_duplicate_id_not_silently_chosen(self):
        m={"studies":[{"accessionId":"GCST90104545", "diseaseTrait":"Ischemic stroke"},
                      {"accessionId":"GCST90104545", "diseaseTrait":"Stroke"}]}
        x=mod.inspect("GCST90104545",m,[{"ancestral_group":"East Asian"}])
        self.assertFalse(x["exact_study_id"])
        self.assertEqual(x["metadata_gate"],"HOLD_METADATA_INCOMPLETE_OR_MISMATCH")


if __name__=="__main__":
    unittest.main()
