"""Regression tests: metadata audits must never promote biological claims."""
import copy
import importlib.util
import json
import pathlib
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
PATH = ROOT / "scripts/research/audit_eight_tracks.py"
spec = importlib.util.spec_from_file_location("audit_eight_tracks", PATH)
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)

class ReadinessTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.data = json.loads((ROOT / "research/operations/daily_20261101_readiness.json").read_text())

    def test_all_eight_source_gates(self):
        tracks = audit.validate_all(self.data)
        self.assertEqual(len(tracks), 8)
        for k in audit.IDS:
            self.assertEqual(audit.validate_track(tracks[k])["audit_check"], "PASS")

    def test_no_biological_pass_from_metadata(self):
        bad = copy.deepcopy(self.data)
        bad["tracks"][0]["scientific_gate"] = "PASS"
        with self.assertRaises(ValueError):
            audit.validate_all(bad)

    def test_duplicate_track_rejected(self):
        bad = copy.deepcopy(self.data)
        bad["tracks"][1]["id"] = "is"
        with self.assertRaises(ValueError):
            audit.validate_all(bad)

    def test_ckd_not_reactivated(self):
        bad = copy.deepcopy(self.data)
        bad["ckd"] = "RUN"
        with self.assertRaises(ValueError):
            audit.validate_all(bad)

    def test_glaucoma_partition(self):
        bad = copy.deepcopy(self.data)
        next(x for x in bad["tracks"] if x["id"] == "glaucoma")["observations"]["remaining_others"] = 83
        with self.assertRaises(ValueError):
            audit.validate_track(next(x for x in bad["tracks"] if x["id"] == "glaucoma"))

    def test_hcc_case_aliquot_not_conflated(self):
        bad = copy.deepcopy(self.data)
        next(x for x in bad["tracks"] if x["id"] == "hcc")["observations"]["cases"] = 331
        with self.assertRaises(ValueError):
            audit.validate_track(next(x for x in bad["tracks"] if x["id"] == "hcc"))

    def test_radiomics_patient_unit(self):
        bad = copy.deepcopy(self.data)
        next(x for x in bad["tracks"] if x["id"] == "radiomics")["observations"]["patients"] = 28
        with self.assertRaises(ValueError):
            audit.validate_track(next(x for x in bad["tracks"] if x["id"] == "radiomics"))

    def test_ocular_library_pair_not_donor_count(self):
        bad = copy.deepcopy(self.data)
        next(x for x in bad["tracks"] if x["id"] == "ocular")["observations"]["library_pairs"] = 46
        with self.assertRaises(ValueError):
            audit.validate_track(next(x for x in bad["tracks"] if x["id"] == "ocular"))

    def test_skin_download_claim(self):
        bad = copy.deepcopy(self.data)
        next(x for x in bad["tracks"] if x["id"] == "skin")["observations"]["full_sumstats_downloaded"] = True
        with self.assertRaises(ValueError):
            audit.validate_track(next(x for x in bad["tracks"] if x["id"] == "skin"))

if __name__ == "__main__":
    unittest.main()
