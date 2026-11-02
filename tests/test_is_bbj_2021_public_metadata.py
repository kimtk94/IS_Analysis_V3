import importlib.util
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PATH = ROOT / "scripts/is/audit_bbj_2021_public_metadata.py"
spec = importlib.util.spec_from_file_location("bbj_public_meta", PATH)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)

PHEWEB = "<html><h1>IS: Ischemic stroke</h1><p>22664 cases, 152022 controls</p><p>Sakaue, S. & Kanai (2021)</p></html>"
OPENGWAS = "<html><h1>Ischemic stroke</h1>ebi-a-GCST90018644 Population East Asian ncase 22,664 ncontrol 152,022 Build HG19/GRCh37</html>"


class BbjPublicMetadataTests(unittest.TestCase):
    def test_valid_both(self):
        x=mod.check_sources(lambda name: PHEWEB if name == "pheweb" else OPENGWAS)
        self.assertEqual(x["public_source_metadata_agreement"], "PASS")
        self.assertEqual(x["scientific_status"], "HOLD")
        self.assertFalse(x["EAS5_archive_sha256_verified"])
        self.assertFalse(x["case_control_N_adopted_into_coloc"])

    def test_study_mismatch_fails_closed(self):
        x=mod.check_sources(lambda name: PHEWEB if name=="pheweb" else OPENGWAS.replace("GCST90018644", "GCST90018864"))
        self.assertEqual(x["public_source_metadata_agreement"], "HOLD")

    def test_phenotype_mismatch_fails_closed(self):
        x=mod.check_sources(lambda name: PHEWEB.replace("Ischemic stroke","Intracerebral hemorrhage") if name=="pheweb" else OPENGWAS)
        self.assertEqual(x["public_source_metadata_agreement"], "HOLD")

    def test_sample_count_mismatch_fails_closed(self):
        x=mod.check_sources(lambda name: PHEWEB if name=="pheweb" else OPENGWAS.replace("152,022","152,023"))
        self.assertEqual(x["public_source_metadata_agreement"], "HOLD")

    def test_unreachable_source_remains_hold(self):
        def fetch(name):
            if name=="opengwas":
                raise OSError("timeout")
            return PHEWEB
        x=mod.check_sources(fetch)
        self.assertEqual(x["public_source_metadata_agreement"], "HOLD")
        self.assertEqual(x["source_results"]["opengwas"]["failure_class"], "OSError")

    def test_reject_unknown_source(self):
        with self.assertRaises(ValueError):
            mod.inspect_source_html("https://untrusted.example",PHEWEB)

    def test_secondary_semantics_not_promoted(self):
        x=mod.check_sources(lambda name: PHEWEB if name=="pheweb" else OPENGWAS)
        self.assertFalse(x["effect_allele_primary_file_manifest_verified"])
        self.assertFalse(x["independent_GIGASTROKE_replication_verified"])
        self.assertFalse(x["EAS5_original_GWAS_source_reapproved"])


if __name__ == "__main__":
    unittest.main()
