import importlib.util
import json
import tempfile
import unittest
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "muscle_stage1_gwas_audit.py"
SEED = ROOT / "config" / "muscle" / "stage1_seed_loci.tsv"
HIIT = ROOT / "config" / "muscle" / "stage1_hiit_comparator.tsv"

spec = importlib.util.spec_from_file_location("muscle_stage1", SCRIPT)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


class MuscleStage1AuditTests(unittest.TestCase):
    def test_primary_has_twenty_rt_loci_from_two_studies(self):
        rows = mod.read_tsv(SEED)
        self.assertEqual(len(rows), 20)
        self.assertTrue(all(r["cohort_type"] == "RT" for r in rows))
        counts = Counter(r["study_id"] for r in rows)
        self.assertEqual(counts["Yang_2024_PhysiolGenomics"], 11)
        self.assertEqual(counts["Gu_2026_JCSM"], 9)

    def test_hiit_comparator_has_eight_loci(self):
        rows = mod.read_tsv(HIIT)
        self.assertEqual(len(rows), 8)
        self.assertTrue(all(r["cohort_type"] == "HIIT" for r in rows))

    def test_pvalue_tiering(self):
        self.assertEqual(mod.tier_for_p(4.9e-8), "A_genomewide")
        self.assertEqual(mod.tier_for_p(5e-8), "B_suggestive")
        self.assertEqual(mod.tier_for_p(9e-6), "B_suggestive")
        self.assertEqual(mod.tier_for_p(1e-5), "C_subthreshold")

    def test_all_primary_loci_are_suggestive_by_conventional_threshold(self):
        audited = [mod.audit_row(r) for r in mod.read_tsv(SEED)]
        self.assertTrue(all(r["recalculated_tier"] == "B_suggestive" for r in audited))
        self.assertTrue(all(r["methods_threshold_pass"] == "1" for r in audited))

    def test_2026_abstract_methods_discordance_is_explicit(self):
        audited = [
            mod.audit_row(r)
            for r in mod.read_tsv(SEED)
            if r["study_id"] == "Gu_2026_JCSM"
        ]
        self.assertEqual(len(audited), 9)
        self.assertTrue(all(r["paper_internal_threshold_discordance"] == "1" for r in audited))
        self.assertTrue(all(r["abstract_threshold_pass"] == "0" for r in audited))

    def test_2024_has_no_abstract_threshold_field_but_nonstandard_methods_label(self):
        audited = [
            mod.audit_row(r)
            for r in mod.read_tsv(SEED)
            if r["study_id"] == "Yang_2024_PhysiolGenomics"
        ]
        self.assertEqual(len(audited), 11)
        self.assertTrue(all(r["paper_internal_threshold_discordance"] == "0" for r in audited))
        self.assertTrue(all(r["abstract_threshold_pass"] == "" for r in audited))
        self.assertTrue(all(r["nonstandard_genomewide_threshold"] == "1" for r in audited))

    def test_functional_triage_contains_expected_pilot_loci(self):
        audited = {r["rsid"]: r for r in (mod.audit_row(x) for x in mod.read_tsv(SEED))}
        expected_medium = {
            "rs10212396",  # ROBO2: CADD/RegulomeDB
            "rs74038095",  # SV2B: skeletal-muscle eQTL
            "rs4665972",   # SNX17: skeletal-muscle GTEx
            "rs7924637",   # CCDC15: skeletal-muscle GTEx
            "rs62149957",  # LIMS1: skeletal-muscle GTEx
        }
        for rsid in expected_medium:
            self.assertEqual(
                audited[rsid]["stage1_priority"],
                "MEDIUM_FUNCTIONAL_PRIORITY",
                rsid,
            )

    def test_end_to_end_outputs_with_comparator(self):
        with tempfile.TemporaryDirectory() as tmp:
            outdir = Path(tmp) / "out"
            summary = mod.run(SEED, outdir, HIIT)
            self.assertEqual(summary["version"], "1.1")
            self.assertEqual(summary["n_loci"], 20)
            self.assertEqual(summary["tier_counts"]["B_suggestive"], 20)
            self.assertEqual(summary["n_methods_threshold_pass"], 20)
            self.assertEqual(summary["n_paper_internal_threshold_discordant"], 9)
            self.assertEqual(summary["n_fail_abstract_claimed_threshold"], 9)
            self.assertEqual(summary["comparator"]["n_loci"], 8)
            self.assertEqual(summary["comparator"]["tier_counts"]["B_suggestive"], 8)

            expected = [
                "MUSCLE_STAGE1_GWAS_AUDIT.tsv",
                "MUSCLE_STAGE1_SUMMARY.json",
                "MUSCLE_STAGE1_REPORT.md",
                "MUSCLE_STAGE1_HIIT_COMPARATOR_AUDIT.tsv",
                "MUSCLE_STAGE1_HIIT_COMPARATOR_SUMMARY.json",
            ]
            for name in expected:
                self.assertTrue((outdir / name).exists(), name)

            saved = json.loads((outdir / "MUSCLE_STAGE1_SUMMARY.json").read_text())
            self.assertEqual(saved["n_loci"], 20)


if __name__ == "__main__":
    unittest.main()
