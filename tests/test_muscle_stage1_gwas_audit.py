import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "muscle_stage1_gwas_audit.py"
SEED = ROOT / "config" / "muscle" / "stage1_seed_loci.tsv"

spec = importlib.util.spec_from_file_location("muscle_stage1", SCRIPT)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


class MuscleStage1AuditTests(unittest.TestCase):
    def test_seed_has_nine_unique_loci(self):
        rows = mod.read_tsv(SEED)
        self.assertEqual(len(rows), 9)
        self.assertEqual(len({r["rsid"] for r in rows}), 9)

    def test_pvalue_tiering(self):
        self.assertEqual(mod.tier_for_p(4.9e-8), "A_genomewide")
        self.assertEqual(mod.tier_for_p(5e-8), "B_suggestive")
        self.assertEqual(mod.tier_for_p(9e-6), "B_suggestive")
        self.assertEqual(mod.tier_for_p(1e-5), "C_subthreshold")

    def test_2026_seed_is_reclassified_as_suggestive(self):
        audited = [mod.audit_row(r) for r in mod.read_tsv(SEED)]
        self.assertTrue(all(r["recalculated_tier"] == "B_suggestive" for r in audited))
        self.assertTrue(all(r["claim_threshold_discordance"] == "1" for r in audited))

    def test_functional_priority_flags_robo2_and_sv2b(self):
        audited = {r["rsid"]: r for r in (mod.audit_row(x) for x in mod.read_tsv(SEED))}
        self.assertEqual(audited["rs10212396"]["cadd_ge_12_37"], "1")
        self.assertEqual(audited["rs10212396"]["regulomedb_1_or_2"], "1")
        self.assertEqual(audited["rs74038095"]["skeletal_muscle_eqtl_reported"], "1")

    def test_end_to_end_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            outdir = Path(tmp) / "out"
            summary = mod.run(SEED, outdir)
            self.assertEqual(summary["n_loci"], 9)
            self.assertEqual(summary["tier_counts"]["B_suggestive"], 9)
            self.assertEqual(summary["n_claim_threshold_discordant"], 9)
            self.assertTrue((outdir / "MUSCLE_STAGE1_GWAS_AUDIT.tsv").exists())
            saved = json.loads((outdir / "MUSCLE_STAGE1_SUMMARY.json").read_text())
            self.assertEqual(saved["n_loci"], 9)


if __name__ == "__main__":
    unittest.main()
