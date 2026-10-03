import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "muscle" / "run_muscle_stage1_rt_gwas_audit.py"
spec = importlib.util.spec_from_file_location("muscle_stage1", SCRIPT)
mod = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = mod
spec.loader.exec_module(mod)

FIXTURE_XML = b'''<?xml version="1.0"?>
<article><body><table-wrap><table>
<tr><th>SNP</th><th>Gene</th><th>P value</th></tr>
<tr><td>rs10212396</td><td>ROBO2</td><td>8.20 x 10-6</td></tr>
<tr><td>rs12625907</td><td>ISM1</td><td>2.00 x 10-8</td></tr>
</table></table-wrap></body></article>'''

class MuscleStage1Tests(unittest.TestCase):
    def test_thresholds(self):
        self.assertEqual(mod.classify_p(4.9e-8)[0], "Tier_A_genome_wide")
        self.assertEqual(mod.classify_p(9e-7)[0], "Tier_B_suggestive")
        self.assertEqual(mod.classify_p(2e-4)[0], "Tier_C_other")
        self.assertEqual(mod.classify_p(None)[0], "unresolved")

    def test_seed_unique(self):
        self.assertEqual(len(mod.SEED_2026_RSIDS), 9)
        self.assertEqual(len(set(mod.SEED_2026_RSIDS)), 9)

    def test_xml_parser_and_outputs(self):
        rows, raw = mod.extract_2026_variants(FIXTURE_XML)
        by_id = {r.rsid: r for r in rows}
        self.assertEqual(by_id["rs10212396"].tier, "Tier_B_suggestive")
        self.assertEqual(by_id["rs12625907"].tier, "Tier_A_genome_wide")
        self.assertTrue(raw)
        with tempfile.TemporaryDirectory() as td:
            xml = Path(td) / "fixture.xml"
            xml.write_bytes(FIXTURE_XML)
            out = Path(td) / "out"
            summary = mod.run(out, xml_input=xml)
            self.assertTrue((out / "MUSCLE_STAGE1_SUMMARY.json").exists())
            self.assertTrue((out / "MUSCLE_STAGE1_QC_FLAGS.tsv").exists())
            self.assertEqual(summary["seed_2026_rsids"], 9)
            self.assertEqual(summary["parsed_numeric_p_values"], 2)

if __name__ == "__main__":
    unittest.main()
