import io
import math
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from ocular_paired_panel_validation import evaluate, exact_sign, gene_counts, summarize_scores


class PairedPanelTest(unittest.TestCase):
    def test_exact_two_sided_test_ties_and_two_panel_family(self):
        self.assertEqual(exact_sign([1] * 5 + [0])["two_sided_p"], 0.0625)
        self.assertEqual(exact_sign([-1] * 4 + [1])["two_sided_p"], 0.375)
        self.assertEqual(exact_sign([0, 0])["two_sided_p"], 1)
        self.assertIsNone(exact_sign([])["two_sided_p"])
        result = summarize_scores({"all_241": [1] * 5, "ECM_16": [1] * 4})
        self.assertEqual(result["all_241"]["holm_p"], 0.125)
        self.assertEqual(result["ECM_16"]["holm_p"], 0.125)

    def test_normalization_uses_all_genes_and_run_identity(self):
        text = "GeneID\tveh\tdex\nENSG1\t0\t2\nENSG2\t8\t2\n"
        matrix, genes = gene_counts(io.StringIO(text), {"ENSG1"}, {"dex", "veh"})
        self.assertEqual(genes, 2)
        self.assertEqual(matrix["ENSG1"]["veh"], 0)
        self.assertAlmostEqual(matrix["ENSG1"]["dex"], math.log2(1 + 500000))
        with self.assertRaisesRegex(ValueError, "run metadata"):
            gene_counts(io.StringIO(text), {"ENSG1"}, {"veh", "unrelated"})

    def test_pair_exclusion_is_shared_across_panels(self):
        panel = {"g": {"acute_log2fc": -2, "ecm_member": True}}
        design = {d: {a: {"run_accession": d + a, "experiment_title": "cohort"} for a in ("DEX", "ETH")}
                  for d in ("a", "b", "c")}
        runs = {"aDEX": {"QC_summary": "WARN"}, "aETH": {"QC_summary": "PASS"},
                "bDEX": {"QC_summary": "FAIL"}, "bETH": {"QC_summary": "PASS"},
                "cDEX": {"QC_summary": "PASS"}}
        expression = {"g": {"aDEX": 1, "aETH": 2, "bDEX": 3, "bETH": 1, "cDEX": 0}}
        result = evaluate(panel, design, runs, expression)
        for key in ("ECM_16", "all_241"):
            self.assertEqual(result["primary"][key]["positive"], 1)
            self.assertEqual(result["primary"][key]["pairs"], 1)
            self.assertEqual(result["all_usable"][key]["negative"], 1)
            self.assertEqual(result["all_usable"][key]["pairs"], 2)
        self.assertEqual(result["pair_coverage"][2]["arm_qc"]["ETH"], "not_in_count_release")


if __name__ == "__main__":
    unittest.main()
