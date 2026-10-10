"""Regression tests for the source-only IS P0 audit; never require raw GWAS."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts/is/audit_is_p0_evidence.py"
spec = importlib.util.spec_from_file_location("audit_is_p0_evidence", SCRIPT)
audit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(audit)


class EvidenceAuditTests(unittest.TestCase):
    def test_abf_counts_and_no_automatic_causal_claims(self):
        def rec(locus, tissue, gene, h4, h3=0.1, status="PASS"):
            return {"locus": locus, "dataset_key": tissue, "gene_base": gene,
                    "gene_symbol": gene, "PP.H4": str(h4), "PP.H3": str(h3),
                    "status": status}
        result = audit.analyze_abf([
            rec("L1", "cerebellum", "FGF5", 0.779, 0.217),
            rec("L2", "cerebellum", "CALHM2", 0.798),
            rec("L2", "artery", "CALHM2", 0.21),
            rec("L2", "artery", "X", 0.99, status="ERROR"),
        ])
        self.assertEqual(result["total_test_rows"], 4)
        self.assertEqual(result["pass_test_rows"], 3)
        self.assertEqual(result["h4_ge_0_5"], 2)
        self.assertEqual(result["h4_ge_0_75"], 2)
        self.assertEqual(result["h4_ge_0_8"], 0)
        self.assertEqual(result["unique_genes_tested"], 2)
        self.assertEqual(result["duplicate_test_keys"], 0)
        self.assertAlmostEqual(result["best_by_gene"]["FGF5"]["h4_h3_ratio"], 3.58986, places=4)

    def test_denominator_repeated_pairs_are_detected(self):
        row = {"locus": "L1", "dataset_key": "t", "gene_base": "G",
               "gene_symbol": "G", "PP.H3": ".05", "PP.H4": ".6", "status": "PASS"}
        self.assertEqual(audit.analyze_abf([row, row])["duplicate_test_keys"], 1)

    def test_external_neurl1_anchor_not_counted_as_positional_gene(self):
        mapped, anchors = audit.partition_candidate_rows([
            {"group_id": "IS_XDATA_1", "gene_id_stable": "ENSG111", "gene_symbol": "X", "source": "BBJ"},
            {"group_id": "", "gene_id_stable": "", "gene_symbol": "NEURL1", "source": "LEGACY_ANCHOR_ONLY"},
        ])
        self.assertEqual(len(mapped), 1)
        self.assertEqual(len(anchors), 1)
        self.assertEqual(len({r["gene_id_stable"] for r in mapped}), 1)

    def test_invalid_posterior_is_not_silently_accepted(self):
        with self.assertRaises(ValueError):
            audit.analyze_abf([{"locus": "L1", "dataset_key": "t",
                                "gene_base": "G", "gene_symbol": "G",
                                "PP.H4": "NA", "status": "PASS"}])

    def test_refuse_nonempty_output_directory(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)
            (path / "existing.txt").write_text("do not replace")
            with self.assertRaises(FileExistsError):
                audit.write_results({}, path)
            self.assertEqual((path / "existing.txt").read_text(), "do not replace")

    def test_missing_required_sources(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaisesRegex(FileNotFoundError, "Missing required source files"):
                audit.audit(Path(tmp))


if __name__ == "__main__":
    unittest.main()
