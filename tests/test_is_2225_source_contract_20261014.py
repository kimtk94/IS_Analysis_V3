"""Fail-closed contract tests for IS 2,225-gene workflow audit.

Synthetic tests intentionally do not require server-controlled source data.
"""
import ast
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "scripts" / "is" / "audit_existing_2225_evidence_priority_queue.py"

class TestISWorkflowSourceContract(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source = SOURCE.read_text(encoding="utf-8")
        cls.tree = ast.parse(cls.source)

    def test_has_run_entrypoint(self):
        functions = {node.name: node for node in ast.walk(self.tree) if isinstance(node, ast.FunctionDef)}
        self.assertIn("run", functions)
        self.assertEqual([a.arg for a in functions["run"].args.args],
                         ["gene_path", "region_path", "manifest_path", "out"])

    def test_explicit_no_causal_claim(self):
        self.assertIn('"direct_causal_gene_status":"NOT_ESTABLISHED"', self.source)
        self.assertIn('"causal_genes_established":0', self.source)
        self.assertIn('"research_wide_hard_candidate_filter_applied":False', self.source)

    def test_source_hashes_recorded(self):
        for name in ("source_manifest_sha256", "source_gene_matrix_sha256",
                     "source_region_matrix_sha256"):
            self.assertIn(name, self.source)

    def test_no_in_place_source_writes(self):
        # Static safety guard: the output directory must be the only write destination.
        self.assertNotIn("gene_path.write_", self.source)
        self.assertNotIn("region_path.write_", self.source)
        self.assertNotIn("manifest_path.write_", self.source)

    def test_tier_count_frozen_snapshot(self):
        self.assertIn('"A_EXISTING_SNP_REPLAY_READY":14', self.source)
        self.assertIn('"B_RECOVER_ARCHIVED_ABF_INPUT":29', self.source)

if __name__ == "__main__":
    unittest.main()
