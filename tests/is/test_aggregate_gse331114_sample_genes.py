#!/usr/bin/env python3
"""No real patient data: synthetic tests for complete 8-GSM CI artifact merge."""
import copy
import json
import math
import sys
import tempfile
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/"scripts"/"is"))
from aggregate_gse331114_sample_genes import (  # noqa: E402
    EXPECTED, TARGETS, IS_CANDIDATES, CELL_CONTEXT,
    reconcile, summarize_gene_rows, run,
)


def fixtures(root):
    for (group, idx), gsm in EXPECTED.items():
        folder = root / f"is-gse331114-{group}-{idx}"
        folder.mkdir(parents=True)
        # Distinct subject-library sizes, plus small group-specific differences.
        total = float(10000 + idx * 1000 + (500 if group == "AIS" else 0))
        record = {
            "status": "SAMPLE_MEX_GENE_QC_PASS", "group": group,
            "group_position_0based": idx, "accession": gsm,
            "matched_genotype": False, "clinical_gwas_patient_matching": False,
            "donor_unique_identity": "NOT_INDEPENDENTLY_VERIFIED",
            "genes": 36601, "cell_barcodes": 1000 + idx,
            "total_umi_count": total,
            "target_gene_rows": [
                {
                    "sample_GSM": gsm, "sample_group": group, "gene_symbol": symbol,
                    "gene_present_in_features": True,
                    "umis_sum": float(10 + idx + (5 if group == "AIS" else 0)),
                    "fraction_cells_nonzero": 0.05,
                }
                for symbol in TARGETS
            ],
        }
        (folder/"gene_counts.json").write_text(json.dumps(record), encoding="utf-8")


class EightSampleIntegrationTests(unittest.TestCase):
    def test_exact_eight_samples_correct_groups(self):
        self.assertEqual(len(EXPECTED), 8)
        self.assertEqual(sum(k[0] == "AIS" for k in EXPECTED), 4)
        self.assertEqual(sum(k[0] == "Control" for k in EXPECTED), 4)
        self.assertEqual(len(TARGETS), len(IS_CANDIDATES) + len(CELL_CONTEXT))
        self.assertEqual(len(TARGETS), len(set(TARGETS)))

    def test_full_integration_and_provenance(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)/"inputs"
            fixtures(root)
            samples, genes, evidence = reconcile(root)
            self.assertEqual(len(samples), 8)
            self.assertEqual(len(genes), 8 * len(TARGETS))
            self.assertEqual(len(evidence), 8)
            self.assertEqual(set(x["sample_GSM"] for x in samples), set(EXPECTED.values()))
            summary = summarize_gene_rows(genes)
            self.assertEqual(len(summary), len(TARGETS))
            self.assertTrue(all(x["p_value"] == "NOT_TESTED" for x in summary))
            self.assertTrue(all(x["fdr"] == "NOT_TESTED" for x in summary))
            self.assertTrue(all(x["donor_independence"] == "NOT_INDEPENDENTLY_VERIFIED" for x in summary))
            first = next(x for x in genes if x["sample_GSM"] == EXPECTED[("AIS", 0)] and x["gene_symbol"] == "ABO")
            self.assertAlmostEqual(first["normalized_CPM"], 15/10500*1e6)
            self.assertAlmostEqual(first["log2_CPM_plus_1"], math.log2(15/10500*1e6+1), places=6)

    def test_missing_control_artifact_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fixtures(root)
            folder = root / "is-gse331114-Control-3"
            (folder/"gene_counts.json").unlink()
            with self.assertRaisesRegex(ValueError, "Missing expected"):
                reconcile(root)

    def test_duplicate_gsm_mislabeled_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            fixtures(root)
            p=root/"is-gse331114-AIS-1"/"gene_counts.json"
            x=json.loads(p.read_text(encoding="utf-8"))
            x["accession"]="GSM9739251"
            p.write_text(json.dumps(x), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "Sample provenance"):
                reconcile(root)

    def test_incomplete_gene_target_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            fixtures(root)
            p=root/"is-gse331114-AIS-0"/"gene_counts.json"
            x=json.loads(p.read_text(encoding="utf-8"))
            x["target_gene_rows"].pop()
            p.write_text(json.dumps(x), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "panel size"):
                reconcile(root)

    def test_unevaluated_source_cannot_pass(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            fixtures(root)
            p=root/"is-gse331114-AIS-0"/"gene_counts.json"
            x=json.loads(p.read_text(encoding="utf-8"))
            x["status"]="SOURCE_UNAVAILABLE"
            p.write_text(json.dumps(x), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "Sample provenance"):
                reconcile(root)

    def test_misleading_genotype_match_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)
            fixtures(root)
            p=root/"is-gse331114-AIS-0"/"gene_counts.json"
            x=json.loads(p.read_text(encoding="utf-8"))
            x["matched_genotype"]=True
            p.write_text(json.dumps(x), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "genotype matching"):
                reconcile(root)

    def test_sample_aggregate_outputs_have_no_pvalues(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp)/"input"
            out=Path(tmp)/"out"
            fixtures(root)
            run(root, out)
            for name in ("samples_qc.tsv", "per_sample_candidate_genes_cpm.tsv",
                         "gene_group_descriptive.tsv", "descriptive_gene_medians.svg",
                         "evidence.json", "summary.md"):
                self.assertTrue((out/name).exists(), name)
            status=json.loads((out/"evidence.json").read_text(encoding="utf-8"))
            self.assertEqual(status["status"], "EIGHT_SAMPLE_DESCRIPTIVE_AGGREGATION_PASS")
            self.assertEqual(status["statistical_testing"], "NOT_PERFORMED")
            self.assertEqual(status["pseudobulk_cell_type"], "NOT_PERFORMED")

if __name__ == "__main__":
    unittest.main()
