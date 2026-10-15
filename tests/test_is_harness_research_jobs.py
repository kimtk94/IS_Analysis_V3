"""Unit tests for the strictly allowlisted IS Harness evidence jobs."""
import importlib.util
import io
import json
import pathlib
import tempfile
import unittest
from unittest import mock

SCRIPT = pathlib.Path(__file__).resolve().parents[1] / "scripts/is/harness_research_jobs.py"
spec = importlib.util.spec_from_file_location("is_harness_research_jobs", SCRIPT)
jobs = importlib.util.module_from_spec(spec)
spec.loader.exec_module(jobs)

WORKER_SCRIPT = pathlib.Path(__file__).resolve().parents[1] / "scripts/is/github_issue_queue_worker.py"
worker_spec = importlib.util.spec_from_file_location("is_github_queue_worker", WORKER_SCRIPT)
worker = importlib.util.module_from_spec(worker_spec)
worker_spec.loader.exec_module(worker)


class ScheduledQueueSecurityTest(unittest.TestCase):
    def test_bot_can_request_only_fixed_daily_refresh(self):
        name = "IS Harness daily research refresh 2026-10-15 KST"
        self.assertTrue(worker.is_trusted_issue("github-actions[bot]", "research-refresh", name))
        self.assertTrue(worker.is_trusted_issue("kimtk94", "headers", "IS arbitrary approved issue"))
        self.assertFalse(worker.is_trusted_issue("github-actions[bot]", "tests", name))
        self.assertFalse(worker.is_trusted_issue("github-actions[bot]", "mkdir-output", name))
        self.assertFalse(worker.is_trusted_issue("github-actions[bot]", "research-refresh", "IS Harness bootstrap"))
        self.assertFalse(worker.is_trusted_issue("random-collaborator", "research-refresh", name))


class ReadinessTest(unittest.TestCase):
    def make_frozen(self, *, duplicate=False):
        headers = [
            "gene_id_stable", "new_nominal_QTL_gene_source_status",
            "ischemic_stroke_causal_gene", "coloc_recomputed",
            "GWAS_and_QTL_alleles_coharmonized", "matched_QTL_LD_verified",
        ]
        lines = ["\t".join(headers)]
        n = 0
        for category, count in jobs.EXPECTED_READY.items():
            for i in range(count):
                n += 1
                gene = f"ENSG{n:011d}"
                if duplicate and n == 2:
                    gene = f"ENSG{1:011d}"
                lines.append("\t".join([gene, category, "NOT_ESTABLISHED", "False", "False", "False"]))
        return ("\n".join(lines) + "\n").encode()

    def test_frozen_snapshot_pass(self):
        data = self.make_frozen()
        report = jobs.audit_candidates(data)
        self.assertEqual(report["unique_ensembl_genes"], 2225)
        self.assertEqual(report["readiness_status_counts"], jobs.EXPECTED_READY)
        self.assertEqual(report["integrity"], "PASS")
        self.assertEqual(report["causal_gene_conclusions"], "NOT_ESTABLISHED_BY_THIS_AUDIT")

    def test_duplicate_gene_fails_closed(self):
        with self.assertRaises(ValueError):
            jobs.audit_candidates(self.make_frozen(duplicate=True))

    def test_short_snapshot_fails_closed(self):
        data = self.make_frozen().decode().splitlines()
        with self.assertRaises(ValueError):
            jobs.audit_candidates(("\n".join(data[:-1])+"\n").encode())


class EvidenceAuditTest(unittest.TestCase):
    def write_table(self, path, columns, rows):
        with open(path, "w", encoding="utf-8") as stream:
            stream.write("\t".join(columns)+"\n")
            for row in rows:
                stream.write("\t".join(str(row.get(c, "")) for c in columns)+"\n")

    def test_five_region_source_audit_is_not_coloc(self):
        with tempfile.TemporaryDirectory() as folder:
            root = pathlib.Path(folder)
            a, b, c = [root / x for x in ("variants.tsv", "ready.tsv", "overlap.tsv")]
            rows = []
            for i in range(5):
                row = dict.fromkeys(jobs.REQUIRED_VARIANT_COLUMNS, "")
                row.update({
                    "region_id":f"G00{i+1}", "QTL_dataset":"QTD1", "gene_id":f"ENSG{i:011d}",
                    "BBJ_exact_SNP":"1:1:A:C", "EAS_AIS_exact_SNP":"1:1:A:C",
                    "both_GWAS_sources_nonpalindromic":"True", "status":"MATCHED",
                    "matched_QTL_LD_verified":"False", "coloc_recomputed":"False",
                })
                rows.append(row)
            self.write_table(a, jobs.REQUIRED_VARIANT_COLUMNS, rows)
            self.write_table(b, jobs.REQUIRED_READINESS_COLUMNS, [
                {"region_id":"G001","QTL_dataset":"QTD1","gene_id":"ENSG1"}
            ])
            self.write_table(c, jobs.REQUIRED_OVERLAP_COLUMNS, [
                {"region_id":"G001","QTL_dataset":"QTD1","gene_id":"ENSG1"}
            ])
            with mock.patch.object(jobs, "EAS5_VARIANTS", a), \
                 mock.patch.object(jobs, "EAS5_READINESS", b), \
                 mock.patch.object(jobs, "EAS5_OVERLAP", c):
                report = jobs.audit_eas5()
            self.assertEqual(report["rows"], 5)
            self.assertEqual(report["regions"], 5)
            self.assertEqual(report["dual_gwas_exact_snp_rows"], 5)
            self.assertEqual(report["matched_qtl_ld_verified_rows"], 0)
            self.assertEqual(report["coloc_recomputed_rows"], 0)
            self.assertEqual(report["causal_gene_conclusions"], "NOT_ESTABLISHED_BY_THIS_AUDIT")
            cand = jobs.audit_candidates(ReadinessTest().make_frozen())
            manifest = jobs.manifest(report, cand)
            self.assertEqual(manifest["schema"], "is-research-evidence-v1")
            self.assertEqual(sum(x["value"] for x in manifest["figures"][0]["bars"]), 2225)
            self.assertEqual(manifest["figures"][1]["bars"][-1]["value"], 0)


if __name__ == "__main__":
    unittest.main()
