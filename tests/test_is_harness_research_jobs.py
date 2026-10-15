"""Unit tests for the strictly allowlisted IS Harness evidence jobs."""
import importlib.util
import datetime as dt
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
    def test_bot_can_request_only_approved_scheduled_refresh(self):
        legacy = "IS Harness daily research refresh 2026-10-15 KST"
        slot_17 = "IS Harness scheduled research refresh 2026-10-15 19:17 KST"
        slot_47 = "IS Harness scheduled research refresh 2026-10-15 19:47 KST"
        slots_5m = [
            "IS Harness scheduled research refresh 2026-10-15 19:00 KST",
            "IS Harness scheduled research refresh 2026-10-15 19:05 KST",
            "IS Harness scheduled research refresh 2026-10-15 19:30 KST",
            "IS Harness scheduled research refresh 2026-10-15 19:55 KST",
        ]
        for title in (legacy, slot_17, slot_47, *slots_5m):
            self.assertTrue(worker.is_trusted_issue("github-actions[bot]", "research-refresh", title))
            self.assertFalse(worker.is_trusted_issue("github-actions[bot]", "tests", title))
            self.assertFalse(worker.is_trusted_issue("github-actions[bot]", "mkdir-output", title))
            self.assertFalse(worker.is_trusted_issue("random-collaborator", "research-refresh", title))
        self.assertTrue(worker.is_trusted_issue("kimtk94", "headers", "IS arbitrary approved issue"))
        self.assertFalse(worker.is_trusted_issue("github-actions[bot]", "research-refresh", "IS Harness bootstrap"))
        self.assertFalse(worker.is_trusted_issue("github-actions[bot]", "research-refresh",
                                                  "IS Harness scheduled research refresh 2026-10-15 19:18 KST"))
        self.assertFalse(worker.is_trusted_issue("github-actions[bot]", "research-refresh",
                                                  "IS Harness scheduled research refresh 2026-10-15 25:17 KST"))
        self.assertFalse(worker.is_trusted_issue("github-actions[bot]", "research-refresh",
                                                  "IS Harness scheduled research refresh 2026-10-15 19:47 KST anything"))
        self.assertFalse(worker.is_trusted_issue("github-actions[bot]", "research-refresh",
                                                  "IS Harness scheduled research refresh 2026-10-15 19:06 KST"))
        self.assertFalse(worker.is_trusted_issue("github-actions[bot]", "research-refresh",
                                                  "IS Harness scheduled research refresh 2026-10-15 19:60 KST"))



class FiveMinuteFallbackTest(unittest.TestCase):
    def now(self, minute=5):
        return dt.datetime(2026,10,15,12,minute,23,tzinfo=dt.timezone.utc)

    def test_fallback_queues_missing_slot_with_allowlisted_body(self):
        calls = []
        def fake_api(method,path,token,data=None):
            calls.append((method,path,data))
            if method=="GET":
                return []
            if method=="POST":
                return {"number":202}
            raise AssertionError(method)
        with mock.patch.object(worker,"api",side_effect=fake_api):
            result=worker.ensure_scheduled_research_issue("TEST_TOKEN",self.now())
        self.assertEqual(result,"FALLBACK_ENQUEUED_ISSUE_202")
        _,_,request=calls[-1]
        self.assertEqual(request["title"],"IS Harness scheduled research refresh 2026-10-15 21:05 KST")
        self.assertEqual(json.loads(request["body"]),{"operation":"research-refresh"})
        self.assertEqual(request["labels"],["is-harness-queued"])
        self.assertTrue(worker.is_trusted_issue("kimtk94","research-refresh",request["title"]))

    def test_gh_cron_same_title_prevents_duplicate(self):
        title=worker.scheduled_slot_title(self.now())
        def fake_api(method,path,token,data=None):
            return [{"number":43,"title":title,"labels":[{"name":"is-harness-done"}]}]
        with mock.patch.object(worker,"api",side_effect=fake_api) as client:
            self.assertEqual(worker.ensure_scheduled_research_issue("TOKEN",self.now()),
                             "SLOT_ALREADY_ENQUEUED")
        client.assert_called_once()

    def test_pending_research_blocks_new_slot(self):
        pending={"number":43,"body":'{"operation":"research-refresh"}',
                 "labels":[{"name":"is-harness-running"}]}
        with mock.patch.object(worker,"api",side_effect=[[],[pending]]) as client:
            self.assertEqual(worker.ensure_scheduled_research_issue("TOKEN",self.now()),
                             "PENDING_REFRESH_EXISTS")
        self.assertEqual(client.call_count,2)

    def test_non_slot_minutes_skip_without_api(self):
        with mock.patch.object(worker,"api") as client:
            self.assertEqual(worker.ensure_scheduled_research_issue("TOKEN",self.now(7)),
                             "OUTSIDE_FIVE_MINUTE_SLOT")
        client.assert_not_called()



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
