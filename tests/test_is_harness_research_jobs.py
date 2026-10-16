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

    def test_late_minute_rolls_into_current_bucket(self):
        self.assertEqual(worker.scheduled_slot_title(self.now(7)),
                         "IS Harness scheduled research refresh 2026-10-15 21:05 KST")
        self.assertEqual(worker.scheduled_slot_title(self.now(9)),
                         "IS Harness scheduled research refresh 2026-10-15 21:05 KST")

    def test_unchanged_source_fingerprint_skips_issue_and_recompute(self):
        with mock.patch.object(worker,"unchanged_inputs_within_hour",return_value=True), \
             mock.patch.object(worker,"api") as client:
            self.assertEqual(worker.ensure_scheduled_research_issue("TOKEN",self.now(7)),
                             "UNCHANGED_SOURCES_FAST_CHECK")
        client.assert_not_called()



class HeartbeatAuditTest(unittest.TestCase):
    def test_watcher_writes_heartbeat_without_science_claim(self):
        with tempfile.TemporaryDirectory() as folder:
            path=pathlib.Path(folder)/"harness_heartbeat.json"
            t=dt.datetime(2026,10,16,0,50,35,tzinfo=dt.timezone.utc)
            self.assertEqual(worker.record_heartbeat(t,"UNCHANGED_SOURCES_FAST_CHECK",path),
                             "HEARTBEAT_WRITTEN")
            saved=json.loads(path.read_text())
            self.assertFalse(saved["is_new_scientific_result"])
            self.assertEqual(saved["watcher_status"],"UNCHANGED_SOURCES_FAST_CHECK")
            self.assertEqual(worker.record_heartbeat(t,"UNCHANGED_SOURCES_FAST_CHECK",path),
                             "UNCHANGED_HEARTBEAT")

class ReadinessTest(unittest.TestCase):
    def make_frozen(self, *, duplicate=False):
        headers = [
            "gene_id_stable", "new_nominal_QTL_gene_source_status",
            "ischemic_stroke_causal_gene", "coloc_recomputed",
            "GWAS_and_QTL_alleles_coharmonized", "matched_QTL_LD_verified",
            "canonical_gene_symbol", "region_id", "source_nominal_QTL_ids",
            "n_new_nominal_QTL_gene_variant_rows", "min_source_nominal_QTL_p",
        ]
        lines = ["\t".join(headers)]
        n = 0
        for category, count in jobs.EXPECTED_READY.items():
            for i in range(count):
                n += 1
                gene = f"ENSG{n:011d}"
                if duplicate and n == 2:
                    gene = f"ENSG{1:011d}"
                nominal = category == "RAW_NOMINAL_QTL_RECORDS_AVAILABLE_NO_COLOC"
                lines.append("\t".join([
                    gene, category, "NOT_ESTABLISHED", "False", "False", "False",
                    "GENE_"+str(n), "R1" if nominal else "", "QTD1" if nominal else "",
                    "5" if nominal else "0", "0.00001" if nominal else "",
                ]))
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


class SourceSnpFieldInterpretationTests(unittest.TestCase):
    def test_false_string_is_not_true_snp_match(self):
        self.assertFalse(jobs.source_exact_match("False"))
        self.assertFalse(jobs.source_exact_match("0"))
        self.assertFalse(jobs.source_exact_match(""))
        self.assertTrue(jobs.source_exact_match("True"))
        self.assertTrue(jobs.source_exact_match("chr1:123:A:C"))
        with self.assertRaises(ValueError):
            jobs.source_exact_match("UNKNOWN_ENCODING")

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
            manifest = jobs.manifest(report, cand, input_signature="UNIT_TEST_FROZEN_METADATA")
            self.assertEqual(manifest["schema"], "is-research-evidence-v1")
            self.assertEqual(sum(x["value"] for x in manifest["figures"][0]["bars"]), 2225)
            self.assertEqual(manifest["figures"][1]["bars"][-1]["value"], 0)


    def test_deep_p2_qc_and_p4_gene_triage(self):
        with tempfile.TemporaryDirectory() as folder:
            root = pathlib.Path(folder)
            variants = root / "variants.tsv"
            cols = jobs.REQUIRED_VARIANT_COLUMNS + ("BBJ_p","EAS_AIS_p")
            rows = []
            for i in range(5):
                data = dict.fromkeys(cols,"")
                data.update({
                    "region_id":f"R{i+1}","QTL_dataset":"QTD1",
                    "gene_id":f"ENSG{i:011d}",
                    "QTL_p":"0.005" if i!=2 else "invalid",
                    "BBJ_p":"0.02","EAS_AIS_p":"1.5" if i==3 else "0.03",
                    "QTL_beta_mapped_GRCh37_ALT_OFFICIAL_DEFINITION":"-0.12",
                    "BBJ_exact_SNP":"1:1:A:C","EAS_AIS_exact_SNP":"1:1:A:C",
                    "both_GWAS_sources_nonpalindromic":"True",
                })
                rows.append(data)
            self.write_table(variants, cols, rows)
            readiness = root/"readiness.tsv"
            self.write_table(readiness,jobs.REQUIRED_READINESS_COLUMNS,[
                {"region_id":"R1","QTL_dataset":"QTD1","gene_id":"ENSG1",
                 "scientific_status":"BLOCKED","QTL_source_ancestry_verified":"False",
                 "molecular_QTL_matched_LD_verified":"False"}
            ])
            overlap = root/"overlap.tsv"
            self.write_table(overlap,jobs.REQUIRED_OVERLAP_COLUMNS,[{
                "region_id":"R1","QTL_dataset":"QTD1","gene_id":"ENSG1"
            }])
            with mock.patch.object(jobs,"EAS5_VARIANTS",variants), \
                 mock.patch.object(jobs,"EAS5_READINESS",readiness), \
                 mock.patch.object(jobs,"EAS5_OVERLAP",overlap):
                e=jobs.audit_eas5()
                p1=jobs.audit_p1_ld_gates()
                meta=jobs.input_metadata_signature()
            self.assertEqual(e["p2_detail"]["group_count"],5)
            self.assertEqual(e["p2_detail"]["qc"]["QTL_p_invalid_nonmissing"],1)
            self.assertEqual(e["p2_detail"]["qc"]["EAS_AIS_p_invalid_nonmissing"],1)
            self.assertEqual(p1["rows"],1)
            self.assertEqual(p1["by_region"]["R1"].get("molecular_qtl_ld_verified",0),0)
            self.assertEqual(len(meta),64)
            c=jobs.audit_candidates(ReadinessTest().make_frozen())
            self.assertEqual(c["p4_detail"]["nominal_qtl_candidate_count"],149)
            self.assertEqual(len(c["p4_detail"]["candidates"]),149)
            summary=jobs.manifest(e,c,p1,input_signature=meta)
            self.assertEqual(summary["research_qc"]["p2_qc"]["invalid_qtl_p"],1)
            self.assertEqual(summary["research_qc"]["p4_candidates"]["nominal_qtl_candidate_count"],149)
            self.assertEqual(summary["research_qc"]["p1_ld"]["matched_ld_verified_rows"],0)



if __name__ == "__main__":
    unittest.main()
