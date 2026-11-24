"""Exercise stale request holds without any network or production job execution."""
import copy
import datetime as dt
import importlib.util
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import textwrap
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests/fixtures/operations/queued_refresh.json"
SPEC = importlib.util.spec_from_file_location("queue_operations_test", ROOT / "scripts/is/github_issue_queue_worker.py")
worker = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(worker)
NOW = dt.datetime(2026, 11, 24, 0, 0, tzinfo=dt.timezone.utc)


def request():
    return json.loads(FIXTURE.read_text())


class QueueOperationsTest(unittest.TestCase):
    def test_queued_age_cannot_be_reset_by_an_unrelated_update(self):
        issue = request()
        issue["updated_at"] = NOW.isoformat()
        hold = worker.stale_request_hold(issue, NOW)
        self.assertIn("HOLD_STALE_REQUEST", hold)
        self.assertIn("issue=#17", hold)

    def test_recently_started_running_request_is_not_classified_as_stale(self):
        issue = request()
        issue["labels"] = [{"name": "is-harness-running"}]
        issue["updated_at"] = (NOW - dt.timedelta(minutes=5)).isoformat()
        self.assertIsNone(worker.stale_request_hold(issue, NOW))
        issue["updated_at"] = (NOW - dt.timedelta(hours=1)).isoformat()
        self.assertIn("state=running", worker.stale_request_hold(issue, NOW))

    def test_stale_running_job_cannot_be_hidden_by_unchanged_inputs(self):
        issue = request()
        issue["labels"] = [{"name": "is-harness-running"}]
        with patch.object(worker, "unchanged_inputs_within_hour", return_value=True) as fast_path, \
             patch.object(worker, "api", side_effect=[[], [issue]]) as api:
            result = worker.ensure_scheduled_research_issue("SYNTHETIC_TOKEN", NOW)
        self.assertIn("HOLD_STALE_REQUEST", result)
        self.assertTrue(all(call.args[0] == "GET" for call in api.call_args_list))
        self.assertIn("labels=is-harness-running", api.call_args_list[-1].args[1])
        fast_path.assert_not_called()

    def test_pending_inventory_limit_holds_without_unbounded_pagination(self):
        batch = [dict(request(), number=number) for number in range(100)]
        with patch.object(worker, "api", return_value=batch) as api:
            with self.assertRaisesRegex(worker.QueueHold, "HOLD_QUEUE_INVENTORY_LIMIT"):
                worker.ensure_scheduled_research_issue("SYNTHETIC_TOKEN", NOW)
        api.assert_called_once()
        self.assertEqual(api.call_args.args[0], "GET")

    def test_unknown_invalid_or_future_request_age_holds(self):
        for timestamp in (None, "bad timestamp", "2026-11-24T00:00:00", "2030-01-01T00:00:00Z"):
            issue = request()
            issue["created_at"] = timestamp
            self.assertIn("HOLD_REQUEST_TIME_UNKNOWN", worker.stale_request_hold(issue, NOW))

    def test_stale_request_is_not_claimed_or_executed(self):
        issue = request()
        original = copy.deepcopy(issue)
        with patch.dict(os.environ, {"IS_QUEUE_GITHUB_TOKEN": "SYNTHETIC_TOKEN"}), \
             patch.object(sys, "argv", ["worker", "--once"]), \
             patch.object(worker, "ensure_scheduled_research_issue", return_value="PENDING_REFRESH_EXISTS"), \
             patch.object(worker, "record_heartbeat", return_value="HEARTBEAT_WRITTEN"), \
             patch.object(worker, "api", return_value=[issue]) as api, \
             patch.object(worker, "run_job") as job, patch("sys.stdout", new_callable=io.StringIO) as output:
            with self.assertRaises(SystemExit) as stopped:
                worker.main()
        self.assertEqual(stopped.exception.code, 2)
        self.assertIn("HOLD_STALE_REQUEST", output.getvalue())
        job.assert_not_called()
        self.assertEqual(api.call_count, 1)
        self.assertEqual(api.call_args.args[0], "GET")
        self.assertEqual(issue, original)

    def test_one_activation_executes_at_most_one_fresh_job(self):
        first = request()
        first["created_at"] = dt.datetime.now(dt.timezone.utc).isoformat()
        second = copy.deepcopy(first)
        second["number"] = 18
        def api_result(method, path, token, data=None):
            return [first, second] if method == "GET" else {}
        with patch.dict(os.environ, {"IS_QUEUE_GITHUB_TOKEN": "SYNTHETIC_TOKEN"}), \
             patch.object(sys, "argv", ["worker", "--once"]), \
             patch.object(worker, "ensure_scheduled_research_issue", return_value="PENDING_REFRESH_EXISTS"), \
             patch.object(worker, "record_heartbeat", return_value="HEARTBEAT_WRITTEN"), \
             patch.object(worker, "api", side_effect=api_result) as api, \
             patch.object(worker, "run_job", return_value="TECHNICAL_ONLY") as job, \
             patch("sys.stdout", new_callable=io.StringIO):
            worker.main()
        job.assert_called_once_with("research-refresh")
        self.assertFalse(any("/issues/18/" in call.args[1] for call in api.call_args_list))

    def test_ambiguous_api_failure_after_claim_is_held_without_retrying_job(self):
        issue = request()
        issue["created_at"] = dt.datetime.now(dt.timezone.utc).isoformat()
        def api_result(method, path, token, data=None):
            if method == "GET":
                return [issue]
            if method == "POST":
                raise OSError("synthetic uncertain API response")
            return {}
        with patch.dict(os.environ, {"IS_QUEUE_GITHUB_TOKEN": "SYNTHETIC_TOKEN"}), \
             patch.object(sys, "argv", ["worker", "--once"]), \
             patch.object(worker, "ensure_scheduled_research_issue", return_value="PENDING_REFRESH_EXISTS"), \
             patch.object(worker, "record_heartbeat", return_value="HEARTBEAT_WRITTEN"), \
             patch.object(worker, "api", side_effect=api_result) as api, \
             patch.object(worker, "run_job") as job, patch("sys.stdout", new_callable=io.StringIO) as output:
            with self.assertRaises(SystemExit):
                worker.main()
        job.assert_not_called()
        self.assertIn("HOLD_QUEUE_STATE_UNCERTAIN", output.getvalue())
        posted = [call for call in api.call_args_list if call.args[0] == "POST"]
        self.assertEqual(len(posted), 1)
        self.assertEqual(posted[0].args[3], {"labels": ["is-harness-running"]})

    def test_tests_operation_uses_the_repository_maintained_entrypoint(self):
        completed = subprocess.CompletedProcess([], 0, "fixture baseline passed", "")
        with patch.object(worker.subprocess, "run", return_value=completed) as child:
            worker._run_job_unlocked("tests")
        self.assertEqual(child.call_args.args[0], ["bash", "scripts/codex_smoke_test.sh"])
        self.assertTrue((ROOT / child.call_args.args[0][1]).is_file())


@unittest.skipUnless(shutil.which("node"), "Node is optional; Python queue hold coverage always runs")
class SchedulerHoldTest(unittest.TestCase):
    def invoke(self, pending, recent=None):
        workflow = (ROOT / ".github/workflows/is-harness-daily-research.yml").read_text()
        source = textwrap.dedent(workflow.split("script: |\n", 1)[1])
        harness = """
const pending = PENDING;
const recent = RECENT;
const output = {failed: [], warnings: [], created: []};
Date.now = () => Date.parse('2026-11-24T00:00:00Z');
const core = {setFailed: value => output.failed.push(value), warning: value => output.warnings.push(value), notice: () => {}};
const github = {
  rest: {issues: {
    listForRepo: async options => ({data: options.labels ? pending.filter(issue => issue.labels.some(label => label.name === options.labels)) : recent}),
    create: async value => {output.created.push(value); return {data: {number: 99, html_url: 'https://example.invalid/synthetic'}};}
  }}
};
await (async () => { SOURCE })();
process.stdout.write(JSON.stringify(output));
""".replace("PENDING", json.dumps(pending)).replace("RECENT", json.dumps(recent or [])).replace("SOURCE", source)
        result = subprocess.run([shutil.which("node"), "--input-type=module", "-"],
                                input=harness, text=True, capture_output=True, check=True, timeout=15)
        return json.loads(result.stdout)

    def test_stale_request_fails_visibly_even_after_a_recent_success(self):
        success = request()
        success["number"] = 99
        success["labels"] = [{"name": "is-harness-done"}]
        success["updated_at"] = NOW.isoformat()
        success["title"] = "IS Harness scheduled research refresh 2026-11-24 09:00 KST"
        output = self.invoke([request()], [success])
        self.assertEqual(len(output["failed"]), 1)
        self.assertIn("HOLD_STALE_REQUEST", output["failed"][0])
        self.assertEqual(output["created"], [])

    def test_fresh_request_waits_without_failure_or_duplicate(self):
        issue = request()
        issue["created_at"] = (NOW - dt.timedelta(minutes=5)).isoformat()
        output = self.invoke([issue])
        self.assertEqual(output["failed"], [])
        self.assertEqual(output["created"], [])
        self.assertEqual(len(output["warnings"]), 1)

    def test_empty_queue_creates_one_refresh(self):
        output = self.invoke([])
        self.assertEqual(output["failed"], [])
        self.assertEqual(len(output["created"]), 1)

    def test_scheduler_inventory_limit_is_an_explicit_hold(self):
        output = self.invoke([dict(request(), number=number) for number in range(100)])
        self.assertIn("HOLD_QUEUE_INVENTORY_LIMIT", output["failed"][0])
        self.assertEqual(output["created"], [])


if __name__ == "__main__":
    unittest.main()
