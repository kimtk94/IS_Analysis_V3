"""Read-only runtime audit against isolated synthetic roots and command outputs."""
import argparse
import builtins
import datetime as dt
import fcntl
import hashlib
import importlib.util
import io
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "tests/fixtures/operations"
SCRIPT = ROOT / "scripts/is/audit_queue_runtime.py"
SPEC = importlib.util.spec_from_file_location("queue_runtime_audit_test", SCRIPT)
runtime = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(runtime)
REVISION = "a" * 40
NOW = dt.datetime(2026, 11, 24, 0, 5, tzinfo=dt.timezone.utc)
SECRET = "SYNTHETIC_SECRET_MUST_NOT_APPEAR"


class RuntimeAuditTest(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.base = Path(self.directory.name)
        self.args = argparse.Namespace(
            repo_root=self.base / "repo", install_root=self.base / "installed",
            unit_root=self.base / "units", runtime_dir=self.base / "runtime",
            results_root=self.base / "results", proc_root=self.base / "proc",
            freshness_seconds=3600, journal_lines=100,
        )
        for name in ("repo_root", "install_root", "unit_root", "runtime_dir", "results_root", "proc_root"):
            getattr(self.args, name).mkdir()
        bundle = FIXTURES / "install_bundle"
        hashes = {}
        for name in runtime.MODULES + runtime.UNITS:
            prefix = "scripts/is" if name in runtime.MODULES else "deploy/is-harness"
            destination = (self.args.install_root if name in runtime.MODULES else self.args.unit_root) / name
            shutil.copyfile(bundle / prefix / name, destination)
            hashes[name] = hashlib.sha256(destination.read_bytes()).hexdigest()
        (self.args.install_root / "installed_revision.json").write_text(json.dumps({
            "source_commit": REVISION, "files_sha256": hashes, "ignored_private_data": SECRET,
        }))
        evidence = self.args.results_root / "evidence"
        evidence.mkdir()
        shutil.copyfile(FIXTURES / "runtime_heartbeat.json", evidence / "harness_heartbeat.json")
        shutil.copyfile(FIXTURES / "runtime_research_manifest.json", evidence / "research_manifest.json")
        (self.args.runtime_dir / "is-harness-queue.lock").touch()
        (self.args.results_root / "is_research_execution.lock").touch()
        # A decoy token file must never be read by the diagnostic.
        (self.args.install_root / "github.env").write_text("IS_QUEUE_GITHUB_TOKEN=" + SECRET)
        for process in json.loads((FIXTURES / "runtime_processes.json").read_text()):
            self.make_process(process["pid"], process["argv"])
        self.outputs = {
            "git": REVISION + "\n",
            "service": (FIXTURES / "runtime_service.properties").read_text(),
            "timer": (FIXTURES / "runtime_timer.properties").read_text(),
            "journal": (FIXTURES / "runtime_journal.jsonl").read_text(),
        }

    def make_process(self, pid, argv):
        path = self.args.proc_root / str(pid)
        path.mkdir()
        (path / "cmdline").write_bytes(b"\0".join(arg.encode() for arg in argv) + b"\0")

    def fake_command(self, arguments, **kwargs):
        if arguments[0] == "git":
            key = "git"
            self.assertIn("--no-optional-locks", arguments)
        elif arguments[0] == "systemctl":
            self.assertEqual(arguments[1:3], ["--user", "show"])
            self.assertNotIn("Environment", " ".join(arguments))
            key = "timer" if arguments[3].endswith(".timer") else "service"
        elif arguments[0] == "journalctl":
            self.assertIn("--since=-24h", arguments)
            self.assertIn("--output-fields=MESSAGE,PRIORITY,__REALTIME_TIMESTAMP", arguments)
            key = "journal"
        else:
            self.fail("Unexpected external command")
        value = self.outputs[key]
        return value if isinstance(value, dict) else {"status": "ok", "returncode": 0, "output": value}

    def report(self):
        with patch.object(runtime, "command", side_effect=self.fake_command), \
             patch.object(runtime.os, "geteuid", return_value=1000):
            return runtime.audit(self.args, NOW)

    def test_inactive_successful_oneshot_with_enabled_timer_is_normal(self):
        report = self.report()
        self.assertEqual(report["collection_status"], "COMPLETE")
        self.assertEqual(report["status"], "OPERATIONAL_OK")
        self.assertEqual(report["service"]["properties"]["ActiveState"], "inactive")
        self.assertEqual(report["service"]["properties"]["TimeoutStartSeconds"], 600)
        self.assertFalse(report["recovery_authorized"])
        self.assertNotIn(SECRET, json.dumps(report))

    def test_missing_or_modified_installation_is_not_healthy(self):
        (self.args.install_root / runtime.MODULES[0]).unlink()
        (self.args.unit_root / runtime.UNITS[0]).write_text("[Service]\nType=simple\n")
        report = self.report()
        self.assertEqual(report["status"], "ATTENTION_REQUIRED")
        self.assertIn("INSTALLATION_UNVERIFIED", report["problems"])
        self.assertEqual(report["installation"]["files"][runtime.MODULES[0]]["status"], "missing")
        self.assertFalse(report["installation"]["files"][runtime.UNITS[0]]["matches_receipt"])

    def test_receipt_commit_and_repository_disagreement_is_visible(self):
        self.outputs["git"] = "b" * 40
        self.assertIn("REPOSITORY_INSTALL_REVISION_MISMATCH", self.report()["problems"])

    def test_missing_receipt_cannot_self_certify_installed_hashes(self):
        (self.args.install_root / "installed_revision.json").unlink()
        report = self.report()
        self.assertIsNone(report["installation"]["receipt_commit"])
        self.assertEqual(report["installation"]["status"], "unverified")

    def test_symlink_to_token_file_is_rejected_without_reading_it(self):
        path = self.args.install_root / runtime.MODULES[0]
        path.unlink()
        path.symlink_to(self.args.install_root / "github.env")
        report = self.report()
        self.assertEqual(report["installation"]["files"][runtime.MODULES[0]]["status"], "symlink")
        self.assertNotIn(SECRET, json.dumps(report))

    def test_fifo_is_rejected_without_waiting_for_a_writer(self):
        path = self.args.install_root / runtime.MODULES[0]
        path.unlink()
        os.mkfifo(path)
        self.assertEqual(runtime.read_file(path)["status"], "not_regular")
        self.assertEqual(runtime.existing_lock(path)["reason"], "not_regular")
        self.assertIn("INSTALLATION_UNVERIFIED", self.report()["problems"])

    def test_unexpected_json_field_types_are_local_unknowns_not_global_failure(self):
        heartbeat_path = self.args.results_root / "evidence/harness_heartbeat.json"
        result_path = self.args.results_root / "evidence/research_manifest.json"
        for value in ([SECRET], {"private": SECRET}):
            heartbeat = json.loads(heartbeat_path.read_text())
            heartbeat["watcher_status"] = value
            heartbeat_path.write_text(json.dumps(heartbeat))
            result = json.loads(result_path.read_text())
            result["integrity"] = value
            result_path.write_text(json.dumps(result))
            report = self.report()
            self.assertEqual(report["heartbeat"]["queue_decision"], "unknown")
            self.assertEqual(report["result_metadata"]["integrity"], "unknown")
            self.assertNotIn("DIAGNOSTIC_INTERNAL_ERROR", report["problems"])
            self.assertNotIn(SECRET, json.dumps(report))

    def test_unavailable_and_failed_service_are_observable(self):
        self.outputs["service"] = {"status": "command_failed", "returncode": 1}
        first = self.report()
        self.assertEqual(first["collection_status"], "PARTIAL")
        self.assertIn("SERVICE_STATE_UNKNOWN", first["problems"])
        self.outputs["service"] = (FIXTURES / "runtime_service.properties").read_text().replace(
            "ActiveState=inactive", "ActiveState=failed").replace("ExecMainStatus=0", "ExecMainStatus=2")
        second = self.report()
        self.assertIn("SERVICE_BLOCKED_OR_FAILED", second["problems"])
        self.assertEqual(second["collection_status"], "COMPLETE")

    def test_disabled_timer_and_never_executed_service_are_not_healthy(self):
        self.outputs["timer"] = self.outputs["timer"].replace("UnitFileState=enabled", "UnitFileState=disabled")
        self.outputs["service"] = self.outputs["service"].replace("ExecMainStartTimestampMonotonic=1000000", "ExecMainStartTimestampMonotonic=0")
        problems = self.report()["problems"]
        self.assertIn("TIMER_NOT_ACTIVE_AND_ENABLED", problems)
        self.assertIn("SERVICE_EXECUTION_NOT_OBSERVED", problems)

    def test_loaded_service_timeout_must_match_reviewed_unit_contract(self):
        self.outputs["service"] = self.outputs["service"].replace("TimeoutStartUSec=10min", "TimeoutStartUSec=90s")
        report = self.report()
        self.assertEqual(report["status"], "ATTENTION_REQUIRED")
        self.assertIn("SERVICE_TIMEOUT_MISMATCH", report["problems"])

    def test_stopping_and_active_exited_service_states_require_attention(self):
        baseline = self.outputs["service"]
        for active, sub in (("deactivating", "stop"), ("active", "exited")):
            self.outputs["service"] = baseline.replace("ActiveState=inactive", "ActiveState=" + active).replace("SubState=dead", "SubState=" + sub)
            report = self.report()
            self.assertEqual(report["status"], "ATTENTION_REQUIRED")
            self.assertIn("SERVICE_STATE_UNSUPPORTED", report["problems"])

    def test_elapsed_timer_without_pending_activation_is_not_healthy(self):
        self.outputs["timer"] = self.outputs["timer"].replace("SubState=waiting", "SubState=elapsed").replace("NextElapseUSecMonotonic=62000000", "NextElapseUSecMonotonic=0")
        report = self.report()
        self.assertEqual(report["status"], "ATTENTION_REQUIRED")
        self.assertIn("TIMER_STATE_INCOMPATIBLE", report["problems"])

    def test_waiting_timer_needs_a_pending_activation_while_service_is_idle(self):
        self.outputs["timer"] = self.outputs["timer"].replace("NextElapseUSecMonotonic=62000000", "NextElapseUSecMonotonic=0")
        report = self.report()
        self.assertIn("TIMER_NO_PENDING_ACTIVATION", report["problems"])
        self.assertEqual(report["status"], "ATTENTION_REQUIRED")

    def test_running_timer_can_await_completion_of_current_oneshot(self):
        self.outputs["timer"] = self.outputs["timer"].replace("SubState=waiting", "SubState=running").replace("NextElapseUSecMonotonic=62000000", "NextElapseUSecMonotonic=0")
        baseline = self.outputs["service"]
        for active, sub in (("activating", "start-pre"), ("activating", "start"), ("activating", "start-post")):
            self.outputs["service"] = baseline.replace("ActiveState=inactive", "ActiveState=" + active).replace("SubState=dead", "SubState=" + sub)
            self.assertEqual(self.report()["status"], "OPERATIONAL_OK")

    def test_active_running_service_does_not_match_expected_oneshot(self):
        self.outputs["service"] = self.outputs["service"].replace("ActiveState=inactive", "ActiveState=active").replace("SubState=dead", "SubState=running")
        report = self.report()
        self.assertEqual(report["status"], "ATTENTION_REQUIRED")
        self.assertIn("SERVICE_STATE_UNSUPPORTED", report["problems"])

    def test_running_timer_without_an_executing_service_requires_attention(self):
        self.outputs["timer"] = self.outputs["timer"].replace("SubState=waiting", "SubState=running")
        self.assertIn("TIMER_STATE_INCOMPATIBLE", self.report()["problems"])

    def test_error_journal_emits_categories_and_counts_without_raw_messages(self):
        self.outputs["journal"] = (FIXTURES / "runtime_error_journal.jsonl").read_text()
        report = self.report()
        counts = report["journal"]["error_counts"]
        self.assertEqual(counts["token_missing"], 1)
        self.assertEqual(counts["stale_request_hold"], 1)
        self.assertEqual(counts["authorization_error"], 1)
        self.assertEqual(counts["execution_failed"], 1)
        self.assertIn("ERRORS_IN_BOUNDED_JOURNAL_WINDOW", report["problems"])
        rendered = json.dumps(report)
        self.assertNotIn(SECRET, rendered)
        self.assertNotIn("MESSAGE", rendered)
        self.assertNotIn("issue #17", rendered)

    def test_command_errors_and_unrecognized_property_values_are_redacted(self):
        self.outputs["service"] = self.outputs["service"].replace("ActiveState=inactive", "ActiveState=" + SECRET)
        self.outputs["git"] = "fatal: " + SECRET
        report = self.report()
        self.assertIn("SERVICE_STATE_UNKNOWN", report["problems"])
        self.assertIn("REPOSITORY_HEAD_UNKNOWN", report["problems"])
        self.assertNotIn(SECRET, json.dumps(report))

    def test_missing_stale_or_held_heartbeat_requires_attention(self):
        path = self.args.results_root / "evidence/harness_heartbeat.json"
        heartbeat = json.loads(path.read_text())
        heartbeat["observed_at"] = "2020-01-01T00:00:00Z"
        path.write_text(json.dumps(heartbeat))
        self.assertIn("HEARTBEAT_STALE_OR_UNKNOWN", self.report()["problems"])
        heartbeat["observed_at"] = NOW.isoformat()
        heartbeat["watcher_status"] = "HOLD_STALE_REQUEST issue=#17 " + SECRET
        path.write_text(json.dumps(heartbeat))
        held = self.report()
        self.assertIn("RECORDED_QUEUE_HOLD", held["problems"])
        self.assertNotIn(SECRET, json.dumps(held))
        path.unlink()
        self.assertIn("HEARTBEAT_STALE_OR_UNKNOWN", self.report()["problems"])

    def test_processes_are_counted_without_rendering_command_arguments(self):
        for pid in (11002, 11003):
            self.make_process(pid, ["/usr/bin/python3", "/home/example/is-harness/github_issue_queue_worker.py", "--token", SECRET])
        report = self.report()
        self.assertEqual(report["processes"]["counts"]["queue_worker"], 2)
        self.assertIn("MULTIPLE_QUEUE_WORKERS", report["problems"])
        self.assertNotIn(SECRET, json.dumps(report))

    def test_existing_lock_states_never_create_a_missing_lock(self):
        path = self.args.runtime_dir / "is-harness-queue.lock"
        self.assertEqual(runtime.existing_lock(path)["status"], "free")
        with path.open("rb") as active:
            fcntl.flock(active.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
            self.assertEqual(runtime.existing_lock(path)["status"], "held")
        path.unlink()
        self.assertEqual(runtime.existing_lock(path)["status"], "missing")
        self.assertFalse(path.exists())

    def test_audit_does_not_write_files_or_read_token_configuration(self):
        original_open = os.open
        original_builtin = builtins.open
        original_io = io.open
        def guarded_open(path, flags, *args, **kwargs):
            self.assertFalse(flags & (os.O_WRONLY | os.O_RDWR | os.O_CREAT | os.O_TRUNC | os.O_APPEND))
            self.assertNotEqual(Path(path).name, "github.env")
            return original_open(path, flags, *args, **kwargs)
        def guard_stream(delegate):
            def guarded(path, mode="r", *args, **kwargs):
                self.assertFalse(any(character in mode for character in "wax+"))
                if not isinstance(path, int):
                    self.assertNotEqual(Path(path).name, "github.env")
                return delegate(path, mode, *args, **kwargs)
            return guarded
        before = {str(path.relative_to(self.base)): (path.read_bytes(), path.stat().st_mtime_ns)
                  for path in self.base.rglob("*") if path.is_file()}
        with patch.object(runtime.os, "open", side_effect=guarded_open), \
             patch("builtins.open", side_effect=guard_stream(original_builtin)), \
             patch("io.open", side_effect=guard_stream(original_io)), \
             patch.object(runtime.os, "mkdir", side_effect=AssertionError("file mutation")), \
             patch.object(runtime.os, "replace", side_effect=AssertionError("file mutation")), \
             patch.object(runtime.os, "unlink", side_effect=AssertionError("file mutation")):
            report = self.report()
        after = {str(path.relative_to(self.base)): (path.read_bytes(), path.stat().st_mtime_ns)
                 for path in self.base.rglob("*") if path.is_file()}
        self.assertEqual(before, after)
        self.assertEqual(report["status"], "OPERATIONAL_OK")

    def test_exit_two_still_prints_sanitized_machine_report(self):
        with patch.object(runtime, "audit", return_value={"status": "ATTENTION_REQUIRED", "problems": ["SERVICE_STATE_UNKNOWN"], "read_only": True}), \
             patch("sys.stdout", new_callable=io.StringIO) as output:
            code = runtime.main([])
        self.assertEqual(code, 2)
        self.assertEqual(json.loads(output.getvalue())["problems"], ["SERVICE_STATE_UNKNOWN"])


class BoundedCommandTest(unittest.TestCase):
    def test_large_output_is_discarded_after_the_limit(self):
        result = runtime.command([sys.executable, "-I", "-B", "-c", "print('x'*10000)"], limit=50)
        self.assertEqual(result["status"], "oversized")
        self.assertNotIn("output", result)

    def test_timeout_does_not_emit_child_stderr(self):
        result = runtime.command([sys.executable, "-I", "-B", "-c", "import sys,time;sys.stderr.write('SYNTHETIC_SECRET_MUST_NOT_APPEAR');time.sleep(2)"], timeout=0.05)
        self.assertEqual(result["status"], "timeout")
        self.assertNotIn(SECRET, json.dumps(result))

    def test_nonzero_command_output_is_not_rendered(self):
        result = runtime.command([sys.executable, "-I", "-B", "-c", "print('SYNTHETIC_SECRET_MUST_NOT_APPEAR');raise SystemExit(1)"])
        self.assertEqual(result["status"], "command_failed")
        self.assertNotIn(SECRET, json.dumps(result))


class WorkflowRouteTest(unittest.TestCase):
    def test_worker_status_streams_reviewed_code_without_install_or_execution(self):
        source = (ROOT / ".github/workflows/is-server-readonly-ssh-audit.yml").read_text()
        step = source.split("- name: Inspect worker metadata without installation or execution", 1)[1].split("- name:", 1)[0]
        self.assertIn("'python3 -I -B -' < scripts/is/audit_queue_runtime.py", step)
        self.assertIn("StrictHostKeyChecking=yes", step)
        self.assertIn("ForwardAgent=no", step)
        for forbidden in ("scp ", "--once", "systemctl", "github.env", "issues.create", "install_worker"):
            self.assertNotIn(forbidden, step)
        self.assertNotIn("issues: write", source)


if __name__ == "__main__":
    unittest.main()
