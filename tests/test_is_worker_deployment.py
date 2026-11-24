"""Execute the actual remote-install shell with a synthetic command boundary."""
import hashlib
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

import yaml


ROOT = Path(__file__).resolve().parents[1]
WORKFLOW = ROOT / ".github/workflows/is-harness-safe-worker-update.yml"
FIXTURES = ROOT / "tests/fixtures/operations"


class WorkerDeploymentTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.commands = self.root / "commands.log"
        self.bin = self.root / "bin"
        shutil.copytree(FIXTURES / "deploy_commands", self.bin)
        for path in self.bin.iterdir():
            path.chmod(0o700)
        self.stage = self.root / "stage"
        bundle = self.stage / "bundle"
        bundle.mkdir(parents=True)
        installer = FIXTURES / "install_bundle/scripts/is/github_issue_queue_worker.py"
        self.payload = bundle / "install_worker.py"
        shutil.copyfile(installer, self.payload)
        self.digest = hashlib.sha256(self.payload.read_bytes()).hexdigest()
        self.workflow = yaml.safe_load(WORKFLOW.read_text())
        install_step = next(step for step in self.workflow["jobs"]["install"]["steps"]
                            if step.get("name") == "Install verified modules and units without activation")
        self.script = install_step["run"].split("<<'REMOTE'\n", 1)[1].rsplit("\nREMOTE", 1)[0]

    def run_install(self, *, digest=None, **environment):
        env = dict(os.environ, PATH=f"{self.bin}:{os.environ['PATH']}",
                   DEPLOY_TEST_COMMANDS=str(self.commands), **environment)
        return subprocess.run(
            ["bash", "-se", "--", str(self.stage), "a" * 40, digest or self.digest],
            input=self.script, capture_output=True, text=True, env=env, timeout=10)

    def calls(self):
        return self.commands.read_text() if self.commands.exists() else ""

    def test_inactive_runtime_installs_then_reloads_without_activation(self):
        result = self.run_install()
        self.assertEqual(result.returncode, 0, result.stderr)
        calls = self.calls().splitlines()
        self.assertEqual(sum(line.startswith("INSTALL ") for line in calls), 1)
        self.assertEqual(calls[-1], "--user daemon-reload")
        self.assertIn("--revision " + "a" * 40, calls[-2])
        self.assertIn("--execution-lock /srv/is-analysis/results/is/github_actions_runs/is_research_execution.lock", calls[-2])
        self.assertIn("WORKER_BUNDLE_INSTALLED", result.stdout)

    def test_active_or_transitioning_units_hold_before_install(self):
        for unit in ("is-harness-queue.service", "is-harness-queue.timer",
                     "is-auto-next.service", "is-auto-next.timer"):
            for state in ("active", "activating", "deactivating", "unknown"):
                with self.subTest(unit=unit, state=state):
                    result = self.run_install(DEPLOY_TEST_ACTIVE_UNIT=unit,
                                              DEPLOY_TEST_ACTIVE_STATE=state)
                    self.assertEqual(result.returncode, 3)
                    self.assertIn("HOLD_ACTIVE_OR_UNKNOWN_UNIT", result.stdout)
        self.assertNotIn("INSTALL ", self.calls())
        self.assertNotIn("daemon-reload", self.calls())

    def test_unavailable_systemd_does_not_become_success(self):
        result = self.run_install(DEPLOY_TEST_SYSTEMD_ERROR="1")
        self.assertEqual(result.returncode, 3)
        self.assertIn("HOLD_SYSTEMD_UNAVAILABLE", result.stdout)
        self.assertNotIn("INSTALL ", self.calls())

    def test_explicitly_absent_auto_next_unit_allows_install(self):
        result = self.run_install(DEPLOY_TEST_ABSENT_UNIT="is-auto-next.timer")
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("INSTALL ", self.calls())

    def test_running_process_or_unknown_inventory_holds(self):
        for status, reason in (("0", "HOLD_RESEARCH_PROCESS_PRESENT"),
                               ("2", "HOLD_PROCESS_INVENTORY_UNAVAILABLE")):
            with self.subTest(status=status):
                result = self.run_install(DEPLOY_TEST_PROCESS_STATUS=status)
                self.assertEqual(result.returncode, 3)
                self.assertIn(reason, result.stdout)
        self.assertNotIn("INSTALL ", self.calls())

    def test_installer_digest_is_verified_before_execution(self):
        self.payload.write_text("tampered transfer\n")
        result = self.run_install()
        self.assertNotEqual(result.returncode, 0)
        self.assertNotIn("INSTALL ", self.calls())

    def test_failed_install_does_not_reload_or_report_success(self):
        result = self.run_install(DEPLOY_TEST_INSTALL_EXIT="23")
        self.assertEqual(result.returncode, 23)
        self.assertNotIn("daemon-reload", self.calls())
        self.assertNotIn("WORKER_BUNDLE_INSTALLED", result.stdout)

    def test_reload_failure_is_not_masked_as_success(self):
        result = self.run_install(DEPLOY_TEST_RELOAD_EXIT="17")
        self.assertEqual(result.returncode, 17)
        self.assertNotIn("WORKER_BUNDLE_INSTALLED", result.stdout)

    def test_root_account_is_rejected_before_install(self):
        result = self.run_install(DEPLOY_TEST_UID="0")
        self.assertEqual(result.returncode, 2)
        self.assertEqual(self.calls(), "")


if __name__ == "__main__":
    unittest.main()
