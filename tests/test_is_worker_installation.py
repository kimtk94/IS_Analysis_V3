"""Install synthetic bundles into a temporary home; never touch user services."""
import fcntl
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests/fixtures/operations/install_bundle"
SPEC = importlib.util.spec_from_file_location("worker_install_test", ROOT / "deploy/is-harness/install_worker.py")
installer = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(installer)


class WorkerInstallationTest(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.base = Path(self.directory.name)
        self.repo = self.base / "source"
        shutil.copytree(FIXTURE, self.repo)
        self.home = self.base / "home"
        self.runtime = self.base / "runtime"
        self.runtime.mkdir()
        self.git("init", "-q")
        self.git("config", "user.name", "Synthetic Fixture")
        self.git("config", "user.email", "fixture@example.invalid")

    def git(self, *arguments):
        return subprocess.run(["git", "-C", str(self.repo), *arguments],
                              check=True, capture_output=True, text=True).stdout.strip()

    def revision(self):
        self.git("add", ".")
        self.git("commit", "-qm", "Synthetic installation fixture")
        return self.git("rev-parse", "HEAD")

    def sentinel(self):
        target = self.home / "is-harness/github_issue_queue_worker.py"
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text("existing installation remains intact\n")
        return target

    def test_all_modules_and_units_come_from_one_commit_not_worktree_edits(self):
        revision = self.revision()
        helper = self.repo / "scripts/is/harness_archived_evidence.py"
        pinned = helper.read_bytes()
        helper.write_text("raise RuntimeError('uncommitted version must not be installed')\n")
        receipt_path = installer.install_worker(self.repo, revision, self.home, self.runtime)
        receipt = json.loads(receipt_path.read_text())
        self.assertEqual(receipt["source_commit"], revision)
        self.assertFalse(receipt["service_started"])
        self.assertEqual((self.home / "is-harness" / helper.name).read_bytes(), pinned)
        for name in installer.MODULES + installer.UNITS:
            destination = (self.home / "is-harness" if name in installer.MODULES else
                           self.home / ".config/systemd/user") / name
            self.assertTrue(destination.is_file())
            self.assertEqual(receipt["files_sha256"][name], hashlib.sha256(destination.read_bytes()).hexdigest())

    def test_missing_companion_cannot_partially_replace_existing_installation(self):
        (self.repo / "scripts/is/harness_archived_evidence.py").unlink()
        revision = self.revision()
        target = self.sentinel()
        before = target.read_bytes()
        with self.assertRaises(subprocess.CalledProcessError):
            installer.install_worker(self.repo, revision, self.home, self.runtime)
        self.assertEqual(target.read_bytes(), before)
        self.assertFalse((target.parent / "installed_revision.json").exists())

    def test_runtime_import_failure_prevents_any_replacement(self):
        (self.repo / "scripts/is/harness_archived_evidence.py").write_text("raise ImportError('synthetic failure')\n")
        revision = self.revision()
        target = self.sentinel()
        before = target.read_bytes()
        with self.assertRaises(subprocess.CalledProcessError):
            installer.install_worker(self.repo, revision, self.home, self.runtime)
        self.assertEqual(target.read_bytes(), before)

    def test_busy_worker_holds_installation_without_overwriting_files(self):
        revision = self.revision()
        target = self.sentinel()
        before = target.read_bytes()
        with (self.runtime / "is-harness-queue.lock").open("a+") as active:
            fcntl.flock(active.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
            with self.assertRaisesRegex(RuntimeError, "Queue worker is active"):
                installer.install_worker(self.repo, revision, self.home, self.runtime)
        self.assertEqual(target.read_bytes(), before)

    def test_branch_names_are_rejected_before_installation(self):
        with self.assertRaisesRegex(ValueError, "full lowercase"):
            installer.install_worker(self.repo, "main", self.home, self.runtime)
        self.assertFalse((self.home / "is-harness").exists())


if __name__ == "__main__":
    unittest.main()
