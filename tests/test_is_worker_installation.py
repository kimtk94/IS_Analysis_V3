"""Install synthetic bundles into a temporary home; never touch user services."""
import fcntl
import hashlib
import importlib.util
import io
import json
import os
from pathlib import Path
import shutil
import stat
import subprocess
import tempfile
import unittest
from unittest import mock

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

    def snapshot(self):
        return {str(path.relative_to(self.home)): (path.read_bytes(), stat.S_IMODE(path.stat().st_mode))
                for path in self.home.rglob("*") if path.is_file()}

    def exported(self):
        revision = self.revision()
        bundle = self.base / "exported"
        installer.export_bundle(self.repo, revision, bundle)
        return revision, bundle

    def installed_bundle(self):
        revision, bundle = self.exported()
        installer.install_bundle(bundle, revision, self.home, self.runtime)
        return revision, bundle

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

    def test_export_includes_pinned_installer_and_exact_inventory(self):
        revision = self.revision()
        expected = {}
        for name in installer.BUNDLE_FILES:
            folder = "scripts/is" if name in installer.MODULES else "deploy/is-harness"
            path = self.repo / folder / name
            expected[name] = path.read_bytes()
            path.write_text("raise RuntimeError('working-tree file must never be imported or exported')\n")
        bundle = self.base / "bundle"
        manifest_path = installer.export_bundle(self.repo, revision, bundle)
        self.assertEqual(set(installer.BUNDLE_FILES) | {installer.MANIFEST},
                         {path.name for path in bundle.iterdir()})
        manifest = json.loads(manifest_path.read_text())
        self.assertEqual(installer.BUNDLE_SCHEMA, manifest["schema"])
        self.assertEqual(revision, manifest["source_commit"])
        for name, pinned in expected.items():
            self.assertEqual(pinned, (bundle / name).read_bytes())
            self.assertEqual(hashlib.sha256(pinned).hexdigest(), manifest["files_sha256"][name])

    def test_export_never_imports_even_the_pinned_modules(self):
        (self.repo / "scripts/is/harness_archived_evidence.py").write_text(
            "raise RuntimeError('export must not execute imports')\n")
        revision, bundle = self.exported()
        self.assertEqual(revision, json.loads((bundle / installer.MANIFEST).read_text())["source_commit"])

    def test_bundle_receipt_identifies_all_six_installed_files_and_installer(self):
        revision, bundle = self.installed_bundle()
        receipt = json.loads((self.home / "is-harness/installed_revision.json").read_text())
        manifest = json.loads((bundle / installer.MANIFEST).read_text())
        self.assertEqual(revision, receipt["source_commit"])
        self.assertEqual(installer.BUNDLE_SCHEMA, receipt["bundle_schema"])
        self.assertFalse(receipt["service_started"])
        self.assertEqual(set(installer.MODULES + installer.UNITS), set(receipt["files_sha256"]))
        self.assertEqual(manifest["files_sha256"][installer.INSTALLER], receipt["installer_sha256"])
        self.assertEqual(hashlib.sha256((bundle / installer.MANIFEST).read_bytes()).hexdigest(),
                         receipt["bundle_manifest_sha256"])
        for name in installer.MODULES + installer.UNITS:
            root = self.home / ("is-harness" if name in installer.MODULES else ".config/systemd/user")
            self.assertEqual((bundle / name).read_bytes(), (root / name).read_bytes())
        self.assertFalse((self.home / "is-harness" / installer.INSTALLER).exists())

    def test_existing_export_directory_is_preserved(self):
        revision = self.revision()
        bundle = self.base / "bundle"
        bundle.mkdir()
        marker = bundle / "marker"
        marker.write_bytes((self.repo / "scripts/is/harness_archived_evidence.py").read_bytes())
        before = marker.read_bytes()
        with self.assertRaisesRegex(ValueError, "must be new"):
            installer.export_bundle(self.repo, revision, bundle)
        self.assertEqual(before, marker.read_bytes())
        self.assertEqual([marker], list(bundle.iterdir()))

    def test_tampered_payload_or_installer_preserves_existing_installation(self):
        revision, original = self.installed_bundle()
        before = self.snapshot()
        for number, name in enumerate((installer.MODULES[0], installer.INSTALLER, installer.UNITS[0])):
            with self.subTest(file=name):
                bundle = self.base / f"tampered-{number}"
                shutil.copytree(original, bundle)
                path = bundle / name
                path.write_bytes(path.read_bytes() + b"\n# changed after export\n")
                with self.assertRaisesRegex(ValueError, "hash mismatch"):
                    installer.install_bundle(bundle, revision, self.home, self.runtime)
                self.assertEqual(before, self.snapshot())

    def test_malformed_incomplete_or_wrong_version_manifest_is_rejected(self):
        revision, original = self.installed_bundle()
        before = self.snapshot()
        mutations = [
            lambda value: value.update(schema="unsupported-v2"),
            lambda value: value.update(source_commit="0" * 40),
            lambda value: value.pop("files_sha256"),
            lambda value: value.update(unexpected=True),
            lambda value: value["files_sha256"].pop(installer.INSTALLER),
            lambda value: value["files_sha256"].update(unexpected="0" * 64),
            lambda value: value["files_sha256"].update({installer.INSTALLER: "not-a-digest"}),
        ]
        for number, mutation in enumerate(mutations):
            with self.subTest(mutation=number):
                bundle = self.base / f"manifest-{number}"
                shutil.copytree(original, bundle)
                path = bundle / installer.MANIFEST
                value = json.loads(path.read_text())
                mutation(value)
                path.write_text(json.dumps(value))
                with self.assertRaises(ValueError):
                    installer.install_bundle(bundle, revision, self.home, self.runtime)
                self.assertEqual(before, self.snapshot())
        with self.assertRaisesRegex(ValueError, "expected revision"):
            installer.install_bundle(original, "f" * 40, self.home, self.runtime)
        self.assertEqual(before, self.snapshot())

    def test_duplicate_json_keys_and_oversized_manifest_are_rejected(self):
        revision, original = self.installed_bundle()
        before = self.snapshot()
        for label, content in (("duplicate", '{"schema":"x","schema":"y"}'),
                               ("oversized", " " * (installer.MAX_MANIFEST_BYTES + 1))):
            bundle = self.base / label
            shutil.copytree(original, bundle)
            (bundle / installer.MANIFEST).write_text(content)
            with self.assertRaises(ValueError):
                installer.install_bundle(bundle, revision, self.home, self.runtime)
            self.assertEqual(before, self.snapshot())

    def test_missing_extra_symlink_hardlink_or_directory_payloads_are_rejected(self):
        revision, original = self.installed_bundle()
        before = self.snapshot()
        for kind in ("missing", "extra", "symlink", "hardlink", "directory", "manifest-link", "root-link"):
            with self.subTest(kind=kind):
                bundle = self.base / kind
                if kind == "root-link":
                    bundle.symlink_to(original, target_is_directory=True)
                else:
                    shutil.copytree(original, bundle)
                    target = bundle / installer.MODULES[0]
                    if kind == "extra":
                        (bundle / "unexpected.txt").write_bytes(target.read_bytes())
                    elif kind == "manifest-link":
                        target = bundle / installer.MANIFEST
                        target.unlink()
                        target.symlink_to(original / installer.MANIFEST)
                    else:
                        target.unlink()
                        if kind == "symlink":
                            target.symlink_to(original / installer.MODULES[0])
                        elif kind == "hardlink":
                            os.link(original / installer.MODULES[0], target)
                        elif kind == "directory":
                            target.mkdir()
                with self.assertRaises((ValueError, OSError)):
                    installer.install_bundle(bundle, revision, self.home, self.runtime)
                self.assertEqual(before, self.snapshot())

    def test_git_symlink_is_rejected_before_export_directory_exists(self):
        path = self.repo / "scripts/is/harness_archived_evidence.py"
        path.unlink()
        path.symlink_to("harness_science_pipeline.py")
        revision = self.revision()
        bundle = self.base / "bundle"
        with self.assertRaisesRegex(ValueError, "regular, not a link"):
            installer.export_bundle(self.repo, revision, bundle)
        self.assertFalse(bundle.exists())

    def test_unit_command_changes_are_rejected_even_with_updated_manifest_hash(self):
        revision, bundle = self.installed_bundle()
        before = self.snapshot()
        path = bundle / installer.UNITS[0]
        path.write_text(path.read_text().replace("ExecStart=", "ExecStartPre=/usr/bin/true\nExecStart="))
        manifest_path = bundle / installer.MANIFEST
        manifest = json.loads(manifest_path.read_text())
        manifest["files_sha256"][path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
        manifest_path.write_text(json.dumps(manifest))
        with self.assertRaisesRegex(ValueError, "unit directives"):
            installer.install_bundle(bundle, revision, self.home, self.runtime)
        self.assertEqual(before, self.snapshot())

    def test_bundle_import_failure_preserves_existing_installation(self):
        revision, bundle = self.installed_bundle()
        before = self.snapshot()
        path = self.repo / "scripts/is/harness_archived_evidence.py"
        path.write_text("raise ImportError('synthetic bundle failure')\n")
        changed_revision = self.revision()
        changed = self.base / "changed"
        installer.export_bundle(self.repo, changed_revision, changed)
        with self.assertRaises(subprocess.CalledProcessError):
            installer.install_bundle(changed, changed_revision, self.home, self.runtime)
        self.assertEqual(before, self.snapshot())

    def test_shared_execution_lock_blocks_install_and_missing_path_is_not_created(self):
        revision, bundle = self.installed_bundle()
        before = self.snapshot()
        execution_lock = self.runtime / "shared-execution.lock"
        with execution_lock.open("a+") as active:
            fcntl.flock(active.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
            with self.assertRaisesRegex(RuntimeError, "Shared research execution lock is busy"):
                installer.install_bundle(bundle, revision, self.home, self.runtime,
                                         execution_lock=execution_lock)
        self.assertEqual(before, self.snapshot())
        missing = self.base / "must-not-create-results/shared.lock"
        with self.assertRaises(FileNotFoundError):
            installer.install_bundle(bundle, revision, self.home, self.runtime, execution_lock=missing)
        self.assertFalse(missing.parent.exists())
        self.assertEqual(before, self.snapshot())

    def test_runtime_is_rechecked_under_both_locks_before_replacement(self):
        revision, bundle = self.installed_bundle()
        execution_lock = self.runtime / "shared.lock"
        execution_lock.touch()
        original_run = subprocess.run
        checked_units = []

        def run(command, **kwargs):
            if command[0] in {"systemctl", "pgrep"}:
                for path in (self.runtime / "is-harness-queue.lock", execution_lock):
                    with path.open("rb") as handle:
                        with self.assertRaises(BlockingIOError):
                            fcntl.flock(handle.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
                if command[0] == "systemctl":
                    checked_units.append(command[3])
                    return subprocess.CompletedProcess(command, 0, "LoadState=loaded\nActiveState=inactive\n", "")
                return subprocess.CompletedProcess(command, 1, "0\n", "")
            return original_run(command, **kwargs)

        with mock.patch.object(installer.subprocess, "run", side_effect=run):
            installer.install_bundle(bundle, revision, self.home, self.runtime,
                                     execution_lock=execution_lock, require_inactive=True)
        self.assertEqual(list(installer.RUNTIME_UNITS), checked_units)

    def test_active_unit_or_orphan_process_holds_existing_installation(self):
        revision, bundle = self.installed_bundle()
        before = self.snapshot()
        execution_lock = self.runtime / "shared.lock"
        execution_lock.touch()
        original_run = subprocess.run
        for active in ("unit", "process"):
            def run(command, **kwargs):
                if command[0] == "systemctl":
                    state = "active" if active == "unit" else "inactive"
                    return subprocess.CompletedProcess(command, 0, f"LoadState=loaded\nActiveState={state}\n", "")
                if command[0] == "pgrep":
                    return subprocess.CompletedProcess(command, 0, "1\n", "")
                return original_run(command, **kwargs)
            with self.subTest(active=active), mock.patch.object(installer.subprocess, "run", side_effect=run):
                with self.assertRaisesRegex(RuntimeError, "active or unverifiable"):
                    installer.install_bundle(bundle, revision, self.home, self.runtime,
                                             execution_lock=execution_lock, require_inactive=True)
            self.assertEqual(before, self.snapshot())

    def test_replace_failure_rolls_back_modules_units_modes_and_receipt(self):
        revision, original = self.installed_bundle()
        old_module = self.home / "is-harness" / installer.MODULES[0]
        old_module.write_bytes(b"")  # Preserve even a previously truncated file on rollback.
        old_module.chmod(0o640)
        before = self.snapshot()
        helper = self.repo / "scripts/is/harness_archived_evidence.py"
        helper.write_bytes(helper.read_bytes() + b"\nREVISION_TWO = True\n")
        changed_revision = self.revision()
        changed = self.base / "changed"
        installer.export_bundle(self.repo, changed_revision, changed)
        original_replace = os.replace
        for name in (installer.MODULES[2], installer.UNITS[0], "installed_revision.json"):
            failed = False

            def replace(source, destination):
                nonlocal failed
                if Path(destination).name == name and not failed:
                    failed = True
                    raise OSError("synthetic replacement failure")
                return original_replace(source, destination)

            with self.subTest(target=name), mock.patch.object(installer.os, "replace", side_effect=replace):
                with self.assertRaisesRegex(OSError, "synthetic replacement failure"):
                    installer.install_bundle(changed, changed_revision, self.home, self.runtime)
            self.assertTrue(failed)
            self.assertEqual(before, self.snapshot())

    def test_export_cli_allows_root_but_install_cli_rejects_root_and_mixed_sources(self):
        revision = self.revision()
        bundle = self.base / "cli-export"
        with mock.patch.object(installer.os, "geteuid", return_value=0), \
                mock.patch("sys.stdout", new_callable=io.StringIO):
            installer.main(["--repo-root", str(self.repo), "--revision", revision,
                            "--export-bundle", str(bundle)])
        self.assertTrue((bundle / installer.MANIFEST).is_file())
        invalid = [
            ["--bundle-dir", str(bundle)],
            ["--repo-root", str(self.repo), "--bundle-dir", str(bundle)],
            ["--bundle-dir", str(bundle), "--export-bundle", str(self.base / "invalid")],
        ]
        for arguments in invalid:
            with mock.patch.object(installer.os, "geteuid", return_value=0), \
                    mock.patch("sys.stderr", new_callable=io.StringIO):
                with self.assertRaises(SystemExit):
                    installer.main(["--revision", revision, *arguments])
        self.assertFalse(self.home.exists())

    def test_cli_rejects_execution_lock_environment_mismatch_before_install(self):
        with mock.patch.object(installer.os, "geteuid", return_value=1000), \
                mock.patch.dict(installer.os.environ, {"IS_HARNESS_EXECUTION_LOCK": "/configured/custom.lock"}), \
                mock.patch.object(installer, "install_bundle") as install, \
                mock.patch("sys.stderr", new_callable=io.StringIO) as errors:
            with self.assertRaises(SystemExit):
                installer.main(["--revision", "a" * 40, "--bundle-dir", str(self.base / "bundle"),
                                "--execution-lock", "/different/default.lock"])
        self.assertIn("does not match", errors.getvalue())
        install.assert_not_called()


if __name__ == "__main__":
    unittest.main()
