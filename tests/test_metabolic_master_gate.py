import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest import mock


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/metabolic_resilience/drive_v2/91_master_analysis_gate.py"
FIXTURE = json.loads((ROOT / "tests/fixtures/metabolic_master_gate/inventory_inputs.json").read_text())
OUTPUT = "results/metabolic_resilience/MASTER_ANALYSIS_GATE.json"


class MasterAnalysisGateTests(unittest.TestCase):
    def populate(self, root, invalid=False):
        paths = []
        for relative in FIXTURE["required_paths"]:
            path = root / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            content = FIXTURE["invalid_input"] if invalid else FIXTURE[
                "json_input" if path.suffix == ".json" else "tsv_input"
            ]
            path.write_text(content, encoding="utf-8")
            paths.append(path)
        return paths

    def run_gate(self, root, use_cli=True, extra_env=None):
        env = os.environ.copy()
        env["IS_ANALYSIS_ROOT"] = str(root)
        env.update(extra_env or {})
        command = [sys.executable, str(SCRIPT)]
        if use_cli:
            command += ["--root", str(root)]
        completed = subprocess.run(command, env=env, capture_output=True, text=True, timeout=10)
        self.assertEqual(completed.returncode, 2, completed.stdout + completed.stderr)
        payload = json.loads(completed.stdout)
        self.assertEqual(payload, json.loads((root / OUTPUT).read_text()))
        self.assertEqual(payload["scientific_status"], "HOLD")
        self.assertIs(payload["scientific_validation_performed"], False)
        self.assertNotIn("PASS_CORE_READY", completed.stdout)
        return payload

    def test_complete_inventory_does_not_certify_science(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self.populate(root)
            payload = self.run_gate(root)
            self.assertEqual(payload["blocking_missing"], 0)
            self.assertEqual(payload["status"], "FILES_PRESENT_NOT_SCIENTIFICALLY_VALIDATED")
            self.assertEqual(sum(row["status"] == "EXISTS" for row in payload["checks"]), 9)

    def test_internal_hold_and_review_remain_unresolved_and_inputs_unchanged(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            paths = self.populate(root)
            next(p for p in paths if p.suffix == ".json").write_text(FIXTURE["hold_json"])
            before = {p: p.read_bytes() for p in paths}
            payload = self.run_gate(root)
            self.assertEqual(payload["status"], "FILES_PRESENT_NOT_SCIENTIFICALLY_VALIDATED")
            self.assertEqual(before, {p: p.read_bytes() for p in paths})

    def test_nonempty_malformed_inputs_are_only_inventory(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            self.populate(root, invalid=True)
            payload = self.run_gate(root)
            self.assertEqual(payload["status"], "FILES_PRESENT_NOT_SCIENTIFICALLY_VALIDATED")
            self.assertEqual(payload["validation_scope"], "nonempty_regular_files_only")

    def test_missing_core_files_block(self):
        with tempfile.TemporaryDirectory() as tmp:
            payload = self.run_gate(Path(tmp))
            self.assertEqual(payload["status"], "HOLD_CORE_INCOMPLETE")
            self.assertEqual(payload["blocking_missing"], 9)

    def test_empty_file_and_directory_do_not_count_as_outputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            paths = self.populate(root)
            paths[0].write_bytes(b"")
            paths[1].unlink()
            paths[1].mkdir()
            payload = self.run_gate(root)
            self.assertEqual(payload["status"], "HOLD_CORE_INCOMPLETE")
            self.assertEqual(payload["blocking_missing"], 2)

    def test_cli_root_overrides_environment_and_environment_remains_supported(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "requested"
            other = Path(tmp) / "environment"
            self.populate(root)
            self.run_gate(root, extra_env={"IS_ANALYSIS_ROOT": str(other)})
            self.assertFalse(other.exists())
            self.run_gate(root, use_cli=False)

    def test_import_does_not_write_or_exit(self):
        with tempfile.TemporaryDirectory() as tmp, mock.patch.dict(os.environ, {"IS_ANALYSIS_ROOT": tmp}):
            spec = importlib.util.spec_from_file_location("metabolic_master_gate", SCRIPT)
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            self.assertTrue(callable(module.main))
            self.assertEqual(list(Path(tmp).iterdir()), [])


if __name__ == "__main__":
    unittest.main()
