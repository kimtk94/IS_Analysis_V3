"""Contract tests for fail-closed AUTO_NEXT decisions; no production operations."""
import importlib.util
import pathlib
import sys
import tempfile
import types
import unittest
from unittest.mock import patch

SCRIPT = pathlib.Path(__file__).resolve().parents[1] / "scripts" / "is" / "auto_next.py"
stub = types.ModuleType("github_issue_queue_worker")
stub.run_job = lambda op: "stubbed audit " + op
sys.modules.setdefault("github_issue_queue_worker", stub)
spec = importlib.util.spec_from_file_location("auto_next_under_test", SCRIPT)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class AutoNextTests(unittest.TestCase):
    def test_starts_dry_run_without_creating_state(self):
        with tempfile.TemporaryDirectory() as root:
            state = pathlib.Path(root) / "state.json"
            result = module.execute(state)
            self.assertEqual(result["status"], "DRY_RUN")
            self.assertEqual(result["next"]["phase"], "P1")
            self.assertFalse(state.exists())

    def test_blocked_step_stops_later_steps(self):
        for stopped in ("RUNNING", "REVIEW_REQUIRED", "BLOCKED"):
            with self.subTest(stopped=stopped):
                state = {"tasks": {"P1-inventory": {"status": stopped}}}
                step, decision = module.choose(state)
                self.assertEqual(step["id"], "P1-inventory")
                self.assertEqual(decision, "BLOCKED_" + stopped)

    def test_dependency_fail_closed(self):
        state = {"tasks": {"P1-inventory": {"status": "PASS"},
                           "P2-headers": {"status": "FAIL"}}}
        step, decision = module.choose(state)
        self.assertEqual(step["id"], "P2-headers")
        self.assertEqual(decision, "READY")

    def test_one_step_per_invocation_and_finish(self):
        with tempfile.TemporaryDirectory() as root:
            state = pathlib.Path(root) / "state.json"
            with patch.object(module, "run_job", return_value="ok") as runner:
                for idx, expected in enumerate(module.STEPS, 1):
                    result = module.execute(state, run=True)
                    self.assertEqual(result["status"], "PASS")
                    self.assertEqual(result["task"], expected["id"])
                    self.assertEqual(runner.call_count, idx)
                self.assertEqual(module.execute(state, run=True)["status"], "ALL_TASKS_COMPLETED")

    def test_error_never_advances(self):
        with tempfile.TemporaryDirectory() as root:
            state = pathlib.Path(root) / "state.json"
            with patch.object(module, "run_job", side_effect=RuntimeError("failure")) as runner:
                self.assertEqual(module.execute(state, run=True)["status"], "REVIEW_REQUIRED")
                self.assertEqual(module.execute(state, run=True)["status"], "BLOCKED_REVIEW_REQUIRED")
                self.assertEqual(runner.call_count, 1)

    def test_running_from_crash_never_repeats(self):
        with tempfile.TemporaryDirectory() as root:
            state = pathlib.Path(root) / "state.json"
            module.write_atomic(state, {"schema": "is-auto-next-v1",
                                        "history": [],
                                        "tasks": {"P1-inventory": {"status": "RUNNING"}}})
            with patch.object(module, "run_job") as runner:
                self.assertEqual(module.execute(state, run=True)["status"], "BLOCKED_RUNNING")
                runner.assert_not_called()

    def test_rejects_corrupt_schema(self):
        with tempfile.TemporaryDirectory() as root:
            state = pathlib.Path(root) / "state.json"
            state.write_text('{"schema":"unknown","tasks":{}}')
            with self.assertRaises(ValueError):
                module.execute(state, run=True)


if __name__ == "__main__":
    unittest.main()
