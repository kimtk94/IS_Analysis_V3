"""Tests for the shared execution lock between issue queue and AUTO_NEXT."""
import fcntl
import importlib.util
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "is" / "github_issue_queue_worker.py"
spec = importlib.util.spec_from_file_location("is_queue_worker_lock_test", SCRIPT)
worker = importlib.util.module_from_spec(spec)
spec.loader.exec_module(worker)


class SharedExecutionLockTests(unittest.TestCase):
    def test_executes_under_uncontended_lock(self):
        with tempfile.TemporaryDirectory() as root:
            path = Path(root) / "execution.lock"
            with patch.dict(os.environ, {"IS_HARNESS_EXECUTION_LOCK": str(path)}), \
                 patch.object(worker, "_run_job_unlocked", return_value="OK") as operation:
                self.assertEqual(worker.run_job("inventory"), "OK")
                operation.assert_called_once_with("inventory")

    def test_rejects_concurrent_job(self):
        with tempfile.TemporaryDirectory() as root:
            path = Path(root) / "execution.lock"
            with path.open("a+") as competing:
                fcntl.flock(competing.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
                with patch.dict(os.environ, {"IS_HARNESS_EXECUTION_LOCK": str(path)}), \
                     patch.object(worker, "_run_job_unlocked") as operation:
                    with self.assertRaisesRegex(RuntimeError, "SHARED_EXECUTION_LOCK_BUSY"):
                        worker.run_job("research-refresh")
                    operation.assert_not_called()


if __name__ == "__main__":
    unittest.main()
